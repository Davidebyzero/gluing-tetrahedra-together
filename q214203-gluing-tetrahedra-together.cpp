#include <stdio.h>
#include <fcntl.h>
#include <iostream>
#include <string.h>
#include <array>
#include <vector>
#include <chrono>

//#define USE_GMP
#define MEMORY_POOL_INITIAL_SIZE (64uLL * 1024)  // in bytes; if the goal is to use more than half of available RAM, this must be preallocated at full expected size
#define MEMORY_POOL_GROW_RATIO 1/16  // what proportion of the memory size to grow it by when more space is needed
#define HASH_TABLE_RATIO 6
#define SHOW_PROGRESS 16  // if defined, show progress starting at this term
//#define PRINT_POLYTETS // requires USE_GMP

#define MAXIMUM_TETCOUNT 17 // 28

#define WRITE_TO_FILES
#define RESUME_FROM_FILE

#define FILE_CHUNK_SIZE (1uLL << 30)  // needs to be less than 1<<31

#ifdef USE_GMP
#   include <gmp.h>
#   if GMP_NUMB_BITS != 64
#       error This is hard-coded for 64-bit limbs
#   endif
#endif

template <typename T, size_t N> size_t countof( T ( & arr )[ N ] ) {return std::extent< T[ N ] >::value;}

auto startTime = std::chrono::steady_clock::now();

void quitMemory()
{
    auto currentTime = std::chrono::steady_clock::now();
    std::cerr << "Out of memory [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms]" << std::endl;
    exit(-1);
}

typedef uint8_t TetIndex;

#if MAXIMUM_TETCOUNT > 17
typedef uint64_t HashIndex;
#else
typedef uint32_t HashIndex;
#endif

#if MAXIMUM_TETCOUNT > 23
typedef unsigned __int128 CompressedPolytetBits; // Can handle up to 44 terms, while a 64-bit size_t can only handle up to about 28 terms
#else
typedef uint64_t          CompressedPolytetBits; // Can handle up to 23 terms
#endif

typedef uint64_t CompressedSubpolytet; // compressed in branchless format; for caching known-overlapping subpolytet paths

#ifdef USE_GMP
class Tetrahedron
{
public:
    mpz_t t[4][3];
    Tetrahedron(const Tetrahedron &_t) : Tetrahedron()
    {
        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_set(t[p][d], _t.t[p][d]);
    }
    Tetrahedron()
    {
        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_init(t[p][d]);
    }
    ~Tetrahedron()
    {
        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_clear(t[p][d]);
    }
};
#else
typedef __int128 Coord;
typedef std::array<Coord, 3> Coord3;
typedef std::array<Coord3, 4> Tetrahedron;
#endif

class Tet
{
    void initFaces()
    {
        faceAttached[0] = NULL; // t[0],t[1],t[2]
        faceAttached[1] = NULL; // t[0],t[1],t[3]
        faceAttached[2] = NULL; // t[0],t[2],t[3]
        faceAttached[3] = NULL; // t[1],t[2],t[3]
    }
    void tagSkipOverlapCheckHelper(int depth, const struct RotationTable *thisRotationTable, int faceRotation = 0, CompressedSubpolytet curCompressedPath = 0, CompressedSubpolytet trit = 1);
public:
    Tetrahedron t;
    Tet    *faceAttached    [4];
    uint8_t faceAttachedFace[4];
    TetIndex index; // 1-based; 0=unassigned
    CompressedSubpolytet compressedPath; // 0 = skip overlap check
    Tet(                    ) : t( ) {initFaces();}
    Tet(const Tetrahedron &t) : t(t) {initFaces();}
    void assignIndex(TetIndex &nextIndex)
    {
        if (index == 0)
            index = nextIndex++;
    }
    void tagSkipOverlapCheck(int depth);
};
class Polytet : public std::vector<Tet>
{
public:
    TetIndex nextIndex;
    void resetIndexing(size_t first) // This function should not be called if the polytet hasn't yet been populated
    {
        
        for (auto t=begin(); t!=end(); ++t)
            t->index = 0;
        (*this)[first].index = 1;
        nextIndex = 2;
    }
};

// vertex indices of faces with identical chirality
static const int tetrahedronFaces[4][4] =
{
    {0, 1, 2, 3},
    {0, 3, 1, 2},
    {0, 2, 3, 1},
    {1, 3, 2, 0},
};

const int tetrahedronEdges_face3 = 3;
static const int tetrahedronEdges[][2] =
{
    {0, 1},
    {0, 2},
    {0, 3},
    // The following edges, on face[3], can be skipped when the new tetrahedron is one whose edges are being checked
    {1, 2},
    {1, 3},
    {2, 3},
};

static const int faceRotateReflect[2][3][3] =
{
    {
        {0, 1, 2},
        {1, 2, 0},
        {2, 0, 1},
    },
    {
        {1, 2, 0},
        {0, 1, 2},
        {2, 0, 1},
    },
};

struct RotationTable
{
    uint8_t faceMap[3];
    uint8_t rotation[3];
};
RotationTable rotationTable[4];

void initLookupTables()
{
#if 0 // initializing "tetrahedronFaces" at runtime results in a slower executable with GCC
    for (int i=0; i<4; i++)
    {
        int oppositeVertex = 3 - i;
        for (int j=0; j<3; j++)
            tetrahedronFaces[i][j] = j + (j >= oppositeVertex ? 1 : 0);
        if (i & 1)
            std::swap(tetrahedronFaces[i][1], tetrahedronFaces[i][2]);
        tetrahedronFaces[i][3] = oppositeVertex;
    }
#endif
    for (int i=0; i<4; i++)
    {
        int vertexMap[4];
        for (int j=0; j<4; j++)
        {
            int vertexNum = tetrahedronFaces[3][j];
            vertexMap[vertexNum] = tetrahedronFaces[i][j];
            if (vertexNum != 0)
                rotationTable[i].faceMap[3 - vertexNum] = 3 - tetrahedronFaces[i][j];
        }
        for (int rotatedFaceNum=0; rotatedFaceNum<3; rotatedFaceNum++)
        {
            int faceNum = rotationTable[i].faceMap[rotatedFaceNum];
            int rotation = 0;
            while (vertexMap[tetrahedronFaces[rotatedFaceNum][rotation]] != tetrahedronFaces[faceNum][0])
                rotation++;
            rotationTable[i].rotation[rotatedFaceNum] = rotation;
        }
    }
}

#ifndef USE_GMP
Coord3 operator+(const Coord3 &a, const Coord3 &b)
{
    Coord3 c;
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
    return c;
}
Coord3 operator-(const Coord3 &a, const Coord3 &b)
{
    Coord3 c;
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
    return c;
}
Coord3 operator*(const Coord3 &a, const Coord b)
{
    Coord3 c;
    c[0] = a[0] * b;
    c[1] = a[1] * b;
    c[2] = a[2] * b;
    return c;
}
Coord dot(const Coord3 &a, const Coord3 &b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
#endif

// This class must be used in exactly the way it is in this program. That is, setA() and setB() must be called first, before operator(),
// and point to tetrahedrons in the same Polytet. The set(A) tetrahedron must be the last one newly attached to the polytet. If operator()
// returns "true", both setA() and setB() must be called before the next operator() call.
class TetrahedronOverlap
{
    const Tet *a;
    const Tet *b;
public:
    void setA(const Tet &x) {a = &x;}
    void setB(const Tet &x) {b = &x;}
#ifdef USE_GMP
private:
    mpz_t intersectNumerator, intersectDenominator, uNumerator, vNumerator, uvDenominator, uvNumeratorSum;
    mpz_t center[3], normal[3], tmp[3], p0p1[3], intersectionPoint[3], delta[3], edge1[3], edge2[3];
    mpz_t triangle[3][3];
    void dot(mpz_t &result, const mpz_t _a[3], const mpz_t _b[3])
    {
        mpz_mul   (result, _a[0], _b[0]);
        mpz_addmul(result, _a[1], _b[1]);
        mpz_addmul(result, _a[2], _b[2]);
    }
    template <void (*MPZ_CALL)(mpz_t x), void (*MPZ_CALLS)(mpz_t x, ...)> void mpz_initOrClear()
    {
        MPZ_CALLS(intersectNumerator, intersectDenominator, uNumerator, vNumerator, uvDenominator, uvNumeratorSum, NULL);
        for (int d=0; d<3; d++)
        {
            MPZ_CALL(center[d]);
            MPZ_CALL(normal[d]);
            MPZ_CALL(tmp[d]);
            MPZ_CALL(p0p1[d]);
            MPZ_CALL(intersectionPoint[d]);
            MPZ_CALL(delta[d]);
            MPZ_CALL(edge1[d]);
            MPZ_CALL(edge2[d]);
        }
        for (int p=0; p<3; p++)
            for (int d=0; d<3; d++)
                MPZ_CALL(triangle[p][d]);
    }
public:
     TetrahedronOverlap() {mpz_initOrClear<mpz_init , mpz_inits >();}
    ~TetrahedronOverlap() {mpz_initOrClear<mpz_clear, mpz_clears>();}
    bool operator()(const mpz_t maximalTouchingSqrDistance)
    {
        // Skip the longer overlap checking algorithm if the two tetrahedrons' centers are sufficiently separated.
        // Get center of tetrahedron "a" by averaging its vertices' coordinates, without dividing by 4.
        for (int d=0; d<3; d++)
            mpz_set(center[d], a->t.t[0][d]);
        for (int p=1; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_add(center[d], center[d], a->t.t[p][d]);
        // Subtract center of tetrahedron "b" by averaging its vertices' coordinates, without dividing by 4.
        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_sub(center[d], center[d], b->t.t[p][d]);
        // Compare the sum of the squares of the orthogonal distances against the threshold squared distance.
        dot(tmp[0], center, center);
        if (mpz_cmp(tmp[0], maximalTouchingSqrDistance) >= 0)
            return false;

        // In some circumstances, two tetrahedrons can intersect such that only the edges of one tetrahedron intersect with the
        // faces of the other, and not the other way around. So we need to check both.
        int maxEdgeNum = tetrahedronEdges_face3; // when "a" is the newly added tetrahedron, we can enable an optimization
        for (int swapped=0; swapped<2; swapped++)
        {
            // Take advantage of the fact that the two tetrahedrons are regular and congruent, and
            // just check if any edge of tetrahedron "a" intersects with any face of tetrahedron "b".
            // Don't count it if only the endpoint of an edge intersects.
            for (int edgeNum=0; edgeNum<maxEdgeNum; edgeNum++)
            {
                const mpz_t *p0 = a->t.t[tetrahedronEdges[edgeNum][0]];
                const mpz_t *p1 = a->t.t[tetrahedronEdges[edgeNum][1]];
                for (int faceNum=0; faceNum<4; faceNum++)
                {
                    // For speed, process only one out of every pair of attached faces (which share the exact same 3 vertices).
                    // This will still process the very first attached face twice, since that is a face[3] attached to another face[3],
                    // but it's probably not worth the extra machinery that would be necessary to special-case that.
                    if (swapped ? faceNum==3
                                : faceNum!=3 && b->faceAttached[faceNum])
                        continue;
                    const mpz_t *normalizedTetrahedron[4][3]; // first 3 points are the face, and the 4th point is for calculating the normal
                    for (int i=0; i<4; i++)
                        for (int d=0; d<3; d++)
                            normalizedTetrahedron[i][d] = &b->t.t[tetrahedronFaces[faceNum][i]][d];
                    // Center coordinates will be multiplied by 3 compared to original coordinates.
                    // Get center of face by averaging its vertices' coordinates; the
                    // division by 3 is implied by omitting the multiplication by 3.
                    for (int d=0; d<3; d++)
                        mpz_set(center[d], *(normalizedTetrahedron[0][d]));
                    for (int p=1; p<3; p++)
                        for (int d=0; d<3; d++)
                            mpz_add(center[d], center[d], *(normalizedTetrahedron[p][d]));
                    for (int d=0; d<3; d++)
                    {
                        mpz_neg(normal[d], center[d]);
                        mpz_addmul_ui(normal[d], *(normalizedTetrahedron[3][d]), 3);
                    }
                    for (int d=0; d<3; d++)
                    {
                        mpz_sub(tmp [d], *(normalizedTetrahedron[0][d]), p0[d]);
                        mpz_sub(p0p1[d],                         p1[d],  p0[d]);
                    }
                    dot(intersectNumerator  , normal, tmp);
                    dot(intersectDenominator, normal, p0p1);
                    int cmp = mpz_cmp_ui(intersectDenominator, 0);
                    if (cmp == 0)
                        continue; // edge is parallel to face, which we don't count as an overlap
                    if (cmp < 0)
                    {
                        mpz_neg(intersectNumerator  , intersectNumerator  );
                        mpz_neg(intersectDenominator, intersectDenominator);
                    }
                    if (mpz_cmp_ui(intersectNumerator, 0) <= 0 || mpz_cmp(intersectNumerator, intersectDenominator) >= 0)
                        continue;
                    // These coordinates are all multiplied by intersectDenominator
                    for (int d=0; d<3; d++)
                    {
                        mpz_mul(intersectionPoint[d], p0[d], intersectDenominator);
                        mpz_addmul(intersectionPoint[d], p0p1[d], intersectNumerator);
                    }
                    for (int i=0; i<3; i++)
                        for (int d=0; d<3; d++)
                            mpz_mul(triangle[i][d], *(normalizedTetrahedron[i][d]), intersectDenominator);
                    // Check if the intersection point is inside the triangle
                    for (int d=0; d<3; d++)
                    {
                        mpz_sub(delta[d], intersectionPoint[d], triangle[0][d]);
                        mpz_sub(edge1[d], triangle[1][d]      , triangle[0][d]);
                        mpz_sub(edge2[d], triangle[2][d]      , triangle[0][d]);
                    }
                    mpz_mul   (uNumerator, delta[1], edge2[0]);
                    mpz_submul(uNumerator, delta[0], edge2[1]);
                    mpz_mul   (vNumerator, delta[0], edge1[1]);
                    mpz_submul(vNumerator, delta[1], edge1[0]);
                    mpz_mul   (uvDenominator, edge1[1], edge2[0]);
                    mpz_submul(uvDenominator, edge1[0], edge2[1]);
                    cmp = mpz_cmp_ui(uvDenominator, 0);
                    if (cmp == 0)
                        continue;
                    if (cmp < 0)
                    {
                        mpz_neg(uNumerator, uNumerator);
                        mpz_neg(vNumerator, vNumerator);
                        mpz_neg(uvDenominator, uvDenominator);
                    }
                    if (mpz_cmp_ui(uNumerator, 0) <= 0 || mpz_cmp_ui(vNumerator, 0) <= 0)
                        continue;
                    mpz_add(uvNumeratorSum, uNumerator, vNumerator);
                    if (mpz_cmp(uvNumeratorSum, uvDenominator) < 0)
                        return true;
                }
            }
            std::swap(a, b);
            maxEdgeNum = countof(tetrahedronEdges); // disable the optimization for when "b" is the newly added tetrahedron
        }
        return false;
    }
#else
public:
    bool operator()(const Coord &maximalTouchingSqrDistance)
    {
        // Skip the longer overlap checking algorithm if the two tetrahedrons' centers are sufficiently separated.
        {
            // Get center of tetrahedron "a" by averaging its vertices' coordinates, without dividing by 4.
            Coord3 center = {0, 0, 0};
            for (int p=0; p<4; p++)
                for (int d=0; d<3; d++)
                    center[d] += a->t[p][d];
            // Subtract center of tetrahedron "b" by averaging its vertices' coordinates, without dividing by 4.
            for (int p=0; p<4; p++)
                for (int d=0; d<3; d++)
                    center[d] -= b->t[p][d];
            // Compare the sum of the squares of the orthogonal distances against the threshold squared distance.
            if (dot(center, center) >= maximalTouchingSqrDistance)
                return false;
        }
        // In some circumstances, two tetrahedrons can intersect such that only the edges of one tetrahedron intersect with the
        // faces of the other, and not the other way around. So we need to check both.
        int maxEdgeNum = tetrahedronEdges_face3; // when "a" is the newly added tetrahedron, we can enable an optimization
        for (int swapped=0; swapped<2; swapped++)
        {
            // Take advantage of the fact that the two tetrahedrons are regular and congruent, and
            // just check if any edge of tetrahedron "a" intersects with any face of tetrahedron "b".
            // Don't count it if only the endpoint of an edge intersects.
            for (int edgeNum=0; edgeNum<maxEdgeNum; edgeNum++)
            {
                Coord3 p0 = a->t[tetrahedronEdges[edgeNum][0]];
                Coord3 p1 = a->t[tetrahedronEdges[edgeNum][1]];
                for (int faceNum=0; faceNum<4; faceNum++)
                {
                    // For speed, process only one out of every pair of attached faces (which share the exact same 3 vertices).
                    // This will still process the very first attached face twice, since that is a face[3] attached to another face[3],
                    // but it's probably not worth the extra machinery that would be necessary to special-case that.
                    if (swapped ? faceNum==3
                                : faceNum!=3 && b->faceAttached[faceNum])
                        continue;
                    Tetrahedron normalizedTetrahedron; // first 3 points are the face, and the 4th point is for calculating the normal
                    for (int i=0; i<4; i++)
                        normalizedTetrahedron[i] = b->t[tetrahedronFaces[faceNum][i]];
                    Coord3 center = {{0, 0, 0}}; // multiplied by 3 compared to original coordinates
                    // Get center of face by averaging its vertices' coordinates; the
                    // division by 3 is implied by omitting the multiplication by 3.
                    for (int p=0; p<3; p++)
                        for (int d=0; d<3; d++)
                            center[d] += normalizedTetrahedron[p][d];
                    Coord3 normal;
                    for (int d=0; d<3; d++)
                        normal[d] = normalizedTetrahedron[3][d] * 3 - center[d];
                    Coord intersectNumerator   = dot(normal, normalizedTetrahedron[0] - p0);
                    Coord intersectDenominator = dot(normal,                       p1 - p0);
                    if (intersectDenominator == 0)
                        continue; // edge is parallel to face, which we don't count as an overlap
                    if (intersectDenominator < 0)
                    {
                        intersectNumerator   = -intersectNumerator;
                        intersectDenominator = -intersectDenominator;
                    }
                    if (intersectNumerator <= 0 || intersectNumerator >= intersectDenominator)
                        continue;
                    // These coordinates are all multiplied by intersectDenominator
                    Coord3 intersectionPoint = p0 * intersectDenominator + (p1 - p0) * intersectNumerator;
                    Coord3 triangle[3];
                    for (int i=0; i<3; i++)
                        triangle[i] = normalizedTetrahedron[i] * intersectDenominator;
                    // Check if the intersection point is inside the triangle
                    Coord3 delta = intersectionPoint - triangle[0];
                    Coord3 edge1 = triangle[1]       - triangle[0];
                    Coord3 edge2 = triangle[2]       - triangle[0];
                    Coord uNumerator = delta[1]*edge2[0] - delta[0]*edge2[1];
                    Coord vNumerator = delta[0]*edge1[1] - delta[1]*edge1[0];
                    Coord uvDenominator = edge1[1]*edge2[0] - edge1[0]*edge2[1];
                    if (uvDenominator == 0)
                        continue;
                    if (uvDenominator < 0)
                    {
                        uNumerator = -uNumerator;
                        vNumerator = -vNumerator;
                        uvDenominator = -uvDenominator;
                    }
                    if (uNumerator <= 0 || vNumerator <= 0)
                        continue;
                    if (uNumerator + vNumerator < uvDenominator)
                        return true;
                }
            }
            std::swap(a, b);
            maxEdgeNum = countof(tetrahedronEdges); // disable the optimization for when "b" is the newly added tetrahedron
        }
        return false;
    }
#endif
};

void attachNewTet(Tet &t, Tet &tetToAttachTo, const int faceNum)
{
#ifdef USE_GMP
    mpz_t *newVertex = t.t.t[0];
    // Get center of face by averaging its vertices' coordinates.
    for (int d=0; d<3; d++)
        mpz_set(newVertex[d], tetToAttachTo.t.t[tetrahedronFaces[faceNum][0]][d]);
    for (int p=1; p<3; p++)
        for (int d=0; d<3; d++)
            mpz_add(newVertex[d], newVertex[d], tetToAttachTo.t.t[tetrahedronFaces[faceNum][p]][d]);
    // Finalize the new vertex
    for (int d=0; d<3; d++)
    {
        mpz_div_ui(newVertex[d], newVertex[d], 3);
        mpz_mul_ui(newVertex[d], newVertex[d], 2);
        mpz_sub   (newVertex[d], newVertex[d], tetToAttachTo.t.t[3 - faceNum][d]);
    }
#else
    Coord3 &newVertex = t.t[0];
    newVertex = {{0, 0, 0}};
    // Get center of face by averaging its vertices' coordinates.
    for (int p=0; p<4; p++)
    {
        if (p == 3 - faceNum)
            continue;
        for (int d=0; d<3; d++)
            newVertex[d] += tetToAttachTo.t[p][d];
    }
    // Finalize the new vertex
    for (int d=0; d<3; d++)
        newVertex[d] = newVertex[d]/3 * 2 - tetToAttachTo.t[3 - faceNum][d];
#endif
    // Copy the other vertices
    for (int p=0; p<3; p++)
    {
        int p1 = tetrahedronFaces[faceNum][p];
        for (int d=0; d<3; d++)
#ifdef USE_GMP
            mpz_set(t.t.t[1+p][d], tetToAttachTo.t.t[p1][d]);
#else
            t.t[1+p][d] = tetToAttachTo.t[p1][d];
#endif
    }
    t.faceAttached    [0] = NULL;
    t.faceAttached    [1] = NULL;
    t.faceAttached    [2] = NULL;
    t.faceAttached    [3] = &tetToAttachTo;
    t.faceAttachedFace[3] = faceNum;
    tetToAttachTo.faceAttached    [faceNum] = &t;
    tetToAttachTo.faceAttachedFace[faceNum] = 3;
}

// First two tetrahedrons are implied. Each element is a subsequent tetrahedron, with the value indicating where
// it's attached. The lower 2 bits indicate which face (can only have 3 different values, because at least 1 face
// will always already be attached). The remaining bits indicate which tetrahedron (which can never be zero,
// because that one is attached implicitly).
class CompressedPolytet
{
    void uncompressHelper(Polytet &polytet, TetIndex index, TetIndex &nextIndex)
    {
        for (int faceNum=0; faceNum<3; faceNum++)
        {
            if (value & ((CompressedPolytetBits)1 << ((index - 1) * 3 + faceNum)))
            {
                TetIndex thisIndex = nextIndex++;
                attachNewTet(polytet[thisIndex], polytet[index], faceNum);
                uncompressHelper(polytet, thisIndex, nextIndex);
            }
        }
    }
public:
    CompressedPolytetBits value;
    CompressedPolytet() : value(0) {}
    void append(Polytet &polytet, Tet &tetToCompress, const RotationTable *thisRotationTable, int faceRotation, int reflect)
    // indices of vertexMap[] are compressed-output vertices; elements of vertexMap[] are the original vertices of tetToCompress
    {
        tetToCompress.assignIndex(polytet.nextIndex);
        for (int _faceNum=0; _faceNum<3; _faceNum++)
        {
#if 0
            int rotatedFaceNum = ((reflect ? _faceNum ^ (_faceNum <= 1) : _faceNum) + faceRotation) % 3;
#else
            int rotatedFaceNum = faceRotateReflect[reflect][_faceNum][faceRotation];
#endif
            int faceNum = thisRotationTable->faceMap[rotatedFaceNum];
            Tet *attachedTet = tetToCompress.faceAttached[faceNum];
            if (!attachedTet)
                continue;
            attachedTet->assignIndex(polytet.nextIndex);
            value |= (CompressedPolytetBits)1 << ((tetToCompress.index - 1 - 1) * 3 + _faceNum);

            int attachedFace = tetToCompress.faceAttachedFace[faceNum];
            int rotation = thisRotationTable->rotation[rotatedFaceNum];
            append(polytet, *attachedTet, &rotationTable[attachedFace], rotation, reflect);
        }
    }
    void uncompress(Polytet &polytet)
    {
        TetIndex index = 2;
        uncompressHelper(polytet, 1, index);
        if (index != polytet.size() - 1)
        {
            std::cerr << "Error! Got " << (unsigned)index << ", expected " << polytet.size() - 1 << std::endl;
            exit(-1);
        }
    }
    size_t hash() const
    {
        std::size_t seed  = std::hash<uint64_t>{}(((uint64_t*)&value)[0]);
#if MAXIMUM_TETCOUNT > 23
        {         } seed ^= std::hash<uint64_t>{}(((uint64_t*)&value)[1]) + (seed << 6) + (seed >> 2);
#endif
        return seed;
    }
};

void Tet::tagSkipOverlapCheckHelper(int depth, const RotationTable *thisRotationTable, int faceRotation/* = 0*/, CompressedSubpolytet curCompressedPath/* = 0*/, CompressedSubpolytet trit/* = 1*/)
{
    depth--;
    compressedPath = depth < 0 ? curCompressedPath : UINT64_MAX;
    for (int _faceNum=0; _faceNum<3; _faceNum++)
    {
        int rotatedFaceNum = faceRotateReflect[0][_faceNum][faceRotation];
        int faceNum = thisRotationTable->faceMap[rotatedFaceNum];
        Tet *attachedTet = faceAttached[faceNum];
        if (!attachedTet)
            continue;
        int attachedFace = faceAttachedFace[faceNum];
        int rotation = thisRotationTable->rotation[rotatedFaceNum];
        attachedTet->tagSkipOverlapCheckHelper(depth, &rotationTable[attachedFace], rotation, curCompressedPath + trit * (_faceNum + 1), trit * 3);
    }
}
void Tet::tagSkipOverlapCheck(int depth)
{
    depth--;
    compressedPath = UINT64_MAX;
    Tet *t = faceAttached[3];
    int attachedFace = faceAttachedFace[3];
    t->tagSkipOverlapCheckHelper(depth, &rotationTable[attachedFace]);
}

#ifdef USE_GMP
void printPolytet(Polytet &polytet)
{
    bool first = true;
    for (auto thisTet=polytet.cbegin(); thisTet!=polytet.cend(); ++thisTet)
    {
        printf(first ? "{" : ",\n");
        first = false;
        gmp_printf("{{%Zd, %Zd, %Zd}, {%Zd, %Zd, %Zd}, {%Zd, %Zd, %Zd}, {%Zd, %Zd, %Zd}}",
            thisTet->t.t[0][0], thisTet->t.t[0][1], thisTet->t.t[0][2],
            thisTet->t.t[1][0], thisTet->t.t[1][1], thisTet->t.t[1][2],
            thisTet->t.t[2][0], thisTet->t.t[2][1], thisTet->t.t[2][2],
            thisTet->t.t[3][0], thisTet->t.t[3][1], thisTet->t.t[3][2]);
    }
    printf("}\n\n");
}
#endif

#if defined(WRITE_TO_FILES) || defined(RESUME_FROM_FILE)
const char *getCompressedPolytetFilename(int tetCount)
{
    static char filename[100];
    sprintf(filename, "polytets_compressed_term_%d.bin", tetCount);
    return filename;
}

void writeFile(const char *filename, uint8_t *ptr, size_t size)
{
    FILE *f = fopen(filename, "wb");
    int fd = fileno(f);
    while (size > FILE_CHUNK_SIZE)
    {
        write(fd, ptr, FILE_CHUNK_SIZE);
        ptr  += FILE_CHUNK_SIZE;
        size -= FILE_CHUNK_SIZE;
    }
    write(fd, ptr, (uint32_t)size);
    fclose(f);
}

bool readAndCloseOpenedFile(FILE *f, uint8_t *ptr, size_t size)
{
    int fd = fileno(f);
    lseek64(fd, 0, SEEK_SET);
    while (size > FILE_CHUNK_SIZE)
    {
        if (read(fd, ptr, FILE_CHUNK_SIZE) != FILE_CHUNK_SIZE)
        {
            fclose(f);
            return false;
        }
        ptr  += FILE_CHUNK_SIZE;
        size -= FILE_CHUNK_SIZE;
    }
    if (read(fd, ptr, size) != size)
    {
        fclose(f);
        return false;
    }
    fclose(f);
    return true;
}
#endif

#ifdef USE_GMP
void mul_start_3(Tetrahedron &start, mpz_t maximalTouchingSqrDistance, CompressedSubpolytet &minUnseenCompressedSubpolytet)
#else
void mul_start_3(Tetrahedron &start, Coord &maximalTouchingSqrDistance, CompressedSubpolytet &minUnseenCompressedSubpolytet)
#endif
{
    for (int p=0; p<4; p++)
        for (int d=0; d<3; d++)
#ifdef USE_GMP
            mpz_mul_ui(start.t[p][d], start.t[p][d], 3);
#else
            start[p][d] *= 3;
#endif
#ifdef USE_GMP
    mpz_mul_ui(maximalTouchingSqrDistance, maximalTouchingSqrDistance, 3*3);
#else
    maximalTouchingSqrDistance *= 3*3;
#endif
    minUnseenCompressedSubpolytet = minUnseenCompressedSubpolytet * 3 + 1;
}

int main(int argc, char *argv[])
{
    if (uint64_t x=1; !*(uint8_t*)&x)
    {
        std::cerr << "Error: This program is hard-coded for little-endian byte order" << std::endl;
        exit(-1);
    }
    // The maximum distance between two touching congruent regular tetrahedrons is twice the radius (distance between the center of a tetrahedron and one of its vertices)
    const unsigned MAXIMAL_TOUCHING_SQR_DISTANCE =
        9*9    // squared coordinate of "start" tetrahedron
        * 3    // number of dimensions
        * 4*4  // squared number of vertices (to avoid dividing by 4 when averaging to get center coordinates)
        * 2*2; // twice the radius, so we square that too
    initLookupTables();
#ifdef USE_GMP
    Tetrahedron start;
    for (int d=0; d<3; d++)
    {
        mpz_set_si(start.t[0][d], -9);
        for (int p=1; p<4; p++)
            mpz_set_si(start.t[p][d], 9);
    }
    for (int p=1; p<4; p++)
        mpz_set_si(start.t[p][p-1], -9);
    mpz_t      maximalTouchingSqrDistance;
    mpz_init  (maximalTouchingSqrDistance);
    mpz_set_ui(maximalTouchingSqrDistance, MAXIMAL_TOUCHING_SQR_DISTANCE);
#else
    static Tetrahedron start =
    {{
        {{-9,-9,-9}},
        {{-9, 9, 9}},
        {{ 9,-9, 9}},
        {{ 9, 9,-9}}
    }};
    static Coord maximalTouchingSqrDistance = MAXIMAL_TOUCHING_SQR_DISTANCE;
#endif

    CompressedSubpolytet minUnseenCompressedSubpolytet = 0;
    bool *overlapBitmap; // An actual bitmap was tried, and was a bit slower; so, we'll use 8 times as much RAM to get slightly better speed
    {
        size_t overlapBitmapSize = 1;
        for (int i=0; i<MAXIMUM_TETCOUNT-2; i++)
            overlapBitmapSize *= 3;
        overlapBitmap = (bool*)malloc(overlapBitmapSize);
        memset(overlapBitmap, 0, overlapBitmapSize);
        std::cout << "Allocated " << overlapBitmapSize << " bytes for overlap caching" << std::endl;
    }

    size_t poolSize;
    void *pool = NULL;

    TetrahedronOverlap overlap;
    int minOverlapDepth = 2;
    bool foundOverlaps = false;

    size_t prevPolytetCount = 0;
    size_t polytetCount = 1;
    size_t memoryUsage = 0;
    
    int tetCount=1;
    bool resumedFromFile = false;
    size_t polytetChiralCount;
#ifdef RESUME_FROM_FILE
    {
        FILE *resumeFile = NULL;
        const char *filename;
        for (int i=3;; i++)
        {
            if (FILE *f = fopen(filename = getCompressedPolytetFilename(i), "rb"))
            {
                if (resumeFile) fclose(resumeFile);
                resumeFile = f;
                tetCount = i;
                mul_start_3(start, maximalTouchingSqrDistance, minUnseenCompressedSubpolytet);
            }
            else
                break;
        }
        if (resumeFile)
        {
            fseeko64(resumeFile, 0, SEEK_END);
            size_t size = ftello64(resumeFile);
            poolSize = size;
            if (poolSize < MEMORY_POOL_INITIAL_SIZE)
                poolSize = MEMORY_POOL_INITIAL_SIZE;
            pool = malloc(poolSize);
            if (!pool)
            {
                fclose(resumeFile);
                quitMemory();
            }
            int polytetsCompressedSize = ((tetCount - 2) * 3 + 8-1) / 8;
            polytetCount = size / polytetsCompressedSize;
            if (!readAndCloseOpenedFile(resumeFile, (uint8_t*)pool, size))
            {
                std::cerr << "Error reading file \"" << filename << "\"" << std::endl;
                goto errorQuit;
            }
            resumedFromFile = true;
        }
    }
    if (!pool)
#endif
    {
        pool = malloc(poolSize = MEMORY_POOL_INITIAL_SIZE);
        if (!pool) quitMemory();
    }

    polytetChiralCount = 0;
    for (;;)
    {
        auto currentTime = std::chrono::steady_clock::now();
        std::cout << tetCount << ": ";
        if (resumedFromFile)
            std::cout << "resumed";
        else
            std::cout << polytetCount + polytetChiralCount;
        std::cout << " (" << polytetCount << ") [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms";
        if (memoryUsage)
            std::cout << ", " << memoryUsage << " bytes";
        std::cout << "]" << std::endl;
        if (prevPolytetCount > polytetCount)
        {
            std::cerr << "Quit due to apparent overflow" << std::endl;
            break;
        }
        prevPolytetCount = polytetCount;
        if (++tetCount <= 2)
            continue;
        if (tetCount > MAXIMUM_TETCOUNT)
            break;

        int basePolytetCompressedSize = ((tetCount - 3) * 3 + 8-1) / 8;
        size_t hashTableSize = polytetCount * HASH_TABLE_RATIO;
        uint8_t *basePolytetTable;
        HashIndex *hashTable;
        void *polytetTable;
        for (;;)
        {
            basePolytetTable = (uint8_t*)pool;
            size_t basePolytetTableSize = basePolytetCompressedSize * polytetCount;
            hashTable = (HashIndex*)(basePolytetTable + basePolytetTableSize);
            polytetTable = hashTable + hashTableSize;
            size_t minSize = (uint8_t*)polytetTable - (uint8_t*)pool;
            if (minSize < poolSize)
                break;
            void *newPool = realloc(pool, basePolytetTableSize); // so that if the below realloc() results in a move, only what memory actually needs to be moved will be moved
            if (!newPool) {free(pool); quitMemory();}
            void *newPool2 = realloc(newPool, poolSize = minSize + minSize * MEMORY_POOL_GROW_RATIO);
            if (!newPool2) {free(newPool); quitMemory();}
            pool = newPool2;
        }
        memset(hashTable, 0, hashTableSize * sizeof(HashIndex));
        const int newPolytetsCompressedSize = ((tetCount - 2) * 3 + 8-1) / 8;
        const int polytetTableElementSize = newPolytetsCompressedSize + sizeof(HashIndex);
        size_t newPolytetCount = 0;
#ifdef SHOW_PROGRESS
        size_t nextProgressOutput = tetCount < SHOW_PROGRESS ? UINT64_MAX : 0;
        size_t progressOutputInterval = polytetCount / 1000;
#endif

        Polytet polytet;
        polytet.reserve(tetCount); // Important, to ensure pointers don't change
        Tet &t0    = polytet.emplace_back(start);
        attachNewTet(polytet.emplace_back(), t0, 3);

        polytet.resize(tetCount);
        polytetChiralCount = 0;
        for (size_t basePolytetI=0; basePolytetI<polytetCount; basePolytetI++)
        {
#ifdef SHOW_PROGRESS
            if (basePolytetI >= nextProgressOutput)
            {
                nextProgressOutput = basePolytetI + progressOutputInterval;
                unsigned perthouProgress = (basePolytetI * 1000) / polytetCount;
                std::cout << perthouProgress / 10 << "." << perthouProgress % 10 << "%\r";
                std::cout.flush();
            }
#endif
            CompressedPolytet *basePolytet = (CompressedPolytet*)(basePolytetTable + basePolytetI * basePolytetCompressedSize);
            {
                CompressedPolytet tmp;
                memcpy(&tmp.value, &basePolytet->value, basePolytetCompressedSize);
                tmp.uncompress(polytet);
            }

            Tet &newTet = polytet[tetCount - 1];
            for (int tetNumToAttachTo = 0; tetNumToAttachTo < tetCount-1; tetNumToAttachTo++)
            {
                Tet &tetToAttachTo = polytet[tetNumToAttachTo];
                for (int faceNum=0; faceNum<3; faceNum++) // skip last face because it's always already attached
                {
                    if (tetToAttachTo.faceAttached[faceNum])
                        continue;
                    attachNewTet(newTet, tetToAttachTo, faceNum);
                    // Canonicalize the rotation of this new polytet in compressed form, so that it can be compared against others
                    bool haveRunningLeast[2] = {false, false};
                    CompressedPolytet runningLeastPolytet[2];

                    Tet *t = &polytet[1];
                    for (int i=0; i<tetCount; i++)
                    {
                        const RotationTable *thisRotationTable;
                        {
                            Tet &singlyAttachedTet = polytet[i];
                            for (int j=0; j<3; j++)
                                if (singlyAttachedTet.faceAttached[j])
                                    goto skipThisTet; // not a singly attached tet
                            t = singlyAttachedTet.faceAttached[3];
                            int attachedFace = singlyAttachedTet.faceAttachedFace[3];
                            thisRotationTable = &rotationTable[attachedFace];
                        }
                        for (int reflect=0; reflect<2; reflect++)
                        for (int rotationStep=0; rotationStep<3; rotationStep++)
                        {
                            polytet.resetIndexing(i);
                            CompressedPolytet newRotatedPolytet;
                            newRotatedPolytet.append(polytet, *t, thisRotationTable, rotationStep, reflect);

                            // Update the running "least" rotation
                            if (!haveRunningLeast[reflect] || newRotatedPolytet.value < runningLeastPolytet[reflect].value)
                            {
                                haveRunningLeast[reflect] = true;
                                runningLeastPolytet[reflect] = newRotatedPolytet;
                            }
                        }
                    skipThisTet:;
                    }

                    bool isChiral = runningLeastPolytet[1].value != runningLeastPolytet[0].value;
                    if (isChiral && runningLeastPolytet[1].value <  runningLeastPolytet[0].value)
                        runningLeastPolytet[0] = runningLeastPolytet[1];

                    HashIndex *index = &hashTable[runningLeastPolytet[0].hash() % hashTableSize];
                    void *entry;
                    for (;;)
                    {
                        if (*index == 0)
                        {
                            entry = (uint8_t*)polytetTable + newPolytetCount * polytetTableElementSize;
                            break; // no duplicate of runningLeastPolytet[0] was found in hash table
                        }
                        entry = (uint8_t*)polytetTable + (size_t)(*index - 1) * polytetTableElementSize;
                        if (memcmp(entry, &runningLeastPolytet[0].value, newPolytetsCompressedSize) == 0)
                            goto skipDuplicate;
                        index = (HashIndex*)((uint8_t*)entry + newPolytetsCompressedSize);
                    }
                    // Check for overlap between this newly attached tetrahedron and the existing ones,
                    // and defer this until after the deduplication, to save a lot of time
                    overlap.setA(newTet);
                    // Set up the "skipOverlapCheck" flags to skip overlap checking up to a depth of 5
                    newTet.tagSkipOverlapCheck(minOverlapDepth);
                    for (auto tetCheckIntersection=polytet.cbegin(); tetCheckIntersection!=polytet.cend(); ++tetCheckIntersection)
                    {
                        if ((*tetCheckIntersection).compressedPath == UINT64_MAX)
                            continue; // skip this check for speed (it'll always be false anyway)
                        CompressedSubpolytet compressedPath = (tetCheckIntersection->compressedPath - 1) / 3;
                        if (compressedPath >= minUnseenCompressedSubpolytet)
                        {
                            overlap.setB(*tetCheckIntersection);
                            if (overlap(maximalTouchingSqrDistance))
                            {
                                overlapBitmap[compressedPath - 1] = true;

                                CompressedSubpolytet compressedPathReflected = 0, trit = 1;
                                while (compressedPath)
                                {
                                    compressedPath--;
                                    unsigned next = compressedPath % 3;
                                    if (next < 2)
                                        next ^= 1;
                                    compressedPathReflected += (next + 1) * trit;
                                    compressedPath /= 3;
                                    trit *= 3;
                                }
                                overlapBitmap[compressedPathReflected - 1] = true;

                                foundOverlaps = true;
                                goto skipDueToOverlap;
                            }
                        }
                        else
                        {
                            if (overlapBitmap[compressedPath - 1])
                            {
                                foundOverlaps = true;
                                goto skipDueToOverlap;
                            }
                        }
                    }
                    // No overlap found, so add runningLeastPolytet[0] to hash table and chiral count
#if defined(USE_GMP) && defined(PRINT_POLYTETS)
                    printPolytet(polytet);
#endif
                    polytetChiralCount += isChiral;
                    if ((uint8_t*)entry + polytetTableElementSize - (uint8_t*)pool > poolSize)
                    {
                        void *newPool = realloc(pool, poolSize += poolSize * MEMORY_POOL_GROW_RATIO);
                        if (!newPool) {free(pool); quitMemory();}
                        ptrdiff_t diff = (uint8_t*)newPool - (uint8_t*)pool;
                        pool = newPool;
                        (uint8_t*&)basePolytetTable += diff;
                        (uint8_t*&)hashTable        += diff;
                        (uint8_t*&)polytetTable     += diff;
                        (uint8_t*&)index            += diff;
                        (uint8_t*&)entry            += diff;
                    }
                    memcpy(entry, &runningLeastPolytet[0].value, newPolytetsCompressedSize);
                    *index = ++newPolytetCount;
                    *(HashIndex*)((uint8_t*)entry + newPolytetsCompressedSize) = 0; // pointer to next hash collision
                skipDuplicate:
                skipDueToOverlap:

                    tetToAttachTo.faceAttached[faceNum] = NULL;
                }
            }

            polytet[1].faceAttached[0] = NULL;
            polytet[1].faceAttached[1] = NULL;
            polytet[1].faceAttached[2] = NULL;
        }

        memoryUsage = (uint8_t*)polytetTable + newPolytetCount * polytetTableElementSize - (uint8_t*)pool;

        polytetCount = newPolytetCount;
        for (size_t i=0; i<polytetCount; i++)
            memcpy(
                basePolytetTable       + i *                          newPolytetsCompressedSize,
                (uint8_t*)polytetTable + i * polytetTableElementSize, newPolytetsCompressedSize);

#ifdef WRITE_TO_FILES
        writeFile(getCompressedPolytetFilename(tetCount), basePolytetTable, polytetCount * newPolytetsCompressedSize);
#endif

        resumedFromFile = false;
        mul_start_3(start, maximalTouchingSqrDistance, minUnseenCompressedSubpolytet);
        if (!foundOverlaps)
            minOverlapDepth++;
    }

errorQuit:
    free(pool);
    free(overlapBitmap);
#ifdef USE_GMP
    mpz_clear(maximalTouchingSqrDistance);
#endif
    return 0;
}
