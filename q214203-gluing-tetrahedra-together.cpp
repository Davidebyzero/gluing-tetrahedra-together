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
//#define KEEP_GOING

#define WRITE_TO_FILES
#define RESUME_FROM_FILE

#define FILE_CHUNK_SIZE (1 << 31)  // needs to be less than 1<<32

#ifdef USE_GMP
#   include <gmp.h>
#   if GMP_NUMB_BITS != 64
#       error This is hard-coded for 64-bit limbs
#   endif
#endif

auto startTime = std::chrono::steady_clock::now();

void quitMemory()
{
    auto currentTime = std::chrono::steady_clock::now();
    std::cerr << "Out of memory [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms]" << std::endl;
    exit(-1);
}

typedef uint8_t TetIndex;

#ifdef KEEP_GOING
typedef uint64_t HashIndex;
#else // just large enough for tetCount==17
typedef uint32_t HashIndex;
#endif

typedef unsigned __int128 CompressedPolytetBits; // Can handle up to 44 terms, while a 64-bit size_t can only handle up to about 28 terms

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
public:
    Tetrahedron t;
    Tet    *faceAttached    [4];
    uint8_t faceAttachedFace[4];
    TetIndex index; // 1-based; 0=unassigned
    bool skipOverlapCheck;
    Tet(                    ) : t( ) {initFaces();}
    Tet(const Tetrahedron &t) : t(t) {initFaces();}
    void assignIndex(TetIndex &nextIndex)
    {
        if (index == 0)
            index = nextIndex++;
    }
    void tagSkipOverlapCheck(int depth)
    {
        if (skipOverlapCheck)
            return;
        skipOverlapCheck = true;
        if (--depth <= 0)
            return;
        for (int faceNum=0; faceNum<4; faceNum++)
        {
            if (auto attached = faceAttached[faceNum])
                attached->tagSkipOverlapCheck(depth);
        }
    }
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

static const uint8_t vertexMapTable[4][4] =
{
    {3, 0, 2, 1},
    {2, 0, 1, 3},
    {1, 0, 3, 2},
    {0, 1, 2, 3},
};

static const int tetrahedronEdges[6][2] =
{
    {0, 1},
    {1, 2},
    {2, 0},
    {0, 3},
    {1, 3},
    {2, 3},
};

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

class TetrahedronOverlap
{
    const Tetrahedron *a;
    const Tet         *b;
public:
    void setA(const Tetrahedron &x) {a = &x;}
    void setB(const Tet         &x) {b = &x;}
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
public:
    TetrahedronOverlap()
    {
        mpz_inits(intersectNumerator, intersectDenominator, uNumerator, vNumerator, uvDenominator, uvNumeratorSum, NULL);
        for (int d=0; d<3; d++)
        {
            mpz_init(center[d]);
            mpz_init(normal[d]);
            mpz_init(tmp[d]);
            mpz_init(p0p1[d]);
            mpz_init(intersectionPoint[d]);
            mpz_init(delta[d]);
            mpz_init(edge1[d]);
            mpz_init(edge2[d]);
        }
        for (int p=0; p<3; p++)
            for (int d=0; d<3; d++)
                mpz_init(triangle[p][d]);
    }
    ~TetrahedronOverlap()
    {
        mpz_clears(intersectNumerator, intersectDenominator, uNumerator, vNumerator, uvDenominator, uvNumeratorSum, NULL);
        for (int d=0; d<3; d++)
        {
            mpz_clear(center[d]);
            mpz_clear(normal[d]);
            mpz_clear(tmp[d]);
            mpz_clear(p0p1[d]);
            mpz_clear(intersectionPoint[d]);
            mpz_clear(delta[d]);
            mpz_clear(edge1[d]);
            mpz_clear(edge2[d]);
        }
        for (int p=0; p<3; p++)
            for (int d=0; d<3; d++)
                mpz_clear(triangle[p][d]);
    }
    bool operator()(const mpz_t maximalTouchingSqrDistance)
    {
        // Skip the longer overlap checking algorithm if the two tetrahedrons' centers are sufficiently separated.
        // Get center of tetrahedron "a" by averaging its vertices' coordinates, without dividing by 4.
        for (int d=0; d<3; d++)
            mpz_set(center[d], a->t[0][d]);
        for (int p=1; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_add(center[d], center[d], a->t[p][d]);
        // Subtract center of tetrahedron "b" by averaging its vertices' coordinates, without dividing by 4.
        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_sub(center[d], center[d], b->t.t[p][d]);
        // Compare the sum of the squares of the orthogonal distances against the threshold squared distance.
        dot(tmp[0], center, center);
        if (mpz_cmp(tmp[0], maximalTouchingSqrDistance) >= 0)
            return false;

        // Take advantage of the fact that the two tetrahedrons are regular and congruent, and
        // just check if any edge of tetrahedron "a" intersects with any face of tetrahedron "b".
        // Don't count it if only the endpoint of an edge intersects.
        for (int edgeNum=0; edgeNum<6; edgeNum++)
        {
            const mpz_t *p0 = a->t[tetrahedronEdges[edgeNum][0]];
            const mpz_t *p1 = a->t[tetrahedronEdges[edgeNum][1]];
            for (int faceNum=0; faceNum<4; faceNum++)
            {
                // For speed, process only one out of every pair of attached faces (which share the exact same 3 vertices).
                // This will still process the very first attached face twice, since that is a face[3] attached to another face[3],
                // but it's probably not worth the extra machinery that would be necessary to special-case that.
                if (faceNum!=3 && b->faceAttached[faceNum])
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
                    center[d] += (*a)[p][d];
            // Subtract center of tetrahedron "b" by averaging its vertices' coordinates, without dividing by 4.
            for (int p=0; p<4; p++)
                for (int d=0; d<3; d++)
                    center[d] -= b->t[p][d];
            // Compare the sum of the squares of the orthogonal distances against the threshold squared distance.
            if (dot(center, center) >= maximalTouchingSqrDistance)
                return false;
        }
        // Take advantage of the fact that the two tetrahedrons are regular and congruent, and
        // just check if any edge of tetrahedron "a" intersects with any face of tetrahedron "b".
        // Don't count it if only the endpoint of an edge intersects.
        for (int edgeNum=0; edgeNum<6; edgeNum++)
        {
            Coord3 p0 = (*a)[tetrahedronEdges[edgeNum][0]];
            Coord3 p1 = (*a)[tetrahedronEdges[edgeNum][1]];
            for (int faceNum=0; faceNum<4; faceNum++)
            {
                // For speed, process only one out of every pair of attached faces (which share the exact same 3 vertices).
                // This will still process the very first attached face twice, since that is a face[3] attached to another face[3],
                // but it's probably not worth the extra machinery that would be necessary to special-case that.
                if (faceNum!=3 && b->faceAttached[faceNum])
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
    void append(Polytet &polytet, Tet &tetToCompress, const uint8_t *vertexMap, int faceRotation, int reflect)
    // indices of vertexMap[] are compressed-output vertices; elements of vertexMap[] are the original vertices of tetToCompress
    {
        tetToCompress.assignIndex(polytet.nextIndex);
        for (int _faceNum=0; _faceNum<3; _faceNum++)
        {
#if 0
            int rotatedFaceNum = ((reflect ? _faceNum ^ (_faceNum <= 1) : _faceNum) + faceRotation) % 3;
#else
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
            int rotatedFaceNum = faceRotateReflect[reflect][_faceNum][faceRotation];
#endif
            int faceNum = 3 - vertexMap[3 - rotatedFaceNum];
            Tet *attachedTet = tetToCompress.faceAttached[faceNum];
            if (!attachedTet)
                continue;
            attachedTet->assignIndex(polytet.nextIndex);
            value |= (CompressedPolytetBits)1 << ((tetToCompress.index - 1 - 1) * 3 + _faceNum);

            int attachedFace = tetToCompress.faceAttachedFace[faceNum];
            int rotation = 0;
            while (vertexMap[tetrahedronFaces[rotatedFaceNum][rotation]] != tetrahedronFaces[faceNum][0])
                rotation++;
            append(polytet, *attachedTet, vertexMapTable[attachedFace], rotation, reflect);
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
        {         } seed ^= std::hash<uint64_t>{}(((uint64_t*)&value)[1]) + (seed << 6) + (seed >> 2);
        return seed;
    }
};

#if defined(WRITE_TO_FILES) || defined(RESUME_FROM_FILE)
const char *getCompressedPolytetFilename(int tetCount)
{
    static char filename[100];
    sprintf(filename, "polytets_compressed_term_%d.bin", tetCount);
    return filename;
}
#endif

#ifdef USE_GMP
void mul_start_3(Tetrahedron &start, mpz_t maximalTouchingSqrDistance)
#else
void mul_start_3(Tetrahedron &start, Coord &maximalTouchingSqrDistance)
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

    size_t poolSize;
    void *pool = NULL;

    TetrahedronOverlap overlap;
    int minOverlapDepth = 2;
    bool foundOverlaps = false;

    size_t prevPolytetCount = 0;
    size_t polytetCount = 1;
    size_t memoryUsage = 0;
    
    int tetCount=1;
#ifdef RESUME_FROM_FILE
    {
        FILE *resumeFile = NULL;
        for (int i=3;; i++)
        {
            if (FILE *f = fopen(getCompressedPolytetFilename(i), "rb"))
            {
                if (resumeFile) fclose(resumeFile);
                resumeFile = f;
                tetCount = i;
                mul_start_3(start, maximalTouchingSqrDistance);
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
            int fd = fileno(resumeFile);
            lseek64(fd, 0, SEEK_SET);
            uint8_t *ptr = (uint8_t*)pool;
            while (size > FILE_CHUNK_SIZE)
            {
                read(fd, ptr, FILE_CHUNK_SIZE);
                ptr  += FILE_CHUNK_SIZE;
                size -= FILE_CHUNK_SIZE;
            }
            read(fd, ptr, size);
            fclose(resumeFile);
        }
    }
    if (!pool)
#endif
    {
        pool = malloc(poolSize = MEMORY_POOL_INITIAL_SIZE);
        if (!pool) quitMemory();
    }

    bool resumedFromFile = tetCount > 1;
    size_t polytetChiralCount = 0;
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
#ifndef KEEP_GOING
        if (tetCount > 17)
            break;
#endif

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
                        const uint8_t *vertexMap;
                        {
                            Tet &singlyAttachedTet = polytet[i];
                            for (int j=0; j<3; j++)
                                if (singlyAttachedTet.faceAttached[j])
                                    goto skipThisTet; // not a singly attached tet
                            t = singlyAttachedTet.faceAttached[3];
                            int attachedFace = singlyAttachedTet.faceAttachedFace[3];
                            vertexMap = vertexMapTable[attachedFace];
                        }
                        for (int reflect=0; reflect<2; reflect++)
                        for (int rotationStep=0; rotationStep<3; rotationStep++)
                        {
                            polytet.resetIndexing(i);
                            CompressedPolytet newRotatedPolytet;
                            newRotatedPolytet.append(polytet, *t, vertexMap, rotationStep, reflect);

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
                    overlap.setA(newTet.t);
                    // Set up the "skipOverlapCheck" flags to skip overlap checking up to a depth of 5
                    for (auto tetCheckIntersection=polytet.begin(); tetCheckIntersection!=polytet.end(); ++tetCheckIntersection)
                        (*tetCheckIntersection).skipOverlapCheck = false;
                    newTet.tagSkipOverlapCheck(minOverlapDepth);
                    for (auto tetCheckIntersection=polytet.cbegin(); tetCheckIntersection!=polytet.cend(); ++tetCheckIntersection)
                    {
                        if ((*tetCheckIntersection).skipOverlapCheck)
                            continue; // skip this check for speed (it'll always be false anyway)
                        overlap.setB(*tetCheckIntersection);
                        if (overlap(maximalTouchingSqrDistance))
                        {
                            foundOverlaps = true;
                            goto skipDueToOverlap;
                        }
                    }
                    // No overlap found, so add runningLeastPolytet[0] to hash table and chiral count
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
        FILE *f = fopen(getCompressedPolytetFilename(tetCount), "wb");
        int fd = fileno(f);
        size_t size = polytetCount * newPolytetsCompressedSize;
        uint8_t *ptr = basePolytetTable;
        while (size > FILE_CHUNK_SIZE)
        {
            write(fd, ptr, FILE_CHUNK_SIZE);
            ptr  += FILE_CHUNK_SIZE;
            size -= FILE_CHUNK_SIZE;
        }
        write(fd, ptr, (uint32_t)size);
        fclose(f);
#endif

        resumedFromFile = false;
        mul_start_3(start, maximalTouchingSqrDistance);
        if (!foundOverlaps)
            minOverlapDepth++;
    }

    free(pool);
#ifdef USE_GMP
    mpz_clear(maximalTouchingSqrDistance);
#endif
    return 0;
}
