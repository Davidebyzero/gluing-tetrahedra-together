// See https://codegolf.stackexchange.com/a/283991/17216 and https://oeis.org/A276272

#include <stdio.h>
#include <fcntl.h>
#include <iostream>
#include <string.h>
#include <array>
#include <inttypes.h>
#include <chrono>
#include "config.h"

#ifdef USE_GMP
#   include <gmp.h>
#   if GMP_NUMB_BITS != 64
#       error This is hard-coded for 64-bit limbs
#   endif
#endif

#ifdef MULTITHREADING
    #include <thread>
    #include <boost/thread/mutex.hpp>
    typedef unsigned THREAD_ID;
#endif // MULTITHREADING

#ifndef _countof
    #define _countof(a) (sizeof(a)/sizeof(*(a)))
#endif

auto startTime = std::chrono::steady_clock::now();

void quitMemory()
{
    auto currentTime = std::chrono::steady_clock::now();
    std::cerr << "Out of memory [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms]" << std::endl;
    exit(-1);
}

typedef int TetIndex;

#define MIN_OVERLAP_DEPTH 5

#if MAXIMUM_TETCOUNT > 17
typedef uint64_t HashIndex;
#else
typedef uint32_t HashIndex;
#endif

#if MAXIMUM_TETCOUNT > 23
typedef unsigned __int128 CompressedPolytetBits; // Can handle up to 44 terms, while a 64-bit size_t can only handle up to about 28 terms
typedef   signed __int128 CompressedPolytetBitsSigned;
#define COMPRESSEDPOLYTETBITS_MAX (unsigned __int128)-1
#else
typedef uint64_t          CompressedPolytetBits; // Can handle up to 23 terms
typedef  int64_t          CompressedPolytetBitsSigned;
#define COMPRESSEDPOLYTETBITS_MAX UINT64_MAX
#endif

typedef uint64_t CompressedSubpolytet; // compressed in branchless format, in bijective trinary; for caching known-overlapping subpolytet paths

#ifdef USE_GMP
class Tetrahedron
{
public:
    mpz_t t[4][3];
    Tetrahedron &operator=(const Tetrahedron &_t)
    {
        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_set(t[p][d], _t.t[p][d]);
        return *this;
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

#if defined(PRINT_POLYTETS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
#   define CALCULATE_COORDINATES_FOR_PRINTING true
#else
#   define CALCULATE_COORDINATES_FOR_PRINTING false
#endif

class Polytet;
class Tet
{
    void tagSkipOverlapCheckHelper(Polytet &polytet, int depth, const struct RotationTable *thisRotationTable, int faceRotation = 0, CompressedSubpolytet curCompressedPath = 0, CompressedSubpolytet trit = 1);
public:
    Tetrahedron t;
    TetIndex faceAttached    [4];
    uint8_t  faceAttachedFace[4];
    TetIndex index; // 1-based; 0=unassigned
    CompressedSubpolytet compressedPath; // 0 = skip overlap check
    bool isLeaf;
    Tet() : t() {}
    void init()
    {
        faceAttached[0] = -1; // t[0],t[1],t[2]
        faceAttached[1] = -1; // t[0],t[1],t[3]
        faceAttached[2] = -1; // t[0],t[2],t[3]
        faceAttached[3] = -1; // t[1],t[2],t[3]
    }
    void assignIndex(TetIndex &nextIndex)
    {
        if (index == 0)
            index = nextIndex++;
    }
    void tagSkipOverlapCheck(Polytet &polytet, int depth);
};
class Polytet : public std::array<Tet, MAXIMUM_TETCOUNT>
{
    int tetCount;
public:
    TetIndex nextIndex;

    template <bool CALCULATE_COORDINATES = CALCULATE_COORDINATES_FOR_PRINTING>
    void init(const Tetrahedron &start)
    {
        (*this)[0].init();
        if (CALCULATE_COORDINATES)
            (*this)[0].t = start;
        (*this)[0].isLeaf = true;
        attachNewTet<CALCULATE_COORDINATES>(1, 0, 3);
    }
    void resetIndexing(size_t first) // This function should not be called if the polytet hasn't yet been populated
    {
        for (int i=0; i<tetCount; i++)
            (*this)[i].index = 0;
        (*this)[first].index = 1;
        nextIndex = 2;
    }
    void setSize(int x) {tetCount = x;}
    int size() {return tetCount;}

    template <bool CALCULATE_COORDINATES = CALCULATE_COORDINATES_FOR_PRINTING>
    void attachNewTet(TetIndex i, TetIndex iAttachTo, const int faceNum)
    {
        Polytet &polytet = *this;
        Tet &t             = polytet[i];
        Tet &tetToAttachTo = polytet[iAttachTo];
        if (CALCULATE_COORDINATES)
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
        }
        t.faceAttached    [0] = -1;
        t.faceAttached    [1] = -1;
        t.faceAttached    [2] = -1;
        t.faceAttached    [3] = iAttachTo;
        t.faceAttachedFace[3] = faceNum;
        t.isLeaf = true;
        tetToAttachTo.faceAttached    [faceNum] = i;
        tetToAttachTo.faceAttachedFace[faceNum] = 3;
        tetToAttachTo.isLeaf &= faceNum == 3;
    }
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
Coord3 &operator+=(Coord3 &a, const Coord3 &b)
{
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
    return a;
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

#ifdef PRINT_POLYTETS_EDGE_CASES
#   define OVERLAP_EDGES_TO_EDGES -1
#else
#   define OVERLAP_EDGES_TO_EDGES 1
#endif

// This class must be used in exactly the way it is in this program. That is, setA() and setB() must be called first, before operator(),
// and point to tetrahedrons in the same Polytet. Both tetrahedrons must be end pieces, singly attached. If operator() returns "true",
// both setA() and setB() must be called before the next operator() call.
class TetrahedronOverlap
{
    const Tet *a;
    const Tet *b;
public:
    void setA(const Tet &x) {a = &x;}
    void setB(const Tet &x) {b = &x;}
#ifdef USE_GMP
private:
    mpz_t intersectNumerator, intersectDenominator, uNumerator, vNumerator, uvDenominator, uvNumeratorSum, pdot, qdot, absDiff, sdot;
    mpz_t aCenter4[3], displacement4[3], normal[3], tmp[3], p0p1[3], intersectionPoint[3], delta[3], edge1[3], edge2[3], p[3], q[3], intersect[3];
    mpz_t face[3][3];
    void dot(mpz_t &result, const mpz_t _a[3], const mpz_t _b[3])
    {
        mpz_mul   (result, _a[0], _b[0]);
        mpz_addmul(result, _a[1], _b[1]);
        mpz_addmul(result, _a[2], _b[2]);
    }
    template <void (*MPZ_CALL)(mpz_t x), void (*MPZ_CALLS)(mpz_t x, ...)> void mpz_initOrClear()
    {
        MPZ_CALLS(intersectNumerator, intersectDenominator, uNumerator, vNumerator, uvDenominator, uvNumeratorSum, pdot, qdot, absDiff, sdot, NULL);
        for (int d=0; d<3; d++)
            MPZ_CALLS(aCenter4[d], displacement4[d], normal[d], tmp[d], p0p1[d], intersectionPoint[d], delta[d], edge1[d], edge2[d], p[d], q[d], intersect[d], NULL);
        for (int p=0; p<3; p++)
            for (int d=0; d<3; d++)
                MPZ_CALL(face[p][d]);
    }
public:
     TetrahedronOverlap() {mpz_initOrClear<mpz_init , mpz_inits >();}
    ~TetrahedronOverlap() {mpz_initOrClear<mpz_clear, mpz_clears>();}
    int8_t operator()(const mpz_t maximalTouchingSqrDistance)
    {
        for (int d=0; d<3; d++)
        {
            // Get center of tetrahedron "a" by averaging its vertices' coordinates, without dividing by 4.
            mpz_add(aCenter4[d], a->t.t[0][d], a->t.t[1][d]);
            for (int p=2; p<4; p++)
                mpz_add(aCenter4[d], aCenter4[d], a->t.t[p][d]);
            // Skip the longer overlap checking algorithm if the two tetrahedrons' centers are sufficiently separated.
            // Subtract center of tetrahedron "b" by averaging its vertices' coordinates, without dividing by 4.
            mpz_sub(displacement4[d], aCenter4[d], b->t.t[0][d]);
            for (int p=1; p<4; p++)
                mpz_sub(displacement4[d], displacement4[d], b->t.t[p][d]);
        }
        // Compare the sum of the squares of the orthogonal distances against the threshold squared distance.
        dot(tmp[0], displacement4, displacement4);
        if (mpz_cmp(tmp[0], maximalTouchingSqrDistance) >= 0)
            return false;

        const unsigned multiplier = 4;
        mpz_div_ui(sdot, maximalTouchingSqrDistance, 3);
        int8_t edgeEdgeIntersectionCount[4] = {}; // number of edges found to intersect with each face
        // In some circumstances, two tetrahedrons can intersect such that only the edges of one tetrahedron intersect with the
        // faces of the other, and not the other way around. So we need to check both.
        for (int swapped=0; swapped<2; swapped++)
        {
            // Take advantage of the fact that the two tetrahedrons are regular and congruent, and
            // just check if any edge of tetrahedron "b" intersects with any face of tetrahedron "a".
            // Don't count it if only the endpoint of an edge intersects.
            // Since it's guaranteed that due to the overlapCache, both tetrahedrons are single-attached end pieces, and all
            // collisions with nearer pieces are already cached, we can safely skip face[3]
            for (int faceNum=0; faceNum<3; faceNum++)
            {
                static const int vertexRemap[4] = {1, 2, 3, 0};
                for (int i=0; i<3; i++)
                {
                    int vertexNum = vertexRemap[(faceNum + i) % 3];
                    for (int d=0; d<3; d++)
                    {
                        mpz_mul_ui(face[i][d], a->t.t[vertexNum][d], multiplier);
                        mpz_sub   (face[i][d], face[i][d], aCenter4[d]);
                    }
                }
                for (int edgeNum=0; edgeNum<3; edgeNum++)
                {
                    for (int d=0; d<3; d++)
                    {
                        mpz_sub(p[d], b->t.t[vertexRemap[edgeNum]][d], a->t.t[vertexRemap[3]][d]);
                        mpz_sub(q[d], b->t.t[vertexRemap[3      ]][d], a->t.t[vertexRemap[3]][d]);
                        mpz_mul_ui(p[d], p[d], multiplier);
                        mpz_mul_ui(q[d], q[d], multiplier);
                    }

                    dot(pdot, face[0], p);
                    dot(qdot, face[0], q);

                    // Check if the signs differ
                    if ((mpz_sgn(pdot) ^ mpz_sgn(qdot)) < 0)
                    {
                        for (int d=0; d<3; d++)
                        {
                            mpz_mul(p[d], p[d], qdot);
                            mpz_mul(q[d], q[d], pdot);
                        }
                        if (mpz_cmp(qdot, pdot) > 0)
                        {
                            for (int d=0; d<3; d++)
                                mpz_sub(intersect[d], p[d], q[d]);
                        }
                        else
                        {
                            for (int d=0; d<3; d++)
                                mpz_sub(intersect[d], q[d], p[d]);
                        }

                        dot(uNumerator, face[1], intersect);
                        dot(vNumerator, face[2], intersect);
                        if (mpz_cmp(pdot, qdot) > 0)
                            mpz_sub(absDiff, pdot, qdot);
                        else
                            mpz_sub(absDiff, qdot, pdot);
                        mpz_add(uvNumeratorSum, uNumerator, vNumerator);
                        mpz_mul(uvDenominator, sdot, absDiff);

                        int uCmp = mpz_cmp_ui(uNumerator, 0);
                        int vCmp = mpz_cmp_ui(vNumerator, 0);
                        int uvCmp = mpz_cmp(uvNumeratorSum, uvDenominator);
                        if (uCmp > 0 && vCmp > 0 && uvCmp < 0)
                            return true;
                        if (uCmp == 0 || vCmp == 0 || uvCmp == 0)
                            edgeEdgeIntersectionCount[faceNum]++;
                    }
                }
            }
            std::swap(a, b);
        }
        // The above will fail to detect an overlap in which 3 edges of one tetrahedron perfectly intersect with 3 edges of the other.
        // But we can detect that case by checking the edgeEdgeIntersectionCount[] values.
        return edgeEdgeIntersectionCount[0] == 4 &&
               edgeEdgeIntersectionCount[1] == 4 &&
               edgeEdgeIntersectionCount[2] == 4 ? OVERLAP_EDGES_TO_EDGES : 0;
    }
#else
public:
    int8_t operator()(const Coord &maximalTouchingSqrDistance)
    {
        // Get center of tetrahedron "a" by averaging its vertices' coordinates, without dividing by 4.
        Coord3 aCenter4 = a->t[0];
        for (int p=1; p<4; p++)
            aCenter4 += a->t[p];
        // Skip the longer overlap checking algorithm if the two tetrahedrons' centers are sufficiently separated.
        {
            Coord3 displacement4 = aCenter4;
            // Subtract center of tetrahedron "b" by averaging its vertices' coordinates, without dividing by 4.
            for (int p=0; p<4; p++)
                for (int d=0; d<3; d++)
                    displacement4[d] -= b->t[p][d];
            // Compare the sum of the squares of the orthogonal distances against the threshold squared distance.
            if (dot(displacement4, displacement4) >= maximalTouchingSqrDistance)
                return false;
        }

        const int multiplier = 4;
        const Coord sdot = maximalTouchingSqrDistance / 3;
        int8_t edgeEdgeIntersectionCount[4] = {}; // number of edges found to intersect with each face
        // In some circumstances, two tetrahedrons can intersect such that only the edges of one tetrahedron intersect with the
        // faces of the other, and not the other way around. So we need to check both.
        for (int swapped=0; swapped<2; swapped++)
        {
            // Take advantage of the fact that the two tetrahedrons are regular and congruent, and
            // just check if any edge of tetrahedron "b" intersects with any face of tetrahedron "a".
            // Don't count it if only the endpoint of an edge intersects.
            // Since it's guaranteed that due to the overlapCache, both tetrahedrons are single-attached end pieces, and all
            // collisions with nearer pieces are already cached, we can safely skip face[3]
            for (int faceNum=0; faceNum<3; faceNum++)
            {
                static const int vertexRemap[4] = {1, 2, 3, 0};
                Coord3 face[3] =
                {
                    a->t[vertexRemap[ faceNum      % 3]] * multiplier - aCenter4,
                    a->t[vertexRemap[(faceNum + 1) % 3]] * multiplier - aCenter4,
                    a->t[vertexRemap[(faceNum + 2) % 3]] * multiplier - aCenter4
                };
                for (int edgeNum=0; edgeNum<3; edgeNum++)
                {
                    Coord3 p = (b->t[vertexRemap[edgeNum]] - a->t[vertexRemap[3]]) * multiplier;
                    Coord3 q = (b->t[vertexRemap[3]]       - a->t[vertexRemap[3]]) * multiplier;

                    Coord pdot = dot(face[0], p);
                    Coord qdot = dot(face[0], q);

                    // Check if the signs differ
                    if ((pdot ^ qdot) < 0)
                    {
                        Coord3 intersect = qdot > pdot ?
                            p * qdot - q * pdot :
                            q * pdot - p * qdot;

                        Coord uNumerator = dot(face[1], intersect);
                        Coord vNumerator = dot(face[2], intersect);
                        Coord absDiff = pdot > qdot ? pdot - qdot : qdot - pdot;
                        Coord uvNumeratorSum = uNumerator + vNumerator, uvDenominator = sdot * absDiff;

                        if (uNumerator > 0 && vNumerator > 0 && uvNumeratorSum < uvDenominator)
                            return true;
                        if (uNumerator == 0 || vNumerator == 0 || uvNumeratorSum == uvDenominator)
                            edgeEdgeIntersectionCount[faceNum]++;
                    }
                }
            }
            std::swap(a, b);
        }
        // The above will fail to detect an overlap in which 3 edges of one tetrahedron perfectly intersect with 3 edges of the other.
        // But we can detect that case by checking the edgeEdgeIntersectionCount[] values.
        return edgeEdgeIntersectionCount[0] == 4 &&
               edgeEdgeIntersectionCount[1] == 4 &&
               edgeEdgeIntersectionCount[2] == 4 ? OVERLAP_EDGES_TO_EDGES : 0;
    }
#endif
};

// First two tetrahedrons are implied. Each element is a subsequent tetrahedron, with the value indicating where
// it's attached. The lower 2 bits indicate which face (can only have 3 different values, because at least 1 face
// will always already be attached). The remaining bits indicate which tetrahedron (which can never be zero,
// because that one is attached implicitly).
class CompressedPolytet
{
    void uncompressHelper(TetIndex index, TetIndex &nextIndex)
    {
        for (int faceNum=0; faceNum<3; faceNum++)
        {
            if (value & ((CompressedPolytetBits)1 << ((index - 1) * 3 + faceNum)))
            {
                TetIndex thisIndex = nextIndex++;
                polytet.attachNewTet(thisIndex, index, faceNum);
                uncompressHelper(thisIndex, nextIndex);
            }
        }
    }
public:
    CompressedPolytetBits value, runningLeastValue;
    Polytet &polytet;
    int reflect;
    CompressedPolytet(Polytet &_polytet) : polytet(_polytet), value(0) {}
    bool append(Tet &tetToCompress, const RotationTable *thisRotationTable, int faceRotation)
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
            if (tetToCompress.faceAttached[faceNum] < 0)
                continue;
            Tet *attachedTet = &polytet[tetToCompress.faceAttached[faceNum]];
            attachedTet->assignIndex(polytet.nextIndex);
            value |= (CompressedPolytetBits)1 << ((tetToCompress.index - 1 - 1) * 3 + _faceNum);
            if (value > runningLeastValue)
                return false;

            int attachedFace = tetToCompress.faceAttachedFace[faceNum];
            int rotation = thisRotationTable->rotation[rotatedFaceNum];
            if (!append(*attachedTet, &rotationTable[attachedFace], rotation))
                return false;
        }
        return true;
    }
    void uncompress()
    {
        TetIndex index = 2;
        uncompressHelper(1, index);
        if (index != polytet.size() - 1)
        {
            std::cerr << "Error! Got " << (unsigned)index << ", expected " << polytet.size() - 1 << std::endl;
            exit(-1);
        }
    }
};
size_t hash(const CompressedPolytetBits &value)
{
    std::size_t seed  = std::hash<uint64_t>{}(((uint64_t*)&value)[0]);
#if MAXIMUM_TETCOUNT > 23
    {         } seed ^= std::hash<uint64_t>{}(((uint64_t*)&value)[1]) + (seed << 6) + (seed >> 2);
#endif
    return seed;
}

void Tet::tagSkipOverlapCheckHelper(Polytet &polytet, int depth, const RotationTable *thisRotationTable, int faceRotation/* = 0*/, CompressedSubpolytet curCompressedPath/* = 0*/, CompressedSubpolytet trit/* = 1*/)
{
    depth--;
    compressedPath = depth < 0 ? curCompressedPath : UINT64_MAX;
    for (int _faceNum=0; _faceNum<3; _faceNum++)
    {
        int rotatedFaceNum = faceRotateReflect[0][_faceNum][faceRotation];
        int faceNum = thisRotationTable->faceMap[rotatedFaceNum];
        if (faceAttached[faceNum] < 0)
            continue;
        Tet *attachedTet = &polytet[faceAttached[faceNum]];
        int attachedFace = faceAttachedFace[faceNum];
        int rotation = thisRotationTable->rotation[rotatedFaceNum];
        attachedTet->tagSkipOverlapCheckHelper(polytet, depth, &rotationTable[attachedFace], rotation, curCompressedPath + trit * (_faceNum + 1), trit * 3);
    }
}
void Tet::tagSkipOverlapCheck(Polytet &polytet, int depth)
{
    depth--;
    compressedPath = UINT64_MAX;
    Tet *t = &polytet[faceAttached[3]];
    int attachedFace = faceAttachedFace[3];
    t->tagSkipOverlapCheckHelper(polytet, depth, &rotationTable[attachedFace]);
}

CompressedSubpolytet reflectCompressedPath(CompressedSubpolytet compressedPath)
{
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
    return compressedPathReflected;
}

#ifdef USE_GMP
#   if defined(DISABLE_OVERLAP_CHECKING) && (defined(PRINT_POLYTETS) || defined(PRINT_POLYTETS_EDGE_CASES) || defined(PRINT_POLYTETS_WITH_SYMMETRY))
#       error Cannot print polytets with overlap checking disabled
#   endif
mpz_t printCenter[3], printTmp[3];
void printPolytet(CompressedPolytetBits value, Polytet &polytet,
    const char *valuePrefix="", const char *valueSuffix="\n", bool nestedParens=true, const char *openParen="{", const char *closeParen="}")
{
    std::cout << valuePrefix;
#if MAXIMUM_TETCOUNT > 23
    printf("0x%" PRIX64 "%08" PRIX64, ((uint64_t*)&value)[1], ((uint64_t*)&value)[0]);
#else
    printf("0x%" PRIX64, value);
#endif
    std::cout << valueSuffix;
    bool first = true;
    for (int d=0; d<3; d++)
        mpz_set_ui(printCenter[d], 0);
    for (int i=0; i<polytet.size(); i++)
        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_add(printCenter[d], printCenter[d], polytet[i].t.t[p][d]);
    for (int d=0; d<3; d++)
        mpz_div_ui(printCenter[d], printCenter[d], polytet.size() * 4);
    for (int i=0; i<polytet.size(); i++)
    {
        fputs(first ? (nestedParens ? openParen : "") : ",\n", stdout);
        first = false;
        if (nestedParens) std::cout << openParen;
        for (int p=0; p<4; p++)
        {
            for (int d=0; d<3; d++)
                mpz_sub(printTmp[d], polytet[i].t.t[p][d], printCenter[d]);
            gmp_printf("%s%s%Zd, %Zd, %Zd%s", p ? ", " : "", openParen, printTmp[0], printTmp[1], printTmp[2], closeParen);
        }
        if (nestedParens) std::cout << closeParen;
    }
    if (nestedParens) std::cout << closeParen;
    std::cout << std::endl << std::endl;
}
#endif

#define OVERLAP_BITMAP_FILENAME "polytets_overlap.bin"

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

bool readFile(const char *filename, uint8_t *ptr, size_t size)
{
    FILE *f = fopen(filename, "rb");
    if (!f)
        return false;
    return readAndCloseOpenedFile(f, ptr, size);
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

#ifdef MULTITHREADING
    static boost::mutex hashTableMutex[HASH_TABLE_SHARDS]; // shard the hash table locking
    static boost::mutex workAssignmentMutex;
    static boost::mutex printPolytetMutex;
    static bool hashTableInUse = false;

    struct WorkerJob
    {
#endif // MULTITHREADING
        Polytet workerPolytet;
#ifdef MULTITHREADING
    };
    static WorkerJob workerJobs[WORKER_THREADS];

    #define LOCAL_STORAGE(member) workerJobs[threadID].member
#else // !MULTITHREADING
    #define LOCAL_STORAGE(member) member
#endif

template <typename T> int sign(T x)
{
    return (T(0) < x) - (x < T(0));
}

enum SymmetryType
{
    SymmetryType_Chiral,
    SymmetryType_AchiralNonmirror,
    SymmetryType_AchiralMirror,
    SymmetryType_AchiralMirrorFace, // tetrahedral face as mirror plane
    SymmetryType_COUNT
};
static const char *const symmetryTypeString[SymmetryType_COUNT] =
{
    "",
    ", achiral",
    ", mirror",
    ", mirror2",
};

struct SymmetryRoot
{
    CompressedPolytetBits value;
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
    TetIndex tetI;
#endif
};

void canonicalizePolytet(Polytet &polytet, int tetCount, SymmetryRoot runningLeastPolytet[2], CompressedPolytet &newRotatedPolytet, bool &isChiral
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
    , SymmetryType &symmetryType
    , int &symmetry  // number of rotations under which the polytet is identical
#endif
    )
{
    runningLeastPolytet[1].value = runningLeastPolytet[0].value = COMPRESSEDPOLYTETBITS_MAX;
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
    TetIndex symmetricRootList[MAXIMUM_TETCOUNT * 3];
    CompressedPolytetBits lastValue = COMPRESSEDPOLYTETBITS_MAX;
#endif

    Tet *t = &polytet[1];
    for (int i=0; i<tetCount; i++)
    {
        static const uint8_t bitRotationTable[1 << 3][2 * 2] =
        {
            {0, 2,   0, 2}, // 000
            {0, 0,   2, 2}, // 001 -> 001
            {1, 1,   0, 0}, // 010 -> 001
            {0, 0,   0, 0}, // 011 -> 011
            {2, 2,   1, 1}, // 100 -> 001
            {2, 2,   2, 2}, // 101 -> 011
            {1, 1,   1, 1}, // 110 -> 011
            {0, 2,   0, 2}, // 111
        };
        const uint8_t *rotationStepRange;

        const RotationTable *thisRotationTable;
        {
            const Tet &singlyAttachedTet = polytet[i];
            if (!singlyAttachedTet.isLeaf)
                goto skipThisTet; // not a singly attached tet
            t = &polytet[singlyAttachedTet.faceAttached[3]];
            int attachedFace = singlyAttachedTet.faceAttachedFace[3];
            thisRotationTable = &rotationTable[attachedFace];
        }
        {
            // Calculate what the top 3 bits of the serialized value will be, in the same way as the top-level "append" call
            unsigned topValue = 0;
            for (int _faceNum=0; _faceNum<3; _faceNum++)
            {
                int rotatedFaceNum = faceRotateReflect[0][_faceNum][0];
                int faceNum = thisRotationTable->faceMap[rotatedFaceNum];
                if (t->faceAttached[faceNum] >= 0)
                    topValue |= 1 << _faceNum;
            }
            // Boost speed by canonicalizing the lower 3 bits to be only "001", "011", or "111".
            // This overrides the "least binary representation" criterion.
            rotationStepRange = bitRotationTable[topValue];
        }
        for (int reflect=0; reflect<2; reflect++, rotationStepRange+=2)
        {
            newRotatedPolytet.runningLeastValue = runningLeastPolytet[reflect].value;
            for (uint8_t rotationStep=rotationStepRange[0]; rotationStep<=rotationStepRange[1]; rotationStep++)
            {
                polytet.resetIndexing(i);
                newRotatedPolytet.value = 0;
                newRotatedPolytet.reflect = reflect;
                if (!newRotatedPolytet.append(*t, thisRotationTable, rotationStep))
                    continue;

                // Update the running "least" rotation
                if (newRotatedPolytet.runningLeastValue > newRotatedPolytet.value)
                {
                    newRotatedPolytet.runningLeastValue = newRotatedPolytet.value;
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                    runningLeastPolytet[reflect].tetI = i;
                    if (!reflect)
                    {
                        symmetricRootList[0] = i;
                        symmetry = 1;
                    }
#endif
                }
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                else
                if (!reflect && newRotatedPolytet.runningLeastValue == newRotatedPolytet.value)
                    symmetricRootList[symmetry++] = i;
#endif
            }
            runningLeastPolytet[reflect].value = newRotatedPolytet.runningLeastValue;
        }
    skipThisTet:;
    }

    isChiral = runningLeastPolytet[0].value != runningLeastPolytet[1].value;
    if (isChiral)
    {
        if (runningLeastPolytet[0].value > runningLeastPolytet[1].value)
            runningLeastPolytet[0].value = runningLeastPolytet[1].value;
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
        symmetryType = SymmetryType_Chiral;
#endif
    }
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
    else
    {
        if (runningLeastPolytet[0].tetI == runningLeastPolytet[1].tetI)
            symmetryType = SymmetryType_AchiralMirror;
        else
        {
            auto checkMirror = [&](TetIndex tetI) -> SymmetryType
            {
                polytet[                       tetI].tagSkipOverlapCheck(polytet, 1);
                CompressedSubpolytet path0to1 =                       (polytet[runningLeastPolytet[1].tetI].compressedPath - 1) / 3;
                polytet[runningLeastPolytet[1].tetI].tagSkipOverlapCheck(polytet, 1);
                CompressedSubpolytet path1to0 = reflectCompressedPath((polytet[                       tetI].compressedPath - 1) / 3);
                if (path1to0 != path0to1)
                    return SymmetryType_AchiralNonmirror;
                while (path0to1 > 3)
                    path0to1 = ((path0to1 - 1) / 3 - 1) / 3;
                return path0to1 ? SymmetryType_AchiralMirrorFace : SymmetryType_AchiralMirror;
            };
            for (int i=0; i<symmetry; i++)
                if (symmetryType = checkMirror(symmetricRootList[i]); symmetryType > SymmetryType_AchiralNonmirror)
                    goto foundIsMirror;
            symmetryType = SymmetryType_AchiralNonmirror;
            return;
        foundIsMirror:;
            if (symmetryType > SymmetryType_AchiralMirror && (symmetry == 2 || symmetry > 3))
                symmetryType = SymmetryType_AchiralMirror; // collapse "mirror2" into "mirror" for symmetries that have multiple mirror planes
        }
    }
#endif
}

void enumerate(
#ifdef MULTITHREADING
    THREAD_ID threadID,
    size_t &nextWorkAssignment,
    size_t workAssignmentLength,
#endif

    const Tetrahedron &start,
#ifdef USE_GMP
    const mpz_t maximalTouchingSqrDistance,
#else
    const Coord maximalTouchingSqrDistance,
#endif

#ifndef DISABLE_OVERLAP_CHECKING
    int8_t *overlapCache,
#endif

    const int tetCount,
    void *&pool,
    size_t &poolSize,
    const uint8_t *&basePolytetTable,
    const int basePolytetCompressedSize,
    const size_t polytetCount,

    HashIndex *&hashTable,
    const size_t hashTableSize,
    void *&polytetTable,
    const int newPolytetsCompressedSize,
    const int polytetTableElementSize,

    size_t &newPolytetCount,
    size_t &polytetChiralCount
#ifdef PRINT_SYMMETRY_TOTALS
    , size_t polytetSymmetryCount[MAXIMUM_TETCOUNT * 3][SymmetryType_COUNT]
#endif
    )
{
#ifdef SHOW_PROGRESS
#   ifdef MULTITHREADING
#       error SHOW_PROGRESS is currently incompatible with MULTITHREADING.
#   endif
    size_t nextProgressOutput = tetCount < SHOW_PROGRESS ? UINT64_MAX : 0;
    size_t progressOutputInterval = polytetCount / 1000;
#endif

    Polytet &polytet = LOCAL_STORAGE(workerPolytet);
    polytet.init(start);
    polytet.setSize(tetCount);

    CompressedPolytet newRotatedPolytet(polytet);

#ifdef MULTITHREADING
    for (;;)
#endif
    {
#ifdef MULTITHREADING
        size_t basePolytet0;
        size_t basePolytet1;
        {
            boost::mutex::scoped_lock lock(workAssignmentMutex);
            if (nextWorkAssignment >= polytetCount)
                break;
            nextWorkAssignment = basePolytet1 = std::min((basePolytet0 = nextWorkAssignment) + workAssignmentLength, polytetCount);
        }
        for (size_t basePolytetI = basePolytet0; basePolytetI < basePolytet1; basePolytetI++)
#else
        for (size_t basePolytetI=0; basePolytetI<polytetCount; basePolytetI++)
#endif
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
            {
                CompressedPolytetBits *basePolytet = (CompressedPolytetBits*)(basePolytetTable + basePolytetI * basePolytetCompressedSize);
                newRotatedPolytet.value = 0;
                memcpy(&newRotatedPolytet.value, basePolytet, basePolytetCompressedSize);
                newRotatedPolytet.uncompress();
            }

            Tet &newTet = polytet[tetCount - 1];
            for (int tetNumToAttachTo = 0; tetNumToAttachTo < tetCount-1; tetNumToAttachTo++)
            {
                Tet &tetToAttachTo = polytet[tetNumToAttachTo];
                for (int faceNum=0; faceNum<3; faceNum++) // skip last face because it's always already attached
                {
                    if (tetToAttachTo.faceAttached[faceNum] >= 0)
                        continue;
                    bool wasLeaf = tetToAttachTo.isLeaf;
                    polytet.attachNewTet(tetCount - 1, tetNumToAttachTo, faceNum);
                    // Canonicalize the rotation of this new polytet in compressed form, so that it can be compared against others
                    SymmetryRoot runningLeastPolytet[2];

                    bool isChiral;
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                    SymmetryType symmetryType;
                    int symmetry;
#endif
                    canonicalizePolytet(polytet, tetCount, runningLeastPolytet, newRotatedPolytet, isChiral
#if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                        , symmetryType
                        , symmetry
#endif
                        );

                    size_t hashIndex = hash(runningLeastPolytet[0].value) % hashTableSize;
                    HashIndex *index = &hashTable[hashIndex];
#ifdef MULTITHREADING
                    int shard = hashIndex % HASH_TABLE_SHARDS;
#endif
                    {
#ifdef MULTITHREADING
                        boost::mutex::scoped_lock lock(hashTableMutex[shard]);
#endif
                        for (;;)
                        {
                            if (*index == 0)
                                break; // no duplicate of runningLeastPolytet[0].value was found in hash table
                            void *entry = (uint8_t*)polytetTable + (size_t)(*index - 1) * polytetTableElementSize;
                            if (memcmp(entry, &runningLeastPolytet[0].value, newPolytetsCompressedSize) == 0)
                                goto skipDuplicate;
                            index = (HashIndex*)((uint8_t*)entry + newPolytetsCompressedSize);
                        }
                    }
#ifndef DISABLE_OVERLAP_CHECKING
                    // Check for overlap between this newly attached tetrahedron and the existing ones,
                    // and defer this until after the deduplication, to save a lot of time
                    newTet.tagSkipOverlapCheck(polytet, MIN_OVERLAP_DEPTH);
                    for (int tetCheckIntersectionI=0; tetCheckIntersectionI<polytet.size(); tetCheckIntersectionI++)
                    {
                        if (polytet[tetCheckIntersectionI].compressedPath == UINT64_MAX)
                            continue; // skip this check for speed (it'll always be false anyway)
                        CompressedSubpolytet compressedPath = (polytet[tetCheckIntersectionI].compressedPath - 1) / 3;
                        if (compressedPath > 0 && overlapCache[compressedPath - 1] > 0)
                            goto skipDueToOverlap;
                    }
#endif // DISABLE_OVERLAP_CHECKING
                    // No overlap found, so add runningLeastPolytet[0].value to hash table and chiral count
#if defined(USE_GMP) && defined(PRINT_POLYTETS) && !defined(PRINT_POLYTETS_WITH_SYMMETRY)
                    {
                        boost::mutex::scoped_lock lock(printPolytetMutex);
                        printPolytet(runningLeastPolytet[0].value, polytet);
                    }
#endif
                    {
#ifdef MULTITHREADING
                        boost::mutex::scoped_lock lock(hashTableMutex[shard]);
                        if (*index != 0)
                        {
                            for (;;)
                            {
                                void *entry = (uint8_t*)polytetTable + (size_t)(*index - 1) * polytetTableElementSize;
                                if (memcmp(entry, &runningLeastPolytet[0].value, newPolytetsCompressedSize) == 0)
                                    goto skipDuplicate;
                                index = (HashIndex*)((uint8_t*)entry + newPolytetsCompressedSize);
                                if (*index == 0)
                                    break; // no duplicate of runningLeastPolytet[0].value was found in hash table
                            }
                        }
                        size_t newPolytetCountFetched = __atomic_fetch_add(&newPolytetCount, 1, __ATOMIC_RELAXED);
                        void *entry = (uint8_t*)polytetTable + newPolytetCountFetched * polytetTableElementSize;
#else
                        void *entry = (uint8_t*)polytetTable + newPolytetCount        * polytetTableElementSize;
#endif
                        polytetChiralCount += isChiral;
                        if ((uint8_t*)entry + polytetTableElementSize - (uint8_t*)pool > poolSize)
                            quitMemory();
                        memcpy(entry, &runningLeastPolytet[0].value, newPolytetsCompressedSize);
#ifdef MULTITHREADING
                        *index = newPolytetCountFetched + 1;
#else
                        *index = ++newPolytetCount;
#endif
                        *(HashIndex*)((uint8_t*)entry + newPolytetsCompressedSize) = 0; // pointer to next hash collision
                    }
#ifdef PRINT_SYMMETRY_TOTALS
                    polytetSymmetryCount[symmetry - 1][symmetryType]++;
#endif
#ifdef PRINT_POLYTETS_WITH_SYMMETRY
    #ifndef PRINT_POLYTETS
                    if (!isChiral || symmetry > 1)
    #endif
                    {
                        boost::mutex::scoped_lock lock(printPolytetMutex);

                        printf("<%d%s>\n", symmetry, symmetryTypeString[symmetryType]);
                        printPolytet(runningLeastPolytet[0].value, polytet);
                    }
#endif
                skipDuplicate:
                skipDueToOverlap:

                    tetToAttachTo.faceAttached[faceNum] = -1;
                    tetToAttachTo.isLeaf = wasLeaf;
                }
            }

            polytet[1].faceAttached[0] = -1;
            polytet[1].faceAttached[1] = -1;
            polytet[1].faceAttached[2] = -1;
            polytet[1].isLeaf = true;
        }
    }
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
    for (int d=0; d<3; d++)
        mpz_inits(printCenter[d], printTmp[d], NULL);
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

#ifndef DISABLE_OVERLAP_CHECKING
    int8_t *overlapCache; // A bitmap (packed booleans) was tried, and was a bit slower; so, we'll use 8 times as much RAM to get slightly better speed.
                          // But now, we're benefitting from using bytes anyway, because during precomputing, a negative value indicates "no overlap".
                          // The array indices are branchless polytet attachment paths in bijective trinary, with 1 subtracted.
    size_t overlapCacheSize = 0;
    for (int i=0; i<MAXIMUM_TETCOUNT-2; i++)
        overlapCacheSize = overlapCacheSize * 3 + 1;
    overlapCache = (int8_t*)calloc(overlapCacheSize, 1);
    std::cout << "Allocated " << overlapCacheSize << " bytes for overlap caching" << std::endl;
#endif

    size_t poolSize;
    void *pool = NULL;

    size_t prevPolytetCount = 0;
    size_t polytetCount = 1;
    size_t memoryUsage = 0;

    int tetCount=1;
    bool resumedFromFile = false;
#ifdef MULTITHREADING
    size_t polytetChiralCount  [WORKER_THREADS];
#   ifdef PRINT_SYMMETRY_TOTALS
    size_t polytetSymmetryCount[WORKER_THREADS][MAXIMUM_TETCOUNT * 3][SymmetryType_COUNT];
#   endif
#else
    size_t polytetChiralCount;
#   ifdef PRINT_SYMMETRY_TOTALS
    size_t polytetSymmetryCount                [MAXIMUM_TETCOUNT * 3][SymmetryType_COUNT];
#   endif
#endif
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
#if defined(PRINT_POLYTETS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                mul_start_3(start, maximalTouchingSqrDistance, minUnseenCompressedSubpolytet);
#endif
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
            if (!readAndCloseOpenedFile(resumeFile, (uint8_t*)pool, size)
#ifndef DISABLE_OVERLAP_CHECKING
                || !readFile(filename = OVERLAP_BITMAP_FILENAME, (uint8_t*)overlapCache, overlapCacheSize)
#endif
                )
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

    decltype(startTime) currentTime;
#ifndef DISABLE_OVERLAP_CHECKING
    // Precompute overlapCache
    {
#if defined(PRINT_POLYTETS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
        auto start2 = start;
        auto maximalTouchingSqrDistance2 = maximalTouchingSqrDistance;
#else
        auto &start2 = start;
        auto &maximalTouchingSqrDistance2 = maximalTouchingSqrDistance;
#endif
        for (int i=5; i<MAXIMUM_TETCOUNT; i+=2)
            mul_start_3(start2, maximalTouchingSqrDistance2);
        Polytet polytet;
        polytet.init<true>(start);
        int tetCount = 2, topCount = 1, bottomCount = 1;
        int8_t faceNumStack[MAXIMUM_TETCOUNT];
        int stackPos = 0;
        faceNumStack[0] = -1;
        TetrahedronOverlap overlap;

        for (;;)
        {
            if (faceNumStack[stackPos] < 0)
            {
                int8_t foundOverlap = 0;
                if (tetCount > MIN_OVERLAP_DEPTH)
                {
                    Tet  &newEndTet = polytet[tetCount - 1];
                    Tet &prevEndTet = polytet[tetCount - 2];
                    // TODO: Autodetect MIN_OVERLAP_DEPTH
                    newEndTet.tagSkipOverlapCheck(polytet, MIN_OVERLAP_DEPTH);
                    CompressedSubpolytet compressedPath = (prevEndTet.compressedPath - 1) / 3;
                    foundOverlap = overlapCache[compressedPath - 1];
                    if (foundOverlap == 0)
                    {
                        auto reflectedCompressedPath = reflectCompressedPath(compressedPath);
                        overlap.setA( newEndTet);
                        overlap.setB(prevEndTet);
                        foundOverlap = overlap(maximalTouchingSqrDistance2);
                        if (foundOverlap)
                        {
#ifdef PRINT_POLYTETS_EDGE_CASES
                            if (foundOverlap < 0)
                            {
                                SymmetryRoot runningLeastPolytet[2];
                                bool isChiral;
    #if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                                SymmetryType symmetryType;
                                int symmetry;
    #endif
                                polytet.setSize(tetCount);
                                CompressedPolytet newRotatedPolytet(polytet);
                                canonicalizePolytet(polytet, tetCount, runningLeastPolytet, newRotatedPolytet, isChiral
    #if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                                    , symmetryType
                                    , symmetry
    #endif
                                    );
                                printf("%d-cell "
    #if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                                    "<%d%s> "
    #endif
                                    , tetCount
    #if defined(PRINT_SYMMETRY_TOTALS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
                                    , symmetry, symmetryTypeString[symmetryType]
    #endif
                                    );
                                printPolytet(runningLeastPolytet[0].value, polytet);
                            }
#endif
                            overlapCache[         compressedPath - 1] = 1;
                            overlapCache[reflectedCompressedPath - 1] = 1;
                        }
                        else
                        {
                            overlapCache[         compressedPath - 1] = -1;
                            overlapCache[reflectedCompressedPath - 1] = -1;
                        }
                    }
                }
                if (foundOverlap)
                    goto backtrackOverlapChecking;
                faceNumStack[stackPos]++;
            }
            if (tetCount == MAXIMUM_TETCOUNT || faceNumStack[stackPos] == 3)
            {
            backtrackOverlapChecking:
                if (--stackPos < 0)
                    break;
                if (bottomCount > topCount) bottomCount--;
                else                           topCount--;
                tetCount--;
                polytet[tetCount - 2].faceAttached[faceNumStack[stackPos]] = -1;
                polytet[tetCount - 2].isLeaf = true;
                faceNumStack[stackPos]++;
            }
            else
            {
                bool attachToTop = bottomCount > topCount;
                int iAttachTo = tetCount - 2;
                if (attachToTop)    topCount++;
                else             bottomCount++;
                polytet.attachNewTet<true>(tetCount++, iAttachTo, faceNumStack[stackPos]);
                faceNumStack[++stackPos] = -1;
            }
        }
    }
#ifdef WRITE_TO_FILES
    writeFile(OVERLAP_BITMAP_FILENAME, (uint8_t*)overlapCache, overlapCacheSize);
#endif
    currentTime = std::chrono::steady_clock::now();
    std::cout << "Precomputed overlap cache [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms]" << std::endl;
#endif // !DISABLE_OVERLAP_CHECKING

#ifdef MULTITHREADING
    std::thread workers[WORKER_THREADS];

    memset(polytetChiralCount, 0, sizeof(polytetChiralCount));
#else
    polytetChiralCount = 0;
#endif
#ifdef PRINT_SYMMETRY_TOTALS
    memset(polytetSymmetryCount, 0, sizeof(polytetSymmetryCount));
#endif
    for (;;)
    {
#ifdef PRINT_POLYTETS_WITH_SYMMETRY
        if (tetCount < 3)
        {
            if (tetCount == 1)
                printf("<12, mirror>\n"
                    "-\n"
                    "{{{-1, -1, -1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}}}\n\n");
            else // tetCount == 2
                printf("<6, mirror>\n"
                    "0x0\n"
                    "{{{-4, -4, -4}, {-4, 2, 2}, {2, -4, 2}, {2, 2, -4}},\n"
                    "{{4, 4, 4}, {-4, 2, 2}, {2, 2, -4}, {2, -4, 2}}}\n\n");
        }
#endif

        currentTime = std::chrono::steady_clock::now();
        std::cout << tetCount << ": ";
        if (resumedFromFile)
            std::cout << "resumed";
        else
        {
            size_t polytetCountOneSided = polytetCount;
#ifdef MULTITHREADING
            for (THREAD_ID threadID=0; threadID < WORKER_THREADS; threadID++)
                polytetCountOneSided += polytetChiralCount[threadID];
#else
            polytetCountOneSided += polytetChiralCount;
#endif
            std::cout << polytetCountOneSided;
        }
        std::cout << " (" << polytetCount << ") [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms";
        if (memoryUsage)
            std::cout << ", " << memoryUsage << " bytes";
        std::cout << "]" << std::endl;
#ifdef PRINT_SYMMETRY_TOTALS
        if (tetCount < 3)
        {
            if (tetCount == 1)
                printf("    1 <12, mirror>\n");
            else // tetCount == 2
                printf("    1 <6, mirror>\n");
        }
        else
        {
            for (int i=0; i < MAXIMUM_TETCOUNT * 3; i++)
            {
                for (int symmetryType=0; symmetryType < SymmetryType_COUNT; symmetryType++)
                {
                    /*if (i == 1-1 && symmetryType == SymmetryType_Chiral)
                        continue;*/
                    size_t symmetryCount = 0;
    #ifdef MULTITHREADING
                    for (THREAD_ID threadID=0; threadID < WORKER_THREADS; threadID++)
    #endif
                    {
                        symmetryCount += polytetSymmetryCount
    #ifdef MULTITHREADING
                            [threadID]
    #endif
                            [i][symmetryType];
                    }
                    if (symmetryCount)
                        printf("    %zu <%d%s>\n", symmetryCount, i+1, symmetryTypeString[symmetryType]);
                }
            }
            fflush(stdout);
        }
#endif
#ifdef PRINT_POLYTETS_WITH_SYMMETRY
        for (int i=0; i<120; i++) putchar('-'); putchar('\n');
#endif
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
        {
            basePolytetTable = (uint8_t*)pool;
            size_t basePolytetTableSize = basePolytetCompressedSize * polytetCount;
            hashTable = (HashIndex*)(basePolytetTable + basePolytetTableSize);
            polytetTable = hashTable + hashTableSize;
            size_t minSize = (uint8_t*)polytetTable - (uint8_t*)pool;
            if (minSize >= poolSize)
                quitMemory();
        }
        memset(hashTable, 0, hashTableSize * sizeof(HashIndex));
        const int newPolytetsCompressedSize = ((tetCount - 2) * 3 + 8-1) / 8;
        const int polytetTableElementSize = newPolytetsCompressedSize + sizeof(HashIndex);
        size_t newPolytetCount = 0;

        size_t nextWorkAssignment = 0;
        size_t workAssignmentLength = polytetCount / WORKER_THREADS;
        if (workAssignmentLength < 1)
            workAssignmentLength = 1;
        else
        if (workAssignmentLength > MAXIMUM_WORK_ASSIGNMENT)
            workAssignmentLength = MAXIMUM_WORK_ASSIGNMENT;
#ifdef PRINT_SYMMETRY_TOTALS
        memset(polytetSymmetryCount, 0, sizeof(polytetSymmetryCount));
#endif
#ifdef MULTITHREADING
        memset(polytetChiralCount, 0, WORKER_THREADS * sizeof(*polytetChiralCount));
        for (THREAD_ID threadID=0; threadID < WORKER_THREADS; threadID++)
        {
            workers[threadID] = std::thread(enumerate,
                threadID, std::ref(nextWorkAssignment), workAssignmentLength,
                std::ref(start), maximalTouchingSqrDistance,
#ifndef DISABLE_OVERLAP_CHECKING
                overlapCache,
#endif
                tetCount, std::ref(pool), std::ref(poolSize), std::ref((const uint8_t *&)basePolytetTable), basePolytetCompressedSize, polytetCount,
                std::ref(hashTable), hashTableSize, std::ref(polytetTable), newPolytetsCompressedSize, polytetTableElementSize,
                std::ref(newPolytetCount), std::ref(polytetChiralCount[threadID])
    #ifdef PRINT_SYMMETRY_TOTALS
                , polytetSymmetryCount[threadID]
    #endif
                );
        }
        for (THREAD_ID threadID=0; threadID < WORKER_THREADS; threadID++)
            workers[threadID].join();
#else
        polytetChiralCount = 0;
        enumerate(
            start, maximalTouchingSqrDistance,
#ifndef DISABLE_OVERLAP_CHECKING
            overlapCache,
#endif
            tetCount, pool, poolSize, (const uint8_t *&)basePolytetTable, basePolytetCompressedSize, polytetCount,
            hashTable, hashTableSize, polytetTable, newPolytetsCompressedSize, polytetTableElementSize,
            newPolytetCount, polytetChiralCount
    #ifdef PRINT_SYMMETRY_TOTALS
            , polytetSymmetryCount
    #endif
            );
#endif

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
#if defined(PRINT_POLYTETS) || defined(PRINT_POLYTETS_WITH_SYMMETRY)
        mul_start_3(start, maximalTouchingSqrDistance, minUnseenCompressedSubpolytet);
#endif
    }

errorQuit:
    free(pool);
#ifndef DISABLE_OVERLAP_CHECKING
    free(overlapCache);
#endif
#ifdef USE_GMP
    mpz_clear(maximalTouchingSqrDistance);
    for (int d=0; d<3; d++)
        mpz_clears(printCenter[d], printTmp[d], NULL);
#endif
    return 0;
}
