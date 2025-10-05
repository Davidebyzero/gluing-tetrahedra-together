#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cstdint>
#include <cinttypes>
#include <array>
#include <algorithm>
#include <vector>
#include <list>
#include <functional>
#include <unordered_set>
#include <chrono>

//#define USE_GMP
//#define DEBUG_PRINT

#ifdef USE_GMP
#include <gmp.h>
#   if GMP_NUMB_BITS != 64
#   error This is hard-coded for 64-bit limbs
#   endif
#endif

auto startTime = std::chrono::steady_clock::now();

typedef uint8_t TetIndex;
typedef uint8_t TetIndexFace; // lowest 2 bits are used for a face index
typedef __int128 Coord;
typedef std::array<Coord, 3> Coord3;
typedef std::array<Coord3, 4> Tetrahedron;
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
    Tet *faceAttached[4];
    TetIndex index; // 1-based; 0=unassigned
    Tet(                    ) : t( ) {initFaces();}
    Tet(const Tetrahedron &t) : t(t) {initFaces();}
    void assignIndex(TetIndex &nextIndex)
    {
        if (index == 0)
            index = nextIndex++;
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
static int tetrahedronFaces[4][4] =
{
    {0, 1, 2, 3}, // 3,0,2,1   0->3, 1->1, 2->0, 3->2
    {0, 3, 1, 2}, // 2,0,1,3   0->0, 1->3, 2->1, 3->2
    {0, 2, 3, 1}, // 1,0,3,2   0->1, 1->0, 2->3, 3->2
    {1, 3, 2, 0}, // 0,1,2,3   0->0, 1->1, 2->2, 3->3
};

static int tetrahedronEdges[6][2] =
{
    {0, 1},
    {1, 2},
    {2, 0},
    {0, 3},
    {1, 3},
    {2, 3},
};

Coord3 operator-(const Coord3 &a)
{
    Coord3 c;
    c[0] = -a[0];
    c[1] = -a[1];
    c[2] = -a[2];
    return c;
}
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
Coord3 cross(const Coord3 a, const Coord3 b)
{
    Coord3 c;
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
}
Coord dot(const Coord3 &a, const Coord3 &b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
bool volumesOverlap(const Tetrahedron &a, const Tetrahedron &b)
{
    // Take advantage of the fact that the two tetrahedrons are regular and congruent, and
    // just check if any edge of tetrahedron "a" intersects with any face of tetrahedron "b".
    // Don't count it if only the endpoint of an edge intersects.
    for (int edgeNum=0; edgeNum<6; edgeNum++)
    {
        Coord3 p0 = a[tetrahedronEdges[edgeNum][0]];
        Coord3 p1 = a[tetrahedronEdges[edgeNum][1]];
        for (int faceNum=0; faceNum<4; faceNum++)
        {
            Tetrahedron normalizedTetrahedron; // first 3 points are the face, and the 4th point is for calculating the normal
            for (int i=0; i<4; i++)
                normalizedTetrahedron[i] = b[tetrahedronFaces[faceNum][i]];
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
            Coord3 d = intersectionPoint - triangle[0];
            Coord3 edge1 = triangle[1] - triangle[0];
            Coord3 edge2 = triangle[2] - triangle[0];
            Coord uNumerator = d[1]*edge2[0] - d[0]*edge2[1];
            Coord vNumerator = d[0]*edge1[1] - d[1]*edge1[0];
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
void attachNewTet(Tet &t, Tet &tetToAttachTo, const int faceNum)
{
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
    // Copy the other vertices
    for (int p=0; p<3; p++)
    {
        int p1 = tetrahedronFaces[faceNum][p];
        for (int d=0; d<3; d++)
            t.t[1+p][d] = tetToAttachTo.t[p1][d];
    }
    t.faceAttached[0] = NULL;
    t.faceAttached[1] = NULL;
    t.faceAttached[2] = NULL;
    t.faceAttached[3] = &tetToAttachTo;
    tetToAttachTo.faceAttached[faceNum] = &t;
}

class NormalizedPolytet : public std::vector<Tetrahedron>
{
public:
    NormalizedPolytet(const Polytet &polytet)
    {
        reserve(polytet.size());
        for (auto tet=polytet.cbegin(); tet!=polytet.cend(); ++tet)
        {
            Tetrahedron &t = emplace_back(tet->t);
            std::sort(t.begin(), t.end());
        }
        std::sort(this->begin(), this->end());
    }
};

bool operator<(const Polytet &_a, const Polytet &_b)
// implicitly assume a.size()==b.size()
{
    NormalizedPolytet a(_a), b(_b);
    for (auto ta=a.cbegin(), tb=b.cbegin(); ta!=a.cend(); ++ta,++tb)
    {
        auto result = *ta <=> *tb;
        if (result != 0)
            return result < 0;
    }
    return false;
}
bool operator>(const Polytet &a, const Polytet &b)
{
    return !(a < b);
}
bool operator==(const Polytet &_a, const Polytet &_b)
// implicitly assume a.size()==b.size()
{
    NormalizedPolytet a(_a), b(_b);
    for (auto ta=a.cbegin(), tb=b.cbegin(); ta!=a.cend(); ++ta,++tb)
    {
        auto result = *ta <=> *tb;
        if (result != 0)
            return false;
    }
    return true;
}

// First two tetrahedrons are implied. Each element is a subsequent tetrahedron, with the value indicating where
// it's attached. The lower 2 bits indicate which face (can only have 3 different values, because at least 1 face
// will always already be attached). The remaining bits indicate which tetrahedron (which can never be zero,
// because that one is attached implicitly).
class CompressedPolytet : public std::vector<TetIndexFace>
{
public:
    void append(Polytet &polytet, Tet &tetToCompress, int vertexMap[4], int faceRotation)
    // indices of vertexMap[] are compressed-output vertices; elements of vertexMap[] are the original vertices of tetToCompress
    {
        tetToCompress.assignIndex(polytet.nextIndex);
        for (int _faceNum=0; _faceNum<3; _faceNum++)
        {
            int rotatedFaceNum = (_faceNum + faceRotation) % 3;
            int faceNum = 3 - vertexMap[3 - rotatedFaceNum];
            if (!tetToCompress.faceAttached[faceNum])
                continue;
            Tet *attachedTet = tetToCompress.faceAttached[faceNum];
            attachedTet->assignIndex(polytet.nextIndex);
            push_back((((TetIndexFace)(tetToCompress.index - 1)) << 2) + _faceNum);
            
            int attachedFace = 0;
            while (attachedTet->faceAttached[attachedFace] != &tetToCompress)
                attachedFace++;
            int vertexMap2[4];
            int rotation = 0;
            while (vertexMap[tetrahedronFaces[rotatedFaceNum][rotation]] != tetrahedronFaces[faceNum][0])
                rotation++;
            vertexMap2[1] = tetrahedronFaces[attachedFace][0];
            vertexMap2[3] = tetrahedronFaces[attachedFace][1];
            vertexMap2[2] = tetrahedronFaces[attachedFace][2];
            vertexMap2[0] = tetrahedronFaces[attachedFace][3];
            if (faceNum == 3)
            {
                int tmp = vertexMap2[3];
                vertexMap2[3] = vertexMap2[2];
                vertexMap2[2] = tmp;
            }
            append(polytet, *attachedTet, vertexMap2, rotation);
        }
    }
};

void verifyTetrahedron(const Tet &t, Coord power9_2)
{
    Coord targetDistSquared = power9_2 * 4;
    for (int edgeNum=0; edgeNum<6; edgeNum++)
    {
        Coord3 p0 = t.t[tetrahedronEdges[edgeNum][0]];
        Coord3 p1 = t.t[tetrahedronEdges[edgeNum][1]];
        Coord3 diff = p1 - p0;
        Coord distSquared = 0;
        for (int d=0; d<3; d++)
            distSquared += diff[d] * diff[d];
        if (distSquared != targetDistSquared)
        {
            auto currentTime = std::chrono::steady_clock::now();
            std::cerr << "Quitting due to detected overflow" << " [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms]" << std::endl;
            exit(-1);
        }
    }
}

namespace std
{
    template<>
    struct hash<Polytet>
    {
        std::size_t operator()(const Polytet &_polytet) const noexcept
        {
            NormalizedPolytet polytet(_polytet);
            std::size_t seed = polytet.size();
            for (auto t=polytet.cbegin(); t!=polytet.cend(); ++t)
            {
                for (auto c=t->cbegin(); c!=t->cend(); ++c)
                {
                    for (auto i=c->cbegin(); i!=c->cend(); ++i)
                    {
                        seed ^= std::hash<uint64_t>{}((uint64_t)(*i      )) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                        seed ^= std::hash<uint64_t>{}((uint64_t)(*i >> 64)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                    }
                }
            }
            return seed;
        }
    };
    template<>
    struct hash<CompressedPolytet>
    {
        std::size_t operator()(const CompressedPolytet &polytet) const noexcept
        {
            std::size_t seed = polytet.size();
            for (auto i=polytet.cbegin(); i!=polytet.cend(); ++i)
                seed ^= std::hash<uint32_t>{}(*i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
};

#ifdef USE_GMP
void mpz_set_int128(mpz_t &dst, const __int128 &src)
{
    __int128 abssrc = src < 0 ? -src : src;
    mpz_import(dst, 2, -1, 4, 0, 0, &abssrc);
    if (src < 0)
        mpz_neg(dst, dst);
}
void mpz_get_int128(__int128 &dst, const mpz_t &src)
{
    dst = 0;
    size_t count;
    uint64_t limbs[2];

    mpz_export(limbs, &count, -1, sizeof(uint64_t), 0, 0, src);

    ((int64_t*)&dst)[0] = limbs[0];
    if (count > 1)
        ((int64_t*)&dst)[1] = limbs[1];

    if (mpz_sgn(src) < 0)
        dst = -dst;
}
#endif

int main(int argc, char *argv[])
{
#if 1 // main
    static Tetrahedron start =
    {{
        {{-9,-9,-9}},
        {{-9, 9, 9}},
        {{ 9,-9, 9}},
        {{ 9, 9,-9}}
    }};

    auto *polytets = new std::unordered_set<CompressedPolytet>;
    polytets->insert(CompressedPolytet()); // add empty vector as the starter polytet (meaning it has two tetrahedrons)
    size_t prevPolytetCount = 0;
    
    size_t blahNum = 0;

    for (int tetCount=1;;)
    {
        auto currentTime = std::chrono::steady_clock::now();
        size_t polytetCount = polytets->size();
        std::cout << tetCount << ": " << polytetCount << " [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms]" << std::endl;
        if (prevPolytetCount > polytetCount)
        {
            std::cerr << "Quit due to apparent overflow" << std::endl;
            break;
        }
        prevPolytetCount = polytetCount;
        if (++tetCount <= 2)
            continue;
        /*if (tetCount > 6)
            break;*/

        Polytet polytet;
        polytet.reserve(tetCount); // Important, to ensure pointers don't change
        Tet &t0    = polytet.emplace_back(start);
        attachNewTet(polytet.emplace_back(), t0, 3);

        auto *newPolytets = new std::unordered_set<CompressedPolytet>;
        polytet.resize(tetCount);
        for (auto basePolytet=polytets->cbegin(); basePolytet!=polytets->cend(); ++basePolytet)
        {
            int tetNumToUncompress = 2;
            for (auto elementToUncompress=basePolytet->cbegin(); elementToUncompress!=basePolytet->cend(); ++elementToUncompress)
            {
                int faceNum          = *elementToUncompress & 3;
                int tetNumToAttachTo = *elementToUncompress >> 2;
                Tet &tetToAttachTo = polytet[tetNumToAttachTo];
                attachNewTet(polytet[tetNumToUncompress++], tetToAttachTo, faceNum);
            }
            if (tetNumToUncompress != tetCount - 1)
            {
                std::cerr << "Error! Got " << tetNumToUncompress << ", expected " << tetCount - 1 << std::endl;
                exit(-1);
            }

            Tet &newTet = polytet[tetCount - 1];
            for (int tetNumToAttachTo = 1; tetNumToAttachTo < tetCount-1; tetNumToAttachTo++)
            {
                Tet &tetToAttachTo = polytet[tetNumToAttachTo];
                for (int faceNum=0; faceNum<3; faceNum++) // skip last face because it's always already attached
                {
                    if (tetToAttachTo.faceAttached[faceNum])
                        continue;
                    attachNewTet(newTet, tetToAttachTo, faceNum);
                    // Check for overlap between this newly attached tetrahedron and the existing ones
                    for (auto tetCheckIntersection=polytet.cbegin(); tetCheckIntersection!=polytet.cend(); ++tetCheckIntersection)
                    {
                        if (&*tetCheckIntersection == &tetToAttachTo || &*tetCheckIntersection == &newTet)
                            continue; // skip this check for speed (it'll always be false anyway)
                        if (volumesOverlap(newTet.t, tetCheckIntersection->t))
                            goto discardThisNewPolytet;
                    }
                    // Canonicalize the rotation of this new polytet in compressed form, so that it can be compared against others
#if 0
                    {
                        polytet.resetIndexing(0);
                        CompressedPolytet newRotatedPolytet;
                        newRotatedPolytet.reserve(tetCount - 2);
                        static int vertexMap[4] = {0, 1, 2, 3};
                        newRotatedPolytet.append(polytet, polytet[1], vertexMap);
                        newPolytets->insert(newRotatedPolytet);
                    }
#else
                    {
                        bool haveRunningLeast = false;
                        CompressedPolytet runningLeastPolytet;

                        int vertexMap[4] = {0, 1, 2, 3};
                        Tet *t = &polytet[1];
                        for (int i=1; i<tetCount; i++)
                        {
                            int attachedFace;
                            if (i > 1)
                            {
                                Tet &singlyAttachedTet = polytet[i];
                                for (int j=0; j<3; j++)
                                    if (singlyAttachedTet.faceAttached[j])
                                        goto skipThisTet; // not a singly attached tet
                                t = singlyAttachedTet.faceAttached[3];
                                attachedFace = 0;
                                while (t->faceAttached[attachedFace] != &singlyAttachedTet)
                                    attachedFace++;
                                static int vertexMapTable[3][4] =
                                {
                                    {3, 0, 2, 1},
                                    {2, 0, 1, 3},
                                    {1, 0, 3, 2},
                                };
                                memcpy(vertexMap, vertexMapTable[attachedFace], sizeof(vertexMap));
                            }
                            for (int rotationStep=0; rotationStep<3; rotationStep++)
                            {
                                polytet.resetIndexing(i > 1 ? i : 0);
                                CompressedPolytet newRotatedPolytet;
                                newRotatedPolytet.reserve(tetCount - 2);
                                newRotatedPolytet.append(polytet, *t, vertexMap, rotationStep);

                                /*if (i==1 && rotationStep==0)
                                {
                                    std::cerr << "[" << blahNum++ << std::endl;
                                    for (int i=0; i<polytet.size(); i++)
                                    {
                                        std::cerr << "  " << i << " (" << (int)polytet[i].index - 1 << "):" << std::endl;
                                        for (int j=0; j<4; j++)
                                            if (polytet[i].faceAttached[j])
                                                std::cerr << "    face" << j << " -> " << (polytet[i].faceAttached[j] - polytet.data()) << std::endl;
                                    }
                                    std::cerr << "]"  << std::endl;
                                }
                                
                                std::cerr << "{ 1->0[3] ";
                                int tetNumToUncompress = 2;
                                for (size_t i=0; i<newRotatedPolytet.size(); i++)
                                {
                                    int faceNum          = newRotatedPolytet[i] & 3;
                                    int tetNumToAttachTo = newRotatedPolytet[i] >> 2;
                                    std::cerr << i+2 << "->" << tetNumToAttachTo << "[" << faceNum << "] ";
                                }
                                std::cerr << "}" << std::endl;*/
                                
                                // Update the running "least" rotation
                                if (!haveRunningLeast ||
                                    std::lexicographical_compare(
                                        runningLeastPolytet.begin(), runningLeastPolytet.end(),
                                        newRotatedPolytet  .begin(), newRotatedPolytet  .end()))
                                {
                                    haveRunningLeast = true;
                                    runningLeastPolytet = newRotatedPolytet;
                                    //std::cerr << "=" << std::endl;
                                }
                                // Switch to the next rotation
                                /*if (i > 1)
                                {
                                    int tmp = vertexMap[tetrahedronFaces[attachedFace][1]];
                                    vertexMap[tetrahedronFaces[attachedFace][1]] = vertexMap[tetrahedronFaces[attachedFace][2]];
                                    vertexMap[tetrahedronFaces[attachedFace][2]] = vertexMap[tetrahedronFaces[attachedFace][0]];
                                    vertexMap[tetrahedronFaces[attachedFace][0]] = tmp;
                                }
                                else
                                {
                                    int tmp = vertexMap[3];
                                    vertexMap[3] = vertexMap[2];
                                    vertexMap[2] = vertexMap[1];
                                    vertexMap[1] = tmp;
                                }*/
                            }
                        skipThisTet:;
                        }
                        newPolytets->insert(runningLeastPolytet);
                        //std::cerr << "+" << std::endl;
                    }
#endif
                discardThisNewPolytet:
                    tetToAttachTo.faceAttached[faceNum] = NULL;
                }
            }

            polytet[1].faceAttached[0] = NULL;
            polytet[1].faceAttached[1] = NULL;
            polytet[1].faceAttached[2] = NULL;
        }

        delete polytets;
        polytets = newPolytets;

        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                start[p][d] *= 3;
    }
#else // main
    static const Tetrahedron start =
    {{
        {{-1,-1,-1}},
        {{-1, 1, 1}},
        {{ 1,-1, 1}},
        {{ 1, 1,-1}}
    }};
    Polytet startPolytet;
    startPolytet.emplace_back(start);

#ifdef USE_GMP
    mpz_t mpz_power3, mpz_power9_2;
    mpz_init(mpz_power3  ); mpz_set_ui(mpz_power3  , 1);
    mpz_init(mpz_power9_2); mpz_set_ui(mpz_power9_2, 2);
    mpz_t
        x , y , z ,
        x0, y0, z0,
        x1, y1, z1,
        x2, y2, z2,
        x3, y3, z3,
        xx_numerator_mpz, xx_denominator_mpz,
        xy_numerator_mpz, xy_denominator_mpz,
        xz_numerator_mpz, xz_denominator_mpz,
        yx_numerator_mpz, yx_denominator_mpz,
        yy_numerator_mpz, yy_denominator_mpz,
        yz_numerator_mpz, yz_denominator_mpz,
        zx_numerator_mpz, zx_denominator_mpz,
        zy_numerator_mpz, zy_denominator_mpz,
        zz_numerator_mpz, zz_denominator_mpz;
    mpz_inits(
        x , y , z ,
        x0, y0, z0,
        x1, y1, z1,
        x2, y2, z2,
        x3, y3, z3,
        xx_numerator_mpz, xx_denominator_mpz,
        xy_numerator_mpz, xy_denominator_mpz,
        xz_numerator_mpz, xz_denominator_mpz,
        yx_numerator_mpz, yx_denominator_mpz,
        yy_numerator_mpz, yy_denominator_mpz,
        yz_numerator_mpz, yz_denominator_mpz,
        zx_numerator_mpz, zx_denominator_mpz,
        zy_numerator_mpz, zy_denominator_mpz,
        zz_numerator_mpz, zz_denominator_mpz, NULL);
#endif
    Coord power3   = 1;
    Coord power9_2 = 2;
    auto *polytets = new std::unordered_set<Polytet>;
    polytets->insert(startPolytet);
    size_t prevPolytetCount = 0;
    for (int tetCount=1;;)
    {
        auto currentTime = std::chrono::steady_clock::now();
        size_t polytetCount = polytets->size();
        std::cout << tetCount << ": " << polytetCount << " [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms]" << std::endl;
        if (prevPolytetCount > polytetCount)
        {
            std::cerr << "Quit due to apparent overflow" << std::endl;
            break;
        }
        prevPolytetCount = polytetCount;
        ++tetCount;
#ifdef DEBUG_PRINT
        if (tetCount > 6)
            break;
#endif
#ifdef USE_GMP
        mpz_mul_ui(mpz_power3  , mpz_power3  , 3);
        mpz_mul_ui(mpz_power9_2, mpz_power9_2, 9);
#endif
        power3   *= 3;
        power9_2 *= 9;
        auto *newPolytets = new std::unordered_set<Polytet>;
        for (auto basePolytet=polytets->cbegin(); basePolytet!=polytets->cend(); ++basePolytet)
        {
            for (auto tetToAttachTo=basePolytet->cbegin(); tetToAttachTo!=basePolytet->cend(); ++tetToAttachTo)
            {
                for (int faceNum=0; faceNum<4; faceNum++)
                {
                    if (!tetToAttachTo->faceAttached[faceNum])
                    {
                        Polytet newPolytet;
                        newPolytet.reserve(tetCount);
                        Tet *tetCopyToAttachTo; // take note of this just for speed
                        for (auto tetToCopy=basePolytet->cbegin(); tetToCopy!=basePolytet->cend(); ++tetToCopy)
                        {
                            Tet &t = newPolytet.emplace_back();
                            for (int p=0; p<4; p++)
                            for (int d=0; d<3; d++)
                                t.t[p][d] = tetToCopy->t[p][d] * 3;
                            memcpy(t.faceAttached, tetToCopy->faceAttached, sizeof(t.faceAttached));
                            if (tetToCopy == tetToAttachTo)
                            {
                                t.faceAttached[faceNum] = true;
                                tetCopyToAttachTo = &t;
                            }
                        }
                        Tet &t = newPolytet.emplace_back();
                        attachNewTet(t, *tetToAttachTo, faceNum);
                        // Check for overlap between this newly attached tetrahedron and the existing ones
                        for (auto tetCheckIntersection=newPolytet.cbegin(); tetCheckIntersection!=newPolytet.cend(); ++tetCheckIntersection)
                        {
                            if (&*tetCheckIntersection == tetCopyToAttachTo || &*tetCheckIntersection == &t)
                                continue; // skip this check for speed (it'll always be false anyway)
                            if (volumesOverlap(t.t, tetCheckIntersection->t))
                                goto discardThisNewPolytet;
                        }
                        // Canonicalize the rotation of this new polytet, so that it can be compared against others
                        bool haveRunningLeast = false;
                        Polytet runningLeastPolytet;
                        for (auto tetToRotateNormalize=newPolytet.cbegin(); tetToRotateNormalize!=newPolytet.cend(); ++tetToRotateNormalize)
                        {
                            // only rotate-normalize to tetrahedrons with exactly 1 attached face; every polytet is guaranteed to have some
                            int faceToRotateNormalize;
                            {
                                int attachCount = 0;
                                for (int i=0; i<4; i++)
                                    if (tetToRotateNormalize->faceAttached[i])
                                    {
                                        faceToRotateNormalize = i;
                                        if (++attachCount > 1)
                                            break;
                                    }
                                if (attachCount > 1)
                                    continue;
                            }
                            Tetrahedron n;
                            for (int i=0; i<4; i++)
                                n[i] = tetToRotateNormalize->t[tetrahedronFaces[faceToRotateNormalize][i]];
                            for (int faceRotation=0; faceRotation<3; faceRotation++)
                            {
                                // Solve for an affine transformation that rotates and translates the polytet so that "tetToRotateNormalize" goes into the same position as the first face of the "start" tetrahedron
                                Coord xx_numerator =   ( n[2][1] - n[3][1])*(n[0][2] - n[1][2])  - (n[0][1] - n[1][1])*(n[2][2] - n[3][2]) ; Coord xx_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord xy_numerator = (-((n[2][0] - n[3][0])*(n[0][2] - n[1][2])) + (n[0][0] - n[1][0])*(n[2][2] - n[3][2])); Coord xy_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord xz_numerator =   ( n[2][0] - n[3][0])*(n[0][1] - n[1][1])  - (n[0][0] - n[1][0])*(n[2][1] - n[3][1]) ; Coord xz_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord yx_numerator =  -((n[1][1] - n[3][1])*(n[0][2] - n[2][2])) + (n[0][1] - n[2][1])*(n[1][2] - n[3][2]) ; Coord yx_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord yy_numerator =   ( n[1][0] - n[3][0])*(n[0][2] - n[2][2])  - (n[0][0] - n[2][0])*(n[1][2] - n[3][2]) ; Coord yy_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord yz_numerator =  -((n[1][0] - n[3][0])*(n[0][1] - n[2][1])) + (n[0][0] - n[2][0])*(n[1][1] - n[3][1]) ; Coord yz_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord zx_numerator =  -((n[0][1] - n[3][1])*(n[1][2] - n[2][2])) + (n[1][1] - n[2][1])*(n[0][2] - n[3][2]) ; Coord zx_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord zy_numerator =   ( n[0][0] - n[3][0])*(n[1][2] - n[2][2])  - (n[1][0] - n[2][0])*(n[0][2] - n[3][2]) ; Coord zy_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord zz_numerator =   ( n[0][0] - n[3][0])*(n[1][1] - n[2][1])  - (n[1][0] - n[2][0])*(n[0][1] - n[3][1]) ; Coord zz_denominator = ( n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] - n[0][0]*n[3][1]*n[2][2] + n[1][0]*n[3][1]*n[2][2] + n[3][0]*(n[1][1]*n[0][2] - n[2][1]*n[0][2] - n[0][1]*n[1][2] + n[2][1]*n[1][2] + n[0][1]*n[2][2] - n[1][1]*n[2][2]) + n[1][0]*n[0][1]*n[3][2] - n[0][0]*n[1][1]*n[3][2] + n[0][0]*n[2][1]*n[3][2] - n[1][0]*n[2][1]*n[3][2] + n[2][0]*(n[3][1]*n[0][2] + n[0][1]*n[1][2] - n[3][1]*n[1][2] - n[0][1]*n[3][2] + n[1][1]*(-n[0][2] + n[3][2])));
                                Coord x0_ = ( n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(n[2][1]*(n[0][2] - n[1][2]) + n[1][1]*(n[0][2] + n[2][2]) - n[0][1]*(n[1][2] + n[2][2])) + n[1][0]*(n[0][1] + n[2][1])*n[3][2] - n[0][0]*(n[1][1] + n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(-n[0][2] + n[1][2]) - n[1][1]*(n[0][2] + n[3][2]) + n[0][1]*(n[1][2] + n[3][2])));
                                Coord y0_ = ( n[1][0]*n[2][1]*n[0][2] + n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] - n[3][0]*((n[1][1] + n[2][1])*n[0][2] + (-n[0][1] + n[2][1])*n[1][2]) - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] + n[3][0]*(n[0][1] + n[1][1])*n[2][2] - n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] + (n[0][0] + n[1][0])*n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(n[0][2] + n[1][2]) + n[0][1]*(n[1][2] - n[3][2]) - n[1][1]*(n[0][2] + n[3][2])));
                                Coord z0_ = (-n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[3][0]*n[2][1]*(n[0][2] + n[1][2]) + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] - n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*n[1][1]*(n[0][2] + n[2][2]) + n[3][0]*n[0][1]*(-n[1][2] + n[2][2]) + n[1][0]*(n[0][1] + n[2][1])*n[3][2] + n[0][0]*(-n[1][1] + n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(n[0][2] + n[1][2]) + n[1][1]*(n[0][2] - n[3][2]) - n[0][1]*(n[1][2] + n[3][2])));
#ifdef USE_GMP
                                mpz_set_int128(xx_numerator_mpz, xx_numerator); mpz_set_int128(xx_denominator_mpz, xx_denominator);
                                mpz_set_int128(xy_numerator_mpz, xy_numerator); mpz_set_int128(xy_denominator_mpz, xy_denominator);
                                mpz_set_int128(xz_numerator_mpz, xz_numerator); mpz_set_int128(xz_denominator_mpz, xz_denominator);
                                mpz_set_int128(yx_numerator_mpz, yx_numerator); mpz_set_int128(yx_denominator_mpz, yx_denominator);
                                mpz_set_int128(yy_numerator_mpz, yy_numerator); mpz_set_int128(yy_denominator_mpz, yy_denominator);
                                mpz_set_int128(yz_numerator_mpz, yz_numerator); mpz_set_int128(yz_denominator_mpz, yz_denominator);
                                mpz_set_int128(zx_numerator_mpz, zx_numerator); mpz_set_int128(zx_denominator_mpz, zx_denominator);
                                mpz_set_int128(zy_numerator_mpz, zy_numerator); mpz_set_int128(zy_denominator_mpz, zy_denominator);
                                mpz_set_int128(zz_numerator_mpz, zz_numerator); mpz_set_int128(zz_denominator_mpz, zz_denominator);
                                mpz_set_int128(x0, x0_);
                                mpz_set_int128(y0, y0_);
                                mpz_set_int128(z0, z0_);
#endif
                                // Apply the affine transformation to the polytet
                                Polytet newRotatedPolytet;
                                newRotatedPolytet.reserve(tetCount);
                                for (auto tetToRotate=newPolytet.cbegin(); tetToRotate!=newPolytet.cend(); ++tetToRotate)
                                {
                                    Tet &nt = newRotatedPolytet.emplace_back();
                                    for (int vertexNum=0; vertexNum<4; vertexNum++)
                                    {
#ifdef USE_GMP
                                        mpz_set_int128(x1, tetToRotate->t[vertexNum][0]);
                                        mpz_set_int128(y1, tetToRotate->t[vertexNum][1]);
                                        mpz_set_int128(z1, tetToRotate->t[vertexNum][2]);

                                        mpz_set(x, x0); mpz_mul(x, x, mpz_power3); mpz_div(x, x, zz_denominator_mpz);
                                        mpz_set(y, y0); mpz_mul(y, y, mpz_power3); mpz_div(y, y, zz_denominator_mpz);
                                        mpz_set(z, z0); mpz_mul(z, z, mpz_power3); mpz_div(z, z, zz_denominator_mpz);

                                        mpz_set(x2, x1); mpz_mul(x2, x2, xx_numerator_mpz); mpz_mul(x2, x2, mpz_power9_2); mpz_div(x2, x2, xx_denominator_mpz);
                                        mpz_set(x3, y1); mpz_mul(x3, x3, xy_numerator_mpz); mpz_mul(x3, x3, mpz_power9_2); mpz_div(x3, x3, xy_denominator_mpz); mpz_add(x2, x2, x3);
                                        mpz_set(x3, z1); mpz_mul(x3, x3, xz_numerator_mpz); mpz_mul(x3, x3, mpz_power9_2); mpz_div(x3, x3, xz_denominator_mpz); mpz_add(x2, x2, x3); mpz_div(x2, x2, mpz_power3); mpz_add(x, x, x2);
                                        mpz_set(y2, x1); mpz_mul(y2, y2, yx_numerator_mpz); mpz_mul(y2, y2, mpz_power9_2); mpz_div(y2, y2, yx_denominator_mpz);
                                        mpz_set(y3, y1); mpz_mul(y3, y3, yy_numerator_mpz); mpz_mul(y3, y3, mpz_power9_2); mpz_div(y3, y3, yy_denominator_mpz); mpz_add(y2, y2, y3);
                                        mpz_set(y3, z1); mpz_mul(y3, y3, yz_numerator_mpz); mpz_mul(y3, y3, mpz_power9_2); mpz_div(y3, y3, yz_denominator_mpz); mpz_add(y2, y2, y3); mpz_div(y2, y2, mpz_power3); mpz_add(y, y, y2);
                                        mpz_set(z2, x1); mpz_mul(z2, z2, zx_numerator_mpz); mpz_mul(z2, z2, mpz_power9_2); mpz_div(z2, z2, zx_denominator_mpz);
                                        mpz_set(z3, y1); mpz_mul(z3, z3, zy_numerator_mpz); mpz_mul(z3, z3, mpz_power9_2); mpz_div(z3, z3, zy_denominator_mpz); mpz_add(z2, z2, z3);
                                        mpz_set(z3, z1); mpz_mul(z3, z3, zz_numerator_mpz); mpz_mul(z3, z3, mpz_power9_2); mpz_div(z3, z3, zz_denominator_mpz); mpz_add(z2, z2, z3); mpz_div(z2, z2, mpz_power3); mpz_add(z, z, z2);

                                        mpz_get_int128(nt.t[vertexNum][0], x);
                                        mpz_get_int128(nt.t[vertexNum][1], y);
                                        mpz_get_int128(nt.t[vertexNum][2], z);

#else
                                        nt.t[vertexNum][0] = x0_*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*xx_numerator*power9_2/xx_denominator + tetToRotate->t[vertexNum][1]*xy_numerator*power9_2/xy_denominator + tetToRotate->t[vertexNum][2]*xz_numerator*power9_2/xz_denominator)/power3;
                                        nt.t[vertexNum][1] = y0_*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*yx_numerator*power9_2/yx_denominator + tetToRotate->t[vertexNum][1]*yy_numerator*power9_2/yy_denominator + tetToRotate->t[vertexNum][2]*yz_numerator*power9_2/yz_denominator)/power3;
                                        nt.t[vertexNum][2] = z0_*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*zx_numerator*power9_2/zx_denominator + tetToRotate->t[vertexNum][1]*zy_numerator*power9_2/zy_denominator + tetToRotate->t[vertexNum][2]*zz_numerator*power9_2/zz_denominator)/power3;
#endif
                                    }
                                    verifyTetrahedron(nt, power9_2);
                                    memcpy(nt.faceAttached, tetToRotate->faceAttached, sizeof(nt.faceAttached));
                                }
                                // Update the running "least" rotation
                                if (!haveRunningLeast ||
                                    runningLeastPolytet > newRotatedPolytet)
                                {
                                    haveRunningLeast = true;
                                    runningLeastPolytet = newRotatedPolytet;
                                }
                                // Switch to the next rotation of this face
                                Coord3 tmpFace = n[2];
                                n[2] = n[1];
                                n[1] = n[0];
                                n[0] = tmpFace;
                            }
                        }
                        newPolytets->insert(runningLeastPolytet);
                    }
                discardThisNewPolytet:;
                }
            }
        }
#ifdef DEBUG_PRINT
        for (auto thisPolytet=newPolytets->cbegin(); thisPolytet!=newPolytets->cend(); ++thisPolytet)
        {
            bool first = true;
            for (auto thisTet=thisPolytet->cbegin(); thisTet!=thisPolytet->cend(); ++thisTet)
            {
                printf(first ? "\n{" : ",\n");
                first = false;
                printf("{{%" PRId64 ", %" PRId64 ", %" PRId64 "}, {%" PRId64 ", %" PRId64 ", %" PRId64 "}, {%" PRId64 ", %" PRId64 ", %" PRId64 "}, {%" PRId64 ", %" PRId64 ", %" PRId64 "}}",
                    (int64_t)thisTet->t[0][0], (int64_t)thisTet->t[0][1], (int64_t)thisTet->t[0][2],
                    (int64_t)thisTet->t[1][0], (int64_t)thisTet->t[1][1], (int64_t)thisTet->t[1][2],
                    (int64_t)thisTet->t[2][0], (int64_t)thisTet->t[2][1], (int64_t)thisTet->t[2][2],
                    (int64_t)thisTet->t[3][0], (int64_t)thisTet->t[3][1], (int64_t)thisTet->t[3][2]);
            }
            printf("}\n");
        }
#endif
        delete polytets;
        polytets = newPolytets;
    }
    delete polytets;
#ifdef USE_GMP
    mpz_clears(
        x , y , z ,
        x1, y1, z1,
        x2, y2, z2,
        x3, y3, z3,
        xx_numerator_mpz, xx_denominator_mpz,
        xy_numerator_mpz, xy_denominator_mpz,
        xz_numerator_mpz, xz_denominator_mpz,
        yx_numerator_mpz, yx_denominator_mpz,
        yy_numerator_mpz, yy_denominator_mpz,
        yz_numerator_mpz, yz_denominator_mpz,
        zx_numerator_mpz, zx_denominator_mpz,
        zy_numerator_mpz, zy_denominator_mpz,
        zz_numerator_mpz, zz_denominator_mpz, NULL);
#endif
#endif // main
	return 0;
}
