#include <stdio.h>
#include <iostream>
#include <string.h>
#include <array>
#include <vector>
#include <unordered_set>
#include <chrono>

//#define USE_GMP

#ifdef USE_GMP
#include <gmp.h>
#   if GMP_NUMB_BITS != 64
#   error This is hard-coded for 64-bit limbs
#   endif
#endif

auto startTime = std::chrono::steady_clock::now();

typedef uint8_t TetIndex;
typedef uint8_t TetIndexFace; // lowest 2 bits are used for a face index

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
    {0, 1, 2, 3},
    {0, 3, 1, 2},
    {0, 2, 3, 1},
    {1, 3, 2, 0},
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
    bool operator()()
    {
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
    bool operator()()
    {
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
    t.faceAttached[0] = NULL;
    t.faceAttached[1] = NULL;
    t.faceAttached[2] = NULL;
    t.faceAttached[3] = &tetToAttachTo;
    tetToAttachTo.faceAttached[faceNum] = &t;
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
            /*if (faceNum == 3)
            {
                int tmp = vertexMap2[3];
                vertexMap2[3] = vertexMap2[2];
                vertexMap2[2] = tmp;
            }*/
            append(polytet, *attachedTet, vertexMap2, rotation);
        }
    }
};

namespace std
{
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

int main(int argc, char *argv[])
{
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
#else
    static Tetrahedron start =
    {{
        {{-9,-9,-9}},
        {{-9, 9, 9}},
        {{ 9,-9, 9}},
        {{ 9, 9,-9}}
    }};
#endif

    TetrahedronOverlap overlap;

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
            for (int tetNumToAttachTo = 0; tetNumToAttachTo < tetCount-1; tetNumToAttachTo++)
            {
                Tet &tetToAttachTo = polytet[tetNumToAttachTo];
                for (int faceNum=0; faceNum<3; faceNum++) // skip last face because it's always already attached
                {
                    if (tetToAttachTo.faceAttached[faceNum])
                        continue;
                    attachNewTet(newTet, tetToAttachTo, faceNum);
                    // Canonicalize the rotation of this new polytet in compressed form, so that it can be compared against others
                    bool haveRunningLeast = false;
                    CompressedPolytet runningLeastPolytet;

                    Tet *t = &polytet[1];
                    for (int i=0; i<tetCount; i++)
                    {
                        int attachedFace;
                        int vertexMap[4];
                        {
                            Tet &singlyAttachedTet = polytet[i];
                            for (int j=0; j<3; j++)
                                if (singlyAttachedTet.faceAttached[j])
                                    goto skipThisTet; // not a singly attached tet
                            t = singlyAttachedTet.faceAttached[3];
                            attachedFace = 0;
                            while (t->faceAttached[attachedFace] != &singlyAttachedTet)
                                attachedFace++;
                            static int vertexMapTable[4][4] =
                            {
                                {3, 0, 2, 1},
                                {2, 0, 1, 3},
                                {1, 0, 3, 2},
                                {0, 1, 2, 3},
                            };
                            memcpy(vertexMap, vertexMapTable[attachedFace], sizeof(vertexMap));
                        }
                        for (int rotationStep=0; rotationStep<3; rotationStep++)
                        {
                            polytet.resetIndexing(i);
                            CompressedPolytet newRotatedPolytet;
                            newRotatedPolytet.reserve(tetCount - 2);
                            newRotatedPolytet.append(polytet, *t, vertexMap, rotationStep);

                            // Update the running "least" rotation
                            if (!haveRunningLeast ||
                                std::lexicographical_compare(
                                    newRotatedPolytet  .begin(), newRotatedPolytet  .end(),
                                    runningLeastPolytet.begin(), runningLeastPolytet.end()))
                            {
                                haveRunningLeast = true;
                                runningLeastPolytet = newRotatedPolytet;
                            }
                        }
                    skipThisTet:;
                    }
                    if (auto [insertedItem, wasInserted] = newPolytets->emplace(runningLeastPolytet); wasInserted)
                    {
                        // Check for overlap between this newly attached tetrahedron and the existing ones,
                        // and defer this until after the deduplication, to save a lot of time
                        overlap.setA(newTet.t);
                        for (auto tetCheckIntersection=polytet.cbegin(); tetCheckIntersection!=polytet.cend(); ++tetCheckIntersection)
                        {
                            if (&*tetCheckIntersection == &tetToAttachTo || &*tetCheckIntersection == &newTet)
                                continue; // skip this check for speed (it'll always be false anyway)
                            overlap.setB(*tetCheckIntersection);
                            if (overlap())
                            {
                                newPolytets->erase(insertedItem);
                                break;
                            }
                        }
                    }
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
#ifdef USE_GMP
                mpz_mul_ui(start.t[p][d], start.t[p][d], 3);
#else
                start[p][d] *= 3;
#endif
    }
	return 0;
}
