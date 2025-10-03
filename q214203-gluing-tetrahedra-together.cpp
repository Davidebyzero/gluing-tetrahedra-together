#include <stdio.h>
#include <string.h>
#include <cstdint>
#include <cinttypes>
#include <array>
#include <algorithm>
#include <vector>
#include <list>
#include <functional>
#include <unordered_set>

//#define DEBUG_PRINT

typedef __int128 Coord;
typedef std::array<Coord, 3> Coord3;
typedef std::array<Coord3, 4> Tetrahedron;
class Tet
{
    void initFaces()
    {
        faceAttached[0] = false; // t[0],t[1],t[2]
        faceAttached[1] = false; // t[0],t[1],t[3]
        faceAttached[2] = false; // t[0],t[2],t[3]
        faceAttached[3] = false; // t[1],t[2],t[3]
    }
public:
    Tetrahedron t;
    bool faceAttached[4];
    Tet(                    ) : t( ) {initFaces();}
    Tet(const Tetrahedron &t) : t(t) {initFaces();}
};
typedef std::vector<Tet> Polytet;

// vertex indices of faces with identical chirality
static int tetrahedronFaces[4][4] =
{
    {0, 1, 2, 3},
    {0, 3, 1, 2},
    {0, 2, 3, 1},
    {1, 3, 2, 0},
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
    static int edges[6][2] =
    {
        {0, 1},
        {1, 2},
        {2, 0},
        {0, 3},
        {1, 3},
        {2, 3},
    };
    for (int edgeNum=0; edgeNum<6; edgeNum++)
    {
        Coord3 p0 = a[edges[edgeNum][0]];
        Coord3 p1 = a[edges[edgeNum][1]];
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

class NormalizedTetrahedron : public Tetrahedron
{
public:
    NormalizedTetrahedron(const Tetrahedron &t)
    {
        *(Tetrahedron*)this = t;
        std::sort(this->begin(), this->end());
    }
};

bool operator<(const Tet &_a, const Tet &_b)
{
    NormalizedTetrahedron a(_a.t), b(_b.t);
    return a < b;
}

bool operator<(const Polytet &_a, const Polytet &_b)
// implicitly assume a.size()==b.size()
{
    Polytet a(_a), b(_b);
    std::sort(begin(a), end(a));
    std::sort(begin(b), end(b));
    for (auto ta=a.cbegin(), tb=b.cbegin(); ta!=a.cend(); ++ta,++tb)
    {
        NormalizedTetrahedron tan(ta->t), tbn(tb->t);
        auto result = tan <=> tbn;
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
    Polytet a(_a), b(_b);
    std::sort(begin(a), end(a));
    std::sort(begin(b), end(b));
    for (auto ta=a.cbegin(), tb=b.cbegin(); ta!=a.cend(); ++ta,++tb)
    {
        NormalizedTetrahedron tan(ta->t), tbn(tb->t);
        auto result = tan <=> tbn;
        if (result != 0)
            return false;
    }
    return true;
}

namespace std
{
    template<>
    struct hash<Polytet>
    {
        std::size_t operator()(const Polytet &_polytet) const noexcept
        {
            Polytet polytet(_polytet);
            std::sort(begin(polytet), end(polytet));
            std::size_t seed = polytet.size();
            for (auto tet=polytet.cbegin(); tet!=polytet.cend(); ++tet)
            {
                NormalizedTetrahedron t(tet->t);
                for (auto c=t.cbegin(); c!=t.cend(); ++c)
                    for (auto i=c->cbegin(); i!=c->cend(); ++i)
                    {
                        seed ^= std::hash<uint64_t>{}((uint64_t)(*i      )) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                        seed ^= std::hash<uint64_t>{}((uint64_t)(*i >> 64)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                    }
            }
            return seed;
        }
    };
};

int main(int argc, char *argv[])
{
    static const Tetrahedron start =
    {{
        {{-1,-1,-1}},
        {{-1, 1, 1}},
        {{ 1,-1, 1}},
        {{ 1, 1,-1}}
    }};
    Polytet startPolytet;
    startPolytet.emplace_back(start);

    Coord power3 = 1;
    auto *polytets = new std::unordered_set<Polytet>;
    polytets->insert(startPolytet);
    for (int outer=1;;)
    {
        printf("%d: %d\n", outer, polytets->size());
        ++outer;
        /*if (outer > 6)
            break;*/
        power3 *= 3;
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
                        Coord3 newVertex = {{0, 0, 0}};
                        // Get center of face by averaging its vertices' coordinates; the
                        // division by 3 is implied by omitting the multiplication by 3.
                        for (int p=0; p<4; p++)
                        {
                            if (p == 3 - faceNum)
                                continue;
                            for (int d=0; d<3; d++)
                                newVertex[d] += tetToAttachTo->t[p][d];
                        }
                        for (int d=0; d<3; d++)
                            newVertex[d] += newVertex[d] - tetToAttachTo->t[3 - faceNum][d] * 3;
                        // Add new tetrahedron
                        Tet &t = newPolytet.emplace_back();
                        for (int p=0; p<3; p++)
                        {
                            int p0 = p;
                            if (!(faceNum & 1))
                                p0 ^=      (p  <  2           ? 1 : 0); // swap first two vertices to preserve chirality
                            int p1  = p0 + (p0 >= 3 - faceNum ? 1 : 0); // skip the opposite vertex
                            for (int d=0; d<3; d++)
                                t.t[p][d] = tetToAttachTo->t[p1][d] * 3;
                        }
                        t.t[3] = newVertex;
                        t.faceAttached[0] = true;
                        t.faceAttached[1] = false;
                        t.faceAttached[2] = false;
                        t.faceAttached[3] = false;
                        // Check for overlap between this newly attached tetrahedron and the existing ones
                        for (auto tetCheckIntersection=newPolytet.cbegin(); tetCheckIntersection!=newPolytet.cend(); ++tetCheckIntersection)
                        {
                            if (&*tetCheckIntersection == tetCopyToAttachTo)
                                continue; // skip this check for speed (it'll always be false anyway)
                            if (volumesOverlap(t.t, tetCheckIntersection->t))
                                goto discardThisNewPolytet;
                        }
                        //
                        std::list<Polytet> rotationsOfThisPolytet;
                        for (auto tetToRotateNormalize=newPolytet.cbegin(); tetToRotateNormalize!=newPolytet.cend(); ++tetToRotateNormalize)
                        {
                            // only rotate-normalize to tetrahedrons with 1 attached face; every polytet is guaranteed to have some
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
                                Coord xx_numerator = power3*(  ( n[2][1] - n[3][1])*(n[0][2] - n[1][2])  - (n[0][1] - n[1][1])*(n[2][2] - n[3][2]) ); Coord xx_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord xy_numerator = power3*((-((n[2][0] - n[3][0])*(n[0][2] - n[1][2])) + (n[0][0] - n[1][0])*(n[2][2] - n[3][2]))); Coord xy_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord xz_numerator = power3*(  ( n[2][0] - n[3][0])*(n[0][1] - n[1][1])  - (n[0][0] - n[1][0])*(n[2][1] - n[3][1]) ); Coord xz_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord yx_numerator = power3*( -((n[1][1] - n[3][1])*(n[0][2] - n[2][2])) + (n[0][1] - n[2][1])*(n[1][2] - n[3][2]) ); Coord yx_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord yy_numerator = power3*(  ( n[1][0] - n[3][0])*(n[0][2] - n[2][2])  - (n[0][0] - n[2][0])*(n[1][2] - n[3][2]) ); Coord yy_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord yz_numerator = power3*( -((n[1][0] - n[3][0])*(n[0][1] - n[2][1])) + (n[0][0] - n[2][0])*(n[1][1] - n[3][1]) ); Coord yz_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord zx_numerator = power3*( -((n[0][1] - n[3][1])*(n[1][2] - n[2][2])) + (n[1][1] - n[2][1])*(n[0][2] - n[3][2]) ); Coord zx_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord zy_numerator = power3*(  ( n[0][0] - n[3][0])*(n[1][2] - n[2][2])  - (n[1][0] - n[2][0])*(n[0][2] - n[3][2]) ); Coord zy_denominator = (-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                Coord zz_numerator = power3*(  ( n[0][0] - n[3][0])*(n[1][1] - n[2][1])  - (n[1][0] - n[2][0])*(n[0][1] - n[3][1]) ); Coord zz_denominator = ( n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] - n[0][0]*n[3][1]*n[2][2] + n[1][0]*n[3][1]*n[2][2] + n[3][0]*(n[1][1]*n[0][2] - n[2][1]*n[0][2] - n[0][1]*n[1][2] + n[2][1]*n[1][2] + n[0][1]*n[2][2] - n[1][1]*n[2][2]) + n[1][0]*n[0][1]*n[3][2] - n[0][0]*n[1][1]*n[3][2] + n[0][0]*n[2][1]*n[3][2] - n[1][0]*n[2][1]*n[3][2] + n[2][0]*(n[3][1]*n[0][2] + n[0][1]*n[1][2] - n[3][1]*n[1][2] - n[0][1]*n[3][2] + n[1][1]*(-n[0][2] + n[3][2])));
                                Coord x0 = ( n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(n[2][1]*(n[0][2] - n[1][2]) + n[1][1]*(n[0][2] + n[2][2]) - n[0][1]*(n[1][2] + n[2][2])) + n[1][0]*(n[0][1] + n[2][1])*n[3][2] - n[0][0]*(n[1][1] + n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(-n[0][2] + n[1][2]) - n[1][1]*(n[0][2] + n[3][2]) + n[0][1]*(n[1][2] + n[3][2])));
                                Coord y0 = ( n[1][0]*n[2][1]*n[0][2] + n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] - n[3][0]*((n[1][1] + n[2][1])*n[0][2] + (-n[0][1] + n[2][1])*n[1][2]) - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] + n[3][0]*(n[0][1] + n[1][1])*n[2][2] - n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] + (n[0][0] + n[1][0])*n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(n[0][2] + n[1][2]) + n[0][1]*(n[1][2] - n[3][2]) - n[1][1]*(n[0][2] + n[3][2])));
                                Coord z0 = (-n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[3][0]*n[2][1]*(n[0][2] + n[1][2]) + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] - n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*n[1][1]*(n[0][2] + n[2][2]) + n[3][0]*n[0][1]*(-n[1][2] + n[2][2]) + n[1][0]*(n[0][1] + n[2][1])*n[3][2] + n[0][0]*(-n[1][1] + n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(n[0][2] + n[1][2]) + n[1][1]*(n[0][2] - n[3][2]) - n[0][1]*(n[1][2] + n[3][2])));
                                // Add this rotation of this polytet
                                Polytet &newRotatedPolytet = rotationsOfThisPolytet.emplace_back();
                                for (auto tetToRotate=newPolytet.cbegin(); tetToRotate!=newPolytet.cend(); ++tetToRotate)
                                {
                                    Tet &nt = newRotatedPolytet.emplace_back();
                                    for (int vertexNum=0; vertexNum<4; vertexNum++)
                                    {
                                        nt.t[vertexNum][0] = x0*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*xx_numerator*2*power3/xx_denominator + tetToRotate->t[vertexNum][1]*xy_numerator*2*power3/xy_denominator + tetToRotate->t[vertexNum][2]*xz_numerator*2*power3/xz_denominator)/power3;
                                        nt.t[vertexNum][1] = y0*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*yx_numerator*2*power3/yx_denominator + tetToRotate->t[vertexNum][1]*yy_numerator*2*power3/yy_denominator + tetToRotate->t[vertexNum][2]*yz_numerator*2*power3/yz_denominator)/power3;
                                        nt.t[vertexNum][2] = z0*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*zx_numerator*2*power3/zx_denominator + tetToRotate->t[vertexNum][1]*zy_numerator*2*power3/zy_denominator + tetToRotate->t[vertexNum][2]*zz_numerator*2*power3/zz_denominator)/power3;
                                    }
                                    memcpy(nt.faceAttached, tetToRotate->faceAttached, sizeof(nt.faceAttached));
                                }
                                Coord3 tmpFace = n[2];
                                n[2] = n[1];
                                n[1] = n[0];
                                n[0] = tmpFace;
                            }
                        }
                        // choose the "least" of these rotations, so it can later be compared for equality by others of its kind
                        auto polytetToCompare = rotationsOfThisPolytet.cbegin();
                        const Polytet *polytetToAdd = &*polytetToCompare;
                        while (polytetToCompare != rotationsOfThisPolytet.cend())
                        {
                            ++polytetToCompare;
                            if (*polytetToAdd >  *polytetToCompare)
                                 polytetToAdd = &*polytetToCompare;
                        }
                        newPolytets->insert(*polytetToAdd);
                    }
                discardThisNewPolytet:;
                }
            }
        }
        
#ifdef DEBUG_PRINT
        for (auto thisPolytet=newPolytets->cbegin(); thisPolytet!=newPolytets->cend(); ++thisPolytet)
        {
            for (auto thisTet=thisPolytet->cbegin(); thisTet!=thisPolytet->cend(); ++thisTet)
            {
                printf("[%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "]\n",
                    thisTet->t[0][0], thisTet->t[0][1], thisTet->t[0][2],
                    thisTet->t[1][0], thisTet->t[1][1], thisTet->t[1][2],
                    thisTet->t[2][0], thisTet->t[2][1], thisTet->t[2][2],
                    thisTet->t[3][0], thisTet->t[3][1], thisTet->t[3][2]);
            }
            printf("---\n");
        }
#endif
        
        delete polytets;
        polytets = newPolytets;
    }
    delete polytets;
	return 0;
}
