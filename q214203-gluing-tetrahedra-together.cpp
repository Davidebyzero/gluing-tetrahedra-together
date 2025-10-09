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

#define USE_GMP
//#define DEBUG_PRINT

#ifdef USE_GMP
#include <gmp.h>
#   if GMP_NUMB_BITS != 64
#   error This is hard-coded for 64-bit limbs
#   endif
#endif

auto startTime = std::chrono::steady_clock::now();

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
    {0, 2, 1, 3},
    {0, 1, 3, 2},
    {0, 3, 2, 1},
    {1, 2, 3, 0},
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
#ifndef USE_GMP
Coord dot(const Coord3 &a, const Coord3 &b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
#endif

#ifdef USE_GMP
void mpz_set_int128(mpz_t &dst, const __int128 &src)
{
    __int128 abssrc = src < 0 ? -src : src;
    mpz_import(dst, 2, -1, 8, 0, 0, &abssrc);
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

class TetrahedronOverlap
{
#ifdef USE_GMP
    mpz_t intersectNumerator, intersectDenominator, uNumerator, vNumerator, uvDenominator, uvNumeratorSum;
    mpz_t center[3], normal[3], tmp[3], p0p1[3], intersectionPoint[3], delta[3], edge1[3], edge2[3];
    mpz_t a[4][3];
    mpz_t b[4][3];
    mpz_t triangle[3][3];
    void setTetrahedron(mpz_t dst[4][3], const Tetrahedron &src)
    {
        for (int p=0; p<4; p++)
            for (int d=0; d<3; d++)
                mpz_set_int128(dst[p][d], src[p][d]);
    }
    void dot(mpz_t &result, const mpz_t a[3], const mpz_t b[3])
    {
        mpz_mul   (result, a[0], b[0]);
        mpz_addmul(result, a[1], b[1]);
        mpz_addmul(result, a[2], b[2]);
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
        for (int p=0; p<4; p++)
        {
            for (int d=0; d<3; d++)
            {
                mpz_init(a[p][d]);
                mpz_init(b[p][d]);
            }
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
        for (int p=0; p<4; p++)
        {
            for (int d=0; d<3; d++)
            {
                mpz_clear(a[p][d]);
                mpz_clear(b[p][d]);
            }
        }
        for (int p=0; p<3; p++)
            for (int d=0; d<3; d++)
                mpz_clear(triangle[p][d]);
    }
    void setA(const Tetrahedron &x) {setTetrahedron(a, x);}
    void setB(const Tetrahedron &x) {setTetrahedron(b, x);}
    bool operator()()
    {
        // Take advantage of the fact that the two tetrahedrons are regular and congruent, and
        // just check if any edge of tetrahedron "a" intersects with any face of tetrahedron "b".
        // Don't count it if only the endpoint of an edge intersects.
        for (int edgeNum=0; edgeNum<6; edgeNum++)
        {
            const mpz_t *p0 = a[tetrahedronEdges[edgeNum][0]];
            const mpz_t *p1 = a[tetrahedronEdges[edgeNum][1]];
            for (int faceNum=0; faceNum<4; faceNum++)
            {
                mpz_t *normalizedTetrahedron[4][3]; // first 3 points are the face, and the 4th point is for calculating the normal
                for (int i=0; i<4; i++)
                    for (int d=0; d<3; d++)
                        normalizedTetrahedron[i][d] = &b[tetrahedronFaces[faceNum][i]][d];
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
    const Tetrahedron *a, *b;
public:
    void setA(const Tetrahedron &x) {a = &x;}
    void setB(const Tetrahedron &x) {b = &x;}
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
                Tetrahedron normalizedTetrahedron; // first 3 points are the face, and the 4th point is for calculating the normal
                for (int i=0; i<4; i++)
                    normalizedTetrahedron[i] = (*b)[tetrahedronFaces[faceNum][i]];
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

void attachNewTet(Tet &t, const Tet &tetToAttachTo, const int faceNum)
{
    Coord3 &newVertex = t.t[3];
    newVertex = {{0, 0, 0}};
    // Get center of face by averaging its vertices' coordinates; the
    // division by 3 is implied by omitting the multiplication by 3.
    for (int p=0; p<4; p++)
    {
        if (p == 3 - faceNum)
            continue;
        for (int d=0; d<3; d++)
            newVertex[d] += tetToAttachTo.t[p][d];
    }
    // Finalize the new vertex
    for (int d=0; d<3; d++)
        newVertex[d] += newVertex[d] - tetToAttachTo.t[3 - faceNum][d] * 3;
    // Copy the other vertices
    for (int p=0; p<3; p++)
    {
        int p1 = tetrahedronFaces[faceNum][p];
        for (int d=0; d<3; d++)
            t.t[p][d] = tetToAttachTo.t[p1][d] * 3;
    }
    t.faceAttached[0] = true;
    t.faceAttached[1] = false;
    t.faceAttached[2] = false;
    t.faceAttached[3] = false;
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

#ifdef USE_GMP
void verifyTetrahedron(mpz_t t[4][3], const mpz_t mpz_power9_8, mpz_t edge[3])
{
    for (int edgeNum=0; edgeNum<6; edgeNum++)
    {
        auto p0 = t[tetrahedronEdges[edgeNum][0]];
        auto p1 = t[tetrahedronEdges[edgeNum][1]];
        for (int d=0; d<3; d++)
        {
            mpz_sub(edge[d], p1[d], p0[d]);
            mpz_mul(edge[d], edge[d], edge[d]);
        }
        for (int d=1; d<3; d++)
            mpz_add(edge[0], edge[0], edge[d]);
        if (mpz_cmp(edge[0], mpz_power9_8) != 0)
        {
            auto currentTime = std::chrono::steady_clock::now();
            std::cerr << "Quitting due to detected overflow" << " [" << std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() << " ms]" << std::endl;
            exit(-1);
        }
    }
}
#else
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
#endif

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

#ifdef USE_GMP
    mpz_t mpz_power3, mpz_power9_2, mpz_power9_8;
    mpz_init(mpz_power3  ); mpz_set_ui(mpz_power3  , 1);
    mpz_init(mpz_power9_2); mpz_set_ui(mpz_power9_2, 2);
    mpz_init(mpz_power9_8); mpz_set_ui(mpz_power9_8, 8);
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
    mpz_t nt_mpz[4][3];
    mpz_t edge[3]; // scratch for verifyTetrahedron()
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
    for (int p=0; p<4; p++)
        for (int d=0; d<3; d++)
            mpz_init(nt_mpz[p][d]);
    for (int d=0; d<3; d++)
        mpz_init(edge[d]);
#endif
    TetrahedronOverlap overlap;
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
        mpz_mul_ui(mpz_power9_8, mpz_power9_8, 9);
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
                        overlap.setA(t.t);
                        for (auto tetCheckIntersection=newPolytet.cbegin(); tetCheckIntersection!=newPolytet.cend(); ++tetCheckIntersection)
                        {
                            if (&*tetCheckIntersection == tetCopyToAttachTo || &*tetCheckIntersection == &t)
                                continue; // skip this check for speed (it'll always be false anyway)
                            overlap.setB(tetCheckIntersection->t);
                            if (overlap())
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

                                        mpz_set(nt_mpz[vertexNum][0], x);
                                        mpz_set(nt_mpz[vertexNum][1], y);
                                        mpz_set(nt_mpz[vertexNum][2], z);
#else
                                        nt.t[vertexNum][0] = x0_*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*xx_numerator*power9_2/xx_denominator + tetToRotate->t[vertexNum][1]*xy_numerator*power9_2/xy_denominator + tetToRotate->t[vertexNum][2]*xz_numerator*power9_2/xz_denominator)/power3;
                                        nt.t[vertexNum][1] = y0_*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*yx_numerator*power9_2/yx_denominator + tetToRotate->t[vertexNum][1]*yy_numerator*power9_2/yy_denominator + tetToRotate->t[vertexNum][2]*yz_numerator*power9_2/yz_denominator)/power3;
                                        nt.t[vertexNum][2] = z0_*power3/zz_denominator + (tetToRotate->t[vertexNum][0]*zx_numerator*power9_2/zx_denominator + tetToRotate->t[vertexNum][1]*zy_numerator*power9_2/zy_denominator + tetToRotate->t[vertexNum][2]*zz_numerator*power9_2/zz_denominator)/power3;
#endif
                                    }
#ifdef USE_GMP
                                    verifyTetrahedron(nt_mpz, mpz_power9_8, edge);
#else
                                    verifyTetrahedron(nt, power9_2);
#endif
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
        mpz_power3, mpz_power9_2, mpz_power9_8,
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
    for (int p=0; p<4; p++)
        for (int d=0; d<3; d++)
            mpz_clear(nt_mpz[p][d]);
    for (int d=0; d<3; d++)
        mpz_clear(edge[d]);
#endif
	return 0;
}
