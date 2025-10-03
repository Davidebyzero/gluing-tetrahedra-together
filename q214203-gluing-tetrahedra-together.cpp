#include <stdio.h>
#include <string.h>
#include <cstdint>
#include <cinttypes>
#include <array>
#include <algorithm>
#include <vector>
#include <list>
#include <functional>

//#define DEBUG_PRINT
//#define DEBUG_PRINT2

typedef int64_t Coord;
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

Coord3 operator-(const Coord3 &a)
{
    Coord3 c;
    c[0] = -a[0];
    c[1] = -a[1];
    c[2] = -a[2];
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
bool pointInTetrahedron(const Coord3 &p, const Tetrahedron &t)
{
    Coord3 n0 = cross(t[1]-t[0], t[2]-t[0]);
    Coord3 n1 = cross(t[1]-t[0], t[3]-t[0]);
    Coord3 n2 = cross(t[2]-t[0], t[3]-t[0]);
    Coord3 n3 = cross(t[2]-t[1], t[3]-t[1]);
    Coord3 center;
    center[0] = (t[0][0] + t[1][0] + t[2][0] + t[3][0]) / 4;
    center[1] = (t[0][1] + t[1][1] + t[2][1] + t[3][1]) / 4;
    center[2] = (t[0][2] + t[1][2] + t[2][2] + t[3][2]) / 4;
    if (dot((center - t[0]), n0) > 0) n0 = -n0;
    if (dot((center - t[0]), n1) > 0) n1 = -n1;
    if (dot((center - t[0]), n2) > 0) n2 = -n2;
    if (dot((center - t[1]), n3) > 0) n3 = -n3;
    return dot(p - t[0], n0) < 0 &&
           dot(p - t[0], n1) < 0 &&
           dot(p - t[0], n2) < 0 &&
           dot(p - t[1], n3) < 0;
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

/*Coord compare(const Coord3 &a, const Coord3 &b)
{
    Coord result = a[0] - b[0]; if (result != 0) return result;
    {   } result = a[1] - b[1]; if (result != 0) return result;
    {   } result = a[2] - b[2];                  return result;
}
Coord compare(const NormalizedTetrahedron &a, const NormalizedTetrahedron &b)
{
    Coord result = compare(a[0], b[0]); if (result != 0) return result;
    {   } result = compare(a[1], b[1]); if (result != 0) return result;
    {   } result = compare(a[2], b[2]); if (result != 0) return result;
    {   } result = compare(a[3], b[3]);                  return result;
}*/

bool operator<(const Polytet &a, const Polytet&b)
// implicitly assume a.size()==b.size()
{
    for (auto ta=a.cbegin(), tb=b.cbegin(); ta!=a.cend(); ++ta,++tb)
    {
        NormalizedTetrahedron tan(ta->t), tbn(tb->t);
        auto result = tan <=> tbn;
        if (result != 0)
            return result < 0;
    }
    return false;
}
bool operator>(const Polytet &a, const Polytet&b)
{
    return !(a < b);
}
bool operator==(const Polytet &a, const Polytet&b)
// implicitly assume a.size()==b.size()
{
    for (auto ta=a.cbegin(), tb=b.cbegin(); ta!=a.cend(); ++ta,++tb)
    {
        NormalizedTetrahedron tan(ta->t), tbn(tb->t);
        auto result = tan <=> tbn;
        if (result != 0)
            return false;
    }
    return true;
}

int main(int argc, char *argv[])
{
    static const Tetrahedron start =
    {{
        {{-1,-1,-1}},
        {{-1, 1, 1}},
        {{ 1,-1, 1}},
        {{ 1, 1,-1}}
    }};
    Coord power3 = 1;
    auto *polytets = new std::list<Polytet>;
    polytets->emplace_back().emplace_back(start);
    for (int outer=1;;)
    {
        printf("%d: %d\n", outer, polytets->size());
        if (++outer > 3)
            break;
        power3 *= 3;
        auto *newPolytets = new std::list<Polytet>;
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
                        for (auto tetCheckIntersection=newPolytet.cbegin(); tetCheckIntersection!=newPolytet.cend(); ++tetCheckIntersection)
                        {
                            if (&*tetCheckIntersection == tetCopyToAttachTo)
                                continue; // skip this check for speed (it'll always be false anyway)
                            if (pointInTetrahedron(newVertex, tetCheckIntersection->t))
                            {
                                //newPolytets->resize(newPolytets->size() - 1);
                                goto discardThisNewPolytet;
                            }
                        }
                        // Add new tetrahedron
                        Tet &t = newPolytet.emplace_back();
                        for (int p=0; p<3; p++)
                        {
                            int p0 = p  ^ (p  <  2           ? 1 : 0); // swap first two vertices to preserve chirality
                            int p1 = p0 + (p0 >= 3 - faceNum ? 1 : 0); // skip the opposite vertex
#ifdef DEBUG_PRINT
                            printf("%d, ", p1);
#endif
                            for (int d=0; d<3; d++)
                                t.t[p][d] = tetToAttachTo->t[p1][d] * 3;
                        }
#ifdef DEBUG_PRINT
                        printf("%d\n", 3 - faceNum);
#endif
                        t.t[3] = newVertex;
                        t.faceAttached[0] = true;
                        t.faceAttached[1] = false;
                        t.faceAttached[2] = false;
                        t.faceAttached[3] = false;
                        //
                        std::list<Polytet> rotationsOfThisPolytet;
                        for (auto tetToRotateNormalize=newPolytet.cbegin(); tetToRotateNormalize!=newPolytet.cend(); ++tetToRotateNormalize)
                        {
                            // only rotate-normalize to tetrahedrons with 1 attached face; every polytet is guaranteed to have some
                            {
                                int attachCount = 0;
                                for (int i=0; i<4; i++)
                                    attachCount += tetToRotateNormalize->faceAttached[i];
                                if (attachCount > 1)
                                    continue;
                            }
#ifdef DEBUG_PRINT
                            printf("[%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "]\n",
                                tetToRotateNormalize->t[0][0], tetToRotateNormalize->t[0][1], tetToRotateNormalize->t[0][2],
                                tetToRotateNormalize->t[1][0], tetToRotateNormalize->t[1][1], tetToRotateNormalize->t[1][2],
                                tetToRotateNormalize->t[2][0], tetToRotateNormalize->t[2][1], tetToRotateNormalize->t[2][2],
                                tetToRotateNormalize->t[3][0], tetToRotateNormalize->t[3][1], tetToRotateNormalize->t[3][2]);
#endif
                            // vertex indices of faces with identical chirality
                            static int transformFaces[4][4] =
                            {
                                {0, 1, 2, 3},
                                {0, 2, 3, 1},
                                {0, 3, 1, 2},
                                {1, 3, 2, 0},
                            };
                            for (int faceToRotateNormalize=0; faceToRotateNormalize<4; faceToRotateNormalize++)
                            {
                                Tetrahedron n;
                                for (int i=0; i<4; i++)
                                    n[i] = tetToRotateNormalize->t[transformFaces[faceToRotateNormalize][i]];
                                for (int faceRotation=0; faceRotation<3; faceRotation++)
                                {
                                    Coord xx_numerator = 3*(  ( n[2][1] - n[3][1])*(n[0][2] - n[1][2])  - (n[0][1] - n[1][1])*(n[2][2] - n[3][2]) ); Coord xx_denominator = 3*(-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                    Coord xy_numerator = 3*((-((n[2][0] - n[3][0])*(n[0][2] - n[1][2])) + (n[0][0] - n[1][0])*(n[2][2] - n[3][2]))); Coord xy_denominator = 3*(-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                    Coord xz_numerator = 3*(  ( n[2][0] - n[3][0])*(n[0][1] - n[1][1])  - (n[0][0] - n[1][0])*(n[2][1] - n[3][1]) ); Coord xz_denominator = 3*(-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                    Coord yx_numerator = 3*( -((n[1][1] - n[3][1])*(n[0][2] - n[2][2])) + (n[0][1] - n[2][1])*(n[1][2] - n[3][2]) ); Coord yx_denominator = 3*(-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                    Coord yy_numerator = 3*(  ( n[1][0] - n[3][0])*(n[0][2] - n[2][2])  - (n[0][0] - n[2][0])*(n[1][2] - n[3][2]) ); Coord yy_denominator = 3*(-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                    Coord yz_numerator = 3*( -((n[1][0] - n[3][0])*(n[0][1] - n[2][1])) + (n[0][0] - n[2][0])*(n[1][1] - n[3][1]) ); Coord yz_denominator = 3*(-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                    Coord zx_numerator = 3*( -((n[0][1] - n[3][1])*(n[1][2] - n[2][2])) + (n[1][1] - n[2][1])*(n[0][2] - n[3][2]) ); Coord zx_denominator = 3*(-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                    Coord zy_numerator = 3*(  ( n[0][0] - n[3][0])*(n[1][2] - n[2][2])  - (n[1][0] - n[2][0])*(n[0][2] - n[3][2]) ); Coord zy_denominator = 3*(-(n[1][0]*n[2][1]*n[0][2]) + n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(-(n[1][1]*n[0][2]) + n[2][1]*n[0][2] + n[0][1]*n[1][2] - n[2][1]*n[1][2] - n[0][1]*n[2][2] + n[1][1]*n[2][2]) + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] - n[0][0]*n[2][1] + n[1][0]*n[2][1])*n[3][2] + n[2][0]*(-(n[3][1]*n[0][2]) - n[0][1]*n[1][2] + n[3][1]*n[1][2] + n[1][1]*(n[0][2] - n[3][2]) + n[0][1]*n[3][2]));
                                    Coord zz_numerator = 3*(  ( n[0][0] - n[3][0])*(n[1][1] - n[2][1])  - (n[1][0] - n[2][0])*(n[0][1] - n[3][1]) ); Coord zz_denominator = 3*( n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] - n[0][0]*n[3][1]*n[2][2] + n[1][0]*n[3][1]*n[2][2] + n[3][0]*(n[1][1]*n[0][2] - n[2][1]*n[0][2] - n[0][1]*n[1][2] + n[2][1]*n[1][2] + n[0][1]*n[2][2] - n[1][1]*n[2][2]) + n[1][0]*n[0][1]*n[3][2] - n[0][0]*n[1][1]*n[3][2] + n[0][0]*n[2][1]*n[3][2] - n[1][0]*n[2][1]*n[3][2] + n[2][0]*(n[3][1]*n[0][2] + n[0][1]*n[1][2] - n[3][1]*n[1][2] - n[0][1]*n[3][2] + n[1][1]*(-n[0][2] + n[3][2])));
                                    Coord x0 = ( n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] + n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*(n[2][1]*(n[0][2] - n[1][2]) + n[1][1]*(n[0][2] + n[2][2]) - n[0][1]*(n[1][2] + n[2][2])) + n[1][0]*(n[0][1] + n[2][1])*n[3][2] - n[0][0]*(n[1][1] + n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(-n[0][2] + n[1][2]) - n[1][1]*(n[0][2] + n[3][2]) + n[0][1]*(n[1][2] + n[3][2])));
                                    Coord y0 = ( n[1][0]*n[2][1]*n[0][2] + n[1][0]*n[3][1]*n[0][2] - n[0][0]*n[2][1]*n[1][2] - n[0][0]*n[3][1]*n[1][2] - n[3][0]*((n[1][1] + n[2][1])*n[0][2] + (-n[0][1] + n[2][1])*n[1][2]) - n[1][0]*n[0][1]*n[2][2] + n[0][0]*n[1][1]*n[2][2] + n[3][0]*(n[0][1] + n[1][1])*n[2][2] - n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + (-(n[1][0]*n[0][1]) + n[0][0]*n[1][1] + (n[0][0] + n[1][0])*n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(n[0][2] + n[1][2]) + n[0][1]*(n[1][2] - n[3][2]) - n[1][1]*(n[0][2] + n[3][2])));
                                    Coord z0 = (-n[1][0]*n[2][1]*n[0][2] - n[1][0]*n[3][1]*n[0][2] + n[0][0]*n[2][1]*n[1][2] + n[0][0]*n[3][1]*n[1][2] - n[3][0]*n[2][1]*(n[0][2] + n[1][2]) + n[1][0]*n[0][1]*n[2][2] - n[0][0]*n[1][1]*n[2][2] - n[0][0]*n[3][1]*n[2][2] - n[1][0]*n[3][1]*n[2][2] + n[3][0]*n[1][1]*(n[0][2] + n[2][2]) + n[3][0]*n[0][1]*(-n[1][2] + n[2][2]) + n[1][0]*(n[0][1] + n[2][1])*n[3][2] + n[0][0]*(-n[1][1] + n[2][1])*n[3][2] + n[2][0]*(n[3][1]*(n[0][2] + n[1][2]) + n[1][1]*(n[0][2] - n[3][2]) - n[0][1]*(n[1][2] + n[3][2])));
                                    /*
                                    Coord xx_numerator =  f[2][1]*(f[1][2] - f[0][2]) - f[1][1]*(f[0][2] + f[2][2]) + f[0][1]*(f[1][2] + f[2][2]); Coord xx_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    Coord xy_numerator =  f[2][0]*(f[0][2] - f[1][2]) + f[1][0]*(f[0][2] + f[2][2]) - f[0][0]*(f[1][2] + f[2][2]); Coord xy_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    Coord xz_numerator =  f[2][0]*(f[1][1] - f[0][1]) - f[1][0]*(f[0][1] + f[2][1]) + f[0][0]*(f[1][1] + f[2][1]); Coord xz_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    Coord yx_numerator =  f[2][1]*(f[0][2] + f[1][2]) + f[1][1]*(f[0][2] - f[2][2]) - f[0][1]*(f[1][2] + f[2][2]); Coord yx_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    Coord yy_numerator = -f[2][0]*(f[0][2] + f[1][2]) + f[1][0]*(f[2][2] - f[0][2]) + f[0][0]*(f[1][2] + f[2][2]); Coord yy_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    Coord yz_numerator =  f[2][0]*(f[0][1] + f[1][1]) + f[1][0]*(f[0][1] - f[2][1]) - f[0][0]*(f[1][1] + f[2][1]); Coord yz_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    Coord zx_numerator =  f[2][1]*(f[0][2] + f[1][2]) + f[0][1]*(f[1][2] - f[2][2]) - f[1][1]*(f[0][2] + f[2][2]); Coord zx_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    Coord zy_numerator = -f[2][0]*(f[0][2] + f[1][2]) + f[1][0]*(f[0][2] + f[2][2]) + f[0][0]*(f[2][2] - f[1][2]); Coord zy_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    Coord zz_numerator =  f[2][0]*(f[0][1] + f[1][1]) + f[0][0]*(f[1][1] - f[2][1]) - f[1][0]*(f[0][1] + f[2][1]); Coord zz_denominator = f[0][2]*(f[1][0]*f[2][1] - f[2][0]*f[1][1]) + f[1][2]*(f[2][0]*f[0][1] - f[0][0]*f[2][1]) + f[2][2]*(f[0][0]*f[1][1] - f[1][0]*f[0][1]);
                                    */
#ifdef DEBUG_PRINT
                                    printf(
                                        "[%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "]  "
                                        "(%" PRId64 ", %" PRId64 ", %" PRId64 ")/%" PRId64 "  "
                                        "%" PRId64 "/%" PRId64 ", %" PRId64 "/%" PRId64 ", %" PRId64 "/%" PRId64 ", %" PRId64 "/%" PRId64 ", %" PRId64 "/%" PRId64 ", %" PRId64 "/%" PRId64 ", %" PRId64 "/%" PRId64 ", %" PRId64 "/%" PRId64 ", %" PRId64 "/%" PRId64 "\n",
                                        n[0][0], n[0][1], n[0][2],
                                        n[1][0], n[1][1], n[1][2],
                                        n[2][0], n[2][1], n[2][2],
                                        n[3][0], n[3][1], n[3][2],
                                        x0, y0, z0, zz_denominator,
                                        xx_numerator, xx_denominator, xy_numerator, xy_denominator, xz_numerator, xz_denominator,
                                        yx_numerator, yx_denominator, yy_numerator, yy_denominator, yz_numerator, yz_denominator,
                                        zx_numerator, zx_denominator, zy_numerator, zy_denominator, zz_numerator, zz_denominator);
                                    fflush(stdout);
#endif
                                    // Add this rotation of this polytet
                                    Polytet &newRotatedPolytet = rotationsOfThisPolytet.emplace_back();
                                    for (auto tetToRotate=newPolytet.cbegin(); tetToRotate!=newPolytet.cend(); ++tetToRotate)
                                    {
#ifdef DEBUG_PRINT
                                        printf("[%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "] ->\n",
                                            tetToRotate->t[0][0], tetToRotate->t[0][1], tetToRotate->t[0][2],
                                            tetToRotate->t[1][0], tetToRotate->t[1][1], tetToRotate->t[1][2],
                                            tetToRotate->t[2][0], tetToRotate->t[2][1], tetToRotate->t[2][2],
                                            tetToRotate->t[3][0], tetToRotate->t[3][1], tetToRotate->t[3][2]);
#endif
                                        Tet &nt = newRotatedPolytet.emplace_back();
                                        for (int vertexNum=0; vertexNum<4; vertexNum++)
                                        {
                                            nt.t[vertexNum][0] = x0*power3/zz_denominator + tetToRotate->t[vertexNum][0]*xx_numerator*2*power3/xx_denominator + tetToRotate->t[vertexNum][1]*xy_numerator*2*power3/xy_denominator + tetToRotate->t[vertexNum][2]*xz_numerator*2*power3/xz_denominator;
                                            nt.t[vertexNum][1] = y0*power3/zz_denominator + tetToRotate->t[vertexNum][0]*yx_numerator*2*power3/yx_denominator + tetToRotate->t[vertexNum][1]*yy_numerator*2*power3/yy_denominator + tetToRotate->t[vertexNum][2]*yz_numerator*2*power3/yz_denominator;
                                            nt.t[vertexNum][2] = z0*power3/zz_denominator + tetToRotate->t[vertexNum][0]*zx_numerator*2*power3/zx_denominator + tetToRotate->t[vertexNum][1]*zy_numerator*2*power3/zy_denominator + tetToRotate->t[vertexNum][2]*zz_numerator*2*power3/zz_denominator;
                                        }
                                        memcpy(nt.faceAttached, tetToRotate->faceAttached, sizeof(nt.faceAttached));
#ifdef DEBUG_PRINT
                                        printf("[%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "], [%" PRId64 ", %" PRId64 ", %" PRId64 "] <--\n",
                                            nt.t[0][0], nt.t[0][1], nt.t[0][2],
                                            nt.t[1][0], nt.t[1][1], nt.t[1][2],
                                            nt.t[2][0], nt.t[2][1], nt.t[2][2],
                                            nt.t[3][0], nt.t[3][1], nt.t[3][2]);
#endif
                                    }
                                    Coord3 tmpFace = n[2];
                                    n[2] = n[1];
                                    n[1] = n[0];
                                    n[0] = tmpFace;
#ifdef DEBUG_PRINT
                                    printf(".\n");
#endif
                                }
                            }
#ifdef DEBUG_PRINT
                            printf("---\n");
#endif
                        }
                        // choose the "least" of these rotations, so it can later be compared for equality by others of its kind
                        auto polytetToCompare = rotationsOfThisPolytet.cbegin();
                        const Polytet *polytetToAdd = &*polytetToCompare;
                    //if (outer<=2)
                        while (polytetToCompare != rotationsOfThisPolytet.cend())
                        {
                            ++polytetToCompare;
                            if (*polytetToAdd >  *polytetToCompare)
                                 polytetToAdd = &*polytetToCompare;
                        }
                        newPolytets->emplace_back(*polytetToAdd);
#ifdef DEBUG_PRINT
                        printf("======\n");
#endif
                    }
                discardThisNewPolytet:;
                }
            }
        }
        
        newPolytets->sort();
        newPolytets->unique();
#ifdef DEBUG_PRINT2
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
