#ifndef COMMON_H
#define COMMON_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <map>
#include <list>
#include <set>

using namespace std;

// 二维坐标点：可用于存储坐标、偏移量
struct POS_2D
{
    float x;
    float y;
    POS_2D() : x(0), y(0) {}
    POS_2D(float _x, float _y) : x(_x), y(_y) {}
    inline void SetZero() { x = y = 0.0f; }
    inline POS_2D operator+(const POS_2D &rhs) const { return POS_2D(x + rhs.x, y + rhs.y); }
    inline POS_2D operator-(const POS_2D &rhs) const { return POS_2D(x - rhs.x, y - rhs.y); }
    inline POS_2D operator*(float c) const { return POS_2D(x * c, y * c); }
    friend inline std::ostream &operator<<(std::ostream &os, const POS_2D &pos)
    {
        os << "(" << pos.x << "," << pos.y << ")";
        return os;
    }
};

// 矩形区域：左下角 ll、右上角 ur
struct CRect
{
    POS_2D ll;
    POS_2D ur;
    CRect() {}
    CRect(const POS_2D &_ll, const POS_2D &_ur) : ll(_ll), ur(_ur) {}
    inline float width() const { return ur.x - ll.x; }
    inline float height() const { return ur.y - ll.y; }
    inline float area() const { return width() * height(); }
    inline POS_2D center() const { return POS_2D((ll.x + ur.x) * 0.5f, (ll.y + ur.y) * 0.5f); }
};

// 三维向量：用于存储梯度
struct VECTOR_3D
{
    float x;
    float y;
    float z;
    VECTOR_3D() : x(0), y(0), z(0) {}
    inline void SetZero() { x = y = z = 0.0f; }
    inline VECTOR_3D operator+(const VECTOR_3D &rhs) const
    {
        VECTOR_3D v;
        v.x = this->x + rhs.x;
        v.y = this->y + rhs.y;
        v.z = this->z + rhs.z;
        return v;
    }
    inline VECTOR_3D operator-(const VECTOR_3D &rhs) const
    {
        VECTOR_3D v;
        v.x = this->x - rhs.x;
        v.y = this->y - rhs.y;
        v.z = this->z - rhs.z;
        return v;
    }
    inline VECTOR_3D operator*(float c) const
    {
        VECTOR_3D v;
        v.x = this->x * c;
        v.y = this->y * c;
        v.z = this->z * c;
        return v;
    }
    inline float operator*(const VECTOR_3D &rhs) const
    {
        return (this->x * rhs.x + this->y * rhs.y + this->z * rhs.z);
    }
    friend inline std::ostream &operator<<(std::ostream &os, const VECTOR_3D &vec)
    {
        os << "[" << vec.x << "," << vec.y << "," << vec.z << "]";
        return os;
    }
};

// 单元放置方向（BookShelf .pl 文件）
enum ORIENT
{
    N = 0,  // North（默认，无旋转）
    S = 1,  // South（翻转）
    W = 2,  // West（旋转）
    E = 3,  // East（旋转）
    FN = 4, // 翻转 + 北
    FS = 5, // 翻转 + 南
    FE = 6, // 翻转 + 东
    FW = 7, // 翻转 + 西
    UNKNOWN_ORIENT = 8
};

// 字符串转方向枚举
static inline ORIENT parseOrientation(const string &s)
{
    if (s == "N") return N;
    if (s == "S") return S;
    if (s == "W") return W;
    if (s == "E") return E;
    if (s == "FN") return FN;
    if (s == "FS") return FS;
    if (s == "FE") return FE;
    if (s == "FW") return FW;
    return UNKNOWN_ORIENT;
}

#endif
