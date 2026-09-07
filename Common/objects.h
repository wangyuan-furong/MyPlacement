#ifndef OBJECTS_H
#define OBJECTS_H
#include "common.h"

class Module;
class SiteRow;
class Pin;
class Net;
class PlaceData;

// 布局单元：可移动的标准单元/宏单元，或固定的终端单元
class Module
{
public:
    int idx;
    string name;
    float width;
    float height;
    float area;
    ORIENT orientation;
    bool isMacro;
    bool isFixed;
    bool isFiller;
    vector<Pin *> modulePins;
    vector<Net *> nets;
    POS_2D ll; // 左下角坐标（BookShelf .pl 文件读出的即左下角）

    Module()
    {
        Init();
    }
    Module(int _index, const string &_name, float _width = 0, float _height = 0)
    {
        Init();
        idx = _index;
        name = _name;
        width = _width;
        height = _height;
        area = width * height;
    }
    void Init()
    {
        idx = -1;
        width = 0;
        height = 0;
        area = 0;
        orientation = N;
        isMacro = false;
        isFiller = false;
        isFixed = false;
        ll.SetZero();
    }
    inline POS_2D center() const { return POS_2D(ll.x + width * 0.5f, ll.y + height * 0.5f); }
    inline POS_2D ur() const { return POS_2D(ll.x + width, ll.y + height); }
    inline POS_2D getLL_2D() const { return ll; }
    inline POS_2D getUR_2D() const { return ur(); }
};

// 布线行（BookShelf .scl 文件的 CoreRow）
class SiteRow
{
public:
    double bottom;      // 行的底部 Y 坐标
    double height;      // 行高
    double step;        // 最小 X 步进
    POS_2D start;       // 行左下角
    POS_2D end;         // 行右下角
    ORIENT orientation; // N(0) 或 S(1)

    SiteRow()
    {
        bottom = 0;
        height = 0;
        step = 0;
        start.SetZero();
        end.SetZero();
        orientation = N;
    }
};

// 引脚
class Pin
{
public:
    int idx;
    Module *module;
    Net *net;
    POS_2D offset; // 引脚相对于所属单元中心的偏移
    int direction; // 0 output, 1 input, -1 not-define

    Pin()
    {
        init();
    }
    void init()
    {
        idx = -1;
        module = NULL;
        net = NULL;
        offset.SetZero();
        direction = -1;
    }
};

// 网络（超边）
class Net
{
public:
    int idx;
    vector<Pin *> netPins;
    double weight;

    Net()
    {
        init();
    }
    void init()
    {
        idx = 0;
        weight = 1.0;
        netPins.clear();
    }
};

#endif
