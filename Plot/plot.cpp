#include "plot.h"
#include <sys/stat.h>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace PLOTTING;

// 将布局坐标映射为图像像素（Y 轴翻转，使图像原点在左上角、芯片左下角在图像下方）
static void toPixel(const PlaceData *db, const POS_2D &p, float unitX, float unitY,
                    int xMargin, int yMargin, int &x, int &y)
{
    x = (int)((p.x - db->chipRegion.ll.x) * unitX + xMargin);
    y = (int)((db->chipRegion.ur.y - p.y) * unitY + yMargin); // 翻转 Y
}

void PLOTTING::plotPlacement(string imageName, PlaceData *db)
{
    string plotPath = "./output/";
    mkdir(plotPath.c_str(), 0777); // 确保输出目录存在（已存在时忽略错误）

    float chipW = db->chipRegion.ur.x - db->chipRegion.ll.x;
    float chipH = db->chipRegion.ur.y - db->chipRegion.ll.y;
    if (chipW <= 0)
        chipW = 1;
    if (chipH <= 0)
        chipH = 1;

    int imageWidth = 1200;
    int imageHeight = (int)(imageWidth * chipH / chipW);
    if (imageHeight < 1)
        imageHeight = 1;
    int xMargin = 30, yMargin = 30;

    float unitX = imageWidth / chipW;
    float unitY = imageHeight / chipH;

    // 白色背景 RGB 图
    CImg<unsigned char> img(imageWidth + 2 * xMargin, imageHeight + 2 * yMargin, 1, 3, 255);

    // 1. 核心布局区域轮廓（灰色边框，仅边框不填充）
    {
        int x1, y1, x2, y2;
        toPixel(db, db->coreRegion.ll, unitX, unitY, xMargin, yMargin, x1, y1);
        toPixel(db, db->coreRegion.ur, unitX, unitY, xMargin, yMargin, x2, y2);
        img.draw_rectangle(x1, y1, x2, y2, Gray, 1.0f, ~0U); // 带 pattern 的重载画边框
    }

    // 2. 可移动单元（宏单元绿 / 标准单元红）
    for (Module *m : db->Nodes)
    {
        int x1, y1, x2, y2;
        toPixel(db, m->getLL_2D(), unitX, unitY, xMargin, yMargin, x1, y1);
        toPixel(db, m->getUR_2D(), unitX, unitY, xMargin, yMargin, x2, y2);
        if (x2 <= x1)
            x2 = x1 + 1; // 保证至少 1 像素可见
        if (y2 <= y1)
            y2 = y1 + 1;
        img.draw_rectangle(x1, y1, x2, y2, m->isMacro ? Green : Red, 0.7f);
    }

    // 3. 固定终端（蓝色）
    for (Module *t : db->Terminals)
    {
        int x1, y1, x2, y2;
        toPixel(db, t->getLL_2D(), unitX, unitY, xMargin, yMargin, x1, y1);
        toPixel(db, t->getUR_2D(), unitX, unitY, xMargin, yMargin, x2, y2);
        if (x2 <= x1)
            x2 = x1 + 1;
        if (y2 <= y1)
            y2 = y1 + 1;
        img.draw_rectangle(x1, y1, x2, y2, Blue, 1.0f);
    }

    img.draw_text(40, 20, imageName.c_str(), Black, NULL, 1, 30);
    img.save_bmp((plotPath + imageName + ".bmp").c_str());
    cout << "bitmap file has been saved: " << plotPath << imageName << ".bmp" << endl;
}

// ======================= GDSII 输出 =======================
// 供 KLayout 打开。坐标放大 1000 倍写入，保留 3 位小数精度。
// 层号约定：1=标准单元  2=宏单元  3=终端

// GDSII 记录类型
enum
{
    GDS_HEADER = 0x00,
    GDS_BGNLIB = 0x01,
    GDS_LIBNAME = 0x02,
    GDS_UNITS = 0x03,
    GDS_ENDLIB = 0x04,
    GDS_BGNSTR = 0x05,
    GDS_STRNAME = 0x06,
    GDS_ENDSTR = 0x07,
    GDS_BOUNDARY = 0x08,
    GDS_LAYER = 0x0D,
    GDS_DATATYPE = 0x0E,
    GDS_XY = 0x10,
    GDS_ENDEL = 0x11
};

// GDSII 8 字节实数：1 字节符号+指数(余 64)，7 字节尾数(大端)
// value = mantissa * 16^(exp-64-14)，尾数归一化到 [16^13, 16^14)
static void gdsReal8(double v, unsigned char b[8])
{
    if (v == 0.0)
    {
        memset(b, 0, 8);
        return;
    }
    bool neg = v < 0;
    if (neg)
        v = -v;
    int e = 64;
    while (v >= 1.0)
    {
        v /= 16.0;
        ++e;
    }
    while (v < 1.0 / 16.0)
    {
        v *= 16.0;
        --e;
    }
    unsigned long long m = (unsigned long long)llround(v * 72057594037927936.0); // 16^14
    if (m >= (1ULL << 56))
    {
        m >>= 4;
        ++e;
    }
    b[0] = (unsigned char)((neg ? 0x80 : 0) | (e & 0x7F));
    for (int i = 0; i < 7; ++i)
        b[1 + i] = (unsigned char)((m >> (8 * (6 - i))) & 0xFF);
}

namespace
{
struct GDSFile
{
    FILE *f;
    GDSFile(const string &path) { f = fopen(path.c_str(), "wb"); }
    ~GDSFile()
    {
        if (f)
            fclose(f);
    }

    void rec(unsigned char type, unsigned char dtype, const vector<unsigned char> &d)
    {
        unsigned int len = 4 + (unsigned int)d.size();
        unsigned char h[2] = {(unsigned char)(len >> 8), (unsigned char)(len & 0xFF)};
        fwrite(h, 1, 2, f);
        unsigned char td[2] = {type, dtype};
        fwrite(td, 1, 2, f);
        if (!d.empty())
            fwrite(d.data(), 1, d.size(), f);
    }
    void recNoData(unsigned char type) { rec(type, 0, {}); }
    // BGNLIB/BGNSTR 的时间戳：12 个 INT2（24 字节）全零
    void recTimestamp(unsigned char type)
    {
        vector<unsigned char> d(24, 0);
        rec(type, 2, d);
    }
    void recInt2(unsigned char type, short v)
    {
        unsigned char d[2] = {(unsigned char)((v >> 8) & 0xFF), (unsigned char)(v & 0xFF)};
        rec(type, 2, vector<unsigned char>(d, d + 2));
    }
    void recStr(unsigned char type, const string &s)
    {
        vector<unsigned char> d(s.begin(), s.end());
        if (d.size() % 2)
            d.push_back(0); // 字符串长度补成偶数
        rec(type, 6, d);
    }
    void recReal8(unsigned char type, const vector<double> &vals)
    {
        vector<unsigned char> d;
        d.reserve(vals.size() * 8);
        for (double v : vals)
        {
            unsigned char b[8];
            gdsReal8(v, b);
            d.insert(d.end(), b, b + 8);
        }
        rec(type, 5, d);
    }
    void boundary(int layer, int x0, int y0, int x1, int y1)
    {
        recNoData(GDS_BOUNDARY);
        recInt2(GDS_LAYER, (short)layer);
        recInt2(GDS_DATATYPE, 0);
        // 闭合多边形（5 个点）
        int pts[10] = {x0, y0, x1, y0, x1, y1, x0, y1, x0, y0};
        vector<unsigned char> d;
        d.reserve(40);
        for (int i = 0; i < 10; ++i)
        {
            int v = pts[i];
            d.push_back((unsigned char)((v >> 24) & 0xFF));
            d.push_back((unsigned char)((v >> 16) & 0xFF));
            d.push_back((unsigned char)((v >> 8) & 0xFF));
            d.push_back((unsigned char)(v & 0xFF));
        }
        rec(GDS_XY, 3, d);
        recNoData(GDS_ENDEL);
    }
};
} // namespace

void PLOTTING::writeGDSII(string fileName, PlaceData *db)
{
    string plotPath = "./output/";
    mkdir(plotPath.c_str(), 0777);
    string path = plotPath + fileName + ".gds";
    GDSFile g(path);
    if (!g.f)
    {
        cerr << "[GDSII] 无法创建文件 " << path << endl;
        return;
    }

    const double SCALE = 1000.0; // 放大 1000 倍，坐标写为 INT4

    g.recInt2(GDS_HEADER, 600);
    g.recTimestamp(GDS_BGNLIB);
    g.recStr(GDS_LIBNAME, "MyPlacement");
    // UNITS：1 用户单位 = 1000 DBU，1 米 = 1e9 DBU
    g.recReal8(GDS_UNITS, {1000.0, 1e9});

    g.recTimestamp(GDS_BGNSTR);
    g.recStr(GDS_STRNAME, "TOP");

    // 标准单元 → 层 1
    for (Module *m : db->Nodes)
    {
        int x0 = (int)llround(m->ll.x * SCALE);
        int y0 = (int)llround(m->ll.y * SCALE);
        int x1 = (int)llround(m->ur().x * SCALE);
        int y1 = (int)llround(m->ur().y * SCALE);
        g.boundary(m->isMacro ? 2 : 1, x0, y0, x1, y1);
    }
    // 终端 → 层 3
    for (Module *t : db->Terminals)
    {
        int x0 = (int)llround(t->ll.x * SCALE);
        int y0 = (int)llround(t->ll.y * SCALE);
        int x1 = (int)llround(t->ur().x * SCALE);
        int y1 = (int)llround(t->ur().y * SCALE);
        g.boundary(3, x0, y0, x1, y1);
    }

    g.recNoData(GDS_ENDSTR);
    g.recNoData(GDS_ENDLIB);
    cout << "GDSII file has been saved: " << path << endl;
}
