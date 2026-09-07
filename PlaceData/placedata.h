#ifndef PLACEDATA_H
#define PLACEDATA_H
#include "objects.h"

// 布局数据库：负责解析 BookShelf 格式，并汇总统计信息
class PlaceData
{
public:
    int moduleCount;   // 单元总数（含终端）
    int MacroCount;    // 宏单元数
    int netCount;      // 网络数
    int pinCount;      // 引脚数
    int terminalCount; // 终端数

    float siteHeight;  // 标准单元行高（来自 .scl）
    float siteWidth;   // site 宽度（来自 .scl）

    vector<Module *> Nodes;      // 可移动单元（标准单元 + 宏）
    vector<Module *> Terminals;  // 固定终端
    vector<Pin *> Pins;
    vector<Net *> Nets;
    vector<SiteRow> SiteRows;

    map<string, Module *> moduleMap; // 单元名 -> 单元指针

    CRect chipRegion;  // 芯片区域（含 core 与终端的最小外包矩形）
    CRect coreRegion;  // 核心布局区域（SiteRow 的外包矩形）

    // BookShelf 各文件路径（aux 目录 + 文件名）
    string basePath;
    string nodesFile, netsFile, wtsFile, plFile, sclFile;

    PlaceData()
    {
        moduleCount = 0;
        MacroCount = 0;
        netCount = 0;
        pinCount = 0;
        terminalCount = 0;
        siteHeight = 0;
        siteWidth = 0;
    }

    // 解析 BookShelf 全套文件（传入 .aux 路径）
    bool ReadBookShelf(const string &auxPath);
    bool ReadAuxFile(const string &auxPath);
    bool ReadSclFile(const string &path);
    bool ReadNodesFile(const string &path);
    bool ReadPlFile(const string &path);
    bool ReadNetsFile(const string &path);
    bool ReadWtsFile(const string &path);

    // 汇总与统计
    void InitChipRegion();
    void PrintSummary();

    // 若 .pl 未指定 terminal 位置（常写成 0 0），将其均匀铺到 core 边界
    void PlaceTerminalsOnBoundary();

    // 辅助：按名称查找单元
    Module *getModuleByName(const string &name);
};

#endif
