#ifndef MYPLACE_H
#define MYPLACE_H

#include "placedata.h"
#include "common.h"
#include <unordered_map>

class Bin_2D;
class MyPlacer;
struct DCTSolver; // 2D 离散余弦变换求解器（定义在 myplace.cpp，自实现零依赖）

// 密度网格单元（全局布局 ePlace 用）
class Bin_2D
{
public:
    POS_2D center;
    POS_2D ll;
    POS_2D ur;
    float width;
    float height;
    float area;

    float nodeDensity;    // 标准单元 + 宏单元密度贡献
    float fillerDensity;  // 填充单元密度贡献
    float terminalDensity; // 终端密度贡献
    float DarkDensity;    // 不可放置区域（暗节点）密度贡献
    float phi;            // 电势

    Bin_2D()
    {
        width = height = area = 0;
        nodeDensity = fillerDensity = terminalDensity = DarkDensity = phi = 0;
        center.SetZero();
        ll.SetZero();
        ur.SetZero();
    }
};

// 布局器：任务4（初始布局）与任务5（全局布局）的核心
class MyPlacer
{
public:
    MyPlacer(PlaceData *_db)
    {
        db = _db;
        dctSolver = nullptr;
        binDimension = 0;
        binWidth = binHeight = 0;
        targetDensity = 0.8f;
        lambda = 1.0f;
        gamma = 1.0f;
        stepSize = 0.0f;
        StdCellArea = MacroArea = 0;
    }

    vector<vector<Bin_2D *>> bins;

    float StdCellArea; // 标准单元总面积
    float MacroArea;   // 宏单元总面积

    PlaceData *db;

    vector<VECTOR_3D> wirelengthGradient; // 线长梯度（节点数）
    vector<VECTOR_3D> densityGradient;    // 密度梯度（节点 + filler）
    vector<VECTOR_3D> totalGradient;      // 总梯度（节点 + filler）
    vector<VECTOR_3D> fillerGradient;     // filler 梯度
    vector<float> pinWeights;             // 每个可移动单元的线网连接度（预条件子用）

    // HPWL（半周长线长）计算
    double calcHPWL();

    // 任务4：初始布局接口
    void RandomInitialPlacement();
    void ClusterInitialPlacement();
    void QuadraticInitialPlacement();

    // 任务5：全局布局（ePlace 电静力学 + Nesterov）
    void GlobalPlacement();

    // 任务7：合法化（Abacus 标准单元行合法化）
    void AbacusLegalization();

    // ---- 任务5 内部状态 ----
    vector<Module *> NodesAndFillers; // 可移动单元（标准+宏）+ filler
    vector<Module *> Fillers;         // filler 单元
    int binDimension;                 // bin 网格维度 m
    float binWidth, binHeight;        // 单个 bin 的宽高
    float targetDensity;              // 目标密度
    float lambda;                     // 惩罚系数
    float gamma;                      // WA 线长模型平滑参数
    float stepSize;                   // Nesterov 步长
    unordered_map<Module *, int> nodeIndex; // 可移动单元 -> 索引
    DCTSolver *dctSolver;             // DCT 求解器（binInit 后创建）

    void fillerInit();
    void binInit();
    void gradientVectorInitialization();
    void updateDensity();
    void updateDensityGradient();
    void updateWirelengthGradient();
    void updateTotalGradient();
    void penaltyFactorInitilization();
    float calcOverflow();
};

#endif
