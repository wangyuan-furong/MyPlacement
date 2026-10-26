#ifndef MYPLACE_H
#define MYPLACE_H

#include "placedata.h"
#include "common.h"
#include <unordered_map>

class Bin_2D;
class MyPlacer;
class DCTSolver;

enum class DensityGradientMethod
{
    ElectricField,  // 资料公式（9）、（10）：在单元中心插值电场，再计算 -q*xi
    DiscreteEnergy // 保留原离散势能导数，供对照及差分测试
};

struct GlobalPlacementOptions
{
    int maxIterations = 1000;
    float stopOverflow = 0.1f;
    float momentum = 0.8f;
    int maxBinDimension = 128; // 自实现 DCT 为 O(m^3)，限制运行成本
    bool usePreconditioner = true; // 可关闭的数值加速扩展
    bool useBBStep = false;        // 仅使用同一 λ、γ 下的梯度差
    bool useStagnationStop = false;
    DensityGradientMethod densityGradientMethod = DensityGradientMethod::ElectricField;
};

struct GlobalPlacementStats
{
    int iterations = 0;
    double initialHPWL = 0, finalHPWL = 0;
    double initialOverflow = 0, finalOverflow = 0;
    string stopReason;
};

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

    double nodeDensity;    // 标准单元 + 加权宏单元的有效面积
    double fillerDensity;  // 加权 filler 有效面积
    double terminalDensity; // 位于布局行内的固定终端并集 × 目标密度
    double DarkDensity;    // 行外暗区面积 × 目标密度
    double phi;            // 电势
    double electricFieldX = 0, electricFieldY = 0; // bin 中心的电场 xi_x、xi_y

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
    ~MyPlacer();
    MyPlacer(const MyPlacer &) = delete;
    MyPlacer &operator=(const MyPlacer &) = delete;

    vector<vector<Bin_2D *>> bins;

    double StdCellArea; // 标准单元总面积
    double MacroArea;   // 宏单元总面积

    PlaceData *db;

    vector<VECTOR_3D> wirelengthGradient; // 线长梯度（节点数）
    vector<VECTOR_3D> densityGradient;    // 密度梯度（节点 + filler）
    vector<VECTOR_3D> totalGradient;      // 原始总梯度；预条件单独用于搜索方向
    vector<VECTOR_3D> fillerGradient;     // filler 梯度
    vector<float> pinWeights;             // 每个可移动单元的线网连接度（预条件子用）

    // HPWL（半周长线长）计算
    double calcHPWL();

    // 任务4：初始布局接口
    void RandomInitialPlacement();
    void ClusterInitialPlacement();
    void QuadraticInitialPlacement();

    // 任务5：全局布局（ePlace 电静力学 + Nesterov）
    bool GlobalPlacement(const GlobalPlacementOptions &options = GlobalPlacementOptions());

    // 任务7：合法化（候选行 + 行内二次优化的标准单元合法化）
    void AbacusLegalization();

    // ---- 任务5 内部状态 ----
    vector<Module *> NodesAndFillers; // 可移动单元（标准+宏）+ filler
    vector<Module *> Fillers;         // filler 单元
    int binDimension;                 // bin 网格维度 m
    float binWidth, binHeight;        // 单个 bin 的宽高
    float targetDensity;              // 目标密度
    double lambda;                    // 惩罚系数
    float gamma;                      // WA 线长模型平滑参数
    float stepSize;                   // Nesterov 步长
    unordered_map<Module *, int> nodeIndex; // 可移动单元 -> 索引
    DCTSolver *dctSolver;             // DCT 求解器（binInit 后创建）
    GlobalPlacementOptions globalOptions;
    GlobalPlacementStats globalStats;
    vector<CRect> placeableRegions, terminalRegions, darkRegions;
    double placeableArea = 0;
    double densityEnergy = 0, wirelengthValue = 0;
    vector<double> preconditioner;

    void preparePlacementRegions();
    void fillerInit();
    void binInit();
    void gradientVectorInitialization();
    void updateDensity();
    void updateDensityGradient();
    void updateWirelengthGradient();
    void updateTotalGradient();
    void penaltyFactorInitilization();
    float calcOverflow();
    void updateSmoothingParameter(double overflow);
    void updatePenaltyFactor(double deltaHPWL);
};

#endif
