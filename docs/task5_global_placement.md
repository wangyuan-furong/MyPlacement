# 任务 5：按参考资料调整全局布局

依据：综合设计III参考资料/任务5参考/全局布局设计.pdf。主要入口为
MyPlace/myplace.cpp 中的 GlobalPlacement；参数和结果结构位于 MyPlace/myplace.h。

## 实现流程

1. 对布局行和固定终端做并集分解。可放置区域为行并集减去终端并集，
   暗区为 core 中行并集之外的区域，终端和暗区互不重复计数。
2. 统计标准单元、宏单元面积。标准单元的有效电荷为真实面积，
   宏单元及 filler 的有效电荷为 targetDensity 乘真实面积。
3. filler 总物理面积由下式计算：
   FillerArea = max(0, PlaceableArea - StdCellArea/targetDensity - MacroArea)。
   因此 StdCellArea + targetDensity*MacroArea + targetDensity*FillerArea
   等于目标密度乘可放置面积（原单元已超出目标容量时不再插入 filler）。
   单个 filler 面积采用中间 80% 单元的平均面积，最后一个补齐余量。
   超过 50 万个的数量预算时增大单个 filler，保留所需总面积。
4. 根据平均单元面积和目标密度估计网格，每 bin 约四个平均单元。
   自实现 DCT 为 O(m^3)，每边默认上限 128，最小 4。
5. 小于 bin 的单元以中心为参考展开密度支持框，面积补偿为
   原宽/展开宽乘原高/展开高，宏单元和 filler 再乘目标密度。
   靠边的支持框平移到 core 内保持电荷量，单元自身物理尺寸不变。
   默认电场模式在真实单元中心取场，不使用支持框平移是否饱和来强制清零梯度。
6. 用正交 DCT-II 解 Poisson 方程，采用 bin 中心采样及物理频率
   wx=πj/(m*binWidth)、wy=πk/(m*binHeight)，零频势置零。
   电势能为 D=1/2 Σ rho*phi*binArea。
   默认按资料公式（10）重建电场：x 方向谱系数为 a*wx/(wx²+wy²)，
   y 方向为 a*wy/(wx²+wy²)，分别做 sin-cos、cos-sin 逆变换。
   电场定义在 bin 中心，用双线性插值得到真实单元中心的 xi_i；
   根据公式（9）计算 densityGradient_i=(-q_i*xi_ix,-q_i*xi_iy)。
   q_i 与第 7 页的密度权重一致：标准单元取面积，宏/filler 取目标密度乘面积。
   core 边界的虚拟格点使用法向电场奇对称、切向电场偶对称的 Neumann 延拓，
   因而墙面法向场为零，近墙合法中心不存在原先整段强制零梯度的区域。
   原离散势能导数保留为 DiscreteEnergy 对照模式。
7. 线长采用资料中的 WA 平滑模型，平滑参数
   gamma = 8*binWidth*10^((tau-0.1)*20/9-1)。
   HPWL 单独计算，用于输出及更新 lambda。
8. 初始 lambda=Σ|WLgradient|/Σ|densityGradient|，去掉原来的 8e-5 基础乘数。
   无线长梯度或密度梯度近零的退化情况使用 lambda=1，避免永久关闭密度力。
9. 每轮保存实际位置 x 和上一轮位移 v，计算并投影前瞻位置 z=x-beta*v，
   在 z 重新统计密度并计算梯度，然后从 z 沿负梯度搜索新位置。
   实际更新后令 v=x-x_new。beta 默认 0.8。
10. 每次回溯固定 lambda、gamma 和前瞻点的预条件因子。
    插值电场力不是 box 沉积离散势能的精确导数，因此默认电场模式采用局部梯度变化回溯：
    alpha <= 0.95*sqrt(sum(p_i*|delta x_i|²)/sum(|delta g_i|²/p_i))。
    delta g 是同一 lambda、gamma 下候选点与前瞻点的原始梯度差；
    p_i 是冻结的前瞻点预条件因子。这是局部 Lipschitz 估计，不保证离散势能每轮下降。
    DiscreteEnergy 对照模式继续使用原来的离散势能二次上界条件。
    单次梯度位移限制为两个 bin 的较短边长；
    这不包含已有动量产生的位移，日志 maxMove 报告实际位置的完整变化。
    数值无效或回溯失败时恢复上轮位置并返回 false，主程序据此停止后续流程。
11. 根据实际位置的 HPWL 变化更新
    lambda_new=lambda_old*1.1^(1-deltaHPWL/350000)，然后更新下一轮 gamma。
    正常数值范围严格采用该公式；lambda 在 log 域 [-60,60]，
    gamma 在 [0.001,20*core最大边长] 内仅作极端值数值保护。
12. 溢出率采用资料定义：
    sum(max(0,nodeDensity+terminalDensity+DarkDensity-targetDensity*binArea))
    / (StdCellArea+targetDensity*MacroArea)。
    filler 参与电势能与密度梯度，但不参与该溢出率的分子和分母。

默认停止条件为 tau<0.1 或执行满 1000 轮；初始状态已达标时允许零轮返回。
日志明确区分目标达标、迭代上限和数值失败，不再默认使用 lambda 冻结或
HPWL 长期无改善就停止的策略。

## 可选扩展

GlobalPlacementOptions 的默认配置：

| 参数 | 默认值 | 含义 |
|---|---:|---|
| maxIterations | 1000 | 最大迭代次数 |
| stopOverflow | 0.1 | 目标溢出率 |
| momentum | 0.8 | Nesterov 动量系数 |
| maxBinDimension | 128 | 每边最大网格数 |
| usePreconditioner | true | 按连接数和有效面积缩放搜索方向 |
| useBBStep | false | 使用同一 lambda、gamma 下的梯度差估计 BB 步长 |
| useStagnationStop | false | 启用额外的 HPWL 停滞停止条件 |
| densityGradientMethod | ElectricField | 默认按资料计算 -q*xi，DiscreteEnergy 为原势能导数对照 |

预条件和回溯是数值实现扩展；不能将整段实现表述为逐字复现 PDF 的全部伪代码。
原始总梯度保持 WLgradient+lambda*densityGradient，预条件单独用于搜索方向。
默认电场模式不再将其当作 box 沉积势能的精确导数来断言能量上界。
BB 开启时只比较本轮前瞻点与候选点在同一参数下
的梯度，不再混用更新前后 lambda、gamma 的梯度差。

代码调用示例：

    GlobalPlacementOptions options;
    // 默认 ElectricField：按资料的电场公式更新。
    // 如需复现旧的离散势能导数模式：
    // options.densityGradientMethod = DensityGradientMethod::DiscreteEnergy;
    bool ok = placer.GlobalPlacement(options);

## 构建和运行

在项目根目录：

    cmake -S . -B build
    cmake --build build -j 2
    ./build/myplace adaptec1/adaptec1.aux

完整命令会先执行任务 4，然后执行任务 5、可视化及任务 7。
--initial-only 的含义保持为只执行任务 4。

启用算法测试：

    cmake -S . -B /tmp/myplacement-task5-build -DCMAKE_BUILD_TYPE=Release -DMYPLACEMENT_BUILD_TESTS=ON
    cmake --build /tmp/myplacement-task5-build -j 2
    ctest --test-dir /tmp/myplacement-task5-build --output-on-failure

## 电场模式测试

当前 CTest 检查电场解析模态及 -q*xi、单元中心与格间插值、标准/宏/filler
电荷权重、近墙非零密度力、墙面法向与切向场、均匀密度零单元力。
原离散势能有限差分测试显式使用 DiscreteEnergy，不将其误用于电场模式。
两种梯度是不同的离散实现，不能宣称在有限网格上严格相等。

## 当前电场模式验证（2026-10-26）

Release 构建与任务 4、5 的 CTest 均通过。实际运行项目完整命令：

    ./myplace /home/lnx/projects/综合设计3/MyPlacement/adaptec1/adaptec1.aux

日志明确显示“密度梯度=电场 -q*xi，步长回溯=局部梯度变化”，
没有启用 DiscreteEnergy 对照模式。

| 指标 | 默认电场模式 |
|---|---:|
| 二次解析起点 HPWL | 42,982,845 |
| 全局布局轮数 | 233 |
| 最终密度溢出率 | 0.0966 |
| 全局 HPWL | 82,000,698 |
| 全局布局耗时 | 43.998 s |
| 合法化后 HPWL | 84,269,479 |
| 合法化 L1 位移 | 3,809,576 |
| 同行重叠/行外/site未对齐/容量回退 | 均为 0 |

停止原因是“密度溢出率达到目标”。完整流程正常退出，
生成三种初始布局图片，以及全局布局、合法化的 BMP 和 GDS。
全流程约 182.66 秒，包含约 127.52 秒的二次解析初始布局，不应与仅全局布局耗时混淆。
下方旧表记录的是更换电场梯度之前的实验，不代表当前默认配置。

## 历史验证：修改为电场模式之前

Release 构建通过，任务 4 和任务 5 的 CTest 均通过。任务 5 检查覆盖：

- 均匀密度产生零电场。
- 已知余弦模态的电势、电场幅度和方向，含非正方形 bin。
- 局部电荷的排斥方向及左右、上下对称性。
- 重叠布局行、重叠终端、终端与暗区的重复覆盖处理。
- 宏单元、filler 和小单元展开的电荷守恒；重复初始化不累积 filler。
- 密度势能和 WA 线长的数值差分梯度。
- lambda 初始化、1.1 反馈公式、gamma 公式及零密度梯度退化处理。
- 独立重建第二轮前瞻位置并检查新坐标，确认求梯度的位置是前瞻位置。
- 最大轮数、初始停止、无预条件的 BB 模式、重复运行及非法动量参数。
- 原 build/myplace 的小型 BookShelf 完整命令入口：
  全局布局 58 轮达到 tau=0.0976，生成初始布局、全局布局和合法化图片及 GDS，
  现有合法化检查四项均为零。

### adaptec1 同起点对照

先用任务 4 的二次解析布局产生并保存浮点坐标，再由新旧全局布局读取同一快照。
初始 HPWL 均为 42,982,845。没有重新随机生成起点。

| 指标 | 修改前 | 修改后 |
|---|---:|---:|
| 全局布局轮数 | 1000 | 241 |
| 全局 HPWL | 80,735,081 | 82,478,516 |
| 全局布局耗时 | 135.448 s | 48.074 s |
| 合法化后 HPWL | 82,244,536 | 84,837,523 |
| 合法化总 L1 位移 | 2,780,008 | 3,811,993 |
| 坐标无效/越核心边界 | 0 | 0 |
| 固定终端移动数 | 0 | 0 |
| 同行重叠/行外/site未对齐/容量回退 | 均为 0 | 均为 0 |

修改后默认使用 128x128 网格、108817 个 filler，可放置面积为 65026488，
初始溢出率约 0.9928、初始 lambda 约 2.3286e-9。
241 轮后以“密度溢出率达到目标”结束，最终 tau=0.0953。

旧版执行满 1000 轮，打印的旧定义溢出率约为 0.1256。
由于新版不计 filler 且修改了密度模型，不能将这两个 tau 直接按同一口径比较。
耗时包含梯度计算及诊断输出，不包含初始布局或合法化；两版网格和模型也不同，
因此不能将全部加速归因于 Nesterov。

本组数据上新版合法化后 HPWL 增加约 3.15%，L1 合法化位移也增加。
此次结果表明参考流程已落实且程序能够完成联调，不代表所有质量指标均改善。
网格密度是离散统计，不能据此宣称全局布局已完全无重叠。

任务 7 的实现没有在本次修改中扩展：其现有检查主要覆盖行内位置约束，
不等同于完整验证所有固定障碍物、跨行或朝向约束。
