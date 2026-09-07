#include "myplace.h"
#include "placedata.h"
#include "plot.h"
#include <iostream>
#include <iomanip>
#include <chrono>
using namespace std;
using namespace std::chrono;

static double elapsedMs(steady_clock::time_point t0)    //计时器，毫秒级，转化为纯数字后返回给count
{
    return duration_cast<duration<double, milli>>(steady_clock::now() - t0).count();
}

int main(int argc, char *argv[])
{
    if (argc < 2)   //检测有无文件
    {
        cerr << "用法: " << argv[0] << " <adaptec1.aux 路径>" << endl;
        cerr << "示例: ./myplace ../综合设计III参考资料/综合设计III参考资料/任务2参考/BookShelf格式的解析/adaptec1/adaptec1.aux" << endl;
        return 1;
    }

    //创建数据库+Myplace操作库实例并读入bookshelf
    PlaceData *db = new PlaceData();
    if (!db->ReadBookShelf(argv[1]))
    {
        cerr << "BookShelf 解析失败" << endl;
        return 1;
    }
    db->PrintSummary();

    MyPlacer *placer = new MyPlacer(db);

    // ===== 任务4：三种初始布局方法对比 =====
    cout << endl
         << "===== 初始布局 HPWL 对比 =====" << endl;

    // 方法1：随机布局
    auto t0 = steady_clock::now();
    placer->RandomInitialPlacement();
    double hpwl1 = placer->calcHPWL();
    double t1 = elapsedMs(t0);
    PLOTTING::plotPlacement("init_random", db);

    // 方法2：聚类驱动布局
    t0 = steady_clock::now();
    placer->ClusterInitialPlacement();
    double hpwl2 = placer->calcHPWL();
    double t2 = elapsedMs(t0);
    PLOTTING::plotPlacement("init_cluster", db);

    // 方法3：二次方程解析布局
    t0 = steady_clock::now();
    placer->QuadraticInitialPlacement();
    double hpwl3 = placer->calcHPWL();
    double t3 = elapsedMs(t0);
    PLOTTING::plotPlacement("init_quadratic", db);

    cout << "  随机布局      HPWL = " << fixed << setprecision(0) << hpwl1 << "   耗时 " << t1 << " ms" << endl;
    cout << "  聚类驱动布局  HPWL = " << hpwl2 << "   耗时 " << t2 << " ms" << endl;
    cout << "  二次解析布局  HPWL = " << hpwl3 << "   耗时 " << t3 << " ms" << endl;

    // ===== 任务5：全局布局（ePlace，从二次解析初值出发） =====
    cout << endl
         << "===== 全局布局 (ePlace) =====" << endl;
    t0 = steady_clock::now();
    placer->GlobalPlacement();
    double hpwl4 = placer->calcHPWL();
    double t4 = elapsedMs(t0);
    cout << "  全局布局      HPWL = " << fixed << setprecision(0) << hpwl4 << "   耗时 " << t4 << " ms" << endl;

    // 任务6：布局可视化（BMP 位图，宏=绿 / 标准单元=红 / 终端=蓝）
    PLOTTING::plotPlacement("gp_result", db);
    // 输出 GDSII 供 KLayout 打开
    PLOTTING::writeGDSII("gp_result", db);

    // ===== 任务7：合法化（Abacus） =====
    cout << endl;
    t0 = steady_clock::now();
    placer->AbacusLegalization();
    double t7 = elapsedMs(t0);
    cout << "  合法化        HPWL = " << fixed << setprecision(0) << placer->calcHPWL()
         << "   耗时 " << t7 << " ms" << endl;
    PLOTTING::plotPlacement("legal_result", db);
    PLOTTING::writeGDSII("legal_result", db);

    return 0;
}
