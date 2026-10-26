#include "myplace.h"
#include <memory>
#include <stdexcept>

// 小算例验证算法行为，不依赖外部 benchmark 或绘图。
struct Fixture
{
    PlaceData db;
    vector<unique_ptr<Module>> modules;
    vector<unique_ptr<Pin>> pins;
    vector<unique_ptr<Net>> nets;

    Fixture() { db.coreRegion = CRect(POS_2D(0, 0), POS_2D(100, 100)); }

    Module *module(float cx, float cy, bool fixed = false, float width = 2, float height = 2)
    {
        modules.emplace_back(new Module((int)modules.size(), "test", width, height));
        Module *m = modules.back().get();
        m->ll = POS_2D(cx - width / 2, cy - height / 2);
        m->isFixed = fixed;
        (fixed ? db.Terminals : db.Nodes).push_back(m);
        return m;
    }

    void net(const vector<Module *> &members, const vector<POS_2D> &offsets = {})
    {
        nets.emplace_back(new Net());
        Net *n = nets.back().get();
        db.Nets.push_back(n);
        for (size_t i = 0; i < members.size(); ++i)
        {
            pins.emplace_back(new Pin());
            Pin *p = pins.back().get();
            p->module = members[i];
            p->net = n;
            if (!offsets.empty())
                p->offset = offsets[i];
            n->netPins.push_back(p);
        }
    }
};

static void require(bool condition, const char *message)
{
    if (!condition)
        throw runtime_error(message);
}

static void checkPositions(const PlaceData &db)
{
    for (Module *m : db.Nodes)
    {
        require(std::isfinite(m->ll.x) && std::isfinite(m->ll.y), "nonfinite coordinate");
        require(m->ll.x >= db.coreRegion.ll.x - 1e-3 && m->ll.y >= db.coreRegion.ll.y - 1e-3 &&
                m->ur().x <= db.coreRegion.ur.x + 1e-3 && m->ur().y <= db.coreRegion.ur.y + 1e-3,
                "cell outside core");
    }
}

static void testDynamicClustering()
{
    Fixture f;
    Module *a = f.module(10, 20, false, 1, 1);
    Module *b = f.module(20, 20, false, 1, 1);
    Module *c = f.module(30, 20, false, 1, 1);
    Module *d = f.module(80, 80, false, 1, 1);
    auto connect = [&](Module *u, Module *v, int count)
    {
        for (int i = 0; i < count; ++i)
            f.net({u, v});
    };
    // 200 个单元 -> 目标 2 簇。D 的星状分量先合并。
    for (int i = 4; i < 200; ++i)
        connect(d, f.module(80, 80, false, 1, 1), 100);
    connect(a, b, 10);
    connect(a, c, 6);
    connect(b, c, 6);
    connect(c, d, 9);
    MyPlacer placer(&f.db);
    placer.ClusterInitialPlacement();
    // AB 合并后到 C 的权重为 12，应先于 CD 的 9 合并，最终得到 ABC 与 D 分量。
    for (Module *m : {a, b, c})
        require(m->center().x < 30 && m->center().y < 30, "cluster weights were not updated");
    for (size_t i = 3; i < f.db.Nodes.size(); ++i)
        require(f.db.Nodes[i]->center().x > 60, "wrong cluster membership");
    require(fabs((a->center().x + b->center().x + c->center().x) / 3 - 20) < 1e-4,
            "wrong area-weighted cluster center");
    require(fabs(a->center().x - b->center().x) > 0.5, "cluster grid collapsed");
    checkPositions(f.db);
}

static void testNetNormalizationAndOffsets()
{
    Fixture f;
    // 固定单元先创建，也验证算法不依赖 Module::idx 恰好等于可移动数组下标。
    Module *left = f.module(10, 50, true);
    Module *right1 = f.module(90, 50, true);
    Module *right2 = f.module(90, 50, true);
    Module *mov = f.module(20, 20);
    f.net({left, mov}, {POS_2D(), POS_2D(-2, 0)});
    f.net({mov, right1, right2}, {POS_2D(2, 0), POS_2D(), POS_2D()});
    MyPlacer placer(&f.db);
    placer.QuadraticInitialPlacement();
    // 二引脚网权重为 1，三引脚网的两条可移动-终端连接各为 1/2。
    // 加上对称的引脚偏移，中心 (50,50) 应保持不变。
    require(fabs(mov->center().x - 50) < 1e-3 && fabs(mov->center().y - 50) < 1e-3,
            "net-degree normalization or pin offsets incorrect");
    require(left->center().x == 10 && right1->center().x == 90, "fixed terminals moved");
    checkPositions(f.db);
}

static void testFloatingComponentsAndBoundary()
{
    Fixture f;
    Module *a = f.module(10, 10);
    Module *b = f.module(20, 20);
    Module *isolated = f.module(30, 30);
    f.net({a, b}, {POS_2D(1, 0), POS_2D(-1, 0)});
    MyPlacer placer(&f.db);
    placer.QuadraticInitialPlacement();
    require(fabs(a->center().x - 50) < 1e-3 && fabs(b->center().x - 52) < 1e-3,
            "unanchored connected component has wrong solution");
    require(fabs(isolated->center().x - 50) < 1e-3, "isolated cell has no reference");
    checkPositions(f.db);

    Fixture edge;
    Module *fixed = edge.module(200, 200, true);
    Module *mov = edge.module(10, 10);
    edge.net({mov, fixed});
    MyPlacer edgePlacer(&edge.db);
    edgePlacer.QuadraticInitialPlacement();
    checkPositions(edge.db);
    require(fabs(mov->ur().x - 100) < 1e-3 && fabs(mov->ur().y - 100) < 1e-3,
            "boundary clipping incorrect");
}

int main()
{
    try
    {
        testDynamicClustering();
        testNetNormalizationAndOffsets();
        testFloatingComponentsAndBoundary();
        Fixture empty;
        MyPlacer placer(&empty.db);
        placer.ClusterInitialPlacement();
        placer.QuadraticInitialPlacement();
        cout << "Task 4 initial placement checks passed" << endl;
    }
    catch (const exception &error)
    {
        cerr << error.what() << endl;
        return 1;
    }
    return 0;
}
