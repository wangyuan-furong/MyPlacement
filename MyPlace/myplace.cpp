#include "myplace.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unordered_map>
#include <queue>
#include <cstdio>

// ===================== 辅助函数 =====================

// pin 绝对位置 = 所属单元中心 + 相对中心偏移
static inline POS_2D pinAbsPos(Pin *pin)
{
    Module *m = pin->module;
    POS_2D c = m->center();
    return POS_2D(c.x + pin->offset.x, c.y + pin->offset.y);
}

// 用两个 int 组成一个 long long 键（i<j）
static inline long long keyOf(int i, int j, int n)
{
    if (i > j)
        swap(i, j);
    return (long long)i * n + j;
}

// ===================== HPWL 计算 =====================

double MyPlacer::calcHPWL()
{
    double hpwl = 0;
    for (size_t e = 0; e < db->Nets.size(); e++)
    {
        Net *net = db->Nets[e];
        int d = (int)net->netPins.size();
        if (d < 2)
            continue;
        float minx = 1e30f, miny = 1e30f, maxx = -1e30f, maxy = -1e30f;
        for (int i = 0; i < d; i++)
        {
            POS_2D p = pinAbsPos(net->netPins[i]);
            minx = min(minx, p.x);
            maxx = max(maxx, p.x);
            miny = min(miny, p.y);
            maxy = max(maxy, p.y);
        }
        hpwl += (double)((maxx - minx) + (maxy - miny));
    }
    return hpwl;
}

// ===================== 方法1：随机布局 =====================

void MyPlacer::RandomInitialPlacement()
{
    srand(12345);
    CRect core = db->coreRegion;
    float cw = core.width();
    float ch = core.height();
    for (size_t i = 0; i < db->Nodes.size(); i++)
    {   //坐标以x举例：左边界+随机权重*（核心区域宽度-单元宽度）
        Module *m = db->Nodes[i];
        float x = core.ll.x + (rand() / (float)RAND_MAX) * (cw - m->width);
        float y = core.ll.y + (rand() / (float)RAND_MAX) * (ch - m->height);
        m->ll = POS_2D(x, y);
    }
}

// ===================== 方法2：聚类驱动布局 =====================

// 并查集（簇）
struct DSU
{
    vector<int> parent, sz;
    DSU(int n)
    {
        parent.resize(n);
        sz.assign(n, 1);
        for (int i = 0; i < n; i++)
            parent[i] = i;
    }
    int find(int x)
    {
        while (parent[x] != x)
        {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void unite(int a, int b)
    {
        a = find(a);
        b = find(b);
        if (a == b)
            return;
        if (sz[a] < sz[b])
            swap(a, b);
        parent[b] = a;
        sz[a] += sz[b];
    }
};

void MyPlacer::ClusterInitialPlacement()
{
    int n = (int)db->Nodes.size();
    if (n == 0)
        return;

    // 1. 遍历网表，累计单元对之间的连接强度（去重）
    unordered_map<long long, double> conn;
    conn.reserve((size_t)8000000);  //分配内存
    for (size_t e = 0; e < db->Nets.size(); e++)
    {
        Net *net = db->Nets[e];
        int d = (int)net->netPins.size();
        if (d < 2)
            continue;
        double w = 1.0 / (d - 1); // 度数归一化，避免大度数网拉近所有单元
        for (int a = 0; a < d; a++)
        {
            Module *ma = net->netPins[a]->module;
            if (!ma || ma->isFixed)
                continue;
            for (int b = a + 1; b < d; b++)
            {
                Module *mb = net->netPins[b]->module;
                if (!mb || mb->isFixed || ma == mb)
                    continue;
                conn[keyOf(ma->idx, mb->idx, n)] += w;  // 累计连接强度
            }
        }
    }

    // 2. 按连接强度降序，用并查集贪心合并到目标簇数 K
    int K = max(1, n / 100); // 目标簇数（约 1% 规模）
    vector<pair<double, pair<int, int>>> edges;
    edges.reserve(conn.size());
    for (auto &kv : conn)
    {
        long long key = kv.first;
        int u = (int)(key / n), v = (int)(key % n);
        edges.push_back(make_pair(kv.second, make_pair(u, v)));
    }
    sort(edges.begin(), edges.end(),    // 按连接强度降序
         [](const pair<double, pair<int, int>> &x, const pair<double, pair<int, int>> &y) { return x.first > y.first; });

    DSU dsu(n);
    int clusters = n;
    for (auto &e : edges)   //贪心合并
    {
        int u = e.second.first, v = e.second.second;
        if (dsu.find(u) != dsu.find(v))
        {
            dsu.unite(u, v);
            clusters--;
            if (clusters <= K)
                break;
        }
    }

    // 3. 收集每个簇的成员与总面积
    unordered_map<int, vector<Module *>> clusterMembers;
    unordered_map<int, double> clusterArea;
    clusterMembers.reserve(K * 2);  //为了性能分配两倍大小
    clusterArea.reserve(K * 2);
    for (int i = 0; i < n; i++)
    {
        int r = dsu.find(i);
        Module *m = db->Nodes[i];
        clusterMembers[r].push_back(m);
        clusterArea[r] += (double)m->width * m->height;
    }

    // 4. 按面积把簇质心放到 core 区域（网格填充），簇内单元围绕质心放置
    CRect core = db->coreRegion;
    srand(12345);
    vector<pair<double, int>> ordered;
    for (auto &kv : clusterArea)
        ordered.push_back(make_pair(kv.second, kv.first));
    sort(ordered.begin(), ordered.end(),    //按面积降序排序
         [](const pair<double, int> &x, const pair<double, int> &y) { return x.first > y.first; });

    int g = (int)ceil(sqrt((double)ordered.size()));    //让形状接近于正方形
    for (size_t ci = 0; ci < ordered.size(); ci++)
    {
        int r = ordered[ci].second;         //先给簇的位置定下
        float gx = core.ll.x + ((ci % g) + 0.5f) * (core.width() / g);
        float gy = core.ll.y + ((ci / g) + 0.5f) * (core.height() / g);
        vector<Module *> &members = clusterMembers[r];
        double a = clusterArea[r];
        float spread = (float)sqrt(a) * 0.6f;
        for (size_t mi = 0; mi < members.size(); mi++)  //簇内成员绕质心随机分布
        {
            Module *m = members[mi];
            float ox = (rand() / (float)RAND_MAX - 0.5f) * spread;
            float oy = (rand() / (float)RAND_MAX - 0.5f) * spread;
            float x = gx + ox - m->width * 0.5f;
            float y = gy + oy - m->height * 0.5f;
            x = max(core.ll.x, min(x, core.ur.x - m->width));
            y = max(core.ll.y, min(y, core.ur.y - m->height));
            m->ll = POS_2D(x, y);
        }
    }
}

// ===================== 方法3：二次方程解析布局（Kraftwerk2 力导向 + BiCGSTAB） =====================

void MyPlacer::QuadraticInitialPlacement()
{
    int n = (int)db->Nodes.size();
    if (n == 0)
        return;

    // 可移动单元 -> 索引映射
    unordered_map<Module *, int> idx;
    idx.reserve(n);
    for (int i = 0; i < n; i++)
        idx[db->Nodes[i]] = i;  // 记录每个可移动单元的索引

    // 初值：把所有可移动单元中心移到 core 中心，以便后续迭代收敛更快
    POS_2D cc = db->coreRegion.center();
    for (int i = 0; i < n; i++)
    {
        Module *m = db->Nodes[i];
        m->ll = POS_2D(cc.x - m->width * 0.5f, cc.y - m->height * 0.5f);
    }

    const float MIN_DIST = 1.0f;
    const int MAX_OUTER = 8;

    // ===== 阶段 A（一次性）：预建 pin 对列表，固化索引与边 ID，迭代阶段零哈希 =====
    struct PinPair
    {
        int eid, i, j;  //eid:边索引，i,j:可移动单元索引
        Pin *p, *q;     //p,q:对应的 pin 指针
    };
    struct TermPair
    {//不存储边索引，因为终端是固定的；不存储终端索引，因为终端不参与迭代
        int k;      // 可移动单元索引
        Pin *mov;   // 可移动端 pin
        Pin *term;  // 终端 pin
    };
    vector<PinPair> pairs;
    vector<TermPair> termPairs; 
    vector<pair<int, int>> edges; // 去重后的单元对 (为什么i<j？因为无向图边是对称的，i<j可以避免重复存储)
    unordered_map<long long, int> edgeId;   // 单元对 -> 边索引（i<j，key=i*n+j）
    edgeId.reserve((size_t)8000000);
    pairs.reserve((size_t)8000000);

    for (size_t e = 0; e < db->Nets.size(); e++)    // 遍历每个网络
    {
        Net *net = db->Nets[e];
        int d = (int)net->netPins.size();
        if (d < 2) 
            continue;
        for (int a = 0; a < d; a++) // 遍历网络中的每个 pin引脚
        {
            Pin *p = net->netPins[a];
            Module *mp = p->module; //引脚所属单元
            if (!mp)
                continue;
            for (int b = a + 1; b < d; b++) // 遍历网络中每个 pin 引脚的后续引脚
            {
                Pin *q = net->netPins[b];
                Module *mq = q->module; 
                if (!mq || mp == mq)
                    continue;
                bool fp = mp->isFixed, fq = mq->isFixed;
                if (!fp && !fq) // 两端都是可移动单元
                {
                    int i = idx[mp], j = idx[mq];
                    if (i == j) //同一个单元的两个引脚，跳过
                        continue;
                    long long key = keyOf(i, j, n); //n是可移动单元总数，key是唯一标识符，key = i*n+j，保证i<j
                    auto it = edgeId.find(key);
                    int eid;
                    if (it == edgeId.end()) // 如果该单元对还没有边索引，则创建一个新的边索引
                    {
                        eid = (int)edges.size();
                        edgeId[key] = eid;
                        edges.push_back(make_pair(min(i, j), max(i, j)));
                    }
                    else
                        eid = it->second;
                    PinPair pp;
                    pp.eid = eid;
                    pp.i = i;   //具体的可移动单元索引
                    pp.j = j;
                    pp.p = p;   //具体的 pin 指针
                    pp.q = q;
                    pairs.push_back(pp);
                }
                else if (!fp && fq)
                {
                    TermPair tp;
                    tp.k = idx[mp];
                    tp.mov = p;
                    tp.term = q;
                    termPairs.push_back(tp);
                }
                else if (fp && !fq)
                {
                    TermPair tp;
                    tp.k = idx[mq];
                    tp.mov = q;
                    tp.term = p;
                    termPairs.push_back(tp);
                }
            }
        }
    }

    int E = (int)edges.size();
    vector<double> diagX(n, 0.0), diagY(n, 0.0), offX(E, 0.0), offY(E, 0.0), bx(n, 0.0), by(n, 0.0);    // diagX/Y:对角线元素，offX/Y:非对角线元素，bx/by:右端项
    CRect core = db->coreRegion;    // 核心布局区域（SiteRow 的外包矩形）

    // ===== 阶段 B（迭代）：纯数组累加 + 求解 =====
    for (int outer = 0; outer < MAX_OUTER; outer++) // 外层迭代次数是8次
    {
        //fill()函数将数组的所有元素设置为指定的值，这里是将所有元素重置为0.0
        fill(diagX.begin(), diagX.end(), 0.0); 
        fill(diagY.begin(), diagY.end(), 0.0);
        fill(offX.begin(), offX.end(), 0.0);
        fill(offY.begin(), offY.end(), 0.0);
        fill(bx.begin(), bx.end(), 0.0);
        fill(by.begin(), by.end(), 0.0);

        // 情况 i：两端都可移动
        for (size_t t = 0; t < pairs.size(); t++)
        {
            PinPair &pp = pairs[t];
            POS_2D p_ = pinAbsPos(pp.p);    //pinAbsPos()函数返回引脚的绝对位置，即所属单元中心 + 相对中心偏移
            POS_2D q_ = pinAbsPos(pp.q);
            float wx = 1.0f / max(fabs(p_.x - q_.x), MIN_DIST); // 计算权重，避免除以零
            float wy = 1.0f / max(fabs(p_.y - q_.y), MIN_DIST);
            int i = pp.i, j = pp.j;
            diagX[i] += wx; //对角线加权重，非对角线负权重，因为在二次方程中，拉力是相互作用的，所以对角线是正的，非对角线是负的
            diagX[j] += wx;
            offX[pp.eid] -= wx;
            diagY[i] += wy;
            diagY[j] += wy;
            offY[pp.eid] -= wy;
            bx[i] += wx * (pp.q->offset.x - pp.p->offset.x);    //offset是引脚相对于所属单元中心的偏移，bx[i]是右端项，表示第i个可移动单元的x方向的偏移量
            bx[j] += wx * (pp.p->offset.x - pp.q->offset.x);    //bx作用是将引脚的偏移量转化为可移动单元的偏移量，从而在求解线性方程组时考虑引脚的影响
            by[i] += wy * (pp.q->offset.y - pp.p->offset.y);
            by[j] += wy * (pp.p->offset.y - pp.q->offset.y);
        }
        // 情况 ii/iii：一端是终端
        for (size_t t = 0; t < termPairs.size(); t++)
        {
            TermPair &tp = termPairs[t];
            POS_2D pm = pinAbsPos(tp.mov);
            POS_2D pt = pinAbsPos(tp.term);
            float wx = 1.0f / max(fabs(pm.x - pt.x), MIN_DIST);
            float wy = 1.0f / max(fabs(pm.y - pt.y), MIN_DIST);
            POS_2D xt = tp.term->module->center();
            diagX[tp.k] += wx;
            diagY[tp.k] += wy;
            bx[tp.k] += wx * (xt.x + tp.term->offset.x - tp.mov->offset.x);
            by[tp.k] += wy * (xt.y + tp.term->offset.y - tp.mov->offset.y);
        }

        // 构建稀疏矩阵（结构固定，值随迭代变化）
        vector<Eigen::Triplet<float>> tx, ty;
        tx.reserve((size_t)n + E * 2);  //n个对角，E条边，每条边有两个非对角元素
        ty.reserve((size_t)n + E * 2);
        for (int i = 0; i < n; i++)
        {
            tx.push_back(Eigen::Triplet<float>(i, i, (float)diagX[i]));
            ty.push_back(Eigen::Triplet<float>(i, i, (float)diagY[i]));
        }
        for (int e = 0; e < E; e++)
        {
            int i = edges[e].first, j = edges[e].second;
            tx.push_back(Eigen::Triplet<float>(i, j, (float)offX[e]));
            tx.push_back(Eigen::Triplet<float>(j, i, (float)offX[e]));
            ty.push_back(Eigen::Triplet<float>(i, j, (float)offY[e]));
            ty.push_back(Eigen::Triplet<float>(j, i, (float)offY[e]));
        }

        Eigen::SparseMatrix<float, Eigen::RowMajor> Ax(n, n), Ay(n, n);
        Ax.setFromTriplets(tx.begin(), tx.end());
        Ay.setFromTriplets(ty.begin(), ty.end());

        Eigen::VectorXf Bx(n), By(n), Xx(n), Xy(n);
        for (int i = 0; i < n; i++)
        {
            Bx[i] = (float)bx[i];
            By[i] = (float)by[i];
            Xx[i] = db->Nodes[i]->center().x;
            Xy[i] = db->Nodes[i]->center().y;
        }

        Eigen::ConjugateGradient<Eigen::SparseMatrix<float, Eigen::RowMajor>> solverX, solverY;
        solverX.setMaxIterations(500);
        solverX.setTolerance(1e-6f);
        solverY.setMaxIterations(500);
        solverY.setTolerance(1e-6f);
        solverX.compute(Ax);
        solverY.compute(Ay);
        Eigen::VectorXf nx = solverX.solveWithGuess(Bx, Xx);
        Eigen::VectorXf ny = solverY.solveWithGuess(By, Xy);

        // 更新位置并 clamp 到 core 内
        for (int i = 0; i < n; i++)
        {
            Module *m = db->Nodes[i];
            float cx = nx[i], cy = ny[i];
            if (cx - m->width * 0.5f < core.ll.x) cx = core.ll.x + m->width * 0.5f;
            if (cx + m->width * 0.5f > core.ur.x) cx = core.ur.x - m->width * 0.5f;
            if (cy - m->height * 0.5f < core.ll.y) cy = core.ll.y + m->height * 0.5f;
            if (cy + m->height * 0.5f > core.ur.y) cy = core.ur.y - m->height * 0.5f;
            m->ll = POS_2D(cx - m->width * 0.5f, cy - m->height * 0.5f);
        }
    }
}

// ===================== 任务5：全局布局（ePlace 电静力学 + Nesterov） =====================

static const float DCT_PI = 3.14159265358979323846f;

// 2D 离散余弦变换求解器（自实现、零外部依赖）
// 依据任务描述公式：对密度 ρ 做 DCT 得到系数 a_jk，再由 a_jk 加权求和得到电场 ξ_x、ξ_y。
// 采用可分离 2D DCT，复杂度 O(m^3)，m<=256 时足够快。
// 注意：此为自实现版本；后续如需更高精度/性能，可替换为 FFTW 等官方 FFT 库实现。
struct DCTSolver
{
    int m;
    vector<float> cosT, sinT;            // cosT[j*m+x]=cos(jπx/m)，sinT 同理
    vector<float> tmpA, a, G, Sx, Sy;    // 复用缓冲区

    DCTSolver(int _m) : m(_m)
    {
        cosT.resize(m * m);
        sinT.resize(m * m);
        for (int j = 0; j < m; j++)
            for (int x = 0; x < m; x++)
            {
                float ang = DCT_PI * (float)j * x / m;
                cosT[j * m + x] = cosf(ang);
                sinT[j * m + x] = sinf(ang);
            }
        tmpA.resize(m * m);
        a.resize(m * m);
        G.resize(m * m);
        Sx.resize(m * m);
        Sy.resize(m * m);
    }

    // 输入 rho[m*m]（行优先 y*m+x），输出电场 ex、ey[m*m]
    void solve(const vector<float> &rho, vector<float> &ex, vector<float> &ey)
    {
        int m2 = m * m;
        // 1) 列方向 DCT：tmpA[k][x] = Σ_y rho(x,y) cos(kπy/m)
        for (int x = 0; x < m; x++)
            for (int k = 0; k < m; k++)
            {
                float s = 0;
                const float *cosk = &cosT[k * m];
                for (int y = 0; y < m; y++)
                    s += rho[y * m + x] * cosk[y];
                tmpA[k * m + x] = s;
            }
        // 2) 行方向 DCT：a[j][k] = (1/m²) Σ_x tmpA[k][x] cos(jπx/m)
        float inv = 1.0f / (float)m2;
        for (int k = 0; k < m; k++)
            for (int j = 0; j < m; j++)
            {
                float s = 0;
                const float *cosj = &cosT[j * m];
                for (int x = 0; x < m; x++)
                    s += tmpA[k * m + x] * cosj[x];
                a[j * m + k] = s * inv;
            }
        // 3) ξ_x：G[j][k]=a·w_j/(w_j²+w_k²)；Sx[j][y]=Σ_k G cos；ex=Σ_j Sx sin
        for (int j = 0; j < m; j++)
        {
            float wj = DCT_PI * (float)j / m;
            for (int k = 0; k < m; k++)
            {
                float wk = DCT_PI * (float)k / m;
                float denom = wj * wj + wk * wk;
                G[j * m + k] = (denom < 1e-12f) ? 0.0f : (a[j * m + k] * wj / denom);
            }
        }
        for (int j = 0; j < m; j++)
            for (int y = 0; y < m; y++)
            {
                float s = 0;
                for (int k = 0; k < m; k++)
                    s += G[j * m + k] * cosT[k * m + y];
                Sx[j * m + y] = s;
            }
        for (int y = 0; y < m; y++)
            for (int x = 0; x < m; x++)
            {
                float s = 0;
                for (int j = 0; j < m; j++)
                    s += Sx[j * m + y] * sinT[j * m + x];
                ex[y * m + x] = s;
            }
        // 4) ξ_y：G[j][k]=a·w_k/(w_j²+w_k²)；Sy[k][x]=Σ_j G cos(jπx/m)；ey=Σ_k Sy sin(kπy/m)
        for (int j = 0; j < m; j++)
        {
            float wj = DCT_PI * (float)j / m;
            for (int k = 0; k < m; k++)
            {
                float wk = DCT_PI * (float)k / m;
                float denom = wj * wj + wk * wk;
                G[j * m + k] = (denom < 1e-12f) ? 0.0f : (a[j * m + k] * wk / denom);
            }
        }
        for (int k = 0; k < m; k++)
            for (int x = 0; x < m; x++)
            {
                float s = 0;
                for (int j = 0; j < m; j++)
                    s += G[j * m + k] * cosT[j * m + x];
                Sy[k * m + x] = s;
            }
        for (int y = 0; y < m; y++)
            for (int x = 0; x < m; x++)
            {
                float s = 0;
                for (int k = 0; k < m; k++)
                    s += Sy[k * m + x] * sinT[k * m + y];
                ey[y * m + x] = s;
            }
    }
};

// ---- filler 初始化：填充空白到目标密度 ----
void MyPlacer::fillerInit()
{
    CRect core = db->coreRegion;
    StdCellArea = MacroArea = 0;
    for (Module *m : db->Nodes)
    {
        if (m->isMacro)
            MacroArea += m->area;
        else
            StdCellArea += m->area;
    }
    // filler 总面积 = 目标密度 × 可放置面积 − 单元面积（填充空白到目标密度）
    // 可放置面积 = core 面积 − 终端（固定 I/O pad）在 core 内的重叠面积
    float termArea = 0;
    for (Module *t : db->Terminals)
    {
        float ox = max(0.0f, min(core.ur.x, t->ur().x) - max(core.ll.x, t->ll.x));
        float oy = max(0.0f, min(core.ur.y, t->ur().y) - max(core.ll.y, t->ll.y));
        termArea += ox * oy;
    }
    float fillerArea = targetDensity * (core.area() - termArea) - StdCellArea - MacroArea;
    if (fillerArea < 0)
        fillerArea = 0;

    // filler 尺寸 = 中间 80% 模块的平均面积
    vector<float> areas;
    areas.reserve(db->Nodes.size());
    for (Module *m : db->Nodes)
        areas.push_back(m->area);
    sort(areas.begin(), areas.end());
    int lo = (int)(areas.size() * 0.1f), hi = (int)(areas.size() * 0.9f);
    float sum = 0;
    int cnt = 0;
    for (int i = lo; i < hi; i++)
    {
        sum += areas[i];
        cnt++;
    }
    float avgFillerArea = cnt > 0 ? sum / cnt : 1.0f;
    float fillerSize = sqrtf(avgFillerArea);

    int nFiller = (int)(fillerArea / avgFillerArea);
    if (nFiller > 500000)
        nFiller = 500000; // 限制 filler 数量，避免内存过大

    srand(54321);
    for (int i = 0; i < nFiller; i++)
    {
        Module *f = new Module();
        f->name = "filler_" + to_string(i);
        f->width = fillerSize;
        f->height = fillerSize;
        f->area = f->width * f->height;
        f->isFiller = true;
        f->isFixed = false;
        f->isMacro = false;
        f->idx = (int)db->Nodes.size() + i;
        float x = core.ll.x + (rand() / (float)RAND_MAX) * (core.width() - f->width);
        float y = core.ll.y + (rand() / (float)RAND_MAX) * (core.height() - f->height);
        f->ll = POS_2D(x, y);
        Fillers.push_back(f);
    }

    NodesAndFillers.clear();
    NodesAndFillers.reserve(db->Nodes.size() + Fillers.size());
    for (Module *m : db->Nodes)
        NodesAndFillers.push_back(m);
    for (Module *f : Fillers)
        NodesAndFillers.push_back(f);
}

// ---- bin 网格初始化 ----
void MyPlacer::binInit()
{
    int total = (int)NodesAndFillers.size();
    binDimension = (int)round(sqrt((double)total));
    binDimension = max(16, min(256, binDimension));
    int m = binDimension;
    CRect core = db->coreRegion;
    binWidth = core.width() / m;
    binHeight = core.height() / m;

    for (auto &row : bins)
        for (auto b : row)
            delete b;
    bins.assign(m, vector<Bin_2D *>(m, nullptr));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
        {
            Bin_2D *b = new Bin_2D();
            b->ll = POS_2D(core.ll.x + i * binWidth, core.ll.y + j * binHeight);
            b->ur = POS_2D(b->ll.x + binWidth, b->ll.y + binHeight);
            b->center = POS_2D((b->ll.x + b->ur.x) * 0.5f, (b->ll.y + b->ur.y) * 0.5f);
            b->width = binWidth;
            b->height = binHeight;
            b->area = binWidth * binHeight;
            bins[i][j] = b;
        }

    // terminalDensity：终端与 bin 的重叠面积 × targetDensity
    for (Module *t : db->Terminals)
    {
        float tx1 = t->ll.x, ty1 = t->ll.y, tx2 = t->ur().x, ty2 = t->ur().y;
        int i0 = max(0, (int)floor((tx1 - core.ll.x) / binWidth));
        int i1 = min(m - 1, (int)floor((tx2 - core.ll.x) / binWidth));
        int j0 = max(0, (int)floor((ty1 - core.ll.y) / binHeight));
        int j1 = min(m - 1, (int)floor((ty2 - core.ll.y) / binHeight));
        for (int i = i0; i <= i1; i++)
            for (int j = j0; j <= j1; j++)
            {
                float ox = max(0.0f, min(bins[i][j]->ur.x, tx2) - max(bins[i][j]->ll.x, tx1));
                float oy = max(0.0f, min(bins[i][j]->ur.y, ty2) - max(bins[i][j]->ll.y, ty1));
                bins[i][j]->terminalDensity += targetDensity * ox * oy;
            }
    }
    // darkDensity：core 内均可放置（暗区域在 core 外，不纳入网格），保持 0

    if (dctSolver)
    {
        delete dctSolver;
        dctSolver = nullptr;
    }
    dctSolver = new DCTSolver(m);
}

// ---- 梯度向量初始化 ----
void MyPlacer::gradientVectorInitialization()
{
    wirelengthGradient.assign(db->Nodes.size(), VECTOR_3D());
    densityGradient.assign(NodesAndFillers.size(), VECTOR_3D());
    totalGradient.assign(NodesAndFillers.size(), VECTOR_3D());
    fillerGradient.assign(Fillers.size(), VECTOR_3D());
    nodeIndex.clear();
    nodeIndex.reserve(db->Nodes.size());
    for (int i = 0; i < (int)db->Nodes.size(); i++)
        nodeIndex[db->Nodes[i]] = i;

    // 预条件子用：每个可移动单元的线网连接度（入射线网数）
    // 参考 RePlAce/DREAMPlace 的 sum_pin_weights_in_nodes，作为近似的 Hessian 对角项
    pinWeights.assign(db->Nodes.size(), 0.0f);
    for (Net *net : db->Nets)
    {
        if (net->netPins.size() < 2)
            continue;
        for (Pin *p : net->netPins)
        {
            if (p->module->isFixed)
                continue;
            auto it = nodeIndex.find(p->module);
            if (it != nodeIndex.end())
                pinWeights[it->second] += 1.0f;
        }
    }
}

// ---- 更新密度（单元 -> bin 重叠面积）----
void MyPlacer::updateDensity()
{
    int m = binDimension;
    CRect core = db->coreRegion;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
        {
            bins[i][j]->nodeDensity = 0;
            bins[i][j]->fillerDensity = 0;
        }
    for (Module *mod : NodesAndFillers)
    {
        float x1 = mod->ll.x, y1 = mod->ll.y, x2 = mod->ur().x, y2 = mod->ur().y;
        int i0 = max(0, (int)floor((x1 - core.ll.x) / binWidth));
        int i1 = min(m - 1, (int)floor((x2 - core.ll.x) / binWidth));
        int j0 = max(0, (int)floor((y1 - core.ll.y) / binHeight));
        int j1 = min(m - 1, (int)floor((y2 - core.ll.y) / binHeight));
        for (int i = i0; i <= i1; i++)
            for (int j = j0; j <= j1; j++)
            {
                Bin_2D *b = bins[i][j];
                float ox = max(0.0f, min(b->ur.x, x2) - max(b->ll.x, x1));
                float oy = max(0.0f, min(b->ur.y, y2) - max(b->ll.y, y1));
                float overlap = ox * oy;
                // 正电荷 = 单元真实面积（ePlace 密度模型：标准/宏/filler 均按面积计，
                // 目标密度通过 filler 补足、再由 DCT 直流项隐式减去均值实现，
                // 而不是对每个单元乘 targetDensity）。固定终端才乘 targetDensity（见 binInit）。
                if (mod->isFiller)
                    b->fillerDensity += overlap;
                else
                    b->nodeDensity += overlap;
            }
    }
}

// ---- 更新密度梯度（DCT 解 Poisson 求电场）----
void MyPlacer::updateDensityGradient()
{
    int m = binDimension;
    CRect core = db->coreRegion;
    vector<float> rho(m * m), ex(m * m), ey(m * m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
        {
            Bin_2D *b = bins[i][j];
            float total = b->nodeDensity + b->fillerDensity + b->terminalDensity + b->DarkDensity;
            rho[j * m + i] = total / b->area;
        }
    dctSolver->solve(rho, ex, ey);

    int nodeCount = (int)db->Nodes.size();
    // 密度梯度 = -Σ_bin (单元与 bin 的重叠面积) × ξ(bin)
    // 这是本布局器密度模型（updateDensity 的 box overlap）的精确负梯度，
    // 而非对单元中心单点采样。对跨度多个 bin 的宏单元，电场在单元内部变化显著，
    // 单点采样（-q·ξ(center)）会严重低估/错估其受力方向，导致宏单元畸变。
    // 参考 DREAMPlace electric_force.cpp 的面积积分做法。
    for (int idx = 0; idx < (int)NodesAndFillers.size(); idx++)
    {
        Module *mod = NodesAndFillers[idx];
        float x1 = mod->ll.x, y1 = mod->ll.y, x2 = mod->ur().x, y2 = mod->ur().y;
        int i0 = max(0, (int)floor((x1 - core.ll.x) / binWidth));
        int i1 = min(m - 1, (int)floor((x2 - core.ll.x) / binWidth));
        int j0 = max(0, (int)floor((y1 - core.ll.y) / binHeight));
        int j1 = min(m - 1, (int)floor((y2 - core.ll.y) / binHeight));
        float gx = 0, gy = 0;
        for (int i = i0; i <= i1; i++)
            for (int j = j0; j <= j1; j++)
            {
                float ox = max(0.0f, min(bins[i][j]->ur.x, x2) - max(bins[i][j]->ll.x, x1));
                float oy = max(0.0f, min(bins[i][j]->ur.y, y2) - max(bins[i][j]->ll.y, y1));
                float overlap = ox * oy;
                gx += overlap * ex[j * m + i];
                gy += overlap * ey[j * m + i];
            }
        densityGradient[idx].x = -gx;
        densityGradient[idx].y = -gy;
        densityGradient[idx].z = 0;
    }
    for (int f = 0; f < (int)Fillers.size(); f++)
        fillerGradient[f] = densityGradient[nodeCount + f];
}

// ---- 更新线长梯度（WA 平滑线长模型）----
void MyPlacer::updateWirelengthGradient()
{
    for (auto &g : wirelengthGradient)
        g.SetZero();
    float g = gamma;
    vector<float> xs, ys, tp, tm;
    for (Net *net : db->Nets)
    {
        int d = (int)net->netPins.size();
        if (d < 2)
            continue;
        xs.resize(d);
        ys.resize(d);
        tp.resize(d);
        tm.resize(d);
        float xmax = -1e30f, xmin = 1e30f, ymax = -1e30f, ymin = 1e30f;
        for (int i = 0; i < d; i++)
        {
            POS_2D a = pinAbsPos(net->netPins[i]);
            xs[i] = a.x;
            ys[i] = a.y;
            xmax = max(xmax, a.x);
            xmin = min(xmin, a.x);
            ymax = max(ymax, a.y);
            ymin = min(ymin, a.y);
        }
        // x 方向
        float bp = 0, cp = 0, bm = 0, cm = 0;
        for (int i = 0; i < d; i++)
        {
            tp[i] = expf((xs[i] - xmax) / g);
            tm[i] = expf(-(xs[i] - xmin) / g);
            bp += tp[i];
            cp += xs[i] * tp[i];
            bm += tm[i];
            cm += xs[i] * tm[i];
        }
        for (int i = 0; i < d; i++)
        {
            Pin *p = net->netPins[i];
            if (p->module->isFixed)
                continue;
            int idx = nodeIndex[p->module];
            float gx = (bp + xs[i] * bp / g - cp / g) / (bp * bp) * tp[i] - (bm - xs[i] * bm / g + cm / g) / (bm * bm) * tm[i];
            wirelengthGradient[idx].x += gx;
        }
        // y 方向
        bp = cp = bm = cm = 0;
        for (int i = 0; i < d; i++)
        {
            tp[i] = expf((ys[i] - ymax) / g);
            tm[i] = expf(-(ys[i] - ymin) / g);
            bp += tp[i];
            cp += ys[i] * tp[i];
            bm += tm[i];
            cm += ys[i] * tm[i];
        }
        for (int i = 0; i < d; i++)
        {
            Pin *p = net->netPins[i];
            if (p->module->isFixed)
                continue;
            int idx = nodeIndex[p->module];
            float gy = (bp + ys[i] * bp / g - cp / g) / (bp * bp) * tp[i] - (bm - ys[i] * bm / g + cm / g) / (bm * bm) * tm[i];
            wirelengthGradient[idx].y += gy;
        }
    }
}

// ---- 更新总梯度（含预条件子）----
// 参考 RePlAce/DREAMPlace 的 PreconditionOp：
//   precond_i = sum_pin_weights_in_nodes_i + alpha * density_weight * node_area_i
//   grad_i   /= max(precond_i, 1.0)
// 作用：λ 的全局平衡使密度项相对线长项很小，若不预条件，高连接度的宏单元
// （线长梯度 ∝ 引脚数，密度梯度 ∝ 面积）的梯度量级是标准单元的数十上百倍，
// 会让宏单元移动过快、拉扯整个布局发散。除以连接度+面积项后各单元步长趋于一致。
void MyPlacer::updateTotalGradient()
{
    int nodeCount = (int)db->Nodes.size();
    const float alpha = 1.0f; // 面积项权重（DREAMPlace 中自适应，这里取常量 1）
    for (int i = 0; i < nodeCount; i++)
    {
        totalGradient[i].x = wirelengthGradient[i].x + lambda * densityGradient[i].x;
        totalGradient[i].y = wirelengthGradient[i].y + lambda * densityGradient[i].y;
        totalGradient[i].z = 0;
        float precond = pinWeights[i] + alpha * lambda * NodesAndFillers[i]->area;
        if (precond < 1.0f)
            precond = 1.0f;
        totalGradient[i].x /= precond;
        totalGradient[i].y /= precond;
    }
    for (int i = nodeCount; i < (int)NodesAndFillers.size(); i++)
    {
        totalGradient[i].x = lambda * densityGradient[i].x;
        totalGradient[i].y = lambda * densityGradient[i].y;
        totalGradient[i].z = 0;
        float precond = alpha * lambda * NodesAndFillers[i]->area;
        if (precond < 1.0f)
            precond = 1.0f;
        totalGradient[i].x /= precond;
        totalGradient[i].y /= precond;
    }
}

// ---- 惩罚系数 λ 初始化（平衡线长梯度与密度梯度量级）----
// 参考 RePlAce/DREAMPlace：λ0 = density_weight_base × (Σ|WLgrad| / Σ|densityGrad|)。
// 关键：必须乘以一个很小的基础权重 8e-5（DREAMPlace params.json 默认值），
// 让密度项在初始时远弱于线长项，从而保住二次布局已优化出的 HPWL；
// 之后靠 RePlAce 的 μ∈[0.95,1.05] 反馈逐步把 λ 抬高到"恰好能铺开"的水平。
// 若直接用裸比值（Σ|WLgrad|/Σ|densityGrad|），宏单元主导的密度梯度会让 λ 偏大
// 约 4 个量级，密度项一开始就压过线长项，把二次布局的 HPWL 彻底搅乱（相差约 10 倍）。
void MyPlacer::penaltyFactorInitilization()
{
    const float densityWeightBase = 8e-5f; // RePlAce/DREAMPlace 密度基础权重
    double numerator = 0, denominator = 0;
    for (auto &g : wirelengthGradient)
        numerator += fabs((double)g.x) + fabs((double)g.y);
    for (auto &g : densityGradient)
        denominator += fabs((double)g.x) + fabs((double)g.y);
    lambda = (denominator < 1e-12) ? 0.0f : (float)(densityWeightBase * numerator / denominator);
    if (lambda < 1e-12f)
        lambda = 1e-12f;
}

// ---- 密度溢出率 τ ----
float MyPlacer::calcOverflow()
{
    int m = binDimension;
    double totalOverflow = 0;
    double normArea = StdCellArea + MacroArea;
    if (normArea < 1e-12)
        return 0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
        {
            Bin_2D *b = bins[i][j];
            // 只计算可移动单元（标准+宏+filler）的溢出；终端/暗节点是固定的，
            // 其占据面积从该 bin 的容量中扣除（不可放置）。
            float movableArea = b->nodeDensity + b->fillerDensity;
            float blockArea = (b->terminalDensity + b->DarkDensity) / targetDensity;
            float capacity = targetDensity * (b->area - blockArea);
            float over = movableArea - capacity;
            if (over > 0)
                totalOverflow += over;
        }
    return (float)(totalOverflow / normArea);
}

// ---- 全局布局主流程（ePlace）----
void MyPlacer::GlobalPlacement()
{
    fillerInit();
    binInit();
    gradientVectorInitialization();
    stepSize = 1.0f * binWidth; // BB 步长初值（量级匹配预条件后的梯度）

    int totalCount = (int)NodesAndFillers.size();
    int maxIter = 1000;
    vector<POS_2D> vel(totalCount);
    for (auto &v : vel)
        v.SetZero();
    CRect core = db->coreRegion;
    const float momentum = 0.8f;
    const float stop_overflow = 0.1f; // DREAMPlace 停止判定的密度溢出阈值
    const float maxMove = 8.0f * binWidth; // 每步最大位移，防止 overshoot 发散
    bool lambdaFrozen = false;             // 密度解开后冻结 λ，避免其无界增长
    double bestHpwl = 1e300;               // 追踪最优 HPWL，用于“已收敛/过度铺开”判定
    int bestIter = 0;

    // Barzilai-Borwein 自适应步长状态（参考 DREAMPlace Nesterov+BB 的 step_bb）：
    // bb_short = <s,y>/<y,y>，s 为位置增量、y 为梯度增量。在“梯度僵硬”的细化阶段
    // 会自动缩小步长、在铺开阶段自动放大，避免固定步长在细化期震荡/过度铺开。
    vector<float> prevCx(totalCount), prevCy(totalCount); // 上一步中心
    vector<float> prevGx(totalCount), prevGy(totalCount); // 上一步梯度

    // 2. 第一次梯度计算 + λ 初始化
    updateDensity();
    float tau = calcOverflow();
    gamma = 4.0f * (binWidth + binHeight) * powf(10.0f, (tau - 0.1f) * 20.0f / 9.0f - 1.0f);
    if (gamma < 1e-3f)
        gamma = 1e-3f;
    if (gamma > 20.0f * max(core.width(), core.height()))
        gamma = 20.0f * max(core.width(), core.height());
    updateDensityGradient();
    updateWirelengthGradient();
    penaltyFactorInitilization();
    updateTotalGradient();
    float lambdaInit = lambda;
    double hpwlPrev = calcHPWL(); // 用于 λ 自适应更新（步骤 2-3）
    double hpwlInit = hpwlPrev;

    // 初始化 BB 步长状态（首步：s=0 前的位置，g=当前梯度）
    for (int i = 0; i < totalCount; i++)
    {
        POS_2D c = NodesAndFillers[i]->center();
        prevCx[i] = c.x;
        prevCy[i] = c.y;
        prevGx[i] = totalGradient[i].x;
        prevGy[i] = totalGradient[i].y;
    }

    printf("===== 全局布局 (ePlace) =====\n");
    printf("  bin 网格: %d x %d   filler: %zu   targetDensity=%.2f   step=%.3f\n",
           binDimension, binDimension, Fillers.size(), targetDensity, stepSize);
    printf("  iter   HPWL          overflow   lambda      step\n");

    // 3. 迭代（动量加速梯度下降 + Barzilai-Borwein 自适应步长）
    for (int it = 0; it < maxIter; it++)
    {
        // 更新速度与位置（带位移限制，防止发散）
        for (int i = 0; i < totalCount; i++)
        {
            Module *mod = NodesAndFillers[i];
            float vx = momentum * vel[i].x + stepSize * totalGradient[i].x;
            float vy = momentum * vel[i].y + stepSize * totalGradient[i].y;
            float mag = sqrtf(vx * vx + vy * vy);
            if (mag > maxMove)
            {
                vx *= maxMove / mag;
                vy *= maxMove / mag;
            }
            vel[i].x = vx;
            vel[i].y = vy;
            POS_2D c = mod->center();
            float cx = c.x - vx;
            float cy = c.y - vy;
            if (cx - mod->width * 0.5f < core.ll.x)
                cx = core.ll.x + mod->width * 0.5f;
            if (cx + mod->width * 0.5f > core.ur.x)
                cx = core.ur.x - mod->width * 0.5f;
            if (cy - mod->height * 0.5f < core.ll.y)
                cy = core.ll.y + mod->height * 0.5f;
            if (cy + mod->height * 0.5f > core.ur.y)
                cy = core.ur.y - mod->height * 0.5f;
            mod->ll = POS_2D(cx - mod->width * 0.5f, cy - mod->height * 0.5f);
        }

        // 重新计算密度与溢出率
        updateDensity();
        tau = calcOverflow();
        if (tau < stop_overflow)
            lambdaFrozen = true;
        gamma = 4.0f * (binWidth + binHeight) * powf(10.0f, (tau - 0.1f) * 20.0f / 9.0f - 1.0f);
        if (gamma < 1e-3f)
            gamma = 1e-3f;
        if (gamma > 20.0f * max(core.width(), core.height()))
            gamma = 20.0f * max(core.width(), core.height());

        // 更新 λ（参考 RePlAce 的密度权重更新，ePlace 论文所用 μ 反馈机制）：
        //   HPWL 改善 → 缓慢增大 λ（逐步加强密度惩罚）；
        //   HPWL 恶化 → 按恶化幅度回落 λ（防止密度项过度挤压线长）。
        // 关键：μ 只在 [0.95, 1.05] 附近小步变化，且允许 <1 回退，
        //       而非上一版 [1.5, 3.0] 的强行单调增长（那会让 λ 十几轮就封顶、
        //       密度项长期压过线长项，把二次布局已优化的 HPWL 破坏掉）。
        double hpwl = calcHPWL();
        double delta = hpwl - hpwlPrev;
        double refHPWL = max(3.5e5, hpwlInit * 0.01);
        // 只在铺开结束（τ 已明显下降）后追踪最优 HPWL；否则初始二次布局的
        // 低 HPWL 会让 bestIter 卡在 0，铺开阶段一结束就误判为"长期无改善"。
        if (tau < 0.5f && hpwl < bestHpwl)
        {
            bestHpwl = hpwl;
            bestIter = it;
        }
        // 密度已基本解开后冻结 λ：否则 RePlAce 的 μ 反馈在 HPWL 微幅波动时
        // 仍保持 μ≈1.03~1.05，λ 会以 3%/步复利无界增长，把已收敛的布局
        // 继续往外挤、反而抬升 HPWL（过度铺开）。
        if (!lambdaFrozen)
        {
            // 密度已基本解开（τ<0.25）且线长开始恶化 → 冻结 λ（即 DREAMPlace
            // Llambda 停止准则“overflow<stop && hpwl>prev”，此处把 0.1 放宽到 0.25，
            // 因本实现 τ 底部约 0.105、0.1 几乎达不到）。否则 μ 在 HPWL 微幅波动时
            // 仍保持 1.03~1.05，λ 复利无界增长，把已收敛布局继续挤开、白损 HPWL。
            if (tau < 0.25f && delta > 0)
            {
                lambdaFrozen = true;
            }
            else
            {
                const float UPPER_PCOF = 1.05f, LOWER_PCOF = 0.95f;
                float mu;
                if (delta < 0)
                {
                    // HPWL 改善：μ ≈ 1.05 × max(0.9999^it, 0.98)，迭代越靠后增幅越小
                    mu = UPPER_PCOF * (float)max(pow(0.9999, (double)it), 0.98);
                }
                else
                {
                    // HPWL 恶化：μ = 1.05 × 1.05^(-ΔHPWL/ref)，夹到 [0.95, 1.05]
                    mu = UPPER_PCOF * (float)pow((double)UPPER_PCOF, -delta / refHPWL);
                    mu = max(LOWER_PCOF, min(UPPER_PCOF, mu));
                }
                lambda *= mu;
                lambda = min(lambda, lambdaInit * 1e15f);
            }
        }
        hpwlPrev = hpwl;

        // 用新 λ 重新计算梯度
        updateDensityGradient();
        updateWirelengthGradient();
        updateTotalGradient();

        // Barzilai-Borwein 自适应步长（s = 位置增量，y = 梯度增量）
        {
            double sdoty = 0, ydoty = 0, sdots = 0;
            for (int i = 0; i < totalCount; i++)
            {
                POS_2D c = NodesAndFillers[i]->center();
                float sx = c.x - prevCx[i], sy = c.y - prevCy[i];
                float yx = totalGradient[i].x - prevGx[i];
                float yy = totalGradient[i].y - prevGy[i];
                sdoty += (double)sx * yx + (double)sy * yy;
                ydoty += (double)yx * yx + (double)yy * yy;
                sdots += (double)sx * sx + (double)sy * sy;
            }
            float newStep;
            if (sdoty > 0 && ydoty > 0)
                newStep = (float)(sdoty / ydoty); // bb_short：细化阶段自动变小
            else
                newStep = (float)min(sqrt(sdots / ydoty), (double)stepSize); // 回退
            if (!(newStep > 0.0f)) // 含 NaN 保护
                newStep = stepSize;
            stepSize = max(0.1f * binWidth, min(8.0f * binWidth, newStep));
            for (int i = 0; i < totalCount; i++)
            {
                POS_2D c = NodesAndFillers[i]->center();
                prevCx[i] = c.x;
                prevCy[i] = c.y;
                prevGx[i] = totalGradient[i].x;
                prevGy[i] = totalGradient[i].y;
            }
        }

        if (it % 5 == 0 || it == maxIter - 1)
        {
            printf("  %4d   %-13.0f  %.4f   %.4e   %.3f\n",
                   it, hpwl, tau, lambda, stepSize);
        }

        // 停止判定 1：密度溢出率低于目标阈值（DREAMPlace stop_overflow）
        if (it > 10 && tau < stop_overflow)
        {
            printf("  密度溢出率 %.4f < %.2f，提前停止\n", tau, stop_overflow);
            break;
        }
        // 停止判定 2：密度已基本解开后 HPWL 长期无改善 → 已收敛，避免过度铺开
        if (it > 200 && tau < 0.25f && it - bestIter > 100)
        {
            printf("  HPWL 已 %d 轮无改善（best=%.0f @%d），停止过度铺开\n",
                   it - bestIter, bestHpwl, bestIter);
            break;
        }
    }

    printf("  全局布局完成，HPWL = %.0f\n", calcHPWL());
}

// ---- 任务7：Abacus 标准单元合法化 ----
// 核心思想（参考 ISPD'08 Abacus）：按 x 坐标从左到右处理标准单元，
// 每个单元放到离它当前 y 最近的行；行内贪心放置以最小化水平位移；
// 行放满则顺延到下一行（溢出处理）。宏单元保持全局布局的位置不动。
void MyPlacer::AbacusLegalization()
{
    vector<SiteRow> &rows = db->SiteRows;
    if (rows.empty())
    {
        cerr << "[合法化] 没有 SiteRow，无法合法化" << endl;
        return;
    }

    double hpwlBefore = calcHPWL();

    // 行按 y（bottom）从小到大排序，保证“下一行”语义正确
    sort(rows.begin(), rows.end(), [](const SiteRow &a, const SiteRow &b)
         { return a.bottom < b.bottom; });

    // 收集标准单元（非宏），按当前 x 排序
    vector<Module *> cells;
    cells.reserve(db->Nodes.size());
    for (Module *m : db->Nodes)
        if (!m->isMacro)
            cells.push_back(m);
    sort(cells.begin(), cells.end(), [](const Module *a, const Module *b)
         { return a->ll.x < b->ll.x; });

    // 每行已放置到的最右 x（初始为行左边界）
    vector<double> rowX(rows.size());
    for (size_t r = 0; r < rows.size(); r++)
        rowX[r] = rows[r].start.x;

    double totalMove = 0;
    int placed = 0;

    for (Module *c : cells)
    {
        double cy = c->center().y;

        // 找 y 最近的行
        int best = 0;
        double bestDist = 1e30;
        for (size_t r = 0; r < rows.size(); r++)
        {
            double ry = rows[r].bottom + rows[r].height * 0.5;
            double d = fabs(cy - ry);
            if (d < bestDist)
            {
                bestDist = d;
                best = (int)r;
            }
        }

        // 从 best 行开始，向下找第一个能放下的行
        int placedRow = -1;
        for (int r = best; r < (int)rows.size(); r++)
        {
            double xPos = max(rowX[r], (double)c->ll.x);
            if (xPos + c->width <= rows[r].end.x + 1e-6)
            {
                placedRow = r;
                break;
            }
        }
        // 向下放不下则向上扫
        if (placedRow < 0)
        {
            for (int r = best - 1; r >= 0; r--)
            {
                double xPos = max(rowX[r], (double)c->ll.x);
                if (xPos + c->width <= rows[r].end.x + 1e-6)
                {
                    placedRow = r;
                    break;
                }
            }
        }
        if (placedRow < 0)
            placedRow = best; // 极端情况：仍放不下，硬放到 best 行（容忍少量溢出）

        // 水平位置：尽量保持原位，否则贴到行当前最右
        double xPos = max(rowX[placedRow], (double)c->ll.x);
        double maxX = rows[placedRow].end.x - c->width;
        if (xPos > maxX)
            xPos = maxX;
        if (xPos < rows[placedRow].start.x)
            xPos = rows[placedRow].start.x;

        // 垂直位置：行底对齐
        double yPos = rows[placedRow].bottom;

        totalMove += fabs(xPos - c->ll.x) + fabs(yPos - c->ll.y);
        c->ll = POS_2D((float)xPos, (float)yPos);
        rowX[placedRow] = xPos + c->width;
        placed++;
    }

    printf("===== 合法化 (Abacus) =====\n");
    printf("  标准单元 %d 个已放入 %zu 行，总位移 = %.0f\n", placed, rows.size(), totalMove);
    printf("  合法化后 HPWL = %.0f（合法化前 %.0f）\n", calcHPWL(), hpwlBefore);
}


