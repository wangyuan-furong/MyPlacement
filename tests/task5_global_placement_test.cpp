#include "myplace.h"
#include "electrostatics.h"
#include <memory>
#include <stdexcept>

static void require(bool ok,const char* message)
{
    if (!ok) throw runtime_error(message);
}
static void close(double actual,double expected,double tolerance,const char* message)
{
    require(std::isfinite(actual)&&fabs(actual-expected)<=tolerance,message);
}
struct Fixture
{
    PlaceData db;
    vector<unique_ptr<Module>> modules;
    vector<unique_ptr<Pin>> pins;
    vector<unique_ptr<Net>> nets;
    Fixture(float size=100)
    {
        db.coreRegion=CRect(POS_2D(0,0),POS_2D(size,size));
        row(0,size,0,size);
    }
    void row(float left,float right,float bottom,float height)
    {
        SiteRow r;
        r.start=POS_2D(left,bottom); r.end=POS_2D(right,bottom);
        r.bottom=bottom; r.height=height; r.step=1;
        db.SiteRows.push_back(r);
    }
    Module* module(float x,float y,float w,float h,bool macro=false,bool fixed=false)
    {
        modules.emplace_back(new Module(modules.size(),"test",w,h));
        Module* m=modules.back().get();
        m->ll=POS_2D(x-w/2,y-h/2); m->isMacro=macro; m->isFixed=fixed;
        (fixed?db.Terminals:db.Nodes).push_back(m);
        return m;
    }
    void connect(Module* a,Module* b)
    {
        nets.emplace_back(new Net());
        Net* n=nets.back().get(); db.Nets.push_back(n);
        for (Module* m:{a,b})
        {
            pins.emplace_back(new Pin());
            Pin* p=pins.back().get(); p->module=m; p->net=n; n->netPins.push_back(p);
        }
    }
};
static void initialize(MyPlacer& p)
{
    p.globalOptions.maxBinDimension=8;
    p.fillerInit(); p.binInit(); p.gradientVectorInitialization(); p.updateDensity(); p.updateDensityGradient();
}
static void testElectricField()
{
    int m=16;
    DCTSolver solver(m,2,3);
    vector<double> rho(m*m,0.8),ex,ey,phi;
    solver.solve(rho,ex,ey,phi);
    for (int i=0;i<m*m;++i)
    {
        close(ex[i],0,1e-12,"uniform density creates x field");
        close(ey[i],0,1e-12,"uniform density creates y field");
    }
    double wx=acos(-1.0)/(m*2);
    for (int y=0;y<m;++y)
        for (int x=0;x<m;++x) rho[y*m+x]=0.8+0.1*cos(acos(-1.0)*(x+0.5)/m);
    solver.solve(rho,ex,ey,phi);
    for (int y=0;y<m;++y)
        for (int x=0;x<m;++x)
        {
            close(phi[y*m+x],0.1*cos(acos(-1.0)*(x+0.5)/m)/(wx*wx),1e-10,"wrong potential spectrum");
            close(ex[y*m+x],0.1*sin(acos(-1.0)*(x+0.5)/m)/wx,1e-10,"wrong field frequency or sign");
            close(ey[y*m+x],0,1e-12,"spurious y field");
        }
    rho.assign(m*m,0);
    for (int y=7;y<=8;++y) for (int x=7;x<=8;++x) rho[y*m+x]=1;
    solver.solve(rho,ex,ey,phi);
    require(ex[8*m+4]<0&&ex[8*m+11]>0,"charge does not repel horizontally");
    require(ey[4*m+8]<0&&ey[11*m+8]>0,"charge does not repel vertically");
    for (int y=0;y<m;++y)
        for (int x=0;x<m;++x)
        {
            close(ex[y*m+x],-ex[y*m+m-1-x],1e-11,"asymmetric x field");
            close(ey[y*m+x],-ey[(m-1-y)*m+x],1e-11,"asymmetric y field");
        }
}
static void testRegionsChargeAndOverflow()
{
    Fixture f;
    f.db.SiteRows.clear(); f.row(0,100,0,40); f.row(0,60,0,40); // 重复覆盖不能重复计面积
    f.module(20,20,20,20,false,true);
    f.module(30,20,20,20,false,true); // 两个终端并集面积 600
    f.module(5,55,10,10,false,true); // 已在暗区，不能重复扣除
    f.module(50,20,2,2);
    f.module(70,20,10,10,true);
    MyPlacer p(&f.db); initialize(p);
    close(p.placeableArea,3400,1e-5,"incorrect region union");
    double real=0,filler=0,blocked=0;
    for (const auto& row:p.bins)
        for (const Bin_2D* b:row)
        {
            real+=b->nodeDensity; filler+=b->fillerDensity;
            blocked+=b->terminalDensity+b->DarkDensity;
        }
    close(real,4+p.targetDensity*100,1e-5,"macro charge or LengthScale incorrect");
    close(real+filler+blocked,p.targetDensity*10000,0.01,"effective charge not conserved");
    close(blocked,p.targetDensity*6600,0.01,"fixed regions double counted");
    size_t fillers=p.Fillers.size();
    p.fillerInit();
    require(p.Fillers.size()==fillers,"reinitialization duplicated fillers");
    p.binInit(); p.gradientVectorInitialization();
    for (auto& row:p.bins)
        for (Bin_2D* b:row) b->nodeDensity=b->fillerDensity=b->terminalDensity=b->DarkDensity=0;
    Bin_2D* b=p.bins[0][0];
    b->nodeDensity=p.targetDensity*b->area+100;
    b->fillerDensity=10000;
    close(p.calcOverflow(),100/(4+p.targetDensity*100),1e-6,"overflow includes filler or wrong macro normalization");
}
static void testDiscreteDensityGradient()
{
    Fixture f;
    f.module(31.3,41.7,9.1,6.3);
    f.module(58.4,65.2,13.4,15.2,true);
    f.module(10,12,4,4,false,true);
    MyPlacer p(&f.db);
    p.globalOptions.densityGradientMethod = DensityGradientMethod::DiscreteEnergy;
    initialize(p);
    if (!p.Fillers.empty()) p.Fillers[0]->ll=POS_2D(72.4,73.1);
    p.updateDensity(); p.updateDensityGradient();
    vector<size_t> indices{0,1};
    if (!p.Fillers.empty()) indices.push_back(f.db.Nodes.size());
    for (size_t index:indices)
        for (int axis=0;axis<2;++axis)
        {
            p.updateDensity(); p.updateDensityGradient();
            double analytic=axis==0?p.densityGradient[index].x:p.densityGradient[index].y;
            Module* m=p.NodesAndFillers[index];
            float& coordinate=axis==0?m->ll.x:m->ll.y;
            float original=coordinate;
            coordinate=original+0.002f; double plusX=coordinate;
            p.updateDensity(); p.updateDensityGradient(); double plus=p.densityEnergy;
            coordinate=original-0.002f; double minusX=coordinate;
            p.updateDensity(); p.updateDensityGradient(); double minus=p.densityEnergy;
            coordinate=original;
            double numerical=(plus-minus)/(plusX-minusX);
            close(analytic,numerical,0.002*max(1.0,fabs(numerical)),"density gradient disagrees with energy difference");
        }
}
static void testFieldDensityGradient()
{
    Fixture f;
    Module *standard = f.module(31.25, 43.75, 4, 4);
    Module *macro = f.module(56.25, 56.25, 20, 20, true);
    MyPlacer p(&f.db);
    initialize(p);
    require(p.globalOptions.densityGradientMethod == DensityGradientMethod::ElectricField,
            "reference electric-field gradient is not the default");
    require(!p.Fillers.empty(), "field charge test requires a filler");
    Module *filler = p.Fillers[0];
    double hx = p.binWidth, hy = p.binHeight;
    standard->ll = POS_2D(1.5 * hx - standard->width / 2, 1.5 * hy - standard->height / 2);
    macro->ll = POS_2D(2.5 * hx - macro->width / 2, 2.5 * hy - macro->height / 2);
    filler->ll = POS_2D(2.5 * hx - filler->width / 2, 1.5 * hy - filler->height / 2);

    // 制造已知密度模态。解析场为 Ex=0.1*sin(wx*x)/wx、Ey=0.07*sin(wy*y)/wy。
    // 采样位置选在 bin 中心，插值应准确返回谱求解得到的电场。
    double wx = acos(-1.0) / f.db.coreRegion.width();
    double wy = acos(-1.0) / f.db.coreRegion.height();
    for (auto &row : p.bins)
        for (Bin_2D *b : row)
        {
            double rho = 0.8 + 0.1 * cos(wx * b->center.x) + 0.07 * cos(wy * b->center.y);
            b->nodeDensity = rho * b->area;
            b->fillerDensity = b->terminalDensity = b->DarkDensity = 0;
        }
    p.updateDensityGradient();
    for (size_t index : {size_t(0), size_t(1), f.db.Nodes.size()})
    {
        Module *m = p.NodesAndFillers[index];
        double q = (double)m->width * m->height * (m->isMacro || m->isFiller ? p.targetDensity : 1.0);
        double expectedX = -q * 0.1 * sin(wx * m->center().x) / wx;
        double expectedY = -q * 0.07 * sin(wy * m->center().y) / wy;
        close(p.densityGradient[index].x, expectedX, 1e-6 * max(1.0, fabs(expectedX)), "wrong -q*Ex at cell center");
        close(p.densityGradient[index].y, expectedY, 1e-6 * max(1.0, fabs(expectedY)), "wrong -q*Ey at cell center");
    }
    require(macro->isMacro && filler->isFiller, "missing weighted charge cases");

    // 两格中心之间的中点必须取两边电场平均值，避免网格下标错半格。
    standard->ll.x = 2 * hx - standard->width / 2;
    p.updateDensityGradient();
    double midpointField = 0.5 * 0.1 / wx * (sin(wx * 1.5 * hx) + sin(wx * 2.5 * hx));
    close(p.densityGradient[0].x, -16 * midpointField, 1e-5, "incorrect field interpolation");

    // 支持框已靠墙时，合法单元中心仍应感受到非零电场，不能再乘 box.dx=0。
    standard->ll.x = 1; // 中心 x=3，小于半个 bin 的宽度
    p.updateDensityGradient();
    double edgeField = (3.0 / (hx / 2)) * 0.1 * sin(wx * hx / 2) / wx;
    close(p.densityGradient[0].x, -16 * edgeField, 1e-5, "near-boundary density force incorrectly vanished");
    require(p.densityGradient[0].x < -1, "near-boundary field has wrong sign");
    standard->ll.x = -standard->width / 2; // 仅检查场延拓：墙面法向场为零
    p.updateDensityGradient();
    close(p.densityGradient[0].x, 0, 1e-9, "Neumann wall field is not zero");
    require(fabs(p.densityGradient[0].y) > 1, "tangential field was incorrectly cancelled at wall");

    for (auto &row : p.bins)
        for (Bin_2D *b : row)
            b->nodeDensity = 0.8 * b->area;
    p.updateDensityGradient();
    for (const VECTOR_3D &g : p.densityGradient)
    {
        close(g.x, 0, 1e-8, "uniform density has nonzero cell x force");
        close(g.y, 0, 1e-8, "uniform density has nonzero cell y force");
    }
}
static void testReferenceParameters()
{
    Fixture f;
    MyPlacer p(&f.db);
    p.wirelengthGradient={VECTOR_3D()};
    p.wirelengthGradient[0].x=3; p.wirelengthGradient[0].y=4;
    p.densityGradient={VECTOR_3D(),VECTOR_3D()};
    p.densityGradient[0].x=1; p.densityGradient[0].y=1; p.densityGradient[1].x=2;
    p.penaltyFactorInitilization();
    close(p.lambda,1.75,1e-12,"unexpected lambda base weight");
    p.lambda=2; p.updatePenaltyFactor(350000); close(p.lambda,2,1e-12,"wrong lambda update");
    p.updatePenaltyFactor(0); close(p.lambda,2.2,1e-12,"wrong 1.1 feedback");
    p.binWidth=10; p.updateSmoothingParameter(0.1); close(p.gamma,8,1e-6,"wrong gamma formula");
    p.densityGradient.clear(); p.penaltyFactorInitilization();
    require(p.lambda>0&&std::isfinite(p.lambda),"zero density gradient permanently disables density force");
}
static void testWirelengthGradient()
{
    Fixture f;
    Module* a=f.module(21.3,33.1,4,4);
    Module* b=f.module(71.2,62.7,4,4);
    f.connect(a,b);
    MyPlacer p(&f.db);
    initialize(p); p.gamma=12.3f;
    for (int axis=0;axis<2;++axis)
    {
        p.updateWirelengthGradient();
        double analytic=axis==0?p.wirelengthGradient[0].x:p.wirelengthGradient[0].y;
        float& coordinate=axis==0?a->ll.x:a->ll.y;
        float original=coordinate;
        coordinate=original+0.002f; double plusX=coordinate;
        p.updateWirelengthGradient(); double plus=p.wirelengthValue;
        coordinate=original-0.002f; double minusX=coordinate;
        p.updateWirelengthGradient(); double minus=p.wirelengthValue;
        coordinate=original;
        close(analytic,(plus-minus)/(plusX-minusX),0.002,"WA gradient disagrees with finite difference");
    }
}
static void configureNAGFixture(Fixture& f)
{
    Module* movable=f.module(13.3,18.7,4,4);
    Module* fixed=f.module(35,32,2,2,false,true);
    f.connect(movable,fixed);
}
static void testLookaheadAndStopping()
{
    Fixture initial(40),one(40),two(40);
    configureNAGFixture(initial); configureNAGFixture(one); configureNAGFixture(two);
    MyPlacer start(&initial.db),first(&one.db),second(&two.db);
    GlobalPlacementOptions options;
    options.maxIterations=1; options.stopOverflow=0; options.maxBinDimension=8;
    start.globalOptions=options; start.fillerInit();
    vector<POS_2D> x0;
    for (Module* m:start.NodesAndFillers) x0.push_back(m->center());
    require(first.GlobalPlacement(options),"one-step NAG failed");
    vector<POS_2D> x1;
    for (Module* m:first.NodesAndFillers) x1.push_back(m->center());
    double lambda=first.lambda;
    float gamma=first.gamma;
    options.maxIterations=2;
    require(second.GlobalPlacement(options),"two-step NAG failed");
    require(second.globalStats.iterations==2,"iteration limit ignored");
    require(second.globalStats.stopReason.find("最大迭代")!=string::npos,"incorrect stop reason");
    // 独立重建第二轮的前瞻点，用一次运行后的下一轮 λ、γ 求梯度。
    for (size_t i=0;i<x0.size();++i)
    {
        Module* m=first.NodesAndFillers[i];
        double x=x1[i].x-options.momentum*(x0[i].x-x1[i].x);
        double y=x1[i].y-options.momentum*(x0[i].y-x1[i].y);
        x=max(m->width/2.0,min(40-m->width/2.0,x));
        y=max(m->height/2.0,min(40-m->height/2.0,y));
        m->ll=POS_2D(x-m->width/2,y-m->height/2);
    }
    first.lambda=lambda; first.gamma=gamma;
    first.updateDensity(); first.updateDensityGradient(); first.updateWirelengthGradient(); first.updateTotalGradient();
    double acceptedStep=second.stepSize/1.1;
    Module* m=first.NodesAndFillers[0];
    double dx=-acceptedStep*first.totalGradient[0].x/first.preconditioner[0];
    double dy=-acceptedStep*first.totalGradient[0].y/first.preconditioner[0];
    double length=hypot(dx,dy),limit=2.0*min(first.binWidth,first.binHeight);
    if (length>limit) { dx*=limit/length; dy*=limit/length; }
    double expectedX=max(2.0,min(38.0,m->center().x+dx));
    double expectedY=max(2.0,min(38.0,m->center().y+dy));
    close(two.db.Nodes[0]->center().x,expectedX,1e-4,"x gradient was not evaluated at lookahead");
    close(two.db.Nodes[0]->center().y,expectedY,1e-4,"y gradient was not evaluated at lookahead");
    for (Module* cell:second.NodesAndFillers)
        require(std::isfinite(cell->ll.x)&&std::isfinite(cell->ll.y)&&cell->ll.x>=0&&cell->ll.y>=0&&
                cell->ur().x<=40.001&&cell->ur().y<=40.001,"invalid NAG output");
}
int main()
{
    try
    {
        testElectricField();
        testRegionsChargeAndOverflow();
        testDiscreteDensityGradient();
        testFieldDensityGradient();
        testReferenceParameters();
        testWirelengthGradient();
        testLookaheadAndStopping();
        cout<<"Task 5 global placement checks passed"<<endl;
    }
    catch (const exception& e) { cerr<<e.what()<<endl; return 1; }
    return 0;
}
