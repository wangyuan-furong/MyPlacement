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
static void testDensityGradient()
{
    Fixture f;
    f.module(31.3,41.7,9.1,6.3);
    f.module(58.4,65.2,13.4,15.2,true);
    f.module(10,12,4,4,false,true);
    MyPlacer p(&f.db); initialize(p);
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
        testDensityGradient();
        testReferenceParameters();
        testLookaheadAndStopping();
        cout<<"Task 5 global placement checks passed"<<endl;
    }
    catch (const exception& e) { cerr<<e.what()<<endl; return 1; }
    return 0;
}
