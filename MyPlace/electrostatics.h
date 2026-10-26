#ifndef MYPLACEMENT_ELECTROSTATICS_H
#define MYPLACEMENT_ELECTROSTATICS_H

#include <vector>

// 正交 DCT-II：在 bin 中心采样，在物理坐标中求解 -Δphi=rho-mean(rho)。
// 零频势设为 0；ex/ey=-grad(phi)，energy=1/2*Σ rho*phi*binArea。
class DCTSolver
{
public:
    DCTSolver(int dimension, double binWidth, double binHeight);
    void solve(const std::vector<double> &rho, std::vector<double> &ex,
               std::vector<double> &ey, std::vector<double> &phi,
               bool computeField = true);
    double energy = 0;

private:
    int m;
    double hx, hy;
    std::vector<double> cosT, sinT, cn, tmp, a, spectral, work;
};

#endif
