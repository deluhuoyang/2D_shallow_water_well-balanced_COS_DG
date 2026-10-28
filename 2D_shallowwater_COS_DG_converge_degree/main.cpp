#include <iostream>
#include "constants.h"
#include "2D_shallowwater.h"
#include <fstream>
#include <iomanip>
#include <string>

int main()
{
    StateArray  uh;   // 空 vector，getProjection 内部 assign 分配
    BottomArray bh;

    getBasic(phi, phiGLL, phi_xi, phi_eta);
    getInterface(phixL, phixR, phiyL, phiyR);

    getProjection(uh, bh);

    std::cout << "Grid: " << Nx << " x " << Ny
              << "  |  dimPk = " << dimPk
              << "  |  Real type = double\n";

    runRK3(uh, bh, tend);

    auto save = [&](const std::string& tag, std::size_t comp)
    {
        std::string fn = tag + "_cell_" + std::to_string(Nx) + "x" + std::to_string(Ny) + ".csv";
        std::ofstream fout(fn);
        fout << "x,y," << tag << "\n";
        fout << std::scientific << std::setprecision(16);
        for (int i = 0; i < Nx; ++i)
        {
            double xc = xa + (i + 0.5) * hx;
            for (int j = 0; j < Ny; ++j)
            {
                double yc = ya + (j + 0.5) * hy;
                fout << xc << "," << yc << ","
                     << uh[idx_state(i,j,0,comp)] << "\n";
            }
        }
        fout.close();
        std::cout << fn << " saved successfully.\n";
    };

    save("h",  0);
    save("hu", 1);
    save("hv", 2);

    return 0;
}