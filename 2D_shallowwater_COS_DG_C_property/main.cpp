// #include <cmath>
// #include <iostream>
// #include <iomanip>
// #include <algorithm>
// #include <cstddef>

#include <iostream>
#include "constants.h"
#include "2D_shallowwater.h"
#include <array>
#include <fstream> // 保存.csv文件
#include <iomanip> // for setprecision
#include <typeinfo> // for typeid

int main() 
{
    static StateArray  uh{}, uh_init{};
    static BottomArray bh{};

    getBasic(phi, phiGLL, phi_xi, phi_eta);
    getInterface(phixL, phixR, phiyL, phiyR);

    getProjection(uh, bh);

    uh_init = uh;

    std::cout << "Grid: " << Nx << " x " << Ny
              << "  |  dimPk = " << dimPk
              << "  |  Real type = "
              << (sizeof(Real) == 4 ? "float" : "double") << "\n";

    runRK3(uh, bh, tend);

    // ---- 误差 ----
    Real L1_h = 0, L1_hu = 0, L1_hv = 0;
    Real Li_h = 0, Li_hu = 0, Li_hv = 0;

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Real dh  = uh[i][j][0][0] - uh_init[i][j][0][0];
            Real dhu = uh[i][j][0][1] - uh_init[i][j][0][1];
            Real dhv = uh[i][j][0][2] - uh_init[i][j][0][2];

            L1_h  += std::abs(dh);
            L1_hu += std::abs(dhu);
            L1_hv += std::abs(dhv);

            Li_h  = std::max(Li_h,  std::abs(dh));
            Li_hu = std::max(Li_hu, std::abs(dhu));
            Li_hv = std::max(Li_hv, std::abs(dhv));
        }
    }

    L1_h  *= hx * hy / ((xb - xa) * (yb - ya));
    L1_hu *= hx * hy / ((xb - xa) * (yb - ya));
    L1_hv *= hx * hy / ((xb - xa) * (yb - ya));

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "\n====  L1 error  ====\n";
    std::cout << "h  : " << L1_h  << "\n";
    std::cout << "hu : " << L1_hu << "\n";
    std::cout << "hv : " << L1_hv << "\n";
    // std::cout << "Size of h: " << typeid(L1_h).name() << " bytes\n";

    std::cout << "\n====  Linf error  ====\n";
    std::cout << "h  : " << Li_h  << "\n";
    std::cout << "hu : " << Li_hu << "\n";
    std::cout << "hv : " << Li_hv << "\n";

    // // ---- 保存单元平均值 h̄ = uh[i][j][0][0] 到 CSV ----
    // {
    //     std::ofstream fout("h_cell_average.csv");
    //     fout << "x,y,h\n";
    //     fout << std::scientific << std::setprecision(16);
    //     for (int i = 0; i < Nx; ++i) 
    //     {
    //         Real xc = xa + (i + 0.5) * hx;
    //         for (int j = 0; j < Ny; ++j) 
    //         {
    //             Real yc = ya + (j + 0.5) * hy;
    //             Real hbar = uh[i][j][0][0];   // 正交基下常数项 = 单元平均
    //             fout << xc << "," << yc << "," << hbar << "\n";
    //         }
    //     }
    //     fout.close();
    //     std::cout << "\nSaved h_cell_average.csv\n";
    // }

    return 0;
}