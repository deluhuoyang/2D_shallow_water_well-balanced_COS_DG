#pragma once

#include <array>
#include <cmath>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace Constants
{
    constexpr int  Nx      { 128 };
    constexpr int  Ny      { 128 };
    constexpr int  dimPk   { 6 };
    constexpr int  NumEq   { 3 };

    constexpr int  NumGLP  { 5 };
    constexpr int  NumGLL  { 3 };

    constexpr double g       { 9.812 };

    constexpr double xa { 0.0 }, xb { 1.0 };
    constexpr double ya { 0.0 }, yb { 1.0 };
    constexpr double hx { (xb - xa) / static_cast<double>(Nx) };
    constexpr double hy { (yb - ya) / static_cast<double>(Ny) };
    constexpr double hx1 { 0.5 * hx };
    constexpr double hy1 { 0.5 * hy };
    
    constexpr double tend    { 0.05 };
    constexpr double CFL     { 0.16 };
    constexpr double dry_tol { 1e-5 };

    constexpr double Ck { 0.2 };

    constexpr std::array<double, dimPk> mm { 1.0, 1.0/3.0, 1.0/3.0, 1.0/9.0, 4.0/45.0, 4.0/45.0 };

    constexpr std::array<double, NumGLP> gl_x {
        -0.9061798459386639927976269,
        -0.5384693101056830910363144,
         0.0,
         0.5384693101056830910363144,
         0.9061798459386639927976269 };
    constexpr std::array<double, NumGLP> gl_w {
        0.2369268850561890875142640,
        0.4786286704993664680412915,
        0.5688888888888888888888889,
        0.4786286704993664680412915,
        0.2369268850561890875142640 };

    constexpr std::array<double, NumGLL> gll_x { -1.0, 0.0, 1.0 };
    constexpr std::array<double, NumGLL> gll_w { 1.0/3.0, 4.0/3.0, 1.0/3.0 };

    inline std::vector<double> phi;
    inline std::vector<double> phiGLL;
    inline std::vector<double> phi_xi;
    inline std::vector<double> phi_eta;

    inline std::vector<double> phixL;
    inline std::vector<double> phixR;
    inline std::vector<double> phiyL;
    inline std::vector<double> phiyR;
}

using namespace Constants;

// =====================================================================
//  一维索引函数
// =====================================================================
inline std::size_t idx_state(std::size_t i, std::size_t j, std::size_t d, std::size_t n) {
    return ((i * static_cast<std::size_t>(Ny) + j) * static_cast<std::size_t>(dimPk) + d)
           * static_cast<std::size_t>(NumEq) + n;
}
inline std::size_t idx_bottom(std::size_t i, std::size_t j, std::size_t d) {
    return (i * static_cast<std::size_t>(Ny) + j) * static_cast<std::size_t>(dimPk) + d;
}
inline std::size_t idx_phi(std::size_t d, std::size_t p, std::size_t q) {
    return (d * static_cast<std::size_t>(NumGLP) + p) * static_cast<std::size_t>(NumGLP) + q;
}
inline std::size_t idx_phiGLL(std::size_t d, std::size_t p, std::size_t q) {
    return (d * static_cast<std::size_t>(NumGLL) + p) * static_cast<std::size_t>(NumGLL) + q;
}
inline std::size_t idx_phi1d(std::size_t d, std::size_t q) {
    return d * static_cast<std::size_t>(NumGLP) + q;
}

// =====================================================================
//  数据别名
// =====================================================================
using StateArray  = std::vector<double>;
using BottomArray = std::vector<double>;