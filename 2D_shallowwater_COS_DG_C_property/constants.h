#pragma once

#include <array>
#include <cmath>
#include <typeinfo>

// =====================================================================
//  全局类型别名
// =====================================================================
using Real = double;   // 改成 float 即切换为单精度

template <typename T, std::size_t D1, std::size_t D2>
using Array2d = std::array<std::array<T, D2>, D1>;

template <typename T, std::size_t D1, std::size_t D2, std::size_t D3>
using Array3d = std::array<std::array<std::array<T, D3>, D2>, D1>;

template <typename T, std::size_t D1, std::size_t D2,
                     std::size_t D3, std::size_t D4>
using Array4d = std::array<
                    std::array<
                        std::array<
                            std::array<T, D4>, D3>, D2>, D1>;

// =====================================================================
//  常量命名空间
// =====================================================================
namespace Constants 
{
    constexpr int  Nx      { 100 };
    constexpr int  Ny      { 100 };
    constexpr int  dimPk   { 6 };
    constexpr int  NumEq   { 3 };

    constexpr int  NumGLP  { 5 };
    constexpr int  NumGLL   { 3 };

    constexpr Real g       { 9.812 };

    constexpr Real xa { 0.0 }, xb { 1.0 };
    constexpr Real ya { 0.0  }, yb { 1.0 };
    constexpr Real hx { (xb - xa) / static_cast<Real>(Nx) };
    constexpr Real hy { (yb - ya) / static_cast<Real>(Ny) };
    constexpr Real hx1 { 0.5 * hx };
    constexpr Real hy1 { 0.5 * hy };

    constexpr Real tend    { 0.1 };
    constexpr Real CFL     { 0.16 };
    // constexpr Real M_tvb   { 0.0 };
    constexpr Real dry_tol { 1e-5 };

    constexpr Real Ck { 0.2 }; // COS 过程的常数 Ck

    constexpr std::array<Real, dimPk> mm { 1.0, 1.0/3.0, 1.0/3.0, 1.0/9.0, 4.0/45.0, 4.0/45.0 };

    constexpr std::array<Real, NumGLP> gl_x { 
                                -0.9061798459386639927976269,
                                -0.5384693101056830910363144,
                                0.0,
                                0.5384693101056830910363144,
                                0.9061798459386639927976269 };
    constexpr std::array<Real, NumGLP> gl_w { 
                                0.2369268850561890875142640,
                                0.4786286704993664680412915,
                                0.5688888888888888888888889,
                                0.4786286704993664680412915,
                                0.2369268850561890875142640 };

    constexpr std::array<Real, NumGLL> gll_x { -1.0, 0.0, 1.0 };
    constexpr std::array<Real, NumGLL> gll_w { 1.0/3.0, 4.0/3.0, 1.0/3.0 };

    inline Array3d<Real, dimPk, NumGLP, NumGLP> phi {};    // Guass point
    inline Array3d<Real, dimPk, NumGLL, NumGLL> phiGLL {}; // Guass-Lobatto point
    inline Array3d<Real, dimPk, NumGLP, NumGLP> phi_xi {};
    inline Array3d<Real, dimPk, NumGLP, NumGLP> phi_eta {};

    inline Array2d<Real, dimPk, NumGLP> phixL {};
    inline Array2d<Real, dimPk, NumGLP> phixR {};

    inline Array2d<Real, dimPk, NumGLP> phiyL {};
    inline Array2d<Real, dimPk, NumGLP> phiyR {};

    
}

using namespace Constants;

// =====================================================================
//  数据别名
// =====================================================================
using StateArray  = Array4d<Real, Nx, Ny, dimPk, NumEq>;
using BottomArray = Array3d<Real, Nx, Ny, dimPk>;