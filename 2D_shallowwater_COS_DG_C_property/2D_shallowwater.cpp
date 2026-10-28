// #include <cmath>

// #include <iostream>
// #include <iomanip>
// #include <algorithm>
// #include <cstddef>


#include "constants.h"
#include "2D_shallowwater.h"
#include "COS_produre.h"
#include <array>
#include <iostream>

// =====================================================================
//  初始资料与flux
// =====================================================================
Real b(Real x, Real y) 
{
    return 0.8 * std::exp(-50.0 * ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)));
}
Real h (Real x, Real y) { return 1.0 - b(x, y); }
Real hu(Real /*x*/,   Real /*y*/)   { return 0.0; }
Real hv(Real /*x*/,   Real /*y*/)   { return 0.0; }

Real f1(Real h, Real hu, Real hv) { return hu; }
Real f2(Real h, Real hu, Real hv) { return hu*hu/h + 0.5 * g * h * h; }
Real f3(Real h, Real hu, Real hv) { return hu*hv/h; }
Real g1(Real h, Real hu, Real hv) { return hv; }
Real g2(Real h, Real hu, Real hv) { return hu*hv/h; }
Real g3(Real h, Real hu, Real hv) { return hv*hv/h + 0.5 * g * h * h; }

// =====================================================================
//  获取高斯点在基函数的值、获取高斯点在界面函数的值
// =====================================================================
void getBasic(Array3d<Real, dimPk, NumGLP, NumGLP>& phi,
              Array3d<Real, dimPk, NumGLL, NumGLL>& phiGLL,  
              Array3d<Real, dimPk, NumGLP, NumGLP>& phi_xi,
              Array3d<Real, dimPk, NumGLP, NumGLP>& phi_eta)
{
    // φ0=1, φ1=ξ, φ2=η, φ3=ξη, φ4=ξ²-1/3, φ5=η²-1/3
    for (std::size_t p { 0 }; p < NumGLP; ++p)
    {
        for (std::size_t q { 0 }; q < NumGLP; ++q)
        {
            phi[0][p][q] = 1.0;
            phi[1][p][q] = gl_x[p];
            phi[2][p][q] = gl_x[q];
            phi[3][p][q] = gl_x[p] * gl_x[q];
            phi[4][p][q] = gl_x[p] * gl_x[p] - 1.0/3.0;
            phi[5][p][q] = gl_x[q] * gl_x[q] - 1.0/3.0;

            phi_xi[0][p][q] = 0.0;
            phi_xi[1][p][q] = 1.0;
            phi_xi[2][p][q] = 0.0;
            phi_xi[3][p][q] = gl_x[q];
            phi_xi[4][p][q] = 2 * gl_x[p];
            phi_xi[5][p][q] = 0.0;

            phi_eta[0][p][q] = 0.0;
            phi_eta[1][p][q] = 0.0;
            phi_eta[2][p][q] = 1.0;
            phi_eta[3][p][q] = gl_x[p];
            phi_eta[4][p][q] = 0.0;
            phi_eta[5][p][q] = 2 * gl_x[q];
        }
    }

    for (std::size_t p { 0 }; p < NumGLL; ++p)
    {
        for (std::size_t q { 0 }; q < NumGLL; ++q)
        {
            phiGLL[0][p][q] = 1.0;
            phiGLL[1][p][q] = gll_x[p];
            phiGLL[2][p][q] = gll_x[q];
            phiGLL[3][p][q] = gll_x[p] * gll_x[q];
            phiGLL[4][p][q] = gll_x[p] * gll_x[p] - 1.0/3.0;
            phiGLL[5][p][q] = gll_x[q] * gll_x[q] - 1.0/3.0;
        }
    }
}

void getInterface(Array2d<Real, dimPk, NumGLP>& phixL,
                  Array2d<Real, dimPk, NumGLP>& phixR,
                  Array2d<Real, dimPk, NumGLP>& phiyL,
                  Array2d<Real, dimPk, NumGLP>& phiyR)
{
    for (std::size_t q { 0 }; q < NumGLP; ++q)
    {
        phixR[0][q] = 1.0;
        phixR[1][q] = 1.0;
        phixR[2][q] = gl_x[q];
        phixR[3][q] = gl_x[q];
        phixR[4][q] = 2.0/3.0;
        phixR[5][q] = gl_x[q] * gl_x[q] - 1.0/3.0;

        phixL[0][q] = 1.0;
        phixL[1][q] = -1.0;
        phixL[2][q] = gl_x[q];
        phixL[3][q] = -gl_x[q];
        phixL[4][q] = 2.0/3.0;
        phixL[5][q] = gl_x[q] * gl_x[q] - 1.0/3.0;
    }

    for (std::size_t p { 0 }; p < NumGLP; ++p)
    {
        
        phiyR[0][p] = 1.0;
        phiyR[1][p] = gl_x[p];
        phiyR[2][p] = 1.0;
        phiyR[3][p] = gl_x[p];
        phiyR[4][p] = gl_x[p] * gl_x[p] - 1.0/3.0;
        phiyR[5][p] = 2.0/3.0;

        phiyL[0][p] = 1.0;
        phiyL[1][p] = gl_x[p];
        phiyL[2][p] = -1.0;
        phiyL[3][p] = -gl_x[p];
        phiyL[4][p] = gl_x[p] * gl_x[p] - 1.0/3.0;
        phiyL[5][p] = 2.0/3.0;
    }
}            

// =====================================================================
//  L² 投影初始化
// =====================================================================
void getProjection(StateArray& uh, BottomArray& bh)
{
    for (auto& a : uh) for (auto& row : a) for (auto& v : row) v.fill(Real(0));
    for (auto& a : bh) for (auto& v : a) v.fill(Real(0));

    for (std::size_t i { 0 }; i < Nx; ++i) 
    {
        Real x_i { xa + static_cast<Real>(i + 0.5) * hx };

        for (std::size_t j { 0 }; j < Ny; ++j) 
        {
            Real y_j { ya + static_cast<Real>(j + 0.5) * hy };

            for (std::size_t d { 0 }; d < dimPk; ++d) 
            {
                // 计算单元I_ij的系数c_d
                for (std::size_t p { 0 }; p < NumGLP; ++p)
                {
                    Real xi_p { gl_x[p] }; 
                    Real w_p { gl_w[p] };

                    for (std::size_t q { 0 }; q < NumGLP; ++q)
                    {
                        Real eta_q { gl_x[q] }; 
                        Real w_q { gl_w[q] };

                        uh[i][j][d][0] += w_p * w_q * h(hx1 * xi_p + x_i, hy1 * eta_q + y_j) * phi[d][p][q]; 
                        bh[i][j][d] += w_p * w_q * b(hx1 * xi_p + x_i, hy1 * eta_q + y_j) * phi[d][p][q]; 
                    }
                }

                uh[i][j][d][0] /= (4 * mm[d]);
                bh[i][j][d] /= (4 * mm[d]);
            }
        }
    }
}

// =====================================================================
//  RK3
// =====================================================================
void runRK3(StateArray& uh, const BottomArray& bh, Real t_end)
{
    Real t = 0.0;
    int step = 0;

    StateArray wh{}; // 后续用于COS过程，w = h + b

    while (t < t_end) 
    {
        // 1. calculate alpha and dt
        Real alpha_1 = 1e-15, alpha_2 = 1e-15;
        for (std::size_t i { 0 }; i < Nx; ++i) 
        {
            for (std::size_t j { 0 }; j < Ny; ++j) 
            {
                Real h { uh[i][j][0][0] };
                Real hu { uh[i][j][0][1] };
                Real hv { uh[i][j][0][2] };

                Real u_bar { hu / h };
                Real v_bar { hv / h };

                alpha_1 = alpha_1 > std::abs(u_bar) + std::sqrt(g * h) ? 
                            alpha_1 : std::abs(u_bar) + std::sqrt(g * h);
                alpha_2 = alpha_2 > std::abs(v_bar) + std::sqrt(g * h) ? 
                            alpha_2 : std::abs(v_bar) + std::sqrt(g * h);
            }
        }

        Real dt { CFL / (alpha_1 / hx + alpha_2 / hy) };

        // 2. 时间步长截断，保证到达tend
        if (t + dt >= t_end)
        {
            dt = t_end - t; 
            t = t_end;
        }
        else
        {
            t += dt;
        }

        // 3. RK3
        static StateArray u1, u2, du;
        // ---- Stage 1 ----
        Lh(uh, bh, du, alpha_1, alpha_2); // 计算du
        for (std::size_t i { 0 }; i < Nx; ++i)
            for (std::size_t j { 0 }; j < Ny; ++j)
                for (std::size_t d { 0 }; d < dimPk; ++d)
                    for (std::size_t n { 0 }; n < NumEq; ++n)
                        u1[i][j][d][n] = uh[i][j][d][n] + dt * du[i][j][d][n];  
        
        positivityPreservingLimiter(uh);                

        // ---- Stage 2 ----
        Lh(u1, bh, du, alpha_1, alpha_2);
        for (std::size_t i { 0 }; i < Nx; ++i)
            for (std::size_t j { 0 }; j < Ny; ++j)
                for (std::size_t d { 0 }; d < dimPk; ++d)
                    for (std::size_t n { 0 }; n < NumEq; ++n)
                        u2[i][j][d][n] = 
                            0.75 * uh[i][j][d][n] + 
                            0.25 * u1[i][j][d][n] + 
                            0.25 * dt * du[i][j][d][n];

        positivityPreservingLimiter(uh); 

        // ---- Stage 3 ----
        Lh(u2, bh, du, alpha_1, alpha_2);
        for (std::size_t i { 0 }; i < Nx; ++i)
            for (std::size_t j { 0 }; j < Ny; ++j)
                for (std::size_t d { 0 }; d < dimPk; ++d)
                    for (std::size_t n { 0 }; n < NumEq; ++n)
                        uh[i][j][d][n] = 
                            (1.0 / 3.0) * uh[i][j][d][n] + 
                            (2.0 / 3.0) * u2[i][j][d][n] + 
                            (2.0 / 3.0) * dt * du[i][j][d][n];

        positivityPreservingLimiter(uh);                     

        // 4. COS 过程                    
        for (std::size_t i { 0 }; i < Nx; ++i)
            for (std::size_t j { 0 }; j < Ny; ++j)
                for (std::size_t d { 0 }; d < dimPk; ++d)
                {
                    wh[i][j][d][0] = uh[i][j][d][0] + bh[i][j][d];  // w = h + b
                    wh[i][j][d][1] = uh[i][j][d][1];                // hu 不变
                    wh[i][j][d][2] = uh[i][j][d][2];                // hv 不变
                }
                        
        COS_produre(wh, bh, dt);

        for (std::size_t i { 0 }; i < Nx; ++i)
            for (std::size_t j { 0 }; j < Ny; ++j)
                for (std::size_t d { 0 }; d < dimPk; ++d)
                {
                    uh[i][j][d][0] = wh[i][j][d][0] - bh[i][j][d];  // h = w - b
                    uh[i][j][d][1] = wh[i][j][d][1];
                    uh[i][j][d][2] = wh[i][j][d][2];
                }
                        
        // 调试信息
        if (step % 10 == 0)
        {
            std::cout << "Step: " << step << ", t = " << t << ", dt = " << dt << "\n"; 
        }            

        ++step;
        // if (step == 1) break;               
    }
}

// =====================================================================
//  空间算子 Lh —— 结果写入 du
// =====================================================================
void Lh(const StateArray&  uh,
        const BottomArray& bh,
        StateArray&        du,
        Real               alpha_1,
        Real               alpha_2)
{
    // 清零
    for (auto& a : du) for (auto& b1 : a) for (auto& v : b1) v.fill(Real(0));

    // ---------- Step 1: 体积分 和 源项 ----------
    for (std::size_t i { 0 }; i < Nx; ++i)
    {
        for (std::size_t j { 0 }; j < Ny; ++j)
        {
            for (std::size_t d { 0 }; d < dimPk; ++d)
            {
                // 计算单元I_ij的系数c_d
                for (std::size_t p { 0 }; p < NumGLP; ++p)
                {
                    // Real xi_p { gl_x[p] }; 
                    Real w_p { gl_w[p] };

                    for (std::size_t q { 0 }; q < NumGLP; ++q)
                    {
                        // Real eta_q { gl_x[q] }; 
                        Real w_q { gl_w[q] };

                        Real h{}, hu{}, hv{};
                        calculateUh(uh, i, j, p, q, h, hu, hv);

                        // 1. 体积分
                        du[i][j][d][0] += w_p * w_q * (
                                            hy1 * f1(h, hu, hv) * phi_xi[d][p][q] + 
                                            hx1 * g1(h, hu, hv) * phi_eta[d][p][q]
                                            );
                        
                        du[i][j][d][1] += w_p * w_q * (
                                            hy1 * f2(h, hu, hv) * phi_xi[d][p][q] + 
                                            hx1 * g2(h, hu, hv) * phi_eta[d][p][q]
                                            );    

                        du[i][j][d][2] += w_p * w_q * (
                                            hy1 * f3(h, hu, hv) * phi_xi[d][p][q] + 
                                            hx1 * g3(h, hu, hv) * phi_eta[d][p][q]
                                            ); 
                                            
                        // 2. 源项
                        Real b_xi{}, b_eta{};
                        calculateBhDeriv(bh, i, j, p, q, b_xi, b_eta);

                        du[i][j][d][1] += hy1 * (-g) * w_p * w_q * h * b_xi * phi[d][p][q];
                        du[i][j][d][2] += hx1 * (-g) * w_p * w_q * h * b_eta * phi[d][p][q];                   
                    }
                }
            }
        }
    }

    // ---------- Step 2: x 方向界面 ----------
    for (std::size_t i { 0 }; i < Nx; ++i)
    {
        for (std::size_t j { 0 }; j < Ny; ++j)
        {
            for (std::size_t d { 0 }; d < dimPk; ++d)
            {
                for (std::size_t q { 0 }; q < NumGLP; ++q)
                {
                    Real w_q { gl_w[q] };

                    // 分别计算单元两侧的左值、右值
                    std::array<Real, NumEq> uR_left{}, uR_right{};
                    std::array<Real, NumEq> uL_left{}, uL_right{};

                    for (std::size_t n { 0 }; n < NumEq; ++n)
                    {
                        calculateUxInterface(uh, i, j, n, q, uL_left[n], uL_right[n], uR_left[n], uR_right[n]);
                    }

                    // b 在单元I_ij沿x方向两侧的左值、右值
                    Real bL_left{0.0}, bL_right{0.0}, bR_left{0.0}, bR_right{0.0};

                    // 计算单元 I_ij 左右相邻的2个单元的索引，周期边界
                    std::size_t i_left{ (i == 0 ? Nx-1 : i-1) }, i_right{ (i == Nx-1 ? 0 : i+1) }; 
                    std::size_t j_left{ j }, j_right{ j };
                    for (int d = 0; d < dimPk; ++d) 
                    {
                        bL_left += bh[i_left][j_left][d] * phixR[d][q];
                        bL_right += bh[i][j][d] * phixL[d][q];

                        bR_left += bh[i][j][d] * phixR[d][q];
                        bR_right += bh[i_right][j_right][d] * phixL[d][q];
                    }

                    // x 方向的静水重构
                    // 1. 单元左侧的左值、右值
                    Real hL_left { uL_left[0] }, hL_right { uL_right[0] };
                    Real hL_left_state { std::max(static_cast<Real>(0.0), hL_left + bL_left - std::max(bL_left, bL_right)) };
                    Real hL_right_state { std::max(static_cast<Real>(0.0), hL_right + bL_right - std::max(bL_left, bL_right)) };

                    uL_left[2] = hL_left_state * (uL_left[2] / uL_left[0]);
                    uL_left[1] = hL_left_state * (uL_left[1] / uL_left[0]);
                    uL_left[0] = hL_left_state;

                    uL_right[2] = hL_right_state * (uL_right[2] / uL_right[0]);
                    uL_right[1] = hL_right_state * (uL_right[1] / uL_right[0]);
                    uL_right[0] = hL_right_state;

                    // 2. 单元右侧的左值、右值
                    Real hR_left { uR_left[0] }, hR_right { uR_right[0] };
                    Real hR_left_state { std::max(static_cast<Real>(0.0), hR_left + bR_left - std::max(bR_left, bR_right)) };
                    Real hR_right_state { std::max(static_cast<Real>(0.0), hR_right + bR_right - std::max(bR_left, bR_right)) };

                    uR_left[2] = hR_left_state * (uR_left[2] / uR_left[0]);
                    uR_left[1] = hR_left_state * (uR_left[1] / uR_left[0]);
                    uR_left[0] = hR_left_state; 

                    uR_right[2] = hR_right_state * (uR_right[2] / uR_right[0]);
                    uR_right[1] = hR_right_state * (uR_right[1] / uR_right[0]);
                    uR_right[0] = hR_right_state;

                    std::array<Real, NumEq> fhatL_state {}, fhatR_state {};
                    for (std::size_t n { 0 }; n < NumEq; ++n)
                    {
                        fhatL_state[n] = LF_Flux_f(n, uL_left, uL_right, alpha_1);
                        fhatR_state[n] = LF_Flux_f(n, uR_left, uR_right, alpha_1);

                        if (n == 1) // 只有 x 方向动量方程 (hu) 需要静水修正
                        {
                            fhatL_state[n] += 0.5 * g * (hL_right * hL_right - hL_right_state * hL_right_state);
                            fhatR_state[n] += 0.5 * g * (hR_left * hR_left - hR_left_state * hR_left_state);
                        }
                    }
                    
                    for (std::size_t n { 0 }; n < NumEq; ++n)
                    {
                        du[i][j][d][n] -= hy1 * w_q * (fhatR_state[n] * phixR[d][q] - fhatL_state[n] * phixL[d][q]);
                    }  
                }                    
            }
        }
    }

    // ---------- Step 3: y 方向界面 ----------
    for (std::size_t i { 0 }; i < Nx; ++i)
    {
        for (std::size_t j { 0 }; j < Ny; ++j)
        {
            for (std::size_t d { 0 }; d < dimPk; ++d)
            {
                for (std::size_t p { 0 }; p < NumGLP; ++p)
                {
                    Real w_p { gl_w[p] };

                    // 分别计算单元两侧的值（D代表下侧 Down，U代表上侧 Up）
                    std::array<Real, NumEq> uU_down{}, uU_up{};
                    std::array<Real, NumEq> uD_down{}, uD_up{};

                    for (std::size_t n { 0 }; n < NumEq; ++n)
                    {
                        calculateUyInterface(uh, i, j, n, p, uD_down[n], uD_up[n], uU_down[n], uU_up[n]);
                    }

                    // b 在单元I_ij沿y方向两侧的值
                    Real bD_down{0.0}, bD_up{0.0}, bU_down{0.0}, bU_up{0.0};

                    std::size_t i_down{ i }, i_up{ i };
                    std::size_t j_down{ (j == 0 ? Ny-1 : j-1) }, j_up{ (j == Ny-1 ? 0 : j+1) };
                    for (int d = 0; d < dimPk; ++d) 
                    {
                        bD_down += bh[i_down][j_down][d] * phiyR[d][p];
                        bD_up += bh[i][j][d] * phiyL[d][p];

                        bU_down += bh[i][j][d] * phiyR[d][p];
                        bU_up += bh[i_up][j_up][d] * phiyL[d][p];
                    }

                    // y 方向的静水重构
                    // 1. 单元左侧的左值、右值
                    Real hD_down { uD_down[0] }, hD_up { uD_up[0] };
                    Real hD_down_state { std::max(static_cast<Real>(0.0), hD_down + bD_down - std::max(bD_down, bD_up)) };
                    Real hD_up_state { std::max(static_cast<Real>(0.0), hD_up + bD_up - std::max(bD_down, bD_up)) };

                    uD_down[2] = hD_down_state * (uD_down[2] / uD_down[0]);
                    uD_down[1] = hD_down_state * (uD_down[1] / uD_down[0]);
                    uD_down[0] = hD_down_state;

                    uD_up[2] = hD_up_state * (uD_up[2] / uD_up[0]);
                    uD_up[1] = hD_up_state * (uD_up[1] / uD_up[0]);
                    uD_up[0] = hD_up_state;

                    // 2. 上界面 (U) 的静水重构
                    Real hU_down { uU_down[0] }, hU_up { uU_up[0] };
                    Real hU_down_state { std::max(static_cast<Real>(0.0), hU_down + bU_down - std::max(bU_down, bU_up)) };
                    Real hU_up_state { std::max(static_cast<Real>(0.0), hU_up + bU_up - std::max(bU_down, bU_up)) };

                    uU_down[2] = hU_down_state * (uU_down[2] / uU_down[0]);
                    uU_down[1] = hU_down_state * (uU_down[1] / uU_down[0]);
                    uU_down[0] = hU_down_state; 

                    uU_up[2] = hU_up_state * (uU_up[2] / uU_up[0]);
                    uU_up[1] = hU_up_state * (uU_up[1] / uU_up[0]);
                    uU_up[0] = hU_up_state;

                    std::array<Real, NumEq> ghatL_state {}, ghatR_state {};
                    for (std::size_t n { 0 }; n < NumEq; ++n)
                    {
                        // 下界面通量
                        ghatL_state[n] = LF_Flux_g(n, uD_down, uD_up, alpha_2);
                        // 上界面通量
                        ghatR_state[n] = LF_Flux_g(n, uU_down, uU_up, alpha_2);

                        if (n == 2) // 只有 y 方向动量方程 (hv) 需要静水修正
                        {
                            ghatL_state[n] += 0.5 * g * (hD_up * hD_up - hD_up_state * hD_up_state);
                            ghatR_state[n] += 0.5 * g * (hU_down * hU_down - hU_down_state * hU_down_state);
                        }
                    }

                    for (std::size_t n { 0 }; n < NumEq; ++n)
                    {
                        du[i][j][d][n] -= hx1 * w_p * (ghatR_state[n] * phiyR[d][p] - ghatL_state[n] * phiyL[d][p]);
                    }
                }
            }
        }
    }

    // ---------- Step 5: 除以 hx·hy·mm[d]----------
    for (std::size_t i { 0 }; i < Nx; ++i)
    {
        for (std::size_t j { 0 }; j < Ny; ++j)
        {
            for (std::size_t d { 0 }; d < dimPk; ++d)
            {
                for (std::size_t n { 0 }; n < NumEq; ++n)
                {
                    du[i][j][d][n] /= hx * hy * mm[d];
                }
            }
        }
    }
}

// =====================================================================
//  Positivity-preserving limiter
// =====================================================================
void positivityPreservingLimiter(StateArray& uh)
{
    for (std::size_t i { 0 }; i < Nx; ++i)
    {
        for (std::size_t j { 0 }; j < Ny; ++j)
        {
            // 1. 计算 m_ij
            Real m_ij { 0.0 };
            for (std::size_t d { 0 }; d < dimPk; ++d) 
            {
                m_ij += uh[i][j][d][0] * phiGLL[d][0][0];
            }

            // 公式 m_ij = min { h_ij }
            for (std::size_t p { 0 }; p < NumGLL; ++p)
            {
                for (std::size_t q { 0 }; q < NumGLL; ++q)
                {
                    // 计算 h_ij 在 Guass-Lobatto 点的值
                    Real h_ij { 0.0 };
                    for (std::size_t d { 0 }; d < dimPk; ++d) 
                    {
                        h_ij += uh[i][j][d][0] * phiGLL[d][p][q];
                    }

                    if (h_ij < m_ij) m_ij = h_ij;
                }
            }

            // 2. 计算 theta
            Real theta{};
            Real h_bar = uh[i][j][0][0];
            theta = std::min(static_cast<Real>(1.0), h_bar / (h_bar - m_ij));

            // 3. 线性放缩
            for (std::size_t d { 1 }; d < dimPk; ++d)
                for (std::size_t n { 0 }; n < NumEq; ++n)
                        uh[i][j][d][n] = theta * uh[i][j][d][n];
        }
    }              
}

// =====================================================================
//  Lax–Friedrichs flux函数
// =====================================================================
Real LF_Flux_f(std::size_t n, std::array<Real, NumEq>& uR_left, 
              std::array<Real, NumEq>& uR_right, Real alpha)
{
    switch (n)
    {
        case 0: return 0.5 * (f1(uR_right[0], uR_right[1], uR_right[2]) + 
                            f1(uR_left[0], uR_left[1], uR_left[2]))
                            - 0.5 * alpha * (uR_right[0] - uR_left[0]);

        case 1: return 0.5 * (f2(uR_right[0], uR_right[1], uR_right[2]) + 
                            f2(uR_left[0], uR_left[1], uR_left[2])) 
                            - 0.5 * alpha * (uR_right[1] - uR_left[1]);

        case 2: return 0.5 * (f3(uR_right[0], uR_right[1], uR_right[2]) + 
                            f3(uR_left[0], uR_left[1], uR_left[2]))
                            - 0.5 * alpha * (uR_right[2] - uR_left[2]);                

        default: return 666;
    }
}

Real LF_Flux_g(std::size_t n, std::array<Real, NumEq>& uR_left, 
              std::array<Real, NumEq>& uR_right, Real alpha)
{
    switch (n)
    {
        case 0: return 0.5 * (g1(uR_right[0], uR_right[1], uR_right[2]) + 
                            g1(uR_left[0], uR_left[1], uR_left[2]))
                            - 0.5 * alpha * (uR_right[0] - uR_left[0]);

        case 1: return 0.5 * (g2(uR_right[0], uR_right[1], uR_right[2]) + 
                            g2(uR_left[0], uR_left[1], uR_left[2])) 
                            - 0.5 * alpha * (uR_right[1] - uR_left[1]);
                        
        case 2: return 0.5 * (g3(uR_right[0], uR_right[1], uR_right[2]) + 
                            g3(uR_left[0], uR_left[1], uR_left[2]))
                            - 0.5 * alpha * (uR_right[2] - uR_left[2]);                

        default: return 666;
    }
}

// =====================================================================
//  求值辅助: h, hu, hv, b_xi, b_eta, uh的界面值
// =====================================================================
void calculateUh(const StateArray& uh, std::size_t i, std::size_t j,
                       std::size_t p, std::size_t q, Real& h, Real& hu, Real& hv)
{
    h = hu = hv = 0.0;
    for (std::size_t d { 0 }; d < dimPk; ++d) 
    {
        h  += uh[i][j][d][0] * phi[d][p][q];
        hu += uh[i][j][d][1] * phi[d][p][q];
        hv += uh[i][j][d][2] * phi[d][p][q];
    }
}

void calculateBhDeriv(const BottomArray& bh, std::size_t i, std::size_t j,
                       std::size_t p, std::size_t q, Real& b_xi, Real& b_eta)
{
    b_xi = b_eta = 0.0;
    for (int d = 0; d < dimPk; ++d) 
    {
        b_xi += bh[i][j][d] * phi_xi[d][p][q];
        b_eta += bh[i][j][d] * phi_eta[d][p][q];
    }
}

void calculateUxInterface(const StateArray& uh, std::size_t i, std::size_t j,
                        std::size_t n, std::size_t q, Real& uL_left, Real& uL_right, 
                        Real& uR_left, Real& uR_right)
{
    uL_left = uL_right = uR_left = uR_right = 0.0;

    // 计算单元 I_ij 左右相邻的2个单元的索引，周期边界
    std::size_t i_left{ (i == 0 ? Nx-1 : i-1) }, i_right{ (i == Nx-1 ? 0 : i+1) }; 
    std::size_t j_left{ j }, j_right{ j };

    for (int d = 0; d < dimPk; ++d) 
    {
            uL_left += uh[i_left][j_left][d][n] * phixR[d][q];
            uL_right += uh[i][j][d][n] * phixL[d][q];

            uR_left += uh[i][j][d][n] * phixR[d][q];
            uR_right += uh[i_right][j_right][d][n] * phixL[d][q];
    }
}

void calculateUyInterface(const StateArray& uh, std::size_t i, std::size_t j,
                        std::size_t n, std::size_t q, Real& uD_down, Real& uD_up, 
                        Real& uU_down, Real& uU_up)
{
    uD_down = uD_up = uU_down = uU_up = 0.0;

    std::size_t i_down{ i }, i_up{ i };
    std::size_t j_down{ (j == 0 ? Ny-1 : j-1) }, j_up{ (j == Ny-1 ? 0 : j+1) };

    for (int d = 0; d < dimPk; ++d) 
    {
        // 下界面 (j-1/2)
        uD_down += uh[i_down][j_down][d][n] * phiyR[d][q];
        uD_up += uh[i][j][d][n] * phiyL[d][q];

        // 上界面 (j+1/2)
        uU_down += uh[i][j][d][n] * phiyR[d][q];
        uU_up += uh[i_up][j_up][d][n] * phiyL[d][q];
    }
}


