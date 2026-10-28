#include "constants.h"
#include "2D_shallowwater.h"
#include "COS_produre.h"
#include <array>
#include <iostream>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =====================================================================
//  初始资料与 flux
// =====================================================================
double b(double x, double y)
{
    double sp = std::sin(M_PI * x);
    double sq = std::sin(M_PI * y);
    return sp * sp * sq * sq;
}
double h(double x, double y)
{
    return 5.0 + std::exp(std::cos(2.0 * M_PI * x)) * std::exp(std::cos(2.0 * M_PI * y));
}
double hu(double x, double) { return std::sin(std::cos(2.0 * M_PI * x)); }
double hv(double, double y) { return std::sin(std::cos(2.0 * M_PI * y)); }

double f1(double, double hu, double) { return hu; }
double f2(double h, double hu, double) { return hu*hu/h + 0.5*g*h*h; }
double f3(double h, double hu, double hv) { return hu*hv/h; }
double g1(double, double, double hv) { return hv; }
double g2(double h, double hu, double hv) { return hu*hv/h; }
double g3(double h, double, double hv) { return hv*hv/h + 0.5*g*h*h; }

// =====================================================================
//  高斯点 / 界面值
// =====================================================================
void getBasic(std::vector<double>& phi,
              std::vector<double>& phiGLL,
              std::vector<double>& phi_xi,
              std::vector<double>& phi_eta)
{
    phi.assign(dimPk * NumGLP * NumGLP, 0.0);
    phiGLL.assign(dimPk * NumGLL * NumGLL, 0.0);
    phi_xi.assign(dimPk * NumGLP * NumGLP, 0.0);
    phi_eta.assign(dimPk * NumGLP * NumGLP, 0.0);

    for (std::size_t p = 0; p < NumGLP; ++p)
        for (std::size_t q = 0; q < NumGLP; ++q)
        {
            phi[idx_phi(0,p,q)] = 1.0;
            phi[idx_phi(1,p,q)] = gl_x[p];
            phi[idx_phi(2,p,q)] = gl_x[q];
            phi[idx_phi(3,p,q)] = gl_x[p] * gl_x[q];
            phi[idx_phi(4,p,q)] = gl_x[p] * gl_x[p] - 1.0/3.0;
            phi[idx_phi(5,p,q)] = gl_x[q] * gl_x[q] - 1.0/3.0;

            phi_xi[idx_phi(0,p,q)] = 0.0;
            phi_xi[idx_phi(1,p,q)] = 1.0;
            phi_xi[idx_phi(2,p,q)] = 0.0;
            phi_xi[idx_phi(3,p,q)] = gl_x[q];
            phi_xi[idx_phi(4,p,q)] = 2 * gl_x[p];
            phi_xi[idx_phi(5,p,q)] = 0.0;

            phi_eta[idx_phi(0,p,q)] = 0.0;
            phi_eta[idx_phi(1,p,q)] = 0.0;
            phi_eta[idx_phi(2,p,q)] = 1.0;
            phi_eta[idx_phi(3,p,q)] = gl_x[p];
            phi_eta[idx_phi(4,p,q)] = 0.0;
            phi_eta[idx_phi(5,p,q)] = 2 * gl_x[q];
        }

    for (std::size_t p = 0; p < NumGLL; ++p)
        for (std::size_t q = 0; q < NumGLL; ++q)
        {
            phiGLL[idx_phiGLL(0,p,q)] = 1.0;
            phiGLL[idx_phiGLL(1,p,q)] = gll_x[p];
            phiGLL[idx_phiGLL(2,p,q)] = gll_x[q];
            phiGLL[idx_phiGLL(3,p,q)] = gll_x[p] * gll_x[q];
            phiGLL[idx_phiGLL(4,p,q)] = gll_x[p] * gll_x[p] - 1.0/3.0;
            phiGLL[idx_phiGLL(5,p,q)] = gll_x[q] * gll_x[q] - 1.0/3.0;
        }
}

void getInterface(std::vector<double>& phixL,
                  std::vector<double>& phixR,
                  std::vector<double>& phiyL,
                  std::vector<double>& phiyR)
{
    phixL.assign(dimPk * NumGLP, 0.0);
    phixR.assign(dimPk * NumGLP, 0.0);
    phiyL.assign(dimPk * NumGLP, 0.0);
    phiyR.assign(dimPk * NumGLP, 0.0);

    for (std::size_t q = 0; q < NumGLP; ++q)
    {
        phixR[idx_phi1d(0,q)] =  1.0;
        phixR[idx_phi1d(1,q)] =  1.0;
        phixR[idx_phi1d(2,q)] =  gl_x[q];
        phixR[idx_phi1d(3,q)] =  gl_x[q];
        phixR[idx_phi1d(4,q)] =  2.0/3.0;
        phixR[idx_phi1d(5,q)] =  gl_x[q]*gl_x[q] - 1.0/3.0;

        phixL[idx_phi1d(0,q)] =  1.0;
        phixL[idx_phi1d(1,q)] = -1.0;
        phixL[idx_phi1d(2,q)] =  gl_x[q];
        phixL[idx_phi1d(3,q)] = -gl_x[q];
        phixL[idx_phi1d(4,q)] =  2.0/3.0;
        phixL[idx_phi1d(5,q)] =  gl_x[q]*gl_x[q] - 1.0/3.0;
    }

    for (std::size_t p = 0; p < NumGLP; ++p)
    {
        phiyR[idx_phi1d(0,p)] =  1.0;
        phiyR[idx_phi1d(1,p)] =  gl_x[p];
        phiyR[idx_phi1d(2,p)] =  1.0;
        phiyR[idx_phi1d(3,p)] =  gl_x[p];
        phiyR[idx_phi1d(4,p)] =  gl_x[p]*gl_x[p] - 1.0/3.0;
        phiyR[idx_phi1d(5,p)] =  2.0/3.0;

        phiyL[idx_phi1d(0,p)] =  1.0;
        phiyL[idx_phi1d(1,p)] =  gl_x[p];
        phiyL[idx_phi1d(2,p)] = -1.0;
        phiyL[idx_phi1d(3,p)] = -gl_x[p];
        phiyL[idx_phi1d(4,p)] =  gl_x[p]*gl_x[p] - 1.0/3.0;
        phiyL[idx_phi1d(5,p)] =  2.0/3.0;
    }
}

// =====================================================================
//  L² 投影
// =====================================================================
void getProjection(StateArray& uh, BottomArray& bh)
{
    uh.assign(Nx * Ny * dimPk * NumEq, 0.0);
    bh.assign(Nx * Ny * dimPk, 0.0);

    for (std::size_t i = 0; i < Nx; ++i)
    {
        double x_i = xa + (static_cast<double>(i) + 0.5) * hx;
        for (std::size_t j = 0; j < Ny; ++j)
        {
            double y_j = ya + (static_cast<double>(j) + 0.5) * hy;
            for (std::size_t d = 0; d < dimPk; ++d)
            {
                for (std::size_t p = 0; p < NumGLP; ++p)
                {
                    double xi_p = gl_x[p];
                    double w_p  = gl_w[p];
                    for (std::size_t q = 0; q < NumGLP; ++q)
                    {
                        double eta_q = gl_x[q];
                        double w_q   = gl_w[q];

                        uh[idx_state(i,j,d,0)] +=
                            w_p * w_q * h(hx1*xi_p + x_i, hy1*eta_q + y_j)
                            * phi[idx_phi(d,p,q)];
                        uh[idx_state(i,j,d,1)] +=
                            w_p * w_q * hu(hx1*xi_p + x_i, hy1*eta_q + y_j)
                            * phi[idx_phi(d,p,q)];
                        uh[idx_state(i,j,d,2)] +=
                            w_p * w_q * hv(hx1*xi_p + x_i, hy1*eta_q + y_j)
                            * phi[idx_phi(d,p,q)];
                            
                        bh[idx_bottom(i,j,d)] +=
                            w_p * w_q * b(hx1*xi_p + x_i, hy1*eta_q + y_j)
                            * phi[idx_phi(d,p,q)];
                    }
                }
                uh[idx_state(i,j,d,0)] /= (4 * mm[d]);
                uh[idx_state(i,j,d,1)] /= (4 * mm[d]);
                uh[idx_state(i,j,d,2)] /= (4 * mm[d]);
                bh[idx_bottom(i,j,d)]  /= (4 * mm[d]);
            }
        }
    }
}

// =====================================================================
//  RK3
// =====================================================================
void runRK3(StateArray& uh, const BottomArray& bh, double t_end)
{
    double t = 0.0;
    int step = 0;

    StateArray wh(Nx * Ny * dimPk * NumEq, 0.0);
    StateArray u1(Nx * Ny * dimPk * NumEq, 0.0);
    StateArray u2(Nx * Ny * dimPk * NumEq, 0.0);
    StateArray du(Nx * Ny * dimPk * NumEq, 0.0);

    while (t < t_end)
    {
        double alpha_1 = 1e-15, alpha_2 = 1e-15;
        for (std::size_t i = 0; i < Nx; ++i)
            for (std::size_t j = 0; j < Ny; ++j)
            {
                double h_  = uh[idx_state(i,j,0,0)];
                double hu_ = uh[idx_state(i,j,0,1)];
                double hv_ = uh[idx_state(i,j,0,2)];
                double u_bar = hu_ / h_;
                double v_bar = hv_ / h_;
                alpha_1 = std::max(alpha_1, std::abs(u_bar) + std::sqrt(g * h_));
                alpha_2 = std::max(alpha_2, std::abs(v_bar) + std::sqrt(g * h_));
            }

        double dt = CFL / (alpha_1 / hx + alpha_2 / hy);
        if (t + dt >= t_end) { dt = t_end - t; t = t_end; }
        else                 { t += dt; }

        // ---- Stage 1 ----
        Lh(uh, bh, du, alpha_1, alpha_2);
        for (std::size_t i = 0; i < Nx; ++i)
            for (std::size_t j = 0; j < Ny; ++j)
                for (std::size_t d = 0; d < dimPk; ++d)
                    for (std::size_t n = 0; n < NumEq; ++n)
                        u1[idx_state(i,j,d,n)] =
                            uh[idx_state(i,j,d,n)] + dt * du[idx_state(i,j,d,n)];
        positivityPreservingLimiter(uh);

        // ---- Stage 2 ----
        Lh(u1, bh, du, alpha_1, alpha_2);
        for (std::size_t i = 0; i < Nx; ++i)
            for (std::size_t j = 0; j < Ny; ++j)
                for (std::size_t d = 0; d < dimPk; ++d)
                    for (std::size_t n = 0; n < NumEq; ++n)
                        u2[idx_state(i,j,d,n)] =
                            0.75 * uh[idx_state(i,j,d,n)] +
                            0.25 * u1[idx_state(i,j,d,n)] +
                            0.25 * dt * du[idx_state(i,j,d,n)];
        positivityPreservingLimiter(uh);

        // ---- Stage 3 ----
        Lh(u2, bh, du, alpha_1, alpha_2);
        for (std::size_t i = 0; i < Nx; ++i)
            for (std::size_t j = 0; j < Ny; ++j)
                for (std::size_t d = 0; d < dimPk; ++d)
                    for (std::size_t n = 0; n < NumEq; ++n)
                        uh[idx_state(i,j,d,n)] =
                            (1.0/3.0) * uh[idx_state(i,j,d,n)] +
                            (2.0/3.0) * u2[idx_state(i,j,d,n)] +
                            (2.0/3.0) * dt * du[idx_state(i,j,d,n)];
        positivityPreservingLimiter(uh);

        // ---- COS ----
        for (std::size_t i = 0; i < Nx; ++i)
            for (std::size_t j = 0; j < Ny; ++j)
                for (std::size_t d = 0; d < dimPk; ++d)
                {
                    wh[idx_state(i,j,d,0)] = uh[idx_state(i,j,d,0)] + bh[idx_bottom(i,j,d)];
                    wh[idx_state(i,j,d,1)] = uh[idx_state(i,j,d,1)];
                    wh[idx_state(i,j,d,2)] = uh[idx_state(i,j,d,2)];
                }
        COS_produre(wh, bh, dt);
        for (std::size_t i = 0; i < Nx; ++i)
            for (std::size_t j = 0; j < Ny; ++j)
                for (std::size_t d = 0; d < dimPk; ++d)
                {
                    uh[idx_state(i,j,d,0)] = wh[idx_state(i,j,d,0)] - bh[idx_bottom(i,j,d)];
                    uh[idx_state(i,j,d,1)] = wh[idx_state(i,j,d,1)];
                    uh[idx_state(i,j,d,2)] = wh[idx_state(i,j,d,2)];
                }

        if (step % 10 == 0)
            std::cout << "Step: " << step
                      << ", t = " << t << ", dt = " << dt << "\n";
        ++step;
        // if (step == 1) break;
    }
}

// =====================================================================
//  空间算子
// =====================================================================
void Lh(const StateArray&  uh,
        const BottomArray& bh,
        StateArray&        du,
        double             alpha_1,
        double             alpha_2)
{
    std::fill(du.begin(), du.end(), 0.0);

    // ---------- Step 1: 体积分 + 源项 ----------
    for (std::size_t i = 0; i < Nx; ++i)
        for (std::size_t j = 0; j < Ny; ++j)
            for (std::size_t d = 0; d < dimPk; ++d)
            {
                for (std::size_t p = 0; p < NumGLP; ++p)
                {
                    double w_p = gl_w[p];
                    for (std::size_t q = 0; q < NumGLP; ++q)
                    {
                        double w_q = gl_w[q];

                        double h_, hu_, hv_;
                        calculateUh(uh, i, j, p, q, h_, hu_, hv_);

                        du[idx_state(i,j,d,0)] += w_p * w_q * (
                            hy1 * f1(h_, hu_, hv_) * phi_xi [idx_phi(d,p,q)] +
                            hx1 * g1(h_, hu_, hv_) * phi_eta[idx_phi(d,p,q)]);

                        du[idx_state(i,j,d,1)] += w_p * w_q * (
                            hy1 * f2(h_, hu_, hv_) * phi_xi [idx_phi(d,p,q)] +
                            hx1 * g2(h_, hu_, hv_) * phi_eta[idx_phi(d,p,q)]);

                        du[idx_state(i,j,d,2)] += w_p * w_q * (
                            hy1 * f3(h_, hu_, hv_) * phi_xi [idx_phi(d,p,q)] +
                            hx1 * g3(h_, hu_, hv_) * phi_eta[idx_phi(d,p,q)]);

                        double b_xi, b_eta;
                        calculateBhDeriv(bh, i, j, p, q, b_xi, b_eta);

                        du[idx_state(i,j,d,1)] += hy1 * (-g) * w_p * w_q * h_ * b_xi  * phi[idx_phi(d,p,q)];
                        du[idx_state(i,j,d,2)] += hx1 * (-g) * w_p * w_q * h_ * b_eta * phi[idx_phi(d,p,q)];
                    }
                }
            }

    // ---------- Step 2: x 方向界面 ----------
    for (std::size_t i = 0; i < Nx; ++i)
        for (std::size_t j = 0; j < Ny; ++j)
            for (std::size_t d = 0; d < dimPk; ++d)
            {
                for (std::size_t q = 0; q < NumGLP; ++q)
                {
                    double w_q = gl_w[q];

                    std::array<double, NumEq> uR_left{}, uR_right{};
                    std::array<double, NumEq> uL_left{}, uL_right{};
                    for (std::size_t n = 0; n < NumEq; ++n)
                        calculateUxInterface(uh, i, j, n, q,
                                             uL_left[n], uL_right[n],
                                             uR_left[n], uR_right[n]);

                    double bL_left = 0.0, bL_right = 0.0, bR_left = 0.0, bR_right = 0.0;
                    std::size_t i_left  = (i == 0     ? Nx-1 : i-1);
                    std::size_t i_right = (i == Nx-1  ? 0    : i+1);
                    for (int dd = 0; dd < dimPk; ++dd)
                    {
                        bL_left  += bh[idx_bottom(i_left, j, dd)]  * phixR[idx_phi1d(dd,q)];
                        bL_right += bh[idx_bottom(i,      j, dd)]  * phixL[idx_phi1d(dd,q)];
                        bR_left  += bh[idx_bottom(i,      j, dd)]  * phixR[idx_phi1d(dd,q)];
                        bR_right += bh[idx_bottom(i_right,j, dd)]  * phixL[idx_phi1d(dd,q)];
                    }

                    double hL_left  = uL_left[0],  hL_right = uL_right[0];
                    double hL_left_state  = std::max(0.0, hL_left  + bL_left  - std::max(bL_left,  bL_right));
                    double hL_right_state = std::max(0.0, hL_right + bL_right - std::max(bL_left,  bL_right));

                    uL_left[2]  = hL_left_state  * (uL_left[2]  / uL_left[0]);
                    uL_left[1]  = hL_left_state  * (uL_left[1]  / uL_left[0]);
                    uL_left[0]  = hL_left_state;
                    uL_right[2] = hL_right_state * (uL_right[2] / uL_right[0]);
                    uL_right[1] = hL_right_state * (uL_right[1] / uL_right[0]);
                    uL_right[0] = hL_right_state;

                    double hR_left  = uR_left[0],  hR_right = uR_right[0];
                    double hR_left_state  = std::max(0.0, hR_left  + bR_left  - std::max(bR_left,  bR_right));
                    double hR_right_state = std::max(0.0, hR_right + bR_right - std::max(bR_left,  bR_right));

                    uR_left[2]  = hR_left_state  * (uR_left[2]  / uR_left[0]);
                    uR_left[1]  = hR_left_state  * (uR_left[1]  / uR_left[0]);
                    uR_left[0]  = hR_left_state;
                    uR_right[2] = hR_right_state * (uR_right[2] / uR_right[0]);
                    uR_right[1] = hR_right_state * (uR_right[1] / uR_right[0]);
                    uR_right[0] = hR_right_state;

                    std::array<double, NumEq> fhatL_state{}, fhatR_state{};
                    for (std::size_t n = 0; n < NumEq; ++n)
                    {
                        fhatL_state[n] = LF_Flux_f(n, uL_left,  uL_right, alpha_1);
                        fhatR_state[n] = LF_Flux_f(n, uR_left,  uR_right, alpha_1);
                        if (n == 1)
                        {
                            fhatL_state[n] += 0.5 * g * (hL_right*hL_right - hL_right_state*hL_right_state);
                            fhatR_state[n] += 0.5 * g * (hR_left *hR_left  - hR_left_state *hR_left_state);
                        }
                    }
                    for (std::size_t n = 0; n < NumEq; ++n)
                        du[idx_state(i,j,d,n)] -= hy1 * w_q *
                            (fhatR_state[n] * phixR[idx_phi1d(d,q)] -
                             fhatL_state[n] * phixL[idx_phi1d(d,q)]);
                }
            }

    // ---------- Step 3: y 方向界面 ----------
    for (std::size_t i = 0; i < Nx; ++i)
        for (std::size_t j = 0; j < Ny; ++j)
            for (std::size_t d = 0; d < dimPk; ++d)
            {
                for (std::size_t p = 0; p < NumGLP; ++p)
                {
                    double w_p = gl_w[p];

                    std::array<double, NumEq> uU_down{}, uU_up{};
                    std::array<double, NumEq> uD_down{}, uD_up{};
                    for (std::size_t n = 0; n < NumEq; ++n)
                        calculateUyInterface(uh, i, j, n, p,
                                             uD_down[n], uD_up[n],
                                             uU_down[n], uU_up[n]);

                    double bD_down = 0.0, bD_up = 0.0, bU_down = 0.0, bU_up = 0.0;
                    std::size_t j_down = (j == 0    ? Ny-1 : j-1);
                    std::size_t j_up   = (j == Ny-1 ? 0    : j+1);
                    for (int dd = 0; dd < dimPk; ++dd)
                    {
                        bD_down += bh[idx_bottom(i, j_down, dd)] * phiyR[idx_phi1d(dd,p)];
                        bD_up   += bh[idx_bottom(i, j,      dd)] * phiyL[idx_phi1d(dd,p)];
                        bU_down += bh[idx_bottom(i, j,      dd)] * phiyR[idx_phi1d(dd,p)];
                        bU_up   += bh[idx_bottom(i, j_up,   dd)] * phiyL[idx_phi1d(dd,p)];
                    }

                    double hD_down = uD_down[0], hD_up = uD_up[0];
                    double hD_down_state = std::max(0.0, hD_down + bD_down - std::max(bD_down, bD_up));
                    double hD_up_state   = std::max(0.0, hD_up   + bD_up   - std::max(bD_down, bD_up));

                    uD_down[2] = hD_down_state * (uD_down[2] / uD_down[0]);
                    uD_down[1] = hD_down_state * (uD_down[1] / uD_down[0]);
                    uD_down[0] = hD_down_state;
                    uD_up[2]   = hD_up_state   * (uD_up[2]   / uD_up[0]);
                    uD_up[1]   = hD_up_state   * (uD_up[1]   / uD_up[0]);
                    uD_up[0]   = hD_up_state;

                    double hU_down = uU_down[0], hU_up = uU_up[0];
                    double hU_down_state = std::max(0.0, hU_down + bU_down - std::max(bU_down, bU_up));
                    double hU_up_state   = std::max(0.0, hU_up   + bU_up   - std::max(bU_down, bU_up));

                    uU_down[2] = hU_down_state * (uU_down[2] / uU_down[0]);
                    uU_down[1] = hU_down_state * (uU_down[1] / uU_down[0]);
                    uU_down[0] = hU_down_state;
                    uU_up[2]   = hU_up_state   * (uU_up[2]   / uU_up[0]);
                    uU_up[1]   = hU_up_state   * (uU_up[1]   / uU_up[0]);
                    uU_up[0]   = hU_up_state;

                    std::array<double, NumEq> ghatL_state{}, ghatR_state{};
                    for (std::size_t n = 0; n < NumEq; ++n)
                    {
                        ghatL_state[n] = LF_Flux_g(n, uD_down, uD_up, alpha_2);
                        ghatR_state[n] = LF_Flux_g(n, uU_down, uU_up, alpha_2);
                        if (n == 2)
                        {
                            ghatL_state[n] += 0.5 * g * (hD_up  *hD_up   - hD_up_state  *hD_up_state);
                            ghatR_state[n] += 0.5 * g * (hU_down*hU_down - hU_down_state*hU_down_state);
                        }
                    }
                    for (std::size_t n = 0; n < NumEq; ++n)
                        du[idx_state(i,j,d,n)] -= hx1 * w_p *
                            (ghatR_state[n] * phiyR[idx_phi1d(d,p)] -
                             ghatL_state[n] * phiyL[idx_phi1d(d,p)]);
                }
            }

    // ---------- 除以 hx·hy·mm[d] ----------
    for (std::size_t i = 0; i < Nx; ++i)
        for (std::size_t j = 0; j < Ny; ++j)
            for (std::size_t d = 0; d < dimPk; ++d)
                for (std::size_t n = 0; n < NumEq; ++n)
                    du[idx_state(i,j,d,n)] /= hx * hy * mm[d];
}

// =====================================================================
//  Positivity-preserving limiter
// =====================================================================
void positivityPreservingLimiter(StateArray& uh)
{
    for (std::size_t i = 0; i < Nx; ++i)
        for (std::size_t j = 0; j < Ny; ++j)
        {
            double m_ij = 0.0;
            for (std::size_t d = 0; d < dimPk; ++d)
                m_ij += uh[idx_state(i,j,d,0)] * phiGLL[idx_phiGLL(d,0,0)];

            for (std::size_t p = 0; p < NumGLL; ++p)
                for (std::size_t q = 0; q < NumGLL; ++q)
                {
                    double h_ij = 0.0;
                    for (std::size_t d = 0; d < dimPk; ++d)
                        h_ij += uh[idx_state(i,j,d,0)] * phiGLL[idx_phiGLL(d,p,q)];
                    if (h_ij < m_ij) m_ij = h_ij;
                }

            double h_bar = uh[idx_state(i,j,0,0)];
            double theta = std::min(1.0, h_bar / (h_bar - m_ij));

            for (std::size_t d = 1; d < dimPk; ++d)
                for (std::size_t n = 0; n < NumEq; ++n)
                    uh[idx_state(i,j,d,n)] = theta * uh[idx_state(i,j,d,n)];
        }
}

// =====================================================================
//  LF flux
// =====================================================================
double LF_Flux_f(std::size_t n, std::array<double, NumEq>& L,
                                   std::array<double, NumEq>& R, double alpha)
{
    switch (n)
    {
        case 0: return 0.5 * (f1(R[0],R[1],R[2]) + f1(L[0],L[1],L[2])) - 0.5*alpha*(R[0]-L[0]);
        case 1: return 0.5 * (f2(R[0],R[1],R[2]) + f2(L[0],L[1],L[2])) - 0.5*alpha*(R[1]-L[1]);
        case 2: return 0.5 * (f3(R[0],R[1],R[2]) + f3(L[0],L[1],L[2])) - 0.5*alpha*(R[2]-L[2]);
        default: return 666.0;
    }
}
double LF_Flux_g(std::size_t n, std::array<double, NumEq>& L,
                                   std::array<double, NumEq>& R, double alpha)
{
    switch (n)
    {
        case 0: return 0.5 * (g1(R[0],R[1],R[2]) + g1(L[0],L[1],L[2])) - 0.5*alpha*(R[0]-L[0]);
        case 1: return 0.5 * (g2(R[0],R[1],R[2]) + g2(L[0],L[1],L[2])) - 0.5*alpha*(R[1]-L[1]);
        case 2: return 0.5 * (g3(R[0],R[1],R[2]) + g3(L[0],L[1],L[2])) - 0.5*alpha*(R[2]-L[2]);
        default: return 666.0;
    }
}

// =====================================================================
//  求值辅助
// =====================================================================
void calculateUh(const StateArray& uh, std::size_t i, std::size_t j,
                 std::size_t p, std::size_t q, double& h_, double& hu_, double& hv_)
{
    h_ = hu_ = hv_ = 0.0;
    for (std::size_t d = 0; d < dimPk; ++d)
    {
        double v = phi[idx_phi(d,p,q)];
        h_  += uh[idx_state(i,j,d,0)] * v;
        hu_ += uh[idx_state(i,j,d,1)] * v;
        hv_ += uh[idx_state(i,j,d,2)] * v;
    }
}

void calculateBhDeriv(const BottomArray& bh, std::size_t i, std::size_t j,
                      std::size_t p, std::size_t q, double& b_xi, double& b_eta)
{
    b_xi = b_eta = 0.0;
    for (int d = 0; d < dimPk; ++d)
    {
        double bval = bh[idx_bottom(i,j,d)];
        b_xi  += bval * phi_xi [idx_phi(d,p,q)];
        b_eta += bval * phi_eta[idx_phi(d,p,q)];
    }
}

void calculateUxInterface(const StateArray& uh, std::size_t i, std::size_t j,
                          std::size_t n, std::size_t q,
                          double& uL_left, double& uL_right,
                          double& uR_left, double& uR_right)
{
    uL_left = uL_right = uR_left = uR_right = 0.0;
    std::size_t i_left  = (i == 0    ? Nx-1 : i-1);
    std::size_t i_right = (i == Nx-1 ? 0    : i+1);

    for (int d = 0; d < dimPk; ++d)
    {
        uL_left  += uh[idx_state(i_left, j, d, n)] * phixR[idx_phi1d(d,q)];
        uL_right += uh[idx_state(i,      j, d, n)] * phixL[idx_phi1d(d,q)];
        uR_left  += uh[idx_state(i,      j, d, n)] * phixR[idx_phi1d(d,q)];
        uR_right += uh[idx_state(i_right,j, d, n)] * phixL[idx_phi1d(d,q)];
    }
}

void calculateUyInterface(const StateArray& uh, std::size_t i, std::size_t j,
                          std::size_t n, std::size_t q,
                          double& uD_down, double& uD_up,
                          double& uU_down, double& uU_up)
{
    uD_down = uD_up = uU_down = uU_up = 0.0;
    std::size_t j_down = (j == 0    ? Ny-1 : j-1);
    std::size_t j_up   = (j == Ny-1 ? 0    : j+1);

    for (int d = 0; d < dimPk; ++d)
    {
        uD_down += uh[idx_state(i, j_down, d, n)] * phiyR[idx_phi1d(d,q)];
        uD_up   += uh[idx_state(i, j,      d, n)] * phiyL[idx_phi1d(d,q)];
        uU_down += uh[idx_state(i, j,      d, n)] * phiyR[idx_phi1d(d,q)];
        uU_up   += uh[idx_state(i, j_up,   d, n)] * phiyL[idx_phi1d(d,q)];
    }
}