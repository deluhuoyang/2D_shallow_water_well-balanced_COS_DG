#include "constants.h"
#include "2D_shallowwater.h"
#include "COS_produre.h"
#include <cmath>
#include <array>

// =====================================================================
// 2D COS 过程
// =====================================================================
void COS_produre(StateArray& wh, const BottomArray& bh, Real dt)
{
    for (std::size_t i { 0 }; i < Nx; ++i) 
    {
        for (std::size_t j { 0 }; j < Ny; ++j) 
        {
            
            // 1. 计算单元I_ij四个边界上的 beta_left, beta_right, beta_down, beta_up
            Real beta_left{}, beta_right{}, beta_down{}, beta_up{};
            calculateBeta(wh, bh, i, j, beta_left, beta_right, beta_down, beta_up);                    

            // 2. 计算单元I_ij四个边界上的 sigma_left, sigma_right, sigma_down, sigma_up
            Real sigma_left{}, sigma_right{}, sigma_down{}, sigma_up{};
            calculateSigma(wh, bh, i, j, sigma_left, sigma_right, sigma_down, sigma_up);

            // 3. 计算 lambda_ij
            Real lambda_ij = std::exp(- dt * 
                ((beta_left * sigma_left + beta_right * sigma_right) / hx
                + (beta_down * sigma_down + beta_up * sigma_up) / hy)
                );

            // 4. 更新 uh_ij
            for (std::size_t d { 1 }; d < dimPk; ++d)
                for (std::size_t n { 0 }; n < NumEq; ++n)
                        wh[i][j][d][n] = lambda_ij * wh[i][j][d][n];
        
        }
    }
}        

// =====================================================================
// 计算 COS 过程的 beta, sigma
// =====================================================================
void calculateBeta(const StateArray& wh, const BottomArray& bh, 
                std::size_t i, std::size_t j, Real& beta_left, Real& beta_right, Real& beta_down, Real& beta_up)
{
    // 1. 计算单元I_ij相邻的4个单元的索引，周期边界
    std::size_t i_left{ (i == 0 ? Nx-1 : i-1) }, i_right{ (i == Nx-1 ? 0 : i+1) }; 
    std::size_t j_left{ j }, j_right{ j };

    std::size_t i_down{ i }, i_up{ i };
    std::size_t j_down{ (j == 0 ? Ny-1 : j-1) }, j_up{ (j == Ny-1 ? 0 : j+1) }; 
    
    // 2. 左右
    Real wave_speed_u_central { std::abs((wh[i][j][0][1]-bh[i][j][0])/ wh[i][j][0][0]) 
                                + std::sqrt(g * wh[i][j][0][0]) }; 
    Real wave_speed_u_left { std::abs((wh[i_left][j_left][0][1]-bh[i_left][j_left][0])/ wh[i_left][j_left][0][0]) 
                                + std::sqrt(g * wh[i_left][j_left][0][0]) }; 
    Real wave_speed_u_right { std::abs((wh[i_right][j_right][0][1]-bh[i_right][j_right][0])/ wh[i_right][j_right][0][0]) 
                                + std::sqrt(g * wh[i_right][j_right][0][0]) };                     

    beta_left = std::max(wave_speed_u_left, wave_speed_u_central);
    beta_right = std::max(wave_speed_u_right, wave_speed_u_central);

    // 3. 上下
    Real wave_speed_v_central { std::abs((wh[i][j][0][2]-bh[i][j][0]) / wh[i][j][0][0]) 
                        + std::sqrt(g * wh[i][j][0][0]) }; 
    Real wave_speed_v_down { std::abs((wh[i_down][j_down][0][2]-bh[i_down][j_down][0]) / wh[i_down][j_down][0][0]) 
                        + std::sqrt(g * wh[i_down][j_down][0][0]) }; 
    Real wave_speed_v_up { std::abs((wh[i_up][j_up][0][2]-bh[i_up][j_up][0]) / wh[i_up][j_up][0][0]) 
                        + std::sqrt(g * wh[i_up][j_up][0][0]) };                     

    beta_down = std::max(wave_speed_v_down, wave_speed_v_central);
    beta_up = std::max(wave_speed_v_up, wave_speed_v_central);                    
}

void calculateSigma(const StateArray& wh, const BottomArray& bh, std::size_t i, std::size_t j, 
                Real& sigma_left, Real& sigma_right, Real& sigma_down, Real& sigma_up)
{
    // 1. 计算单元 I_ij 相邻的4个单元的索引，周期边界
    std::size_t i_left{ (i == 0 ? Nx-1 : i-1) }, i_right{ (i == Nx-1 ? 0 : i+1) }; 
    std::size_t j_left{ j }, j_right{ j };

    std::size_t i_down{ i }, i_up{ i };
    std::size_t j_down{ (j == 0 ? Ny-1 : j-1) }, j_up{ (j == Ny-1 ? 0 : j+1) };

    // 2. 计算单元 I_ij 的平均值
    Real h_bar = wh[i][j][0][0] - bh[i][j][0];   // 真实的 \bar h
    Real hu_bar = wh[i][j][0][1];
    Real hv_bar = wh[i][j][0][2];
    Real u_bar = (h_bar > dry_tol) ? hu_bar / h_bar : 0.0;
    Real v_bar = (h_bar > dry_tol) ? hv_bar / h_bar : 0.0;

    // 3. 计算 epsilon
    Real epsilon {
        static_cast<Real>(1e-8) * std::max(static_cast<Real>(1.0), 
        static_cast<Real>(0.5) * (hu_bar*u_bar + hv_bar*v_bar + g*h_bar*h_bar))
    };
 
    // 4. 计算 theta
    Real h{}, hu{}, hv{}; // 记录左右、上下单元的值
    Real delta_h{}, delta_hu{}, delta_hv{}; // 记录差值

    Real theta_left{ 0.0 }, theta_right{ 0.0 }, theta_down{ 0.0 }, theta_up{ 0.0 };
    for (std::size_t p { 0 }; p < NumGLP; ++p)
    {
        Real w_p { gl_w[p] };

        for (std::size_t q { 0 }; q < NumGLP; ++q)
        {
            Real w_q { gl_w[q] };
            // 左
            calculateUh(wh, i_left, j_left, p, q, h, hu, hv);
            delta_h = h-h_bar;
            delta_hu = hu-hu_bar;
            delta_hv = hv-hv_bar;
            theta_left += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, delta_h, delta_hu, delta_hv);
            // 右
            calculateUh(wh, i_right, j_right, p, q, h, hu, hv);
            delta_h = h-h_bar;
            delta_hu = hu-hu_bar;
            delta_hv = hv-hv_bar;
            theta_right += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, delta_h, delta_hu, delta_hv);
            // 下
            calculateUh(wh, i_down, j_down, p, q, h, hu, hv);
            delta_h = h-h_bar;
            delta_hu = hu-hu_bar;
            delta_hv = hv-hv_bar;
            theta_down += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, delta_h, delta_hu, delta_hv);
            // 上
            calculateUh(wh, i_up, j_up, p, q, h, hu, hv);
            delta_h = h-h_bar;
            delta_hu = hu-hu_bar;
            delta_hv = hv-hv_bar;
            theta_up += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, delta_h, delta_hu, delta_hv);
        }
    }
    theta_left *= 0.25 * hx * hy; theta_down *= 0.25 * hx * hy; 
    theta_right *= 0.25 * hx * hy; theta_up *= 0.25 * hx * hy;

    // 5. sigma 的分子项范数
    Real norm_left{ 0.0 }, norm_right{ 0.0 }, norm_down{ 0.0 }, norm_up{ 0.0 };

    Real h_center{}, hu_center{}, hv_center{}; // 记录左右、上下单元的值
    std::array<Real, dimPk> phi_extrap{}; // 记录外推函数在高斯点的值

    for (std::size_t p { 0 }; p < NumGLP; ++p)
    {
        Real xi_p { gl_x[p] }, w_p { gl_w[p] };

        for (std::size_t q { 0 }; q < NumGLP; ++q)
        {
            Real eta_q { gl_x[q] }, w_q { gl_w[q] };

            // 1. 左
            if (theta_left >= epsilon)
            {
                calculateUh(wh, i_left, j_left, p, q, h, hu, hv);

                // 计算单元 I_ij 的外推值，
                // static<Real>() 用来处理 double 到 float 的窄化转换警告
                phi_extrap = {1.0, xi_p-static_cast<Real>(2.0), eta_q, 
                    (xi_p-static_cast<Real>(2.0))*eta_q, 
                    (xi_p-static_cast<Real>(2.0))*(xi_p-static_cast<Real>(2.0))-static_cast<Real>(1.0/3.0), 
                    eta_q*eta_q-static_cast<Real>(1.0/3.0)};

                h_center = hu_center = hv_center = 0.0;
                for (std::size_t d { 0 }; d < dimPk; ++d) 
                {
                    h_center  += wh[i][j][d][0] * phi_extrap[d];
                    hu_center += wh[i][j][d][1] * phi_extrap[d];
                    hv_center += wh[i][j][d][2] * phi_extrap[d];
                }

                // 计算差值
                delta_h = h-h_center;
                delta_hu = hu-hu_center;
                delta_hv = hv-hv_center;
                norm_left += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, delta_h, delta_hu, delta_hv);
            }
            
            // 2. 右
            if (theta_right >= epsilon)
            {
                calculateUh(wh, i_right, j_right, p, q, h, hu, hv);

                phi_extrap = {1.0, xi_p+static_cast<Real>(2.0), 
                    eta_q, (xi_p+static_cast<Real>(2.0))*eta_q, 
                    (xi_p+static_cast<Real>(2.0))*(xi_p+static_cast<Real>(2.0))-static_cast<Real>(1.0/3.0), 
                    eta_q*eta_q-static_cast<Real>(1.0/3.0)};

                h_center = hu_center = hv_center = 0.0;
                for (std::size_t d { 0 }; d < dimPk; ++d) 
                {
                    h_center  += wh[i][j][d][0] * phi_extrap[d];
                    hu_center += wh[i][j][d][1] * phi_extrap[d];
                    hv_center += wh[i][j][d][2] * phi_extrap[d];
                }

                delta_h = h-h_center;
                delta_hu = hu-hu_center;
                delta_hv = hv-hv_center;
                norm_right += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, delta_h, delta_hu, delta_hv);
            }

            // 3. 下
            if (theta_down > epsilon)
            {
                calculateUh(wh, i_down, j_down, p, q, h, hu, hv);

                phi_extrap = {1.0, xi_p, eta_q-static_cast<Real>(2.0), 
                    xi_p*(eta_q-static_cast<Real>(2.0)), 
                    xi_p*xi_p-static_cast<Real>(1.0/3.0), 
                    (eta_q-static_cast<Real>(2.0))*(eta_q-static_cast<Real>(2.0))-static_cast<Real>(1.0/3.0)};

                h_center = hu_center = hv_center = 0.0;
                for (std::size_t d { 0 }; d < dimPk; ++d) 
                {
                    h_center  += wh[i][j][d][0] * phi_extrap[d];
                    hu_center += wh[i][j][d][1] * phi_extrap[d];
                    hv_center += wh[i][j][d][2] * phi_extrap[d];
                }

                delta_h = h-h_center;
                delta_hu = hu-hu_center;
                delta_hv = hv-hv_center;
                norm_down += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, delta_h, delta_hu, delta_hv);
            }
            
            // 4. 上
            if (theta_up > epsilon)
            {
                calculateUh(wh, i_up, j_up, p, q, h, hu, hv);
            
                phi_extrap = {1.0, xi_p, eta_q+static_cast<Real>(2.0), 
                    xi_p*(eta_q+static_cast<Real>(2.0)), 
                    xi_p*xi_p-static_cast<Real>(1.0/3.0), 
                    (eta_q+static_cast<Real>(2.0))*(eta_q+static_cast<Real>(2.0))-static_cast<Real>(1.0/3.0)};

                h_center = hu_center = hv_center = 0.0;
                for (std::size_t d { 0 }; d < dimPk; ++d) 
                {
                    h_center  += wh[i][j][d][0] * phi_extrap[d];
                    hu_center += wh[i][j][d][1] * phi_extrap[d];
                    hv_center += wh[i][j][d][2] * phi_extrap[d];
                }

                delta_h = h-h_center;
                delta_hu = hu-hu_center;
                delta_hv = hv-hv_center;
                norm_up += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, delta_h, delta_hu, delta_hv);
            }
        }
    }

    norm_left *= 0.25 * hx * hy; norm_down *= 0.25 * hx * hy; 
    norm_right *= 0.25 * hx * hy; norm_up *= 0.25 * hx * hy;

    // 6. 计算 sigma
    sigma_left  = (theta_left  >= epsilon) ? Ck * norm_left  / theta_left  : 0.0;
    sigma_right = (theta_right >= epsilon) ? Ck * norm_right / theta_right : 0.0;
    sigma_down  = (theta_down  >= epsilon) ? Ck * norm_down  / theta_down  : 0.0;
    sigma_up    = (theta_up    >= epsilon) ? Ck * norm_up    / theta_up    : 0.0;
}

Real HessianNorm(Real h_center, Real u_center, Real v_center, 
                Real delta_h, Real delta_hu, Real delta_hv)
{
    return g * h_center * delta_h * delta_h +
            (delta_hu - u_center*delta_h) * (delta_hu - u_center*delta_h) +
            (delta_hv - v_center*delta_h) * (delta_hv - v_center*delta_h);
}