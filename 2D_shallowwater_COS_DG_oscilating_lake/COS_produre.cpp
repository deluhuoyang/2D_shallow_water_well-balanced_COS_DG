#include "constants.h"
#include "2D_shallowwater.h"
#include "COS_produre.h"
#include <cmath>
#include <array>

// =====================================================================
//  COS 滤波器：直接对 uh = (h, hu, hv) 作用
// =====================================================================
void COS_produre(StateArray& uh, const BottomArray& bh, double dt)
{
    for (std::size_t i = 0; i < Nx; ++i)
        for (std::size_t j = 0; j < Ny; ++j)
        {
            double beta_left, beta_right, beta_down, beta_up;
            calculateBeta(uh, i, j, beta_left, beta_right, beta_down, beta_up);

            double sigma_left, sigma_right, sigma_down, sigma_up;
            calculateSigma(uh, bh, i, j, sigma_left, sigma_right, sigma_down, sigma_up);

            double lambda_ij = std::exp(-dt * (
                (beta_left*sigma_left + beta_right*sigma_right) / hx +
                (beta_down*sigma_down + beta_up  *sigma_up)    / hy));

            if (!std::isfinite(lambda_ij) || lambda_ij < 0.0) lambda_ij = 0.0;
            if (lambda_ij > 1.0)                              lambda_ij = 1.0;

            for (std::size_t d = 1; d < dimPk; ++d)
                for (std::size_t n = 0; n < NumEq; ++n)
                    uh[idx_state(i,j,d,n)] = lambda_ij * uh[idx_state(i,j,d,n)];
        }
}

// =====================================================================
//  calculateBeta
// =====================================================================
void calculateBeta(const StateArray& uh,
                   std::size_t i, std::size_t j,
                   double& beta_left, double& beta_right,
                   double& beta_down, double& beta_up)
{
    std::size_t i_left  = (i == 0    ? Nx-1 : i-1);
    std::size_t i_right = (i == Nx-1 ? 0    : i+1);
    std::size_t j_down  = (j == 0    ? Ny-1 : j-1);
    std::size_t j_up    = (j == Ny-1 ? 0    : j+1);

    // 辅助函数：提取 x 方向特征波速 |u| + c
    auto waveSpeed_x = [&](std::size_t ii, std::size_t jj) -> double
    {
        double h = uh[idx_state(ii,jj,0,0)];
        double hu = uh[idx_state(ii,jj,0,1)];
        
        if (h <= dry_tol) return h > 0.0 ? std::sqrt(g * h) : 0.0; // 干区波速为0，完美防除零
        
        double u = hu / h;
        double c = std::sqrt(g * h);
        return std::abs(u) + c;
    };

    // 辅助函数：提取 y 方向特征波速 |v| + c
    auto waveSpeed_y = [&](std::size_t ii, std::size_t jj) -> double
    {
        double h = uh[idx_state(ii,jj,0,0)];
        double hv = uh[idx_state(ii,jj,0,2)];
        
        if (h <= dry_tol) return h > 0.0 ? std::sqrt(g * h) : 0.0;
        
        double v = hv / h;
        double c = std::sqrt(g * h);
        return std::abs(v) + c;
    };

    // 按照论文公式求 max
    beta_left  = std::max(waveSpeed_x(i_left, j),  waveSpeed_x(i, j));
    beta_right = std::max(waveSpeed_x(i_right, j), waveSpeed_x(i, j));
    
    beta_down  = std::max(waveSpeed_y(i, j_down),  waveSpeed_y(i, j));
    beta_up    = std::max(waveSpeed_y(i, j_up),    waveSpeed_y(i, j));
}

// =====================================================================
//  calculateSigma
// =====================================================================
void calculateSigma(const StateArray& uh, const BottomArray& bh,
                    std::size_t i, std::size_t j,
                    double& sigma_left, double& sigma_right,
                    double& sigma_down, double& sigma_up)
{
    (void)bh; // 本函数不再需要 bh

    std::size_t i_left  = (i == 0    ? Nx-1 : i-1);
    std::size_t i_right = (i == Nx-1 ? 0    : i+1);
    std::size_t j_down  = (j == 0    ? Ny-1 : j-1);
    std::size_t j_up    = (j == Ny-1 ? 0    : j+1);

    double h_bar  = uh[idx_state(i,j,0,0)];   // ← 不再减 b
    double hu_bar = uh[idx_state(i,j,0,1)];
    double hv_bar = uh[idx_state(i,j,0,2)];
    double u_bar  = (h_bar > dry_tol) ? hu_bar / h_bar : 0.0;
    double v_bar  = (h_bar > dry_tol) ? hv_bar / h_bar : 0.0;

    double epsilon = 1e-8 * std::max(1.0,
        0.5 * (hu_bar*u_bar + hv_bar*v_bar + g*h_bar*h_bar));

    double h_ = 0.0, hu_ = 0.0, hv_ = 0.0;
    double dh = 0.0, dhu = 0.0, dhv = 0.0;
    double theta_left = 0.0, theta_right = 0.0, theta_down = 0.0, theta_up = 0.0;

    for (std::size_t p = 0; p < NumGLP; ++p)
    {
        double w_p = gl_w[p];
        for (std::size_t q = 0; q < NumGLP; ++q)
        {
            double w_q = gl_w[q];

            calculateUh(uh, i_left,  j,       p, q, h_, hu_, hv_);
            dh = h_ - h_bar; dhu = hu_ - hu_bar; dhv = hv_ - hv_bar;
            theta_left += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, dh, dhu, dhv);

            calculateUh(uh, i_right, j,       p, q, h_, hu_, hv_);
            dh = h_ - h_bar; dhu = hu_ - hu_bar; dhv = hv_ - hv_bar;
            theta_right += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, dh, dhu, dhv);

            calculateUh(uh, i,       j_down,  p, q, h_, hu_, hv_);
            dh = h_ - h_bar; dhu = hu_ - hu_bar; dhv = hv_ - hv_bar;
            theta_down += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, dh, dhu, dhv);

            calculateUh(uh, i,       j_up,    p, q, h_, hu_, hv_);
            dh = h_ - h_bar; dhu = hu_ - hu_bar; dhv = hv_ - hv_bar;
            theta_up += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, dh, dhu, dhv);
        }
    }
    theta_left  *= 0.25 * hx * hy;
    theta_right *= 0.25 * hx * hy;
    theta_down  *= 0.25 * hx * hy;
    theta_up    *= 0.25 * hx * hy;

    double norm_left = 0.0, norm_right = 0.0, norm_down = 0.0, norm_up = 0.0;

    double h_c = 0.0, hu_c = 0.0, hv_c = 0.0;
    std::array<double, dimPk> pe{};

    for (std::size_t p = 0; p < NumGLP; ++p)
    {
        double xi_p = gl_x[p], w_p = gl_w[p];
        for (std::size_t q = 0; q < NumGLP; ++q)
        {
            double eta_q = gl_x[q], w_q = gl_w[q];

            if (theta_left >= epsilon)
            {
                calculateUh(uh, i_left, j, p, q, h_, hu_, hv_);
                pe = {1.0, xi_p-2.0, eta_q, (xi_p-2.0)*eta_q,
                      (xi_p-2.0)*(xi_p-2.0) - 1.0/3.0,
                      eta_q*eta_q - 1.0/3.0};
                h_c = hu_c = hv_c = 0.0;
                for (std::size_t d = 0; d < dimPk; ++d)
                {
                    h_c  += uh[idx_state(i,j,d,0)] * pe[d];
                    hu_c += uh[idx_state(i,j,d,1)] * pe[d];
                    hv_c += uh[idx_state(i,j,d,2)] * pe[d];
                }
                dh = h_ - h_c; dhu = hu_ - hu_c; dhv = hv_ - hv_c;
                norm_left += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, dh, dhu, dhv);
            }

            if (theta_right >= epsilon)
            {
                calculateUh(uh, i_right, j, p, q, h_, hu_, hv_);
                pe = {1.0, xi_p+2.0, eta_q, (xi_p+2.0)*eta_q,
                      (xi_p+2.0)*(xi_p+2.0) - 1.0/3.0,
                      eta_q*eta_q - 1.0/3.0};
                h_c = hu_c = hv_c = 0.0;
                for (std::size_t d = 0; d < dimPk; ++d)
                {
                    h_c  += uh[idx_state(i,j,d,0)] * pe[d];
                    hu_c += uh[idx_state(i,j,d,1)] * pe[d];
                    hv_c += uh[idx_state(i,j,d,2)] * pe[d];
                }
                dh = h_ - h_c; dhu = hu_ - hu_c; dhv = hv_ - hv_c;
                norm_right += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, dh, dhu, dhv);
            }

            if (theta_down > epsilon)
            {
                calculateUh(uh, i, j_down, p, q, h_, hu_, hv_);
                pe = {1.0, xi_p, eta_q-2.0, xi_p*(eta_q-2.0),
                      xi_p*xi_p - 1.0/3.0,
                      (eta_q-2.0)*(eta_q-2.0) - 1.0/3.0};
                h_c = hu_c = hv_c = 0.0;
                for (std::size_t d = 0; d < dimPk; ++d)
                {
                    h_c  += uh[idx_state(i,j,d,0)] * pe[d];
                    hu_c += uh[idx_state(i,j,d,1)] * pe[d];
                    hv_c += uh[idx_state(i,j,d,2)] * pe[d];
                }
                dh = h_ - h_c; dhu = hu_ - hu_c; dhv = hv_ - hv_c;
                norm_down += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, dh, dhu, dhv);
            }

            if (theta_up > epsilon)
            {
                calculateUh(uh, i, j_up, p, q, h_, hu_, hv_);
                pe = {1.0, xi_p, eta_q+2.0, xi_p*(eta_q+2.0),
                      xi_p*xi_p - 1.0/3.0,
                      (eta_q+2.0)*(eta_q+2.0) - 1.0/3.0};
                h_c = hu_c = hv_c = 0.0;
                for (std::size_t d = 0; d < dimPk; ++d)
                {
                    h_c  += uh[idx_state(i,j,d,0)] * pe[d];
                    hu_c += uh[idx_state(i,j,d,1)] * pe[d];
                    hv_c += uh[idx_state(i,j,d,2)] * pe[d];
                }
                dh = h_ - h_c; dhu = hu_ - hu_c; dhv = hv_ - hv_c;
                norm_up += w_p * w_q * HessianNorm(h_bar, u_bar, v_bar, dh, dhu, dhv);
            }
        }
    }
    norm_left  *= 0.25 * hx * hy;
    norm_right *= 0.25 * hx * hy;
    norm_down  *= 0.25 * hx * hy;
    norm_up    *= 0.25 * hx * hy;

    sigma_left  = (theta_left  >= epsilon) ? Ck * norm_left  / theta_left  : 0.0;
    sigma_right = (theta_right >= epsilon) ? Ck * norm_right / theta_right : 0.0;
    sigma_down  = (theta_down  >= epsilon) ? Ck * norm_down  / theta_down  : 0.0;
    sigma_up    = (theta_up    >= epsilon) ? Ck * norm_up    / theta_up    : 0.0;
}

double HessianNorm(double h_center, double u_center, double v_center,
                   double delta_h, double delta_hu, double delta_hv)
{
    return g * h_center * delta_h * delta_h +
           (delta_hu - u_center*delta_h) * (delta_hu - u_center*delta_h) +
           (delta_hv - v_center*delta_h) * (delta_hv - v_center*delta_h);
}