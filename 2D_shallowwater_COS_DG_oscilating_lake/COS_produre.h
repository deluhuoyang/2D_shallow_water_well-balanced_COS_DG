#pragma once

#include "constants.h"

void COS_produre(StateArray& uh, const BottomArray& bh, double dt);

void calculateBeta(const StateArray& uh,
                   std::size_t i, std::size_t j,
                   double& beta_left, double& beta_right,
                   double& beta_down, double& beta_up);

void calculateSigma(const StateArray& uh, const BottomArray& bh,
                    std::size_t i, std::size_t j,
                    double& sigma_left, double& sigma_right,
                    double& sigma_down, double& sigma_up);

double HessianNorm(double h_center, double u_center, double v_center,
                   double delta_h, double delta_hu, double delta_hv);