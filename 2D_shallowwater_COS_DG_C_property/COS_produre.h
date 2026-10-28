#pragma once

#include "constants.h"

void COS_produre(StateArray& wh, const BottomArray& bh, Real dt);

void calculateBeta(const StateArray& wh, const BottomArray& bh, 
                std::size_t i, std::size_t j, Real& beta_left, Real& beta_right, Real& beta_down, Real& beta_up);

void calculateSigma(const StateArray& wh, const BottomArray& bh, std::size_t i, std::size_t j, 
                Real& sigma_left, Real& sigma_right, Real& sigma_down, Real& sigma_up);

Real HessianNorm(Real h_center, Real u_center, Real v_center, 
                Real delta_h, Real delta_hu, Real delta_hv);