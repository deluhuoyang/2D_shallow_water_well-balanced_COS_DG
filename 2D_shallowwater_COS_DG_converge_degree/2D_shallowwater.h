#pragma once

#include "constants.h"
#include <array>
#include <vector>

double b (double x, double y);
double h (double x, double y);
double hu(double,    double);
double hv(double,    double);

double f1(double h, double hu, double hv);
double f2(double h, double hu, double hv);
double f3(double h, double hu, double hv);
double g1(double h, double hu, double hv);
double g2(double h, double hu, double hv);
double g3(double h, double hu, double hv);

void getBasic(std::vector<double>& phi,
              std::vector<double>& phiGLL,
              std::vector<double>& phi_xi,
              std::vector<double>& phi_eta);

void getInterface(std::vector<double>& phixL,
                  std::vector<double>& phixR,
                  std::vector<double>& phiyL,
                  std::vector<double>& phiyR);

void getProjection(StateArray& uh, BottomArray& bh);

void runRK3(StateArray& uh, const BottomArray& bh, double t_end);

void Lh(const StateArray&  uh,
        const BottomArray& bh,
        StateArray&        du,
        double             alpha_1,
        double             alpha_2);

void positivityPreservingLimiter(StateArray& uh);

double LF_Flux_f(std::size_t n, std::array<double, NumEq>& uR_left,
                 std::array<double, NumEq>& uR_right, double alpha);
double LF_Flux_g(std::size_t n, std::array<double, NumEq>& uR_left,
                 std::array<double, NumEq>& uR_right, double alpha);

void calculateUh(const StateArray& uh, std::size_t i, std::size_t j,
                 std::size_t p, std::size_t q, double& h, double& hu, double& hv);

void calculateBhDeriv(const BottomArray& bh, std::size_t i, std::size_t j,
                      std::size_t p, std::size_t q, double& b_xi, double& b_eta);

void calculateUxInterface(const StateArray& uh, std::size_t i, std::size_t j,
                          std::size_t n, std::size_t q,
                          double& uL_left, double& uL_right,
                          double& uR_left, double& uR_right);

void calculateUyInterface(const StateArray& uh, std::size_t i, std::size_t j,
                          std::size_t n, std::size_t q,
                          double& uD_down, double& uD_up,
                          double& uU_down, double& uU_up);