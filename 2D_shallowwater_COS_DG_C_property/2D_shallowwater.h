#pragma once

Real b(Real x, Real y);
Real h (Real x, Real y);
Real hu(Real,   Real);
Real hv(Real,   Real);

Real f1(Real h, Real hu, Real hv);
Real f2(Real h, Real hu, Real hv);  
Real f3(Real h, Real hu, Real hv);
Real g1(Real h, Real hu, Real hv);
Real g2(Real h, Real hu, Real hv);
Real g3(Real h, Real hu, Real hv);

void getBasic(Array3d<Real, dimPk, NumGLP, NumGLP>& phi,
              Array3d<Real, dimPk, NumGLL, NumGLL>& phiGLL, 
              Array3d<Real, dimPk, NumGLP, NumGLP>& phi_xi,
              Array3d<Real, dimPk, NumGLP, NumGLP>& phi_eta);

void getInterface(Array2d<Real, dimPk, NumGLP>& phixL,
                  Array2d<Real, dimPk, NumGLP>& phixR,
                  Array2d<Real, dimPk, NumGLP>& phiyL,
                  Array2d<Real, dimPk, NumGLP>& phiyR);

void getProjection(StateArray& uh, BottomArray& bh);

void runRK3(StateArray& uh, const BottomArray& bh, Real t_end);

void Lh(const StateArray& uh, const BottomArray& bh, StateArray& du, 
        Real alpha_1, Real alpha_2);

void positivityPreservingLimiter(StateArray& uh);        

Real LF_Flux_f(std::size_t n, std::array<Real, NumEq>& uR_left, 
              std::array<Real, NumEq>& uR_right, Real alpha);

Real LF_Flux_g(std::size_t n, std::array<Real, NumEq>& uR_left, 
              std::array<Real, NumEq>& uR_right, Real alpha);

void calculateUh(const StateArray& uh, std::size_t i, std::size_t j,
                       std::size_t p, std::size_t q, Real& h, Real& hu, Real& hv);

void calculateBhDeriv(const BottomArray& bh, std::size_t i, std::size_t j,
                       std::size_t p, std::size_t q, Real& b_xi, Real& b_eta);

void calculateUxInterface(const StateArray& uh, std::size_t i, std::size_t j,
                        std::size_t n, std::size_t q, Real& uL_left, Real& uL_right, 
                        Real& uR_left, Real& uR_right);

void calculateUyInterface(const StateArray& uh, std::size_t i, std::size_t j,
                        std::size_t n, std::size_t q, Real& uD_down, Real& uD_up, 
                        Real& uU_down, Real& uU_up);
 