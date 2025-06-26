#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <random>
#include <cmath>
#include <vector>
#include <functional>
#include <map>
#include <algorithm>
#include "barrier_option.h"
#include <iomanip> // Required for std::fixed and std::setprecision
//#include "util.h"
#include <numeric>
#include <Eigen/SparseLU>

// pi=3.141592653589793;

const double pi = std::atan(1)*4.0;


// Generate correlated standard normal random variables
std::pair<double, double> BarrierOption::generateCorrelatedNormals() const {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double z1 = std::sqrt(-2.0 * std::log(dist(*rng))) * std::cos(2.0 * pi * dist(*rng));
    double z2 = rho * z1 + std::sqrt(1.0 - rho * rho) * 
                std::sqrt(-2.0 * std::log(dist(*rng))) * std::cos(2.0 * pi * dist(*rng));
    return {z1, z2};
}

// Simulate one price and variance path and check barrier condition
double BarrierOption::simulatePath() const {
    double dt = T / numSteps;
    double S = S0;
    double v = v0;
    bool knockedOut = false;

    // Simulate price and variance paths
    for (int i = 0; i < numSteps; ++i) {
        double t = i * dt;
        std::pair<double,double> normals = generateCorrelatedNormals();
        double z_s = normals.first;
        double z_v = normals.second;

//            auto [z_s, z_v] = generateCorrelatedNormals();
        
        // Update variance (Euler discretization)
        v += mu_v(v, t) * dt + sigma_v(v, t) * std::sqrt(dt) * z_v;
        v = std::max(v, 0.0); // Truncate at 0 to avoid negative variance
        
        // Update asset price with local volatility
        double vol = std::sqrt(std::max(v, 0.0) * localVol(S, t) * localVol(S, t));
        S *= std::exp((r - 0.5 * vol * vol) * dt + vol * std::sqrt(dt) * z_s);

        // Check barrier condition
        if (barrierType == BarrierType::UpAndOut && S >= B) {
            knockedOut = true;
            break;
        } else if (barrierType == BarrierType::DownAndOut && S <= B) {
            knockedOut = true;
            break;
        }
    }

    // Calculate payoff if not knocked out
    if (!knockedOut) {
        if (optionType == OptionType::Call) {
            return std::max(S - K, 0.0);
        } else { // Put
            return std::max(K - S, 0.0);
        }
    }
    return 0.0; // Knocked out, payoff is 0
}

// Constructor with config file
BarrierOption::BarrierOption(const std::string& configFile) : rng(new std::mt19937(std::random_device{}())) {
    std::ifstream file(configFile);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open config file: " + configFile);
    }

    std::string line;
//    std::map<std::string, std::string> params;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip empty or comment lines
        std::istringstream iss(line);
        std::string key, value, eq;
        iss >> key >> eq >> value;
        if (eq != "=") throw std::runtime_error("Invalid config format: " + line);
        params[key] = value;
    }
    file.close();

    // Parse parameters
    try {
        S0 = std::stod(params.at("S0"));
        K = std::stod(params.at("K"));
        B = std::stod(params.at("B"));
        T = std::stod(params.at("T"));
        r = std::stod(params.at("r"));
        v0 = std::stod(params.at("v0"));
        kappa = std::stod(params.at("kappa"));
        theta = std::stod(params.at("theta"));
        xi = std::stod(params.at("xi"));
        rho = std::stod(params.at("rho"));
        beta = std::stod(params.at("beta"));
        numSimulations = std::stoi(params.at("numSimulations"));
        numSteps = std::stoi(params.at("numSteps"));

        if(params.find("seed")!=params.end()){
            seed = std::stoul(params.at("seed"));
            rng.reset(new std::mt19937(seed)); 
        }

        if (params.at("optionType") == "Call") {
            optionType = OptionType::Call;
        } else if (params.at("optionType") == "Put") {
            optionType = OptionType::Put;
        } else {
            throw std::runtime_error("Invalid optionType: " + params.at("optionType"));
        }

        if (params.at("barrierType") == "UpAndOut") {
            barrierType = BarrierType::UpAndOut;
        } else if (params.at("barrierType") == "DownAndOut") {
            barrierType = BarrierType::DownAndOut;
        } else {
            throw std::runtime_error("Invalid barrierType: " + params.at("barrierType"));
        }

        if (params.at("exerciseStyle") == "European") {
            exerciseStyle = ExerciseStyle::European;
        } else if (params.at("exerciseStyle") == "American") {
            exerciseStyle = ExerciseStyle::American;
        } else {
            throw std::runtime_error("Invalid exerciseStyle");
        }

        if (params.find("pricingMethod") != params.end()) {
            if (params.at("pricingMethod") == "PDE") {
                pricingMethod = PricingMethod::PDE;
            } else if (params.at("pricingMethod") == "MonteCarlo") {
                pricingMethod = PricingMethod::MonteCarlo;
            } else {
                throw std::runtime_error("Invalid pricingMethod");
            }
        } else {
            pricingMethod = PricingMethod::MonteCarlo; // Default
        }

    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Missing config parameter");
    }

    // Input validation
    if (S0 <= 0 || K <= 0 || B <= 0 || T <= 0 || v0 < 0 || kappa <= 0 || theta <= 0 || xi <= 0) {
        throw std::invalid_argument("Parameters must be positive (except rho)");
    }
    if (std::abs(rho) > 1) {
        throw std::invalid_argument("Correlation must be between -1 and 1");
    }
    if (barrierType == BarrierType::UpAndOut && B < S0) {
        std::cout<<"current S0: "<< S0 << " while barrier is B: "<< B << std::endl;
        throw std::invalid_argument("Up-and-out barrier must be above initial price");
    }
    if (barrierType == BarrierType::DownAndOut && B > S0) {
        throw std::invalid_argument("Down-and-out barrier must be below initial price");
    }

    // Define variance drift and volatility (Heston-like for demonstration)
    mu_v = [this](double v, double t) { return kappa * (theta - v); };
    sigma_v = [this](double v, double t) { return xi * std::sqrt(std::max(v, 0.0)); };

    // Define local volatility (CEV-like for demonstration)
    localVol = [this](double S, double t) { return 0.2 * std::pow(S, beta); };
}

// Price the barrier option using Monte Carlo
double BarrierOption::price() const {

    if (pricingMethod == PricingMethod::PDE)
        {
            std::cout<< "using PDE to solve price (ADI) ..."<< std::endl;
            return priceADI();
        } // isolate PDE solver
    else{
        std::cout<< "using MC to solve price ..."<< std::endl;
        double sumPayoffs = 0.0;

        // Run Monte Carlo simulations
        for (size_t i = 0; i < numSimulations; ++i) {
            sumPayoffs += simulatePath();
        }

        // Calculate average payoff and discount to present value
        return std::exp(-r * T) * (sumPayoffs / numSimulations);
    }
}

// Getter for option type
std::string BarrierOption::getOptionType() const {
    return (optionType == OptionType::Call) ? "Call" : "Put";
}

// Getter for barrier type
std::string BarrierOption::getBarrierType() const {
    return (barrierType == BarrierType::UpAndOut) ? "Up-and-Out" : "Down-and-Out";
}

const std::map<std::string, std::string>& BarrierOption::getConfigParams() const {
    return params;
}

// american option  simulate path overload  -------------------
std::pair<double, std::vector<std::pair<double, double>>> BarrierOption::simulatePath(bool storePath) const {
    double dt = T / numSteps;
    double S = S0;
    double v = v0;
    bool knockedOut = false;
    std::vector<std::pair<double, double>> path(storePath ? numSteps + 1 : 0);
    if (storePath) {
        path[0] = {S, v};
    }

    for (int i = 0; i < numSteps; ++i) {
        double t = i * dt;
        std::pair<double,double> normals = generateCorrelatedNormals();
        double z_s = normals.first;
        double z_v = normals.second;
//        auto [z_s, z_v] = generateCorrelatedNormals();
        
        v += mu_v(v, t) * dt + sigma_v(v, t) * std::sqrt(dt) * z_v;
        v = std::max(v, 0.0);
        
        double vol = std::sqrt(std::max(v, 0.0) * localVol(S, t) * localVol(S, t));
        S *= std::exp((r - 0.5 * vol * vol) * dt + vol * std::sqrt(dt) * z_s);

        if (barrierType == BarrierType::UpAndOut && S >= B) {
            knockedOut = true;
            break;
        } else if (barrierType == BarrierType::DownAndOut && S <= B) {
            knockedOut = true;
            break;
        }

        if (storePath) {
            path[i + 1] = {S, v};
        }
    }

    double payoff = 0.0;
    if (!knockedOut) {
        if (optionType == OptionType::Call) {
            payoff = std::max(S - K, 0.0);
        } else {
            payoff = std::max(K - S, 0.0);
        }
    }
    return {payoff, path};
}


// american option   langstaff schwartz monte carlo  --------------------
double BarrierOption::priceAmerican() const {
    double dt = T / numSteps;
    std::vector<std::vector<std::pair<double, double>>> paths(numSimulations);
    std::vector<double> cashFlows(numSimulations, 0.0);

    // Simulate paths
    for (int i = 0; i < numSimulations; ++i) {
        auto result = simulatePath(true);
        paths[i] = result.second;
        cashFlows[i] = result.first;
//        auto [payoff, path] = simulatePath();
//        paths[i] = path;
//        cashFlows[i] = payoff;
    }

    // Backward induction using Longstaff-Schwartz
    for (int t = numSteps - 1; t >= 1; --t) {
        std::vector<double> X, Y;
        std::vector<int> indices;
        for (int i = 0; i < numSimulations; ++i) {
            double S = paths[i][t].first;
            if (barrierType == BarrierType::UpAndOut && S >= B) continue;
            if (barrierType == BarrierType::DownAndOut && S <= B) continue;
            double exercise = (optionType == OptionType::Call) ? std::max(S - K, 0.0) : std::max(K - S, 0.0);
            if (exercise > 0) {
                X.push_back(S);
                Y.push_back(cashFlows[i] * std::exp(-r * dt * (numSteps - t)));
                indices.push_back(i);
            }
        }

        if (X.empty()) continue;

        // Linear regression: E[Y] = a + b*S + c*S^2
        double sumX = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0, sumY = 0, sumXY = 0, sumX2Y = 0;
        int n = X.size();
        for (int i = 0; i < n; ++i) {
            double x = X[i], y = Y[i];
            sumX += x;
            sumX2 += x * x;
            sumX3 += x * x * x;
            sumX4 += x * x * x * x;
            sumY += y;
            sumXY += x * y;
            sumX2Y += x * x * y;
        }

        // Solve [n, sumX, sumX2; sumX, sumX2, sumX3; sumX2, sumX3, sumX4][a, b, c] = [sumY, sumXY, sumX2Y]
        double det = n * (sumX2 * sumX4 - sumX3 * sumX3) - sumX * (sumX * sumX4 - sumX2 * sumX3) +
                     sumX2 * (sumX * sumX3 - sumX2 * sumX2);
        if (std::abs(det) < 1e-10) continue;

        double a = (sumY * (sumX2 * sumX4 - sumX3 * sumX3) - sumX * (sumXY * sumX4 - sumX2Y * sumX3) +
                    sumX2 * (sumXY * sumX3 - sumX2Y * sumX2)) / det;
        double b = (n * (sumXY * sumX4 - sumX2Y * sumX3) - sumY * (sumX * sumX4 - sumX2 * sumX3) +
                    sumX2 * (sumX * sumX2Y - sumX2 * sumXY)) / det;
        double c = (n * (sumX2 * sumX2Y - sumX3 * sumXY) - sumX * (sumX * sumX2Y - sumX2 * sumXY) +
                    sumY * (sumX * sumX3 - sumX2 * sumX2)) / det;

        // Update cash flows
        for (int i : indices) {
            double S = paths[i][t].first;
            double exercise = (optionType == OptionType::Call) ? std::max(S - K, 0.0) : std::max(K - S, 0.0);
            double continuation = a + b * S + c * S * S;
            if (exercise > continuation) {
                cashFlows[i] = exercise * std::exp(r * dt * (numSteps - t));
            }
        }
    }

    double sumCashFlows = std::accumulate(cashFlows.begin(), cashFlows.end(), 0.0);
    return sumCashFlows / numSimulations;
}
// --- pde solver ADI
double BarrierOption::priceADI() const {
    double epsiln = 1e-8;
    std::cout<<std::fixed << std::setprecision(8);
    // Grid parameters
    double S_max = 1.5 * B; // 180 for B=120
    double v_max = 0.4;     // 10x theta
    int Ns = 200;
    int Nv = 200;
    int Nt = 100;
    double ds = S_max / Ns;
    double dv = v_max / Nv;
    double dt = T / Nt;
//    double theta = 0.5; // Douglas-Rachford parameter
//    Nt = 1;

    // Grid
    std::vector<double> S(Ns + 1);
    std::vector<double> v(Nv + 1);
    for (int i = 0; i <= Ns; ++i) S[i] = i * ds;
    for (int j = 0; j <= Nv; ++j) v[j] = j * dv;

    // Value grid
    Eigen::VectorXd V((Ns + 1) * (Nv + 1)), V_star((Ns + 1) * (Nv + 1));
    std::vector<Eigen::Triplet<double>> triplets;

    // Terminal condition (t = T)
    for (int i = 0; i <= Ns; ++i) {
        for (int j = 0; j <= Nv; ++j) {
            int k = i * (Nv + 1) + j;
            if (S[i] < B) {
                V(k) = optionType == OptionType::Call ? std::max(S[i] - K, epsiln) : std::max(K - S[i], epsiln);
            } else {
                V(k) = 0.0;
            }
        }
    }

    // Time stepping (backward from T to 0)
    for (int n = 0; n < Nt; ++n) {
        double t_n = T - n * dt; // Corrected time, tau = T -t

        // First half-step: implicit in S
//        Eigen::SparseMatrix<double> A_S(Ns - 1, Ns - 1);
//        Eigen::VectorXd rhs_S(Ns - 1);

        Eigen::SparseMatrix<double> A_S(Ns+1, Ns+1);
        Eigen::VectorXd rhs_S(Ns+1);

        for (int j = 1; j < Nv; ++j) {
            triplets.clear();
            rhs_S.setZero();
//            std::cout<< S[0]<<", "<<S[1]<<", "<<S[Ns]<<std::endl;

            double aa = 0.0; //0.5 * v[j] * std::pow(localVol(S[0], t_n), 2) * S[0] * S[0] / (ds * ds);
            double bb = r * S[0] / (2.0 * ds);
            triplets.push_back(Eigen::Triplet<double>(0, 0, 1.0 + 0.5* dt*(2.0*aa + r)));
            triplets.push_back(Eigen::Triplet<double>(0, 1, -0.5 * dt * (aa + bb)));

//            std::cout<< "diag 0" <<", "<< 1.0 + 0.5* dt*(2.0*aa + r) <<std::endl;
//            std::cout<< "diag 01" <<", "<< -0.5* dt*(aa + bb) <<std::endl;

            for (int i = 1; i < Ns; ++i) {
                double aa = 0.5 * v[j] * std::pow(localVol(S[i], t_n), 2) * S[i] * S[i] / (ds * ds);
                double bb = r * S[i] / (2.0 * ds);
                double cc = 0.5 * xi * xi * v[j]; // / (dv * dv);
                double dd = kappa * (theta - v[j]); // / (2 * dv);
                double cross = rho * xi * v[j] * localVol(S[i], t_n) * S[i]; // / (4 * ds * dv);
                int k = i * (Nv + 1) + j;

//                std::cout<< "diag i,i" <<", "<< i <<", "<<1.0 + 0.5* dt*(2.0*aa + r) <<std::endl;
  
                // Implicit S terms: I + 0.5 * dt * (-L_S + r)
//                triplets.push_back(Eigen::Triplet<double>(i - 1, i - 1, 1.0 + 0.5* dt*(2.0*aa + r)));
//                triplets.push_back(Eigen::Triplet<double>(i - 1, i, -0.5 * dt * (aa + bb)));
//                triplets.push_back(Eigen::Triplet<double>(i - 1, i - 2, 0.5 * dt * (-aa + bb)));

                triplets.push_back(Eigen::Triplet<double>(i, i, 1.0 + 0.5* dt*(2.0*aa + r)));
                triplets.push_back(Eigen::Triplet<double>(i, i+1, -0.5 * dt * (aa + bb))); // upper diag
                triplets.push_back(Eigen::Triplet<double>(i, i-1,  0.5 * dt * (-aa + bb)));  // lower diag

                // Explicit terms rhs
                double v1 = (V(k+1) - V(k-1))/(2.0*dv);
                double v2 = (V(k+1) - 2.0*V(k) + V(k-1))/dv/dv;
                double dvcross = (V(k+ Nv + 2) - V(k+ Nv) - V(k-Nv) + V(k-Nv-2))/(4.0*ds*dv);

                rhs_S(i) = V(k) + dt/2.0*(dd*v1 + cc*v2 + cross*dvcross);
            }

//            int i = Ns;
            aa = 0.5 * v[j] * std::pow(localVol(S[Ns], t_n), 2) * S[Ns] * S[Ns] / (ds * ds);
            bb = r * S[Ns] / (2.0 * ds);
//            std::cout<< aa <<", "<<bb<<", "<<std::endl;

            triplets.push_back(Eigen::Triplet<double>(Ns, Ns, 1.0 + 0.5* dt*(2.0*aa + r)));
            triplets.push_back(Eigen::Triplet<double>(Ns, Ns-1,  0.5 * dt * (-aa + bb)));

//            rhs_S(0) = 0.0;
            rhs_S(Ns) = S_max - K * std::exp(r * (T - t_n + 0.5 * dt));

//            for (int i = 0; i <= Ns; ++i) {
//                std::cout <<i<<", rhs: "<< rhs_S(i)<<std::endl;
//            }


            // Boundary conditions for S
            for (int i = 0; i <= Ns; ++i) {
                int k = i * (Nv + 1) + j;
                if (i == 0 || (barrierType == BarrierType::UpAndOut && S[i] -B>= epsiln)) {
                    V_star(k) = 0.0;
                }
            }

            // Solve tridiagonal system
            A_S.setFromTriplets(triplets.begin(), triplets.end());
/*
            for (int i = 1; i <Ns; ++i) {
                std::cout << "(i,j) " <<i<<", "<<j<<", "<< A_S.coeff(i-1,i)<<", "<< A_S.coeff(i,i)<<", "<< A_S.coeff(i,i+1)<<", "<<std::endl;
            }
            std::cout<<"A(0,0)"<< A_S.coeff(0,0)<<", "<< "A(0,1)"<< A_S.coeff(0,1)<<std::endl;
*/

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_S;
            solver_S.compute(A_S);
//            std::cout<< "called S step solver, appears ok ..." << std::endl;
            if (solver_S.info() != Eigen::Success) {
                throw std::runtime_error("Sparse LU decomposition failed in S-step");
            }
            Eigen::VectorXd V_S_new = solver_S.solve(rhs_S);
            for (int i = 0; i <= Ns; ++i) {
//                std::cout<< V_S_new(i)<<std::endl;
                V_star(i * (Nv + 1) + j) = V_S_new(i);
            }
        }

        // Apply v=0 and v=v_max boundaries
        for (int i = 0; i <= Ns; ++i) {
            int k0 = i * (Nv + 1);
            int kNv = i * (Nv + 1) + Nv;
            double S_future = S[i] * std::exp(r * (T - t_n + 0.5 * dt)); // t_{n+1/2} = t_n - dt/2
            V_star(k0) = optionType == OptionType::Call ? 
                         (S_future < B ? std::max(S_future - K, 0.0) : 0.0) :
                         (S_future < B ? std::max(K - S_future, 0.0) : 0.0);
            V_star(kNv) = V_star(kNv - 1); // Neumann
        }

//--------------------------------------------------------------------------------------------------------------
        // Second half-step: implicit in v
        Eigen::SparseMatrix<double> A_v(Nv+1, Nv+1);
        Eigen::VectorXd rhs_v(Nv+1);
        for (int i = 1; i < Ns; ++i) {
            triplets.clear();
            rhs_v.setZero();

            double cc = 0.5 * xi * xi * v[0] / (dv * dv);
            double dd = kappa * (theta - v[0]) / (2.0 * dv);

            triplets.push_back(Eigen::Triplet<double>(0, 0, 1.0 + 0.5* dt*(2.0*cc + r)));
            triplets.push_back(Eigen::Triplet<double>(0, 1, -0.5 * dt * (cc + dd)));

            for (int j = 1; j < Nv; ++j) {
                double aa = 0.5 * v[j] * std::pow(localVol(S[i], t_n), 2) * S[i] * S[i]; // / (ds * ds);
                double bb = r * S[i]; // / (2 * ds);
                double cc = 0.5 * xi * xi * v[j] / (dv * dv);
                double dd = kappa * (theta - v[j]) / (2.0 * dv);
                double cross = rho * xi * v[j] * localVol(S[i], t_n) * S[i]; /// (4.0 * ds * dv);
                int k = i * (Nv + 1) + j;

                // Implicit v terms: I + theta * dt * (L_v - r)
                triplets.push_back(Eigen::Triplet<double>(j, j, 1.0 + 0.5*dt * (2.0*cc + r)));
                triplets.push_back(Eigen::Triplet<double>(j, j+1, -0.5* dt * (cc + dd)));
                triplets.push_back(Eigen::Triplet<double>(j, j-1, -0.5* dt * (cc - dd)));

                // Explicit terms rhs
                double s1 = (V_star(k+Nv+1) - V_star(k-Nv-1))/(2.0*ds);
                double s2 = (V_star(k+Nv+1) - 2.0*V_star(k) + V_star(k-Nv-1))/ds/ds;
                double dvcross = (V_star(k+ Nv + 2) - V_star(k+ Nv) - V_star(k-Nv) + V_star(k-Nv-2))/(4.0*ds*dv);

                rhs_v(j) = V(k) + dt/2.0*(bb*s1 + aa*s2 + cross*dvcross);
            }

            cc = 0.5*xi*xi * v[Nv] / (dv * dv);
            dd = kappa * (theta - v[Nv]) / (2.0 * dv);
//            std::cout<< aa <<", "<<bb<<", "<<std::endl;

            triplets.push_back(Eigen::Triplet<double>(Nv, Nv, 1.0 + 0.5* dt*(2.0*cc + r)));
            triplets.push_back(Eigen::Triplet<double>(Nv, Nv-1, -0.5* dt * (cc -dd)));

            rhs_v[0] = V_star(i*(Nv+1));
            rhs_v(Nv) = V_star(i*(Nv+1)+Nv);

            // Solve tridiagonal system
            A_v.setFromTriplets(triplets.begin(), triplets.end());
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_v;
            solver_v.compute(A_v);

//            std::cout<< "called S step solver, appears ok ..." << std::endl;

            if (solver_v.info() != Eigen::Success) {
                throw std::runtime_error("Sparse LU decomposition failed in v-step");
            }
            Eigen::VectorXd V_v_new = solver_v.solve(rhs_v);
            for (int j = 1; j < Nv; ++j) {
                V(i * (Nv + 1) + j) = V_v_new(j - 1);
            }
        }

        // Apply boundaries after v-step
        for (int i = 0; i <= Ns; ++i) {
            int k0 = i * (Nv + 1);
            int kNv = i * (Nv + 1) + Nv;
            double S_future = S[i] * std::exp(r * (T - t_n)); // t_n
            V(k0) = optionType == OptionType::Call ? 
                    (S_future -B < epsiln ? std::max(S_future - K, 0.0) : 0.0) :
                    (S_future - B<epsiln ? std::max(K - S_future, 0.0) : 0.0);
            V(kNv) = V(kNv - 1); // Neumann
            if (i == 0 || (barrierType == BarrierType::UpAndOut && S[i] >= B)) {
                for (int j = 0; j <= Nv; ++j) {
                    V(i * (Nv + 1) + j) = 0.0;
                }
            }
        }
    }

    // Interpolate at (S0, v0, t=0)
    int i = static_cast<int>(S0 / ds);
    int j = static_cast<int>(v0 / dv);
    if (i >= Ns || j >= Nv) return 0.0;
    double s_frac = (S0 - S[i]) / ds;
    double v_frac = (v0 - v[j]) / dv;
    int k00 = i * (Nv + 1) + j;
    int k10 = (i + 1) * (Nv + 1) + j;
    int k01 = i * (Nv + 1) + (j + 1);
    int k11 = (i + 1) * (Nv + 1) + (j + 1);
    double V00 = V[k00], V10 = V[k10], V01 = V[k01], V11 = V[k11];
    return (1 - s_frac) * (1 - v_frac) * V00 + s_frac * (1 - v_frac) * V10 + 
           (1 - s_frac) * v_frac * V01 + s_frac * v_frac * V11;
}
