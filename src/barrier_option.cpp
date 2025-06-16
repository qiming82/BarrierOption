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
//#include "util.h"

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
    localVol = [this](double S, double t) { return 0.2 * std::pow(S / 100.0, beta); };
}

// Price the barrier option using Monte Carlo
double BarrierOption::price() const {
    double sumPayoffs = 0.0;

    // Run Monte Carlo simulations
    for (size_t i = 0; i < numSimulations; ++i) {
        sumPayoffs += simulatePath();
    }

    // Calculate average payoff and discount to present value
    return std::exp(-r * T) * (sumPayoffs / numSimulations);
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