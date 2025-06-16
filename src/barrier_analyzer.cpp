#include <iostream>
#include <fstream>
//#include <sstream>
//#include <random>
#include <cmath>
#include <vector>
//#include <stdexcept>
#include <map>
#include <algorithm>
#include "barrier_option.h"
#include "util.h"


// Helper to create a config map with modified parameters
std::map<std::string, std::string> getModifiedConfig(const BarrierOption& option, double newS0, double newV0) {
    auto params = option.getConfigParams(); // Copy existing config
    params["S0"] = std::to_string(newS0);
    params["v0"] = std::to_string(newV0);
    return params;
}

// Write config map to a temporary file
std::string writeTempConfig(const std::map<std::string, std::string>& params) {
    std::string tempFile = "temp_config.txt";
    std::ofstream out(tempFile);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open temp config file");
    }
    for (const auto& p : params) {
        out << p.first << " = " << p.second << "\n";
    }
    out.close();
    return tempFile;
}

double BarrierOption::priceWithParams(double newS0, double newV0) const {
    auto params = getModifiedConfig(*this, newS0, newV0);
    std::string tempFile = writeTempConfig(params);
    BarrierOption temp(tempFile);
    return temp.price();
}

BarrierOption::Greeks BarrierOption::calcGreeks(double epsilonS, double epsilonV) const {
    Greeks g;
    double priceBase = price();
    
    double priceUp = priceWithParams(S0 + epsilonS, v0);
    double priceDown = priceWithParams(S0 - epsilonS, v0);
    g.delta = (priceUp - priceDown) / (2 * epsilonS);
    g.gamma = (priceUp - 2 * priceBase + priceDown) / (epsilonS * epsilonS);
    
    double priceVegaUp = priceWithParams(S0, v0 + epsilonV);
    double priceVegaDown = priceWithParams(S0, v0 - epsilonV);
    g.vega = (priceVegaUp - priceVegaDown) / (2 * epsilonV);
    
    return g;
}

void BarrierOption::calcGreeksRange(double S_min, double S_max, int numPoints, const std::string& outputFile) const {
    std::ofstream out(outputFile);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + outputFile);
    }
    out << "StockPrice,Price,Delta,Gamma,Vega\n";
    
    double v = std::stod(params.at("v0"));

    double dS = (S_max - S_min) / (numPoints - 1);
    for (size_t i = 0; i < numPoints-1; ++i) {
        double S = S_min + i * dS;
        auto params = getModifiedConfig(*this,S,v);
        std::string tempFile = writeTempConfig(params);
        BarrierOption temp(tempFile);
//        BarrierOption temp = *this;
//        temp.S0 = S;
        double price = temp.price();
        Greeks g = temp.calcGreeks();
        out << S << "," << price << "," << g.delta << "," << g.gamma << "," << g.vega << "\n";
        std::cout << "computing range ... step: " << i << std::endl;
    }
    out.close();
}


