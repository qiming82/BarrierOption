#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <chrono>
#include "barrier_option.h"


int main() {

    auto start = std::chrono::high_resolution_clock::now(); // timer
    try {
        // Read from config file
        BarrierOption option("config.txt");

        double price;
        const auto& params = option.getConfigParams();
        if (params.at("exerciseStyle") == "American") {
            std::cout<< "American option:"<< std::endl;
            price = option.priceAmerican();
        } else if (params.at("exerciseStyle") == "European") {
            std::cout<< "European option:"<< std::endl;
            price = option.price();
        } else {
            throw std::runtime_error("Invalid exerciseStyle in config");
        }
        std::cout << option.getBarrierType() << " " << option.getOptionType()
                  << " Option Price: " << price << std::endl;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
        std::cout<<"Time cost (one-time calc) in secs:      "<< duration*1e-6 << " seconds" << std::endl;


// ----------- investigate over a range of S0  --------------------------------
        // Compute Greeks for S0 range and save to CSV
        start = std::chrono::high_resolution_clock::now(); // timer
        std::cout<< "start calculating range ..." << std::endl;
        option.calcGreeksRange(80.0, 120.0, 41, "greeks.csv");
        std::cout << "Greeks data saved to greeks.csv\n";
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
        std::cout<<"Time cost for range test in secs:      "<< duration*1e-6 << " seconds" << std::endl;
// ----------- investigate over a range of S0  --------------------------------

/*
// ---------- simple price, greeks test ------------------------
        double price = option.price();
//        std::cout << option.getBarrierType() << " " << option.getOptionType()
//                  << " Option Price: " << price << std::endl;
        
        BarrierOption::Greeks g = option.calcGreeks();
        std::cout << option.getBarrierType() << " " << option.getOptionType()
                  << " Option Price: " << price << std::endl;
        std::cout << "Delta: " << g.delta << "\nGamma: " << g.gamma << "\nVega: " << g.vega << std::endl;
 
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
        std::cout<<"Time cost (one-time calc) in secs:      "<< duration*1e-6 << " seconds" << std::endl;
// ---------- simple price, greeks test ------------------------


// ----------- convergence test  --------------------------------
        start = std::chrono::high_resolution_clock::now(); // timer
        std::cout << "convergence test: ..." << std::endl;
        option.calcConvgTest(200000, 200000, "convergence.csv");
        std::cout << "Convergence test data saved to convergence.csv\n";
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
        std::cout<<"Time cost for range test in secs:      "<< duration*1e-6 << " seconds" << std::endl;
// ----------- convergence test  --------------------------------
*/

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    

//    std::cout<<"Time cost in microsecs: "<< duration << " microseconds" << std::endl;
//    std::cout<<"Time cost in secs:      "<< duration*1e-6 << " seconds" << std::endl;
    return 0;
}
