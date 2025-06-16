#include <iostream>
#include <fstream>
#include <sstream>
//#include <memory>
#include <random>
#include <cmath>
#include <vector>
//#include <functional>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <chrono>
#include "barrier.h"


int main() {

    auto start = std::chrono::high_resolution_clock::now();
    try {
        // Read from config file
        BarrierOption option("config.txt");
        double price = option.price();
        std::cout << option.getBarrierType() << " " << option.getOptionType()
                  << " Option Price: " << price << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
//    std::cout<<"Time cost in microsecs: "<< duration << " microseconds" << std::endl;
    std::cout<<"Time cost in secs:      "<< duration*1e-6 << " seconds" << std::endl;
    return 0;
}
