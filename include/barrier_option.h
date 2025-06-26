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
#include <Eigen/Sparse>

// pi=3.141592653589793;
typedef unsigned long u32;
// Enum for option type
enum class OptionType { Call, Put };

// Enum for barrier type
enum class BarrierType { UpAndOut, DownAndOut };

// Enum for exercise style
enum class ExerciseStyle { European, American };

enum class PricingMethod { PDE, MonteCarlo }; // New enum
//const double pi = std::atan(1)*4.0;

// Class to encapsulate barrier option pricing with general SLV
class BarrierOption {
public:
    // add calculation of greeks 
    struct Greeks {
        double delta;
        double gamma;
        double vega;
    };

    // Constructor with config file
    BarrierOption(const std::string& configFile);

    // Price the barrier option using Monte Carlo
    double price() const;

    // american option
    double priceAmerican() const;

    // pde solver, crank nicolson
//    double pricePDE() const;
    double priceADI() const;

    // Getter for option type
    std::string getOptionType() const;

    // Getter for barrier type
    std::string getBarrierType() const;

    const std::map<std::string, std::string>& getConfigParams() const; // Access config

    // greeks related
    Greeks calcGreeks(double epsilonS = 0.01, double epsilonV = 0.0001) const;
    void calcGreeksRange(double S_min, double S_max, int numPts, const std::string& outputFile) const;

    // convergence
    void calcConvgTest(int minSims, int stepSize, const std::string& outputFile) const;

private:
    double S0; // Initial stock price
    double K;  // Strike price
    double B;  // Barrier level
    double T;  // Time to maturity (years)
    double r;  // Risk-free rate
    double v0; // Initial variance
    double kappa; // Mean reversion rate (for Heston-like variance)
    double theta; // Long-term variance
    double xi;    // Volatility of variance
    double rho;   // Correlation between asset and variance
    double beta;  // exponent in local volatility model
    u32 seed;
    int numSimulations; // Number of Monte Carlo simulations
    int numSteps;      // Number of time steps

    OptionType optionType; // Call or Put
    BarrierType barrierType; // Up-and-out or Down-and-out
    ExerciseStyle exerciseStyle; // European or American
    PricingMethod pricingMethod; // New field
    std::unique_ptr<std::mt19937> rng; // Smart pointer for random number generator
//    std::mt19937* rng;

    std::map<std::string, std::string> params; // Store config parameters

    // Function types for variance drift and volatility
    using VarianceDrift = std::function<double(double, double)>;
    using VarianceVol = std::function<double(double, double)>;
    using LocalVol = std::function<double(double, double)>;

    VarianceDrift mu_v;   // Variance drift function
    VarianceVol sigma_v;  // Variance volatility function
    LocalVol localVol;    // Local volatility function

    // Generate correlated standard normal random variables
    std::pair<double, double> generateCorrelatedNormals() const;

    // Simulate one price and variance path and check barrier condition
    double simulatePath() const;

    // Simulate one price and variance path, return path (American)
    std::pair<double, std::vector<std::pair<double, double>>> simulatePath(bool storePath) const;

    // auxi fun to compute price w changed parameters
    double priceWithParams(double newS0, double newV0) const;

    // auxi fun to compute monte carlon convergence test
    double priceWithSimulations(int sims) const;


};
