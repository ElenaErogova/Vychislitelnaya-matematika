#ifndef VYCHMATY_NEWTON_HPP
#define VYCHMATY_NEWTON_HPP

#include <array>
#include <iostream>
#include <cmath>

double keplerSolver(double ecc, double meanAnomaly, const unsigned int maxIter, double tol){
    double eccAnomaly_temp = meanAnomaly;
    double eccAnomaly = 0;

    for(int i = 1; i < maxIter+1; i++) {
        eccAnomaly = eccAnomaly_temp - (eccAnomaly_temp - ecc * sin(eccAnomaly_temp) - meanAnomaly) /
                                            (1 - ecc * cos(eccAnomaly_temp));
        std::cout << std::log10(std::abs(eccAnomaly - 1.5853138613542081536)) << std::endl;
        if (std::abs(eccAnomaly - eccAnomaly_temp) < tol) {
            return eccAnomaly;
        }
        if (i == maxIter) {
            if (std::abs(eccAnomaly - eccAnomaly_temp) > tol) {
                std::cout << "The number of iterations is not enough" << std::endl;
            }
        }
        eccAnomaly_temp = eccAnomaly;
    }
}

#endif //VYCHMATY_NEWTON_HPP
