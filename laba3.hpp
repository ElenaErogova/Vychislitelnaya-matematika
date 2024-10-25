#ifndef UNTITLED_LABA3_HPP
#define UNTITLED_LABA3_HPP

#include <array>
#include <cmath>

template<typename RealType, unsigned int N>
std::array<RealType, N> Transfer(std::array<RealType, N>& nodes, RealType start, RealType end) {
    std::array<RealType, N> points;

    for(int i = 0; i < N; i++){
        points[i] = (start + end)/2 + (end - start) * nodes[i]/2;
    }

    return points;
}

template<typename RealType, unsigned int N>
RealType Integrator(std::array<RealType, N>& points, std::array<RealType, N>& weight, RealType start, RealType end){
    RealType result = 0;

    for(int j = 0; j < N; j++){
        result += (end - start) * weight[j] * std::sin(points[j])/2;
    }

    return result;
}

#endif //UNTITLED_LABA3_HPP
