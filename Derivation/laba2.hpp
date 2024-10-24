#ifndef UNTITLED_LABA2_HPP
#define UNTITLED_LABA2_HPP

#include <array>
#include "eigen-3.4.0/eigen-3.4.0/Eigen/Dense"

template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N, unsigned int L>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {
    Eigen::VectorX<RealType> b(N + 1);
    b.setZero();
    b[L] = 1;

    Eigen::MatrixX<RealType> A(N + 1, N + 1);

    A.setZero();

    for(int i = 0; i < N+1; i++){
        A(0, i) = 1;
    }

    for(int i = 1; i < N+1; i++){
        for(int j = 1; j < N+1; j++){
            A(j, i) = A(j-1, i) * points[i-1]/j;
        }
    }

    Eigen::VectorX<RealType> x = A.lu().solve(b);

    DerivativeCoef<double,N> coef;

    coef.centralCoef =  x[0];
    for(int i = 0; i < N; i++){
        coef.otherCoefs[i] = x[i+1];
    }

    return coef;
}

#endif //UNTITLED_LABA2_HPP
