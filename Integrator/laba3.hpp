#ifndef UNTITLED_LABA3_HPP
#define UNTITLED_LABA3_HPP

#include <array>
#include <cmath>


template<typename RealType, unsigned int N>
RealType integrate(const RealType& start, const RealType& end){

    std::array<double, N> nodes, weight;

//    if(N == 3){
//        nodes = { std::sqrt(3.0 / 5), 0, -std::sqrt(3.0 / 5) };
//        weight = { 5.0 / 9, 8.0 / 9, 5.0 / 9 };
//    }

    if(N == 4){
        nodes = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
        weight = {0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539};
    }


    RealType result = 0;

    for(int j = 0; j < N; j++){
        result += weight[j] * std::sin((start + end)/2 + (end - start)/2 * nodes[j]);
    }

    return (end - start) / 2 * result;
}

template<typename RealType, unsigned int N>
RealType Integrate(const RealType& start, const RealType& end, const RealType& dx){
    RealType result = 0;

    for(RealType x = start; x <= end; x += dx){
        if (x + dx >= end) {
            result += integrate<RealType, N>(x, end);
            break;
        }
        result += integrate<RealType, N>(x, x + dx);
    }

    return result;
}

#endif //UNTITLED_LABA3_HPP
