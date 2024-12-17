#ifndef VYCHMATY_MULTI_STEP_HPP
#define VYCHMATY_MULTI_STEP_HPP
#include <array>
#include <vector>
#include "../Runge-Kutta methods/Runge-Kutta methods.hpp"

struct BDF4{
    static constexpr unsigned int size = 4;
    static constexpr std::array<double, size> alpha = {11./6, -3, 1.5, -1./3};
};

struct IntegrationParameters{
    double step;  // шаг интегрирования
    double epsilon;  // точность решения нелинейного уравнения
    int maxIter;  // максимальное количество итераций для решения нелинейного уравнения
};

template<typename BDF,  typename RKTable, typename RHS>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
        const typename RHS::StateAndArg& initialState,
        const typename RHS::Argument& endTime,
        const IntegrationParameters& parameters,
        const RHS& rhs
){
    std::vector<typename RHS::StateAndArg> values = integrate<RKTable>(initialState, 2 * parameters.step, parameters.step, rhs);

    typename RHS::StateAndArg initialResult;
    typename RHS::StateAndArg Result;
    std::vector<typename RHS::StateAndArg> RESULT = values;

    for(double time = initialState.arg + 3 * parameters.step; time < endTime; time += parameters.step) {
        initialResult.state = 3 * (parameters.step * rhs.calc(values[2]) - 0.5 * values[2].state
                            + values[1].state - 1./6 * values[0].state);
        initialResult.arg = values[2].arg + parameters.step;

        Result = initialResult;

        for (int i = 0; i < parameters.maxIter; i++) {
            Result.state = (parameters.step * rhs.calc(Result) - BDF::alpha[1] * values[2].state - BDF::alpha[2] * values[1].state - BDF::alpha[3] * values[0].state) / BDF::alpha[0];

            if ((Result.state - initialResult.state).norm() < parameters.epsilon) {
                initialResult.state = Result.state;
                break;
            }

            initialResult.state = Result.state;
        }

        for(int i = 0; i < 2; i++){
            values[i].state = values[i+1].state;
            values[i].arg += parameters.step;
        }
        values[2].state = Result.state;
        values[2].arg = Result.arg;

        RESULT.push_back(Result);
    }

    return RESULT;
}

#endif //VYCHMATY_MULTI_STEP_HPP
