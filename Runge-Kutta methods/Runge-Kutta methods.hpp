#ifndef VYCHMATY_RUNGE_KUTTA_METHODS_HPP
#define VYCHMATY_RUNGE_KUTTA_METHODS_HPP
#include <array>
#include <vector>
#include <iostream>


struct RK4Table{
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = {{{0, 0, 0, 0},
                                                                              {0.5, 0, 0, 0},
                                                                              {0, 0.5, 0, 0},
                                                                              {0, 0, 1, 0}}};
    static constexpr std::array<double, stages> cColumn = {0, 0.5, 0.5, 1};
    static constexpr std::array<double, stages> bString = {1./6, 1./3, 1./3, 1./6};
};

template<typename Table, typename RHS>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
        const typename RHS::StateAndArg& initialState,
        const typename RHS::Argument& endTime,
        double step,
        const RHS& rhs
){
    std::vector<typename RHS::StateAndArg> result;
    std::vector<typename RHS::StateAndArg> coef(Table::stages);

    coef[0].state = rhs.calc(initialState);

    for(int i = 1; i < Table::cColumn.size(); i++) {
        coef[i].state = rhs.calc({initialState.state + step * Table::table[i][i-1] * coef[i-1].state, initialState.arg + step * Table::cColumn[i]});
    }

    result.push_back(initialState);
    result.push_back({initialState.state + step * (Table::bString[0] * coef[0].state + Table::bString[1] * coef[1].state + Table::bString[2] * coef[2].state + Table::bString[3] * coef[3].state), initialState.arg + step});

    for(double time = initialState.arg + step; time < endTime; time += step){
        coef[0].state = rhs.calc({result[result.size()-1].state, result[result.size()-1].arg});

        for(int i = 1; i < Table::cColumn.size(); i++) {
            coef[i].state = rhs.calc({result[result.size()-1].state + step * Table::table[i][i-1] * coef[i-1].state, result[result.size()-1].arg + step * Table::cColumn[i]});
        }

        result.push_back({result[result.size()-1].state + step * (Table::bString[0] * coef[0].state + Table::bString[1] * coef[1].state + Table::bString[2] * coef[2].state + Table::bString[3] * coef[3].state), time + step});
    }

    return result;
}

#endif //VYCHMATY_RUNGE_KUTTA_METHODS_HPP
