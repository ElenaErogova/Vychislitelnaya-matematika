#ifndef VYCHMATY_MPI_HPP
#define VYCHMATY_MPI_HPP

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename Callable, typename RealType>
decltype(auto) solve(
        const Callable& func,                                             // функция F
        const RealType& tau,                                              // шаг тау
        const typename ArgumentGetter<Callable>::Argument& initialGuess,  // начальное приближение
        const unsigned int nIteration                                     // количество итераций
){
    RealType result = initialGuess;

    for(int i = 0; i < nIteration; i++){
        result -= tau * func(result);
    }

    return result;
}


#endif //VYCHMATY_MPI_HPP
