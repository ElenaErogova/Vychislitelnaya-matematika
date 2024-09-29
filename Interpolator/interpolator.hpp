#ifndef VYCHISLITELNAYA_MATEMATIKA_INTERPOLATOR_HPP
#define VYCHISLITELNAYA_MATEMATIKA_INTERPOLATOR_HPP

#include <array>

template <typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
private:
    std::array<yType, N> difference {}; // разделённые разности
    std::array<xType, N> result {};
    double resultt = 0;
public:

    NewtonInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values) noexcept{
        for(int i = 0; i <N; i++){
            difference[i] = values[i];
        }

        double timevalue1 = 0;
        double timevalue2 = 0;

        for(int i = 1; i < N; i++){
            for(int j = i; j < N ; j++){
                if(j==i) {
                    timevalue1 = difference[j];
                    difference[j] = (difference[j] - difference[j - 1]) / (points[j] - points[j - i]);
                }
                else{
                    timevalue2 = difference[j];
                    difference[j] = (difference[j] - timevalue1) / (points[j] - points[j - i]);
                    timevalue1 = timevalue2;
                }
            }
        }
    }
    yType interpolate(const xType& x, const std::array<xType, N> &points) noexcept{
        resultt = difference[0];
        result[0] = difference[0];

        for(int i=1; i<N; i++){
            result[i] = (result[i-1]/difference[i-1]) * difference[i] * (x - points[i-1]);
            resultt += result[i];
        }

        return resultt;
    }
};

#endif //VYCHISLITELNAYA_MATEMATIKA_INTERPOLATOR_HPP
