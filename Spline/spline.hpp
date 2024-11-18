#ifndef UNTITLED_SPLAIN_HPP
#define UNTITLED_SPLAIN_HPP

#include <vector>
#include <type_traits>
#include <cmath>

/** класс для работы с трехдиагональной матрицей **/
template<typename T>
class ThreeDiagonalMatrix {
private:
    std::vector<T> a_;
    std::vector<T> b_;
    std::vector<T> c_;
public:
    ThreeDiagonalMatrix(const std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c): a_(a), b_(b), c_(c){}

    const std::vector<T>& get_a() const{return a_;}
    const std::vector<T>& get_b() const{return b_;}
    const std::vector<T>& get_c() const{return c_;}

    const T& a(const std::size_t& i) const{return a_[i];}
    const T& b(const std::size_t& i) const{return b_[i];}
    const T& c(const std::size_t& i) const{return c_[i];}
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

/** Функция для решения методм  прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve( const ThreeDiagonalMatrix<mType>& A,
                                            const std::vector<cType>& d){
    std::vector<mType> p(d.size());
    std::vector<DivisType<cType, mType>> q(d.size());
    std::vector<DivisType<cType, mType>> x(d.size());

    p[0] = -A.c(0) / A.b(0);
    q[0] = d[0] / A.b(0);

    for(std::size_t i = 1; i < d.size() - 1; i++){
        p[i] = -(A.c(i) / (A.a(i - 1) * p[i - 1] + A.b(i)));
        q[i] = (d[i] - A.a(i - 1) * q[i - 1]) / (A.a(i - 1) * p[i - 1] + A.b(i));
    }

    q[d.size() - 1] = (d[d.size() - 1] - A.a(d.size() - 2) * q[d.size() - 2]) / (A.a(d.size() - 2) * p[d.size() - 2] + A.b(d.size() - 1));

    x[d.size() - 1] = q[d.size() - 1];

    for(std::size_t i = 1; i < d.size(); i++){
        x[d.size() - 1 - i] = p[d.size() - 1 - i] * x[d.size() - i] + q[d.size() - 1 - i];
    }

    return x;
}

/**
* xType - тип аргумента x.
* yType - тип значения функции y
*/
template<typename xType, typename yType>
class CubicSpline {
private:
    std::vector<yType> a;
    std::vector<yType> b;
    std::vector<yType> c;
    std::vector<yType> d;
    std::vector<yType> Points;
public:
    CubicSpline(const std::vector<xType> &points, const std::vector<yType>& values) {
        std::vector<yType> a_(values.size() - 2);
        std::vector<yType> b_(values.size() - 1);
        std::vector<yType> c_(values.size() - 2);
        std::vector<yType> temp_c(values.size() - 2);
        std::vector<yType> u(values.size() - 1); //разделённые разности
        std::vector<xType> h_(points.size());

        Points.push_back(points[0]);
        for(int i = 1; i < points.size(); i ++){
            h_[i] = points[i] - points[i-1];
            Points.push_back(points[i]);
        }

        c.push_back(1);

        for(int i = 0; i < values.size() - 1; i++){
            u[i] = (((values[i+2] - values[i+1]) / (points[i+2] - points[i+1])) - ((values[i+1] - values[i]) / (points[i+1] - points[i]))) / (points[i+2] - points[i]);
        }

        for(int i = 0; i < values.size() - 2; i++){
            a.push_back(values[i]) ;

            a_[i] = h_[i+2] / (h_[i+1] + h_[i+2]);
            b_[i] = 2;
            c_[i] = h_[i+1] / (h_[i+1] + h_[i+2]);
        }
        b_[values.size() - 2] = 2;
        a.push_back(values[values.size() - 2]);
        a.push_back(values[values.size() - 1]);

        ThreeDiagonalMatrix<xType> Matrix(a_, b_, c_);

        temp_c = solve<xType, yType>(Matrix, u);

        for(int i = 0; i < values.size() - 2; i++){
            c.push_back(temp_c[i]);
        }

        c.push_back(exp(10));
        b.push_back(0);
        d.push_back(0);

        for(int i = 1; i < values.size(); i++){
            b.push_back((c[i] * h_[i] / 3) + (c[i-1] * h_[i] / 6) + (values[i] - (values[i-1])) / (points[i] - points[i-1]));
            d.push_back((c[i] - c[i-1]) / h_[i]);
        }
    }

    yType interpolate(const xType& x) const noexcept{
        int index = 0;

        for(int i = 0; i < Points.size(); i++){
            if(x < Points[i]){
                index = i;
                break;
            }
        }

        return a[index] + b[index] * (x - Points[index]) + 0.5 * c[index] * pow(x - Points[index], 2) + d[index] * pow(x - Points[index], 3) / 6;
    }
};

#endif //UNTITLED_SPLAIN_HPP
