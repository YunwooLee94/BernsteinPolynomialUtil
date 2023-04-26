#pragma once
#include "../math_util/math_util.hpp"
#include <type_traits>
class BernsteinPoly{
private:
    float* time_interval; // initial time and terminal time 
    float* bernstein_coeff; // bernstein coefficient
    int degree; // The degree of a polynomial
public:
    BernsteinPoly(float time_interval_[], float bernstein_coeff_[], const int & degree_){time_interval = time_interval_; bernstein_coeff = bernstein_coeff_; degree=degree_;};
    BernsteinPoly(const BernsteinPoly & bern_poly_){this->time_interval = bern_poly_.time_interval; this->bernstein_coeff= bern_poly_.bernstein_coeff; this->degree=bern_poly_.degree;};
    void SetTimeInterval(float time_interval_[]);
    void SetBernsteinCoeff(float bernstein_coeff_[]);
    void SetDegree(int degree_);
    int GetDegree() const;
    bool IsSet(){return ((time_interval != nullptr) and ( bernstein_coeff != nullptr));};
    float* GetTimeInterval(){return time_interval;};
    float* GetBernsteinCoefficient(){return this->bernstein_coeff;};
    float getValue(float t);
    BernsteinPoly ElevateDegree(int m); // Change to higher degree (m)
    BernsteinPoly operator+(const BernsteinPoly &rhs_);
    BernsteinPoly operator-(const BernsteinPoly &rhs_);
    BernsteinPoly operator*(const BernsteinPoly &rhs_);
    BernsteinPoly operator*(const float & scalar_);

};

