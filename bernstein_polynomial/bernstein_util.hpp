#pragma once
#include "../math_util/math_util.hpp"

class BernsteinPoly{
private:
    float* time_interval; // initial time and terminal time 
    float* bernstein_coeff; // bernstein coefficient
public:
    BernsteinPoly(float time_interval_[], float bernstein_coeff_[]){time_interval = time_interval_; bernstein_coeff = bernstein_coeff_;};
    BernsteinPoly(const BernsteinPoly & bern_poly_){this->time_interval = bern_poly_.time_interval; this->bernstein_coeff= bern_poly_.bernstein_coeff;};
    void SetTimeInterval(float time_interval_[]);
    void SetBernsteinCoeff(float bernstein_coeff_[]);
    int GetDegree() const;
    bool IsSet(){return ((time_interval != nullptr) and ( bernstein_coeff != nullptr));};
    float* GetTimeInterval(){return time_interval;};
    float getValue(float t);
    BernsteinPoly ElevateDegree(int m); // Change to higher degree (m)
    BernsteinPoly operator+(const BernsteinPoly &rhs_);
    BernsteinPoly operator-(const BernsteinPoly &rhs_);
    BernsteinPoly operator*(const BernsteinPoly &rhs_);
    BernsteinPoly operator*(const float & scalar_);

};

