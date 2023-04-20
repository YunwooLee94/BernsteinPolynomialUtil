#pragma once
#include <vector>
template<typename T>
class BernsteinPoly{
    private:
        std::vector<T> time_interval; // initial time and terminal time [2]
        std::vector<T> bernstein_coeff; // bernstein coefficient
        
    public:
        BernsteinPoly(){};
        BernsteinPoly(T t0_, T tf_){SetTimeInterval(t0_,tf_);};
        BernsteinPoly(const std::vector<T>& bernstein_coeff_){SetBernsteinCoeff(bernstein_coeff_);};
        BernsteinPoly(const std::vector<T>& bernstein_coeff_,T t0_, T tf_){SetBernsteinCoeff(bernstein_coeff_);{SetTimeInterval(t0_,tf_);}};
        void SetTimeInterval(T t0_, T tf_){time_interval.clear(); time_interval.push_back(t0_);;time_interval.push_back(tf_);};
        void SetBernsteinCoeff(const std::vector<T>& bernstein_coeff_){bernstein_coeff.clear(); bernstein_coeff=bernstein_coeff_;};
        int GetCoefficient(){return bernstein_coeff;};
        bool IsSet(){return not time_interval.empty() and not bernstein_coeff.empty();};
        T* GetTimeInterval(){return time_interval;};
        T GetValue(T t_); // TODO: in cpp
        BernsteinPoly ElevateDegree(const BernsteinPoly<T>& low, int m); // Change to higher degree (m) // TODO: in cpp
        BernsteinPoly operator+(const BernsteinPoly<T> &rhs); // TODO: in cpp
        BernsteinPoly operator-(const BernsteinPoly<T> &rhs); // TODO: in cpp
        BernsteinPoly operator*(const BernsteinPoly<T> &rhs); // TODO: in cpp
};

