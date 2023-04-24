#include "bernstein_util.hpp"

void BernsteinPoly::SetTimeInterval(float time_interval_[])
{
    this->time_interval = time_interval_;
}

void BernsteinPoly::SetBernsteinCoeff(float bernstein_coeff_[])
{
    this->bernstein_coeff = bernstein_coeff_;
}

int BernsteinPoly::GetDegree() const
{
    if(bernstein_coeff==nullptr)
        return 0;
    else
        return sizeof(bernstein_coeff)/sizeof(bernstein_coeff[0])-1;
}

float BernsteinPoly::getValue(float t)
{
    int poly_order = this->GetDegree();
    float value = 0.0f;
    for(int i = 0; i<poly_order+1; i++)
        value+=bernstein_coeff[i]*float(nchoosek(poly_order,i))*(float)pow(t-time_interval[0],i)*pow(time_interval[1]-t,poly_order-i)/(float)pow(time_interval[1]-time_interval[0],poly_order);    
    return value;
}

BernsteinPoly BernsteinPoly::ElevateDegree(int m)
{
    return *this;
}

BernsteinPoly BernsteinPoly::operator+(const BernsteinPoly &rhs_)
{
    BernsteinPoly result(*this);
    int poly_order = result.GetDegree();
    for(int i = 0;i<poly_order+1;i++)
        result.bernstein_coeff[i]+=rhs_.bernstein_coeff[i];
    return result;
}

BernsteinPoly BernsteinPoly::operator-(const BernsteinPoly &rhs_)
{
    BernsteinPoly result(*this);
    int poly_order = result.GetDegree();
    for(int i = 0;i<poly_order+1;i++)
        result.bernstein_coeff[i]-=rhs_.bernstein_coeff[i];
    return result;
}

BernsteinPoly BernsteinPoly::operator*(const BernsteinPoly &rhs_)
{
    int n_lhs = this->GetDegree();
    int n_rhs = rhs_.GetDegree();
    float* dataPtr = new float[n_lhs+n_rhs+1];
    if(n_lhs>=n_rhs){

    }
    else{   //n_lhs<n_rhs

    }

    BernsteinPoly result(this->GetTimeInterval(),dataPtr);
    return result;
}

BernsteinPoly BernsteinPoly::operator*(const float &scalar_)
{
    BernsteinPoly result(*this);
    int poly_order = this->GetDegree();
    for(int i = 0;i<poly_order+1;i++)
        result.bernstein_coeff[i]*=scalar_;
    return result;
}
