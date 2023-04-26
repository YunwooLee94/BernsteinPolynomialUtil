#include "bernstein_util.hpp"

void BernsteinPoly::SetTimeInterval(float time_interval_[])
{
    this->time_interval = time_interval_;
}

void BernsteinPoly::SetBernsteinCoeff(float bernstein_coeff_[])
{
    this->bernstein_coeff = bernstein_coeff_;
}

void BernsteinPoly::SetDegree(int degree_)
{
    this->degree = degree_;
}

int BernsteinPoly::GetDegree() const
{
    return this->degree;
}

float BernsteinPoly::getValue(float t)
{
    int poly_order = this->degree;
    float value = 0.0f;
    for(int i = 0; i<poly_order+1; i++)
        value+=bernstein_coeff[i]*float(nchoosek(poly_order,i))*(float)pow(t-time_interval[0],i)*pow(time_interval[1]-t,poly_order-i)/(float)pow(time_interval[1]-time_interval[0],poly_order);    
    return value;
}

BernsteinPoly BernsteinPoly::ElevateDegree(int m)
{
    int poly_order = this->degree;
    float ** mat;
    mat = new float* [poly_order+1];
    for(int i =0;i<poly_order+1;i++)
        mat[i] = new float [m+1];
    
    for (int i =0;i<poly_order+1;i++){
        for(int j =0;j <= m-poly_order;j++)
            mat[i][i+j] = float(nchoosek(m-poly_order,j))*float(nchoosek(poly_order,i))/float(nchoosek(m,i+j));
    }
    float element;
    float * dataPtr = new float [m+1];
    for(int i =0;i<m+1;i++){
        element = 0.0f;
        for(int j=0;j<poly_order+1;j++){
            element += this->bernstein_coeff[j]*mat[j][i];
        }
        dataPtr[i] = element;
    }
    for(int i =0;i<poly_order+1;i++)
        delete [] mat[i];
    delete []mat;
    BernsteinPoly result(this->GetTimeInterval(),dataPtr,m);
    return result;
}

BernsteinPoly BernsteinPoly::operator+(const BernsteinPoly &rhs_)
{
    int poly_order = this->degree;
    auto *dataPtr = new float[degree+1];
    for(int i = 0;i<poly_order+1;i++)
        dataPtr[i] = this->bernstein_coeff[i]+rhs_.bernstein_coeff[i];
    BernsteinPoly result(this->GetTimeInterval(),dataPtr,poly_order);
    return result;
}

BernsteinPoly BernsteinPoly::operator-(const BernsteinPoly &rhs_)
{
    int poly_order = this->degree;
    auto *dataPtr = new float[degree+1];
    for(int i = 0;i<poly_order+1;i++)
        dataPtr[i] = this->bernstein_coeff[i]-rhs_.bernstein_coeff[i];
    BernsteinPoly result(this->GetTimeInterval(),dataPtr,poly_order);
    return result;
}

BernsteinPoly BernsteinPoly::operator*(const BernsteinPoly &rhs_)
{
    int n_lhs = this->degree;
    int n_rhs = rhs_.GetDegree();
    float* dataPtr = new float[n_lhs+n_rhs+1];
    if(n_lhs>=n_rhs){
        for (int i = 0; i <= n_lhs + n_rhs; i++){
            if (i <= n_rhs){
                for (int j = i; j >= 0; j--)
                    dataPtr[i] += float(nchoosek(n_lhs, j)) *
                                  float(nchoosek(n_rhs, i - j)) / float(nchoosek(n_lhs + n_rhs, i)) * this->bernstein_coeff[j] * rhs_.bernstein_coeff[i - j];
            }
            else if ((i >= n_rhs + 1) and (i <= n_lhs)){
                for (int j = i; j >= i - n_rhs; j--)
                    dataPtr[i] += float(nchoosek(n_lhs, j)) *
                                  float(nchoosek(n_rhs, i - j)) / float(nchoosek(n_lhs + n_rhs, i)) * this->bernstein_coeff[j] * rhs_.bernstein_coeff[i - j];
            }
            else if (i >= n_lhs + 1){
                for (int j = n_lhs; j >= i - n_rhs; j--)
                    dataPtr[i] += float(nchoosek(n_lhs, j)) *
                                  float(nchoosek(n_rhs, i - j)) / float(nchoosek(n_lhs + n_rhs, i)) * this->bernstein_coeff[j] * rhs_.bernstein_coeff[i - j];
            }
        }
    }
    else{   //n_lhs<n_rhs
        for (int i = 0; i <= n_lhs + n_rhs; i++){
            if (i <= n_lhs){
                for (int j = i; j >= 0; j--)
                    dataPtr[i] += float(nchoosek(n_rhs, j)) *
                                  float(nchoosek(n_lhs, i - j)) / float(nchoosek(n_lhs + n_rhs, i)) * rhs_.bernstein_coeff[j] * this->bernstein_coeff[i - j];
            }
            else if ((i >= n_lhs + 1) and (i <= n_rhs)){
                for (int j = i; j >= i - n_rhs; j--)
                    dataPtr[i] += float(nchoosek(n_rhs, j)) *
                                  float(nchoosek(n_lhs, i - j)) / float(nchoosek(n_lhs + n_rhs, i)) * rhs_.bernstein_coeff[j] * this->bernstein_coeff[i - j];
            }
            else if (i >= n_rhs + 1){
                for (int j = n_rhs; j >= i - n_lhs; j--)
                    dataPtr[i] += float(nchoosek(n_rhs, j)) *
                                  float(nchoosek(n_lhs, i - j)) / float(nchoosek(n_lhs + n_rhs, i)) * rhs_.bernstein_coeff[j] * this->bernstein_coeff[i - j];
            }
        }
    }
    BernsteinPoly result(this->GetTimeInterval(),dataPtr,n_lhs+n_rhs);
    return result;
}

BernsteinPoly BernsteinPoly::operator*(const float &scalar_){
    int poly_order = this->degree;
    auto *dataPtr = new float[degree+1];
    for(int i = 0;i<poly_order+1;i++)
        dataPtr[i] = this->bernstein_coeff[i]*scalar_;
    BernsteinPoly result(this->GetTimeInterval(),dataPtr,poly_order);
    return result;
}
