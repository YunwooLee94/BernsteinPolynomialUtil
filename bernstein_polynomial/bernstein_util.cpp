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
    int poly_order = this->GetDegree();
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
    BernsteinPoly result(this->GetTimeInterval(),dataPtr);
    return result;
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
