# BernsteinPolynomial

This package gives utility functions related to **"Bernstein Polynomial"**. 
The package is implemented in C++ and can be used as a library. 



## 1. Functions
Class **"BernsteinPoly "** is defined with following core functions.
### Fundamental Functions:
* `GetDegree()`: Get degree of the polynomial
* `GetTimeInterval()`: Get time interval of the polynomial (return 2-element-array)
* `GetBernsteinCoefficient()`: Get coefficients (return (degree+1)-element-array)
### Operator Functions:
* `+`: Addition between two polynomials
* `-`: Subtract between two polynomials
* `*`: Multiplication between 1) two polynomials and 2) polynomial with scalar
### Additional Functions:
* `ElevateDegree(m)`: Change the degree to 'm' ('m' should be equal or greater than the degree of base polynomial).

## 2. Installation of this package
```
git clone https://github.com/YunwooLee94/BernsteinPolynomialUtil.git
cd BernsteinPolynomialUtil
mkdir build && cd build
cmake..
make
sudo make install
```


##Example Code
```c++
#include <iostream>
using namespace std;
int main(int argc, char **argv){
    float time_interval[2] = {0.0,1.0};
    float coefficient_a[4] = {0.0, 1.0, 2.0, 3.0};
    float coefficient_b[4] = {0.0, 2.0, 4.0, 6.0};
    int degree = 3; // equals to size of an array of coefficients -1
    BernsteinPoly a(time_interval,coefficient_a,degree);
    BernsteinPoly b(time_interval,coefficient_b,degree);
    // Fundamental Functions
    cout<<"Degree of a: "<<a.GetDegree()<<endl; // get degree of polynomial "a"
    cout<<"Time interval of a: ["<<a.GetTimeInterval()[0]<<", "<<a.GetTimeInterval()[1]<<"]"<<endl;
    cout<<"Coefficients of a"<<endl; // get coefficients
    for(int i = 0;i<=degree;i++)
        cout<<i<<"-th coefficient: "a.GetCoefficient()<<endl;
    // Operator Functions
    BernsteinPoly result_addition = a+b; // a plus b
    BernsteinPoly result_subtraction = a-b; // a minus b
    BernsteinPoly result_multiplication = a*b // a multiply b
    BernsteinPoly result_const = a*2.0f // a multiply scalar 2.0f
    // Degree Elevation
    BernsteinPoly result_degree_elevation = a.ElevateDegree(5); // Change representation of "a" into 5-th order polynomial
    return 0;
}
```
