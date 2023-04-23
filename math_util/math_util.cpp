#include "math_util.hpp"

int factorial(int num)
{
        if (num<=1)
            return 1;
        return num * factorial(num - 1);
}
int nchoosek(int n, int r)
{
    long long int number = 1;
    if (n-r > r){
        for(int i = n; i>n-r;i--)
        {
            number = number*i;
        }
        return (int)(number/factorial(r));
    }
    else{
        for(int i = n; i>r;i--)
        {
            number = number*i;
        }
        return (int) (number/factorial(n-r));
    }    
}
