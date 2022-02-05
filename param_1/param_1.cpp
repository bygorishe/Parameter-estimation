#include <iostream>
#include <cmath>
#include <vector>

# define pi 3.14159265358979323846

using namespace std;

class MNK {
private:
    double delta = 1., I = 10.;
    vector<int> A = { 0, 0, 0 }, B = { 100, 0, 0 },
        M1 = { 200, 0, 0 }, M2 = { 500, 0, 0 }, M3 = { 1000, 0, 0 },
        N1 = { 300, 0, 0 }, N2 = { 600, 0, 0 }, N3 = { 1100, 0, 0 };
    vector<double> V_sint, w;
    vector<double> V = { potential_V(M1, N1), potential_V(M2, N2), potential_V(M3, N3) };

    double r(vector<int> a, vector<int> b) {
        return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
    }
   
    double potential_V(vector<int> M, vector<int> N) { //Дифференциал по I, а также V при I = 1 
        return ((1 / r(B, M) - 1 / r(A, M)) - (1 / r(B, N) - 1 / r(A, N))) / (2 * pi * delta);
    }

    double F(double I) {
        double sum=0.;
        for (int i = 0; i < 3; i++)
            sum += w[i] * (I * V[i] - V_sint[i]) * w[i] * (I * V[i] - V_sint[i]);
    }

public:
    //генерируем синт данные, допустим при I = 10
    //считаем весовой коэффициент 
    //считаем функционал 
    //дифф и видимо Гаусс Ньютон потом
    //получаем дельтуI
    //делаем новое приближение
    //повторяем до функ = 0 и видимо там I будет равен искомому 
    // я нне уверен и хачу пиццу
    
    
    MNK() {
        V_sint.resize(3); w.resize(3);
        V_sint[0] = V[0] * I; V_sint[1] = V[1] * I; V_sint[2] = V[2] * I;
        w[0] = 1 / V_sint[0]; w[1] = 1 / V_sint[1]; w[2] = 1 / V_sint[2];
        
        
        
        double FFF = F()
    }

};


int main()
{
    
}

