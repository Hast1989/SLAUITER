#include <iostream>
#include <fstream>
#include<cmath>
#include <iomanip>
#include<string>
int iter, kest;
int indikQR = 0;
int indikG = 0;
template<typename T>
T NormVec1(T* x, int n)
{
    T res = 0;
    for (int i = 0; i < n; i++)
    {
        res = res + abs(x[i]);
    }
    return res;
}
template<typename T>
T NormVecinf(T* x, int n)
{
    T max = 0;
    for (int i = 0; i < n; i++)
    {
        if (max < abs(x[i]))
        {
            max = abs(x[i]);
        }
    }
    return max;
}
template<typename T>
void MultM(T** Matrixl, T** Matrixr, T** Mresult, int n)
{
    T res;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res = 0;
            for (int k = 0; k < n; k++)
            {
                res = res + Matrixl[i][k] * Matrixr[k][j];
            }
            Mresult[i][j] = res;
        }
    }
}
template<typename T>
T NormVect(T* x, int n)
{
    T res = 0;
    for (int i = 0; i < n; i++)
    {
        res = res + x[i] * x[i];
    }
    res = sqrt(res);
    return res;
}
template<typename T>
void MultWV(T** Matrix, T* x, T* result, int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = 0;
        for (int j = 0; j < n; j++)
        {
            result[i] = result[i] + Matrix[i][j] * x[j];
        }
    }
}
template<typename T>
int FindLeader(T** Matrix, int n, int j)
{
    int maxint = j;
    T max = 0;
    for (int i = j; i < n; i++)
    {
        if (abs(Matrix[i][j]) > max)
        {
            max = abs(Matrix[i][j]);
            maxint = i;
        }
    }
    if (max == 0)
    {
        indikG = 1;
    }
    return maxint;
}
template<typename T>
void SwichLines(T** Matrix, T* rightb, int i, int j)
{
    T resb = rightb[i];
    T* resline = Matrix[i];
    Matrix[i] = Matrix[j];
    Matrix[j] = resline;
    rightb[i] = rightb[j];
    rightb[j] = resb;

}
template<typename T>
void Gauss(T** Matrix, T* rightb, int n, T* x)
{
    T d, s;
    for (int k = 0; k < n - 1; k++)
    {
        SwichLines(Matrix, rightb, FindLeader(Matrix, n, k), k);
        if (indikG == 1)
        {
            return;
        }
        for (int j = k + 1; j < n; j++)
        {
            if (abs(Matrix[j][k]) != 0)
            {
                d = Matrix[j][k] / Matrix[k][k];
                Matrix[j][k] = 0;
                for (int i = k + 1; i < n; i++)
                {
                    Matrix[j][i] = Matrix[j][i] - d * Matrix[k][i];
                }
                rightb[j] = rightb[j] - d * rightb[k];
            }
        }
    }
    if (abs(Matrix[n - 1][n - 1]) < 0.00000000000001)
    {
        indikG = 1;
        return;
    }
    for (int k = n; k > 0; k--)
    {
        d = 0;
        for (int j = k; j < n; j++)
        {
            s = Matrix[k - 1][j] * x[j];
            d = d + s;
        }
        x[k - 1] = (rightb[k - 1] - d) / Matrix[k - 1][k - 1];
    }
}
template<typename T>
void Invers(T** Matrix, T** Invers, int n)
{
    T* righte;
    T* xe;
    T** ResM;
    ResM = new T * [n];
    xe = new T[n];
    righte = new T[n];
    for (int i = 0; i < n; i++)
    {
        ResM[i] = new T[n];
    }
    for (int i = 0; i < n; i++)
    {
        for (int o = 0; o < n; o++)
        {
            righte[o] = 0;
            for (int p = 0; p < n; p++)
            {
                ResM[o][p] = Matrix[o][p];
            }
        }
        righte[i] = 1;
        Gauss(ResM, righte, n, xe);
        for (int k = 0; k < n; k++)
        {
            Invers[k][i] = xe[k];
        }
    }
    for (int i = 0; i < n; i++) {
        delete[] ResM[i];
    }
    delete[] ResM;
    delete[] righte;
    delete[] xe;

}
template<typename T>
void Transpose(T** Matrix, int n)
{
    T res;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            res = Matrix[i][j];
            Matrix[i][j] = Matrix[j][i];
            Matrix[j][i] = res;
        }
    }
}
template<typename T>
T NormM1(T** Matrix, int n)
{
    T res;
    T max = 0;
    for (int i = 0; i < n; i++)
    {
        res = 0;
        for (int j = 0; j < n; j++)
        {
            res = res + abs(Matrix[j][i]);
        }
        if (max < res)
        {
            max = res;
        }
    }
    return max;
}
template<typename T>
T NormMInf(T** Matrix, int n)
{
    T res;
    T max = 0;
    for (int i = 0; i < n; i++)
    {
        res = 0;
        for (int j = 0; j < n; j++)
        {
            res = res + abs(Matrix[i][j]);
        }
        if (max < res)
        {
            max = res;
        }
    }
    return max;
}
template <typename T>
void Justiter(int n, T** Matrix, T* rightb, T* x, T tau, T eps)
{
    T* xk;
    T* xkl;
    xk = new T[n];
    xkl = new T[n];
    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0;
        xk[i] = tau * rightb[i];
        if (Matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
            {
                Matrix[i][j] = -Matrix[i][j];
            }
            rightb[i] = -rightb[i];
        }
    }
    while (fabs(NormVec1(xk, n) - NormVec1(xkl, n)) > eps)
    {
        iter++;
        if (iter > kest + 1)
        {
            break;
        }
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        for (int i = 0; i < n; i++)
        {
            xk[i] = xkl[i] + tau * rightb[i];
            for (int j = 0; j < n; j++)
            {
                xk[i] = xk[i] - tau * Matrix[i][j] * xkl[j];
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }

    delete[] xk;
    delete[] xkl;
}
template <typename T>
void Justiter2(int n, T** Matrix, T* rightb, T* x, T tau, T eps)
{
    T* xk;
    T* xkl;
    T* res;
    res = new T[n];
    xk = new T[n];
    xkl = new T[n];

    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0;
        xk[i] = tau * rightb[i];
        if (Matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
            {
                Matrix[i][j] = -Matrix[i][j];
            }
            rightb[i] = -rightb[i];
            res[i] = rightb[i];
        }
    }
    while (NormVect(res, n) > eps)
    {
        iter++;
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        for (int i = 0; i < n; i++)
        {
            xk[i] = xkl[i] + tau * rightb[i];
            for (int j = 0; j < n; j++)
            {
                xk[i] = xk[i] - tau * Matrix[i][j] * xkl[j];
            }
        }
        MultWV(Matrix, xk, res, n);
        for (int i = 0; i < n; i++)
        {
            res[i] = res[i] - rightb[i];
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }
    delete[] res;
    delete[] xk;
    delete[] xkl;
}
template <typename T>
void Jacobi(int n, T** Matrix, T* rightb, T* x, T eps)
{
    T* xk;
    T* xkl;
    xk = new T[n];
    xkl = new T[n];
    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0;
        xk[i] = rightb[i];
        if (Matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
            {
                Matrix[i][j] = -Matrix[i][j];
            }
            rightb[i] = -rightb[i];
        }
    }
    while (fabs(NormVec1(xk, n) - NormVec1(xkl, n)) > eps)
    {
        iter++;
        if (iter > kest + 1)
        {
            break;
        }
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        for (int i = 0; i < n; i++)
        {
            xk[i] = rightb[i] / Matrix[i][i];
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                {
                    xk[i] = xk[i] - (Matrix[i][j] * xkl[j]) / Matrix[i][i];
                }
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }

    delete[] xk;
    delete[] xkl;
}
template <typename T>
void Jacobi2(int n, T** Matrix, T* rightb, T* x, T eps)
{
    T* xk;
    T* xkl;
    T* res;
    res = new T[n];
    xk = new T[n];
    xkl = new T[n];
    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0;
        xk[i] = rightb[i];
        if (Matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
            {
                Matrix[i][j] = -Matrix[i][j];
            }
            rightb[i] = -rightb[i];
            res[i] = rightb[i];
        }
    }
    while (NormVect(res, n) > eps)
    {
        iter++;
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        for (int i = 0; i < n; i++)
        {
            xk[i] = rightb[i] / Matrix[i][i];
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                {
                    xk[i] = xk[i] - (Matrix[i][j] * xkl[j]) / Matrix[i][i];
                }
            }
        }
        MultWV(Matrix, xk, res, n);
        for (int i = 0; i < n; i++)
        {
            res[i] = res[i] - rightb[i];
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }
    delete[] res;
    delete[] xk;
    delete[] xkl;
}
template <typename T>
void Zeidel(int n, T** Matrix, T* rightb, T* x, T eps)
{
    T* xk;
    T* xkl;
    xk = new T[n];
    xkl = new T[n];
    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0;
        xk[i] = rightb[i];
        if (Matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
            {
                Matrix[i][j] = -Matrix[i][j];
            }
            rightb[i] = -rightb[i];
        }
    }
    while (fabs(NormVec1(xk, n) - NormVec1(xkl, n)) > eps)
    {
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        for (int i = 0; i < n; i++)
        {
            xk[i] = rightb[i] / Matrix[i][i];
            for (int j = 0; j < i; j++)
            {
                xk[i] = xk[i] - Matrix[i][j] * xkl[j] / Matrix[i][i];
            }
            for (int j = i + 1; j < n; j++)
            {
                xk[i] = xk[i] - Matrix[i][j] * xkl[j] / Matrix[i][i];
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }

    delete[] xk;
    delete[] xkl;
}
template <typename T>
void Zeidel2(int n, T** Matrix, T* rightb, T* x, T eps)
{
    T* xk;
    T* xkl;
    T* res;
    res = new T[n];
    xk = new T[n];
    xkl = new T[n];
    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0;
        xk[i] = rightb[i];
        if (Matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
            {
                Matrix[i][j] = -Matrix[i][j];
            }
            rightb[i] = -rightb[i];
            res[i] = rightb[i];
        }
    }
    while (NormVect(res, n) > eps)
    {
        iter++;
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        for (int i = 0; i < n; i++)
        {
            xk[i] = rightb[i] / Matrix[i][i];
            for (int j = 0; j < i; j++)
            {
                xk[i] = xk[i] - Matrix[i][j] * xkl[j] / Matrix[i][i];
            }
            for (int j = i + 1; j < n; j++)
            {
                xk[i] = xk[i] - Matrix[i][j] * xkl[j] / Matrix[i][i];
            }
        }
        MultWV(Matrix, xk, res, n);
        for (int i = 0; i < n; i++)
        {
            res[i] = res[i] - rightb[i];
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }
    delete[] res;
    delete[] xk;
    delete[] xkl;
}
template <typename T>
void Relaks(int n, T** Matrix, T* rightb, T* x, T omega, T eps)
{
    T* xk;
    T* xkl;
    xk = new T[n];
    xkl = new T[n];
    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0;
        xk[i] = rightb[i];
        if (Matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
            {
                Matrix[i][j] = -Matrix[i][j];
            }
            rightb[i] = -rightb[i];
        }
    }
    while (fabs(NormVec1(xk, n) - NormVec1(xkl, n)) > eps)
    {
        iter++;
        if (iter > kest + 1)
        {
            break;
        }
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        for (int i = 0; i < n; i++)
        {
            xk[i] = omega * rightb[i] / Matrix[i][i] + (1 - omega) * xkl[i];
            for (int j = 0; j < i; j++)
            {
                xk[i] = xk[i] - omega * Matrix[i][j] * xk[j] / Matrix[i][i];
            }
            for (int j = i + 1; j < n; j++)
            {
                xk[i] = xk[i] - omega * Matrix[i][j] * xkl[j] / Matrix[i][i];
            }
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }

    delete[] xk;
    delete[] xkl;
}
template <typename T>
void Relaks2(int n, T** Matrix, T* rightb, T* x, T omega, T eps)
{
    T* xk;
    T* xkl;
    T* res;
    res = new T[n];
    xk = new T[n];
    xkl = new T[n];
    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0;
        xk[i] = rightb[i];
        if (Matrix[i][i] < 0)
        {
            for (int j = 0; j < n; j++)
            {
                Matrix[i][j] = -Matrix[i][j];
            }
            rightb[i] = -rightb[i];
            res[i] = rightb[i];
        }
    }
    while (NormVect(res, n) > eps)
    {
        iter++;
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        for (int i = 0; i < n; i++)
        {
            xk[i] = omega * rightb[i] / Matrix[i][i] + (1 - omega) * xkl[i];
            for (int j = 0; j < i; j++)
            {
                xk[i] = xk[i] - omega * Matrix[i][j] * xk[j] / Matrix[i][i];
            }
            for (int j = i + 1; j < n; j++)
            {
                xk[i] = xk[i] - omega * Matrix[i][j] * xkl[j] / Matrix[i][i];
            }
        }
        MultWV(Matrix, xk, res, n);
        for (int i = 0; i < n; i++)
        {
            res[i] = res[i] - rightb[i];
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }
    delete[] res;
    delete[] xk;
    delete[] xkl;
}
template <typename T>
void test(int set, int neps, int ntest)
{
    std::string adress;
    T** MC;
    T** Res;
    T** Matrix;
    T* x;
    T* rightb;
    T* Resb;
    T** Inv;
    T Nres[4][4] = { {5.,-7.,12.,4.},{10.,-10.,12.,4.},{5.,300.,4.,0.},{3.,5.,7.,-5.} };
    T res, minNC, tau, omega;
    int n, i, j;
    T eps[2] = { 0.0001, 0.0000001 };
    adress = "test";
    adress = adress + ".txt";
    std::ifstream file;
    file.open(adress);
    file >> n;
    MC = new T * [n];
    Inv = new T * [n];
    Res = new T * [n];
    Matrix = new T * [n];
    x = new T[n];
    rightb = new T[n];
    Resb = new T[n];
    kest = 1000;
    for (i = 0; i < n; i++)
    {
        Matrix[i] = new T[n];
        Res[i] = new T[n];
        Inv[i] = new T[n];
        MC[i] = new T[n];
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            file >> Matrix[i][j];
            Res[i][j] = Matrix[i][j];
        }
        file >> rightb[i];
        Resb[i] = rightb[i];
    }
    file.close();
    std::ofstream ans;
    adress = "answerJustiter";
    adress = adress + ".txt";
    ans.open(adress);
    ans << std::setprecision(set);
    iter = 0;
    tau = 1 / (NormM1(Matrix, n));
    Justiter(n, Matrix, rightb, x, tau, eps[neps - 1]);
    ans << "Matrix" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << rightb[i] << std::endl;
    }
    ans << "Tau" << ' ' << tau << std::endl;
    ans << "Iterations" << ' ' << iter << std::endl;
    ans << "X" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                Matrix[i][j] = 1 - tau * Matrix[i][j];
            }
            else
            {
                Matrix[i][j] = tau * Matrix[i][j];
            }
        }
    }
    ans << "Matrix C" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << "Vector y" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << tau * rightb[i] << ' ';
    }
    ans << std::endl;
    ans << "Norm matrix C1" << ' ' << NormM1(Matrix, n) << std::endl;
    ans << "Norm matrix Cinf" << ' ' << NormMInf(Matrix, n) << std::endl;
    if (NormM1(Matrix, n) > NormMInf(Matrix, n))
    {
        minNC = NormMInf(Matrix, n);
    }
    else
    {
        minNC = NormM1(Matrix, n);
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    iter = 0;
    Justiter(n, Matrix, rightb, x, tau, ((1 - minNC) / minNC) * eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov ((1 - NC) / NC)*eps iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    iter = 0;
    Justiter2(n, Matrix, rightb, x, tau, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov Axk-b iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    ans << std::endl;
    kest = (int)(log((1 / eps[neps - 1])) / (log(1 / minNC)));
    ans << "kest" << ' ' << kest << std::endl;
    iter = 0;
    Justiter(n, Matrix, rightb, x, tau, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke after kest" << ' ' << NormVect(x, n) << std::endl;
    kest = 1000;
    iter = 0;
    tau = 1 / (NormMInf(Matrix, n));
    Justiter(n, Matrix, rightb, x, tau, eps[neps - 1]);
    ans << "Matrix" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << rightb[i] << std::endl;
    }
    ans << "Tau" << ' ' << tau << std::endl;
    ans << "Iterations" << ' ' << iter << std::endl;
    ans << "X" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                Matrix[i][j] = 1 - tau * Matrix[i][j];
            }
            else
            {
                Matrix[i][j] = tau * Matrix[i][j];
            }
        }
    }
    ans << "Matrix C" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << "Vector y" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << tau * rightb[i] << ' ';
    }
    ans << std::endl;
    ans << "Norm matrix C1" << ' ' << NormM1(Matrix, n) << std::endl;
    ans << "Norm matrix Cinf" << ' ' << NormMInf(Matrix, n) << std::endl;
    if (NormM1(Matrix, n) > NormMInf(Matrix, n))
    {
        minNC = NormMInf(Matrix, n);
    }
    else
    {
        minNC = NormM1(Matrix, n);
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    iter = 0;
    Justiter(n, Matrix, rightb, x, tau, ((1 - minNC) / minNC) * eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov ((1 - NC) / NC)*eps iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    iter = 0;
    Justiter2(n, Matrix, rightb, x, tau, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov Axk-b iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    ans << std::endl;
    kest = (int)(log((1 / eps[neps - 1])) / (log(1 / minNC)));
    ans << "kest" << ' ' << kest << std::endl;
    iter = 0;
    Justiter(n, Matrix, rightb, x, tau, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke after kest" << ' ' << NormVect(x, n) << std::endl;
    kest = 1000;
    ans.close();
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    adress = "answerJacobi";
    adress = adress + ".txt";
    ans.open(adress);
    ans << std::setprecision(set);
    iter = 0;
    Jacobi(n, Matrix, rightb, x, eps[neps - 1]);
    ans << "Matrix" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << rightb[i] << std::endl;
    }
    ans << "Iterations" << ' ' << iter << std::endl;
    ans << "X" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        res = Matrix[i][i];
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                Matrix[i][j] = 1 - Matrix[i][j] / res;
            }
            else
            {
                Matrix[i][j] = -Matrix[i][j] / res;
            }
        }
    }
    ans << std::endl;
    ans << "Matrix C" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << "Vector y" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << rightb[i] / abs(Res[i][i]) << ' ';
    }
    ans << std::endl;
    ans << "Norm matrix C1" << ' ' << NormM1(Matrix, n) << std::endl;
    ans << "Norm matrix Cinf" << ' ' << NormMInf(Matrix, n) << std::endl;
    if (NormM1(Matrix, n) > NormMInf(Matrix, n))
    {
        minNC = NormMInf(Matrix, n);
    }
    else
    {
        minNC = NormM1(Matrix, n);
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    iter = 0;
    Jacobi(n, Matrix, rightb, x, ((1 - minNC) / minNC) * eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov ((1 - NC) / NC)*eps iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    iter = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    Jacobi2(n, Matrix, rightb, x, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov Axk-b iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    ans << std::endl;
    kest = (int)(log((1 / eps[neps - 1])) / (log(1 / minNC)));
    ans << "kest" << ' ' << kest << std::endl;
    iter = 0;
    Jacobi(n, Matrix, rightb, x, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke after kest" << ' ' << NormVect(x, n) << std::endl;
    kest = 1000;
    ans.close();
    adress = "answerRelaks";
    adress = adress + ".txt";
    ans.open(adress);
    ans << std::setprecision(set);
    iter = 0;
    omega = 0.9;
    Relaks(n, Matrix, rightb, x, omega, eps[neps - 1]);
    ans << "Matrix" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << rightb[i] << std::endl;
    }
    ans << "Omega" << ' ' << omega << std::endl;
    ans << "Iterations" << ' ' << iter << std::endl;
    ans << "X" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i > j)
            {
                MC[i][j] = 0;
            }
            else
            {
                if (i == j)
                {
                    MC[i][j] = (1 - omega) * Matrix[i][j];
                }
                else

                {
                    MC[i][j] = -omega * Matrix[i][j];
                }

            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i < j)
            {
                Matrix[i][j] = 0;
            }
            else
            {
                if (i == j)
                {
                    Matrix[i][j] = Matrix[i][j];
                }
                else

                {
                    Matrix[i][j] = omega * Matrix[i][j];
                }

            }
        }
    }
    Invers(Matrix, Inv, n);
    MultM(Inv, MC, Matrix, n);
    ans << std::endl;
    ans << "Matrix C" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Inv[i][j] = omega * Inv[i][j];
        }
    }
    MultWV(Inv, rightb, x, n);
    ans << "Vector y" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    ans << std::endl;
    ans << "Norm matrix C1" << ' ' << NormM1(Matrix, n) << std::endl;
    ans << "Norm matrix Cinf" << ' ' << NormMInf(Matrix, n) << std::endl;
    if (NormM1(Matrix, n) > NormMInf(Matrix, n))
    {
        minNC = NormMInf(Matrix, n);
    }
    else
    {
        minNC = NormM1(Matrix, n);
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i <= j)
            {
                Inv[i][j] = 0;
            }
            else
            {
                Inv[i][j] = Matrix[i][j];
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                MC[i][j] = 0;
            }
            else
            {
                MC[i][j] = Matrix[i][j];
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i >= j)
            {
                Matrix[i][j] = 0;
            }
            else
            {
                Matrix[i][j] = Matrix[i][j];
            }
        }
    }
    ans << std::endl;
    ans << "Matrix Cl,CD,CU" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Inv[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << std::endl;
    ans << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << MC[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << std::endl;
    ans << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    iter = 0;
    Relaks(n, Matrix, rightb, x, omega, ((1 - minNC) / minNC) * eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov ((1 - NC) / NC)*eps iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    iter = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    Relaks2(n, Matrix, rightb, x, omega, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov Axk-b iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    ans << std::endl;
    kest = (int)(log((1 / eps[neps - 1])) / (log(1 / minNC)));
    ans << "kest" << ' ' << kest << std::endl;
    iter = 0;
    Relaks(n, Matrix, rightb, x, omega, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke after kest" << ' ' << NormVect(x, n) << std::endl;
    kest = 1000;
    iter = 0;
    omega = 1.5;
    Relaks(n, Matrix, rightb, x, omega, eps[neps - 1]);
    ans << "Matrix" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << rightb[i] << std::endl;
    }
    ans << "Omega" << ' ' << omega << std::endl;
    ans << "Iterations" << ' ' << iter << std::endl;
    ans << "X" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i > j)
            {
                MC[i][j] = 0;
            }
            else
            {
                if (i == j)
                {
                    MC[i][j] = (1 - omega) * Matrix[i][j];
                }
                else

                {
                    MC[i][j] = -omega * Matrix[i][j];
                }

            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i < j)
            {
                Matrix[i][j] = 0;
            }
            else
            {
                if (i == j)
                {
                    Matrix[i][j] = Matrix[i][j];
                }
                else

                {
                    Matrix[i][j] = omega * Matrix[i][j];
                }

            }
        }
    }
    Invers(Matrix, Inv, n);
    MultM(Inv, MC, Matrix, n);
    ans << std::endl;
    ans << "Matrix C" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Inv[i][j] = omega * Inv[i][j];
        }
    }
    MultWV(Inv, rightb, x, n);
    ans << "Vector y" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    ans << std::endl;
    ans << "Norm matrix C1" << ' ' << NormM1(Matrix, n) << std::endl;
    ans << "Norm matrix Cinf" << ' ' << NormMInf(Matrix, n) << std::endl;
    if (NormM1(Matrix, n) > NormMInf(Matrix, n))
    {
        minNC = NormMInf(Matrix, n);
    }
    else
    {
        minNC = NormM1(Matrix, n);
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i <= j)
            {
                Inv[i][j] = 0;
            }
            else
            {
                Inv[i][j] = Matrix[i][j];
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                MC[i][j] = 0;
            }
            else
            {
                MC[i][j] = Matrix[i][j];
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i >= j)
            {
                Matrix[i][j] = 0;
            }
            else
            {
                Matrix[i][j] = Matrix[i][j];
            }
        }
    }
    ans << std::endl;
    ans << "Matrix Cl,CD,CU" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Inv[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << std::endl;
    ans << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << MC[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << std::endl;
    ans << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    iter = 0;
    Relaks(n, Matrix, rightb, x, omega, ((1 - minNC) / minNC) * eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov ((1 - NC) / NC)*eps iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    iter = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    Relaks2(n, Matrix, rightb, x, omega, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov Axk-b iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    ans << std::endl;
    kest = (int)(log((1 / eps[neps - 1])) / (log(1 / minNC)));
    ans << "kest" << ' ' << kest << std::endl;
    iter = 0;
    Relaks(n, Matrix, rightb, x, omega, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke after kest" << ' ' << NormVect(x, n) << std::endl;
    kest = 1000;
    ans.close();
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    adress = "answerZeidel";
    adress = adress + ".txt";
    ans.open(adress);
    ans << std::setprecision(set);
    iter = 0;
    Zeidel(n, Matrix, rightb, x, eps[neps - 1]);
    ans << "Matrix" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << rightb[i] << std::endl;
    }
    ans << "Iterations" << ' ' << iter << std::endl;
    ans << "X" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i > j)
            {
                MC[i][j] = 0;
            }
            else
            {
                if (i == j)
                {
                    MC[i][j] = 0;
                }
                else

                {
                    MC[i][j] = Matrix[i][j];
                }

            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i < j)
            {
                Matrix[i][j] = 0;
            }
            else
            {
                if (i == j)
                {
                    Matrix[i][j] = Matrix[i][j];
                }
                else

                {
                    Matrix[i][j] = Matrix[i][j];
                }

            }
        }
    }
    Invers(Matrix, Inv, n);
    MultM(Inv, MC, Matrix, n);
    ans << std::endl;
    ans << "Matrix C" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    MultWV(Inv, rightb, x, n);
    ans << "Vector y" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    ans << std::endl;
    ans << "Norm matrix C1" << ' ' << NormM1(Matrix, n) << std::endl;
    ans << "Norm matrix Cinf" << ' ' << NormMInf(Matrix, n) << std::endl;

    if (NormM1(Matrix, n) > NormMInf(Matrix, n))
    {
        minNC = NormMInf(Matrix, n);
    }
    else
    {
        minNC = NormM1(Matrix, n);
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i <= j)
            {
                Inv[i][j] = 0;
            }
            else
            {
                Inv[i][j] = Matrix[i][j];
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                MC[i][j] = 0;
            }
            else
            {
                MC[i][j] = Matrix[i][j];
            }
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i >= j)
            {
                Matrix[i][j] = 0;
            }
            else
            {
                Matrix[i][j] = Matrix[i][j];
            }
        }
    }
    ans << std::endl;
    ans << "Matrix Cl,CD,CU" << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Inv[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << std::endl;
    ans << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << MC[i][j] << ' ';
        }
        ans << std::endl;
    }
    ans << std::endl;
    ans << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ans << Matrix[i][j] << ' ';
        }
        ans << std::endl;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    iter = 0;
    Zeidel(n, Matrix, rightb, x, ((1 - minNC) / minNC) * eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov ((1 - NC) / NC)*eps iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    iter = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    Zeidel2(n, Matrix, rightb, x, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Ostanov Axk-b iters" << ' ' << iter << std::endl;
    ans << std::endl;
    ans << "Normmisteke" << ' ' << NormVect(x, n) << std::endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Matrix[i][j] = Res[i][j];
        }
        rightb[i] = Resb[i];
    }
    ans << std::endl;
    kest = (int)(log((1 / eps[neps - 1])) / (log(1 / minNC)));
    ans << "kest" << ' ' << kest << std::endl;
    iter = 0;
    Zeidel(n, Matrix, rightb, x, eps[neps - 1]);
    for (i = 0; i < n; i++)
    {
        x[i] = x[i] - Nres[ntest - 1][i];
    }
    ans << std::endl;
    ans << "Normmisteke after kest" << ' ' << NormVect(x, n) << std::endl;
    kest = 1000;
    ans.close();
    for (i = 0; i < n; i++)
    {
        delete[] Matrix[i];
        delete[] Res[i];
        delete[] Inv[i];
        delete[] MC[i];
    }
    delete[] x;
    delete[] Matrix;
    delete[] rightb;
    delete[] Inv;
    delete[] MC;
    delete[] Res;
    delete[] Resb;
}
template <typename T>
void Relaks3diag(int n, T* L, T* D, T* U, T* rightb, T* x, T omega, T eps)
{
    T* xk;
    T* xkl;
    T* res;
    xk = new T[n];
    xkl = new T[n];
    res = new T[n];
    for (int i = 0; i < n; i++)
    {
        xkl[i] = 0.;
        xk[i] = D[i];
        res[i] = 1;
    }
    while (fabs(NormVecinf(xk, n) - NormVecinf(xkl, n)) > eps)
    {
        iter++;
        for (int i = 0; i < n; i++)
        {
            xkl[i] = xk[i];
        }
        xk[0] = omega * rightb[0] / D[0] + (1 - omega) * xkl[0] - omega * U[0] * xkl[1] / D[0];
        for (int i = 1; i < n - 1; i++)
        {
            xk[i] = omega * rightb[i] / D[i] + (1 - omega) * xkl[i] - omega * U[i] * xkl[i + 1] / D[i] - omega * L[i] * xk[i - 1] / D[i];
        }
        xk[n - 1] = omega * rightb[n - 1] / D[n - 1] + (1 - omega) * xkl[n - 1] - omega * L[n - 1] * xk[n - 2] / D[n - 1];
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = xk[i];
    }

    delete[] xk;
    delete[] xkl;
}
template <typename T>
void Progonka(int n, T* L, T* D, T* U, T* rightb, T* x)
{
    T* alpha;
    T* betta;
    T res;
    alpha = new T[n];
    betta = new T[n];
    alpha[0] = -U[0] / D[0];
    betta[0] = rightb[0] / D[0];
    alpha[n - 1] = 0;
    betta[n - 1] = 0;
    for (int i = 1; i < n - 1; i++)
    {
        res = L[i] * alpha[i - 1] + D[i];
        alpha[i] = -U[i] / res;
        betta[i] = (rightb[i] - L[i] * betta[i - 1]) / res;
    }
    x[n - 1] = (rightb[n - 1] - L[n - 1] * betta[n - 2]) / (L[n - 1] * alpha[n - 2] + D[n - 1]);
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = alpha[i] * x[i + 1] + betta[i];
    }
    delete[] alpha;
    delete[] betta;
}
template <typename T>
void test3diag()
{
    T* x;
    T* rightb;
    T* L;
    T* D;
    T* U;
    T  omega = 0.8;
    int n, i;
    n = 221;
    T eps[2] = { 0.0001, 0.0000001 };
    x = new T[n];
    rightb = new T[n];
    L = new T[n];
    D = new T[n];
    U = new T[n];
    for (i = 0; i < n; i++)
    {
        L[i] = 1.;
        U[i] = 1.;
        D[i] = 4.;
        rightb[i] = 10 - 2 * ((i + 1) % 2);
        x[i] = 0;
    }

    rightb[0] = 6;
    rightb[n - 1] = 9 - 3 * (n % 2);
    U[n - 1] = 0.;
    L[0] = 0.;
    std::ofstream ans;

    ans.open("Relaks3diageps1.txt");
    ans << std::setprecision(5);
    iter = 0;
    Relaks3diag(n, L, D, U, rightb, x, omega, eps[1]);
    ans << "Iterations" << ' ' << iter << std::endl;
    ans << "X" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    ans << std::endl;
    ans.close();
    ans.open("Progonkaeps1.txt");
    ans << std::setprecision(5);
    Progonka(n, L, D, U, rightb, x);
    ans << "X" << std::endl;
    for (i = 0; i < n; i++)
    {
        ans << x[i] << ' ';
    }
    ans << std::endl;
    ans.close();


    delete[] D;
    delete[] rightb;
    delete[] L;
    delete[] U;
}
int main()
{
    test<double>(15,1,4);
}
