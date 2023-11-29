using namespace std;
#define oo INT_MAX
#define LD long double
#include <bits/stdc++.h>
#include "matrix.hpp"

vector<int> Index;
int NormType;
LD L;
LD eps;
ifstream fi("A.inp");

ofstream fo("A.csv");
//
template <typename T>
LD NormInf(const matrix<T> &A)
{
    LD tmp = 0;
    LD res = 0;
    for (int i = 0; i < A.row; ++i)
    {
        for (int j = 0; j < A.col; ++j)
            tmp += abs(A[i][j]);
        res = max(res, tmp);
        tmp = 0;
    }
    return res;
}

template <typename T>
LD Norm1(const matrix<T> &A)
{
    LD tmp = 0;
    LD res = 0;
    for (int j = 0; j < A.col; ++j)
    {
        for (int i = 0; i < A.row; ++i)
            tmp += abs(A[i][j]);
        res = max(res, tmp);
        tmp = 0;
    }
    return res;
}

template <typename T>
LD Norm(const matrix<T> &A)
{
    return (NormType == oo) ? NormInf(A) : Norm1(A);
}
//
template <typename T>
bool Check(const matrix<T> &A)
{
    bool Mark[A.col] = {};
    LD Sum;
    for (int i = 0; i < A.row; i++)
    {
        Index[i] = 0;
        Sum = abs(A[i][0]);
        for (int j = 1; j < A.col; j++)
        {
            Sum += abs(A[i][j]);
            if (abs(A[i][j]) > abs(A[i][Index[i]]))
                Index[i] = j;
        }
        if (Mark[Index[i]] || abs(A[i][Index[i]]) * 2 <= Sum)
            return false;
        else
            Mark[Index[i]] = i;
    }
    return true;
}
//
template <typename T>
bool Normalize(matrix<T> &A, matrix<T> &B)
{
    Index.resize(A.row);
    if (!Check(A))
        if (Check(!A))
        {
            LD maxAii = 0, minAii = abs(A[0][0]);
            vector<int> Temp = Index;
            for (int i = 0; i < A.row; i++)
            {
                Index[Temp[i]] = i;
                maxAii = max(maxAii, abs(A[i][Index[i]]));
                minAii = min(minAii, abs(A[i][Index[i]]));
            }
            NormType = 1;
            L = maxAii / minAii;
        }
        else
            return false;
    else
    {
        NormType = oo;
        L = 1.0;
    }

    matrix<T> _A = A;
    matrix<T> _B = B;
    for (int i = 0; i < A.row; i++)
    {
        A[Index[i]] = _A[i];
        B[Index[i]] = _B[i];
    }
    return true;
}

template <typename T>
T Pow(const T &A, int n)
{
    if (n == 1)
        return A;
    T tmp = Pow(A, n >> 1);
    return (n & 1) ? tmp * tmp * A : tmp * tmp;
}
//OUTPUT
template <typename T>
void Print(const T &X)
{
    fo << fixed << setprecision(-log10(abs(eps)) + 1) << X << '\n';
}
//Lap don
template <typename T>
void SingleLoop(const matrix <T> &A, const matrix<T>A1, const matrix<T> &Alpha, const matrix<T> &Beta, const matrix<T> &D)
{
    matrix<T> X0(Alpha.col, Beta.col);
    matrix<T> X1 = Beta;
    matrix<T> _X0 = X0;
    matrix<T> _X1 = X1;

    LD q = Norm(Alpha);

    //Hau nghiem                O(n^2 * m * cnt)
    fo << "So buoc lap theo Hau Nghiem: ";
    int cnt = 1;
    while (Norm(X1 - X0) * L * q / (1 - q) > eps)
    {
        X0 = X1;
        X1 = Alpha * X1 + Beta;
        cnt++;
    }
    fo << cnt << '\n';
    fo <<"Ma tran nghich dao la: "<<'\n';
    Print(D * X1); //neu la chuan inf thi D = I
    fo << "Kiem tra nhan nguoc lai: "<<'\n';
    Print(D * X1 * A1);

    //Tien nghiem             O(n^2 * m * cnt)
    fo << "So buoc lap theo Tien Nghiem: ";
    int _cnt = ceil(log2(eps * (1 - q) / (Norm(_X1 - _X0) * L)) / log2(q));
    for (int i = 2; i <= _cnt; i++)
    {
        _X0 = _X1;
        _X1 = Alpha * _X1 + Beta;
    }
    fo << _cnt << '\n';
    fo << "Ma tran nghich dao la: "<<'\n';
    Print(D * _X1);
    fo << "Kiem tra nhan nguoc lai: "<<'\n';
    Print(D * _X1 * A1);
    // Nhan nguoc lai

}
//Jacobi
template <typename T>
void Jacobi(const matrix<T> &A, const matrix<T>A1, const matrix<T> &B)
{
    matrix<T> D(A.row, A.col); //Chu T minh dat cho typename roi nen gio thay chu D la ma tran diag nhe
    matrix<T> I(A.row, A.col);
    Id(I); //Ma tran don vi
    for (int i = 0; i < A.row; ++i)
        D[i][i] = (T)1 / (A[i][i]);
    if (NormType == oo)
        SingleLoop(A, A1, I - D *A, D * B, I);

    if (NormType == 1)
        SingleLoop(A, A1, I - A * D, B, D);
}

template<typename T>
void Check(matrix<T> &A, matrix<T> &B, LD &eps)//ham kiem tra dau vao
{
    if (A.col != A.row)
    {
        fo << "Ma tran khong vuong!!";
        exit(0);
    }
    if (eps < 0)
    {
        fo << "Sai so khong hop le";
        exit(0);
    }
    if (!Normalize(A, B))
    {
        fo << "INPUT khong hop le";
        exit(0);
    }
}

int main()
{
    //INPUT
    int n, m;
    fi >> n >> m >> eps;
    matrix<LD> A(n, m), B(n, m) ;//thay CLD = LD de dung voi so thuc
    for ( int i= 0; i< n; i++)
    {
    	for ( int j= 0; j< m; j++)
    	{
    		if ( j==i) B[i][j]= 1;
    		else B[i][j]= 0;
		}
	}
    fi >> A;
    matrix<LD>A1 = A;
    fo <<"Ma tran A: ";
    fo <<'\n'<< A << '\n';
    Check(A, B, eps);
    //WORK
    Jacobi(A, A1, B);
    return 0;
}
