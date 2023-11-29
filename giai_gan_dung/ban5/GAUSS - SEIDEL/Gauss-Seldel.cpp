using namespace std;
#define oo INT_MAX
#define LD long double
#define CLD complex<long double>
#include <bits/stdc++.h>
#include "matrix.hpp"
vector<int> Index;
int NormType;
LD eps;
ifstream fi("A.inp");
ofstream fo("A.csv");
//Các loại chuẩn
template <typename T>
LD NormInf(const matrix<T> &A)//tính chuẩn hàng
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
LD Norm1(const matrix<T> &A)//tính chuẩn cột
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

template <typename T>//Hàm trả về loại chuẩn của ma trận
LD Norm(const matrix<T> &A)
{
    return (NormType == oo) ? NormInf(A) : Norm1(A);
}
//Kiểm tra trội hàng
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
//Kiểm tra dạng chéo trội và ưu tiên kiểm tra trội hàng trước
template <typename T>
bool Normalize(matrix<T> &A, matrix<T> &B)
{
    Index.resize(A.row);
    if (!Check(A))     //Nếu không trội hàng
       {
           if (Check(!A)) //Nếu trội cột (!A là chuyển vị của A)
        {
            vector<int> Temp = Index;
            for (int i = 0; i < A.row; i++) Index[Temp[i]] = i;
            NormType = 1;        //Chuẩn 1
        }
        else
            return false; //A không có dạng chéo trội
      }
    else
        NormType = oo; //Chuẩn vô cùng

    matrix<T> _A = A;
     matrix<T> _B = B;
    for (int i = 0; i < A.row; i++) //Đổi chỗ các hàng để đưa về chéo trội với trường hợp các phần tử trội xếp lẫn lộn
    {
        A[Index[i]] = _A[i];
        B[Index[i]] = _B[i];
    }
    return true;
}
//OUTPUT
template <typename T>
void Print(const T &X)//Đưa ra ma trận X
{
    fo << fixed << setprecision(-log10(abs(eps)) + 1) << X << '\n';
}

template <typename T>
matrix<T> Lower(const matrix<T> &A)//ma trận Lower để áp dụng tính S
{
    matrix<T> C(A.row, A.col);
    for (int i = 0; i < A.row; i++)
        for (int j = 0; j <= i; j++)
            C[i][j] = A[i][j];
    return C;
}
//Hàm đưa ra ma trận tiếp theo bằng pp lặp Gauss-Seidel
template <typename T>
void Calc(const matrix<T> &A, const matrix<T> &B, matrix<T> &X)
{
    for (int i = 0; i < A.row; i++)
    {
        X[i] = B[i];
        for (int j = 0; j < A.col; j++)
            if (i != j)
            X[i] = X[i] +  X[j] * A[i][j];
    }
}
//Hậu nghiệm
template <typename T>
void Loop1(const matrix<T> &Alpha, const matrix<T> &Beta, const matrix<T> &_D,const matrix<T>A1, LD q, LD S)
{
    matrix<T> X0(Alpha.col, Beta.col); //X0 = 0
    matrix<T> X1 = Beta;
    fo << "So buoc lap theo Hau Nghiem: ";
    int cnt = 1;
    while (Norm(X1 - X0) * q / ((1 - q)*(1-S)) > eps)
    {
        X0 = X1;
        Calc(Alpha, Beta, X1);
        cnt++;
    }
    fo << cnt << '\n';
    fo<<"Ma tran nghich dao la : "<<'\n';
    Print(_D * X1);
    fo<<"Kiem tra nhan nguoc lai :"<<'\n';
    Print(_D * X1 *A1);
}
//Tiên nghiệm
template <typename T>
void Loop2(const matrix<T> &Alpha, const matrix<T> &Beta, const matrix<T> &_D,const matrix<T>A1, LD q, LD S)
{
    matrix<T> X0(Alpha.col, Beta.col); //X0 = 0
    matrix<T> X1 = Beta;
    fo << "So buoc lap theo Tien Nghiem: ";
    int _cnt = ceil(log2(eps * (1 - q) * (1 - S) / (Norm(X1 - X0))) / log2(q));
    for (int i = 2; i <= _cnt; i++)
    {
        X0 = X1;
        Calc(Alpha, Beta, X1);
    }
    fo << _cnt << '\n';
    fo<<"Ma tran nghich dao la"<<'\n';
    Print(_D * X1);
    fo<<"Kiem tra nhan nguoc :"<<'\n';
    Print(_D * X1 * A1);
}

template <typename T>
double He_So_Co(const matrix<T> &Alpha)  //dpt hàm này là O(n*n)
{
    int n = Alpha.row;
    double P, Q, lambda = 0;
    matrix<double> Sum(n+1, n+1);
    for (int i = 1; i < n+1; i++)
        for (int j = 1; j < n+1; j++)
            Sum[i][j] = Sum[i-1][j] + Sum[i][j-1] - Sum[i-1][j-1] + abs(Alpha[i-1][j-1]);

    if (NormType == oo)
    {
        for (int i = 1; i < n+1; i++)
        {
            P = Sum[i][i-1] - Sum[i-1][i-1];
            Q = Sum[i][n]   - Sum[i-1][n]   - P;
            lambda = max(lambda , Q / (1 - P));
        }
        return lambda;
    }
    else
    {
        for (int j = 1; j < n+1; j++)
        {
            Q = Sum[j][j] - Sum[j][j-1];
            P = Sum[n][j] - Sum[n][j-1] - Q;
            lambda = max(lambda , Q / (1 - P));
        }
        return lambda;
    }
}
//Gaus_Seidel
template <typename T>
void Gauss_Seidel(const matrix<T> &A, const matrix<T> &B, const matrix<T>A1)
{
    matrix<T>  I(A.row, A.col); Id(I);  //Ma trận đơn vị
    matrix<T>  D(A.row, A.col);         //Ma trận diag
    matrix<T> Alpha, Beta;
    LD q,S;

    for (int i = 0; i < A.row; ++i)
         D[i][i] = (T)1 / (A[i][i]);

    if (NormType == oo)//Nếu ma trận chuẩn hàng
    {
        Alpha = I - D * A;
        Beta  = D * B;
        q = He_So_Co(Alpha);

        Loop1(Alpha, Beta, I, A1,q,0);     //hau nghiem
        Loop2(Alpha, Beta, I,A1, q, 0);  //tiên nghiem
    }
    else//ma trận chuẩn cột
    {
        Alpha = I - A * D;
        Beta  = B;
        q = He_So_Co(Alpha);
        S = Norm1(Lower(Alpha));
        Loop1(Alpha, Beta, D,A1,q, S); //hậu nghiệm
        Loop2(Alpha, Beta, D,A1,q, S); //tien nghiem
    }
}
//Hàm kiểm tra đầu vào
template<typename T>
void Checkone(matrix<T> &A,matrix<T> &B, LD &eps)
{
    if(A.col != A.row)
    {
        fo<<"Ma tran khong vuong, khong thoa man";
        exit(0);
    }
    if(eps<0)
    {
        fo<<"Sai so khong hop le";
        exit(0);
    }
    if(!Normalize(A,B))
    {
        fo<<"Input khong hop le ";
        exit(0);
    }
}
int main()
{
    //INPUT
    int n, m;
    fi >> n >> m >> eps;
    matrix<LD> A(n, m), B(n, m);

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            if(j==i) B[i][j]=1;
            else B[i][j]=0;
        }
    }
    fi >> A;
    matrix<LD>A1 = A;
    fo<<"Ma tran A: ";
    fo<<'\n'<<A<<'\n';
    Checkone(A,B,eps);
    //WORK
    Gauss_Seidel(A,B,A1);
    return 0;
}
