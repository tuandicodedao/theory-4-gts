#include <iostream>
#include<iomanip>
#include<math.h>
#include <fstream>
#include <stdlib.h>
#define max 100
using namespace std;

bool docmatran(string file, double a[max][max], int &n){ // Hàm đọc ma trận từ file
    n = -1;
    std::ifstream input(file);
    string line;
    // Bỏ qua tên ma trận
    getline( input, line );
    int rows, cols;
    for(rows = 0; getline( input, line ); ++rows)
    {
        string element = "";
        cols = 0;
        for(int t = 0; t < line.length(); ++t){
            if(line[t] != ' '){
                element += line[t];
            }
            if(line[t] == ' ' || t == line.length() - 1){
                a[rows][cols] = atoi(element.c_str());
                cols++;
                element = "";
            }
        }
        if(n == -1) n = cols;
        if(n != cols) return false;
    }
    if(n != rows) return false;
    return true;
}
void Xuatmt(double a[max][max], int n){
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j<n; j++)
            cout << a[i][j]<<setprecision(20)<<fixed<<"\t";
        cout << endl;
    }
}
void gauss(double a[max][max], double xt[max][max], int n)
{
		 double a2[max][max], x[max], ratio;
		 int i,j,k;
		 for(i = 0; i < n; i++)
            for(j = 0; j < n; j++)
                a2[i][j] = a[i][j];

		 //Tạo ma trận bổ sung
		 for(i = 0;i < n; i++)
		 {
			  for(j = 0; j < n; j++)
			  {
				   if(i==j)
				   {
				    	a2[i][j+n] = 1;
				   }
				   else
				   {
				    	a2[i][j+n] = 0;
				   }
			  }
		 }

		 //Sử dụng khử Gauss Jordan
		 for(i = 0;i < n; i++)
		 {
            if(a2[i][i] == 0.0)
			  {
				   cout<<"Khong co ma tran nghich dao";
				   exit(0);
			  }
			  for(j = 0;j < n; j++)
			  {
				   if(i!=j)
				   {
					    ratio = a2[j][i]/a2[i][i];
					    for(k=0; k<2*n; k++)
					    {
					     	a2[j][k] = a2[j][k] - ratio*a2[i][k];
					    }
				   }
			  }
		 }
		 //Tạo đường chéo chính thành 1
		 for(i = 0; i < n;i++)
		 {
			  for(j = n; j <= 2*n-1; j++)
			  {
			   	a2[i][j] = a2[i][j]/a2[i][i];
			  }
		 }
		 // tạo ma trận nghịch đảo
		 for(i = 0; i < n; i++)
		 {
			  for(j = 0; j < n;j++)
			  {
			   	xt[i][j] = a2[i][j+n];
			  }
		 }
}
void Daomatran(int n,double eps ,double a[max][max],double xt[max][max],double xs[max][max])
 {
    int i, j, k;
    double g[max][max],y[max][max];
    double t;
    // Tính ma trận g0 = E - a*x0
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
        {
            g[i][j] = 0;
            for(k = 0; k < n; k++)
            {
                g[i][j] = g[i][j] + a[i][k]*xt[k][j];
            }
        }
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if (i==j)
                {
                g[i][j] = 1. - g[i][j];
                }
            else
                {
                g[i][j] = - g[i][j];
                }
        }
    }
    //Tính ||X0||
    double cxt = 0;
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            cxt += xt[i][j]*xt[i][j];
    cxt = sqrt(cxt);
    // Kiem tra dieu kien hoi tu cua ma tran lap la ||g0|| <1
    t=0.;
    for(i = 0; i < n; i++)
    for(j = 0; j < n; j++) t += g[i][j]*g[i][j];
    t = sqrt(t);
    if (t >= 1)
    {
    printf("\tPhep lap khong hoi tu! \n");
    return;
    }
    else //Tinh lap ma tran dao
    {
    int dem = 0;

    while (cxt*pow(t,pow(2,dem))> eps) // Điều kiện lặp ||X0||*(q^(2^k) > epsi
    {
        //Tạo ma trận y = 2E-AX
        for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
        {
            y[i][j] = 0;
        for(k = 0; k < n; k++)
        y[i][j] += a[i][k] * xt[k][j];}

        for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
        {
            if(i == j)
                y[i][j] = 2-y[i][j];
            else y[i][j] = -y[i][j];
        }
        // Ma trận lặp X = X*y = X(2E-AX)
        for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
        {
            xs[i][j] = 0;
            for(k = 0; k<n; k++)
                xs[i][j] += xt[i][k] * y[k][j];
        };
        for(i = 0; i < n; i++)
        for(j = 0; j < n; j++) xt[i][j] = xs[i][j];
        dem++;
    }
    cout << "Ma tran ket thuc sau " << dem << " lan lap"<<endl;
    cout << "Ma tran nghich dao la: "<<endl;
    if(dem == 0)
    {Xuatmt(xt,n); cout << endl;}
    else
    {Xuatmt(xs,n); cout << endl;}
    }
}

int main()
{
    int n;
    double a[max][max], xt[max][max], xs[max][max], eps;
    if(docmatran("matrix3x3.txt", a, n))
    {
        cout << "Nhap sai so : "; cin >> eps;
        cout << "Ma tran ban dau: "<<endl;
        Xuatmt(a,n); cout << endl;
        gauss(a,xt,n);
        Daomatran(n,eps,a,xt,xs);
    }
}
