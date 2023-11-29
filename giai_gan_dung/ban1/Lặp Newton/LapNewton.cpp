#include <stdio.h>
#include <math.h>

#define MaxN 30
int k = 0;

void Tru_Ma_Tran(int size, double A[MaxN][MaxN], double B[MaxN][MaxN], double C[MaxN][MaxN]);  
void Chep_ma_tran(int size, double A[MaxN][MaxN], double B[MaxN][MaxN]);                       
void chuyen_vi(int size, double A[MaxN][MaxN], double B[MaxN][MaxN]);                          
void Nhan_1_so_ma_tran(double num, int size, double A[MaxN][MaxN]);                            
void Nhan_Ma_Tran(int size, double A[MaxN][MaxN], double B[MaxN][MaxN], double C[MaxN][MaxN]); 
void In_Ma_Tran(int size, double A[MaxN][MaxN]);
double Chuan_1(int size, double A[MaxN][MaxN]);
double Chuan_2(int size, double A[MaxN][MaxN]);
double Chuan_vo_cung(int size, double A[MaxN][MaxN]);
double det_A(int size, double A[MaxN][MaxN]);
//===================================//
void Xap_xi_dau(int size, double A[MaxN][MaxN], double Xo[MaxN][MaxN]);
void Ma_tran_nghich(int size, double A[MaxN][MaxN], double Xo[MaxN][MaxN], double ep);
void Gauss(int size, double A[MaxN][MaxN], double Xo[MaxN][MaxN]);

int main()
{
    int size;
    printf("Nhap size cua ma tran A : ");
    scanf("%d", &size);
    double A[MaxN][MaxN], Xo[MaxN][MaxN], ep;
    printf("Nhap epxilon : ");
    scanf("%lf", &ep);

    //===================================//
    FILE *fptr, *fptr3;
    fptr = fopen("input.txt", "r");
    if (fptr == NULL)
    {
        printf("ERROR");
        return 1;
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            fscanf(fptr, "%lf", &A[i][j]);
        }
    }
    //===================================//

    // Kiem tra dieu kien co ma tran nghich dao
    if (det_A(size, A) == 0)
    {
        printf("Khong co ma tran nghich dao !");
        return 0;
    }

    // Ma tran dau vao
    printf("Ma tran A : \n");
    In_Ma_Tran(size, A);

    Xap_xi_dau(size, A, Xo);
    Ma_tran_nghich(size, A, Xo, ep);

    fclose(fptr);
}

void In_Ma_Tran(int size, double A[MaxN][MaxN])
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("  %.15lf  ", A[i][j]);
        }
        printf("\n");
    }
}

void Nhan_Ma_Tran(int size, double A[MaxN][MaxN], double B[MaxN][MaxN], double C[MaxN][MaxN])
{
    for (int i = 0; i < size; i++)
    { // size
        for (int j = 0; j < size; j++)
        { // col
            float sum = 0;
            for (int k = 0; k < size; k++)
            { // size
                sum = sum + A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

void Nhan_1_so_ma_tran(double num, int size, double A[MaxN][MaxN])
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            A[i][j] = num * A[i][j];
        }
    }
}

double Chuan_1(int size, double A[MaxN][MaxN])
{
    double max = -10000, temp = 0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            temp = temp + fabs(A[j][i]);
        }
        if (temp > max)
            max = temp;
        temp = 0;
    }
    return max;
}

double Chuan_2(int size, double A[MaxN][MaxN])
{
    double sum = 0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            sum = sum + A[i][j] * A[i][j];
        }
    }
    sum = sqrt(sum);
    return sum;
}

double Chuan_vo_cung(int size, double A[MaxN][MaxN])
{
    double max = -10000, temp = 0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            temp = temp + fabs(A[i][j]);
        }
        if (max < temp)
            max = temp;
        temp = 0;
    }
    return max;
}

void chuyen_vi(int size, double A[MaxN][MaxN], double B[MaxN][MaxN])
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            B[j][i] = A[i][j];
        }
    }
}

void Xap_xi_dau(int size, double A[MaxN][MaxN], double Xo[MaxN][MaxN])
{
    double E[MaxN][MaxN], ma_tran1[MaxN][MaxN], ma_tran2[MaxN][MaxN], ma_tran3[MaxN][MaxN], ma_tran4[MaxN][MaxN], ma_tran_2E[MaxN][MaxN]; // Ma tran don vi E
    double ma_tran5[MaxN][MaxN];
    int attempt = 0, max_attempt = 2;

    // Tao ma tran E
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                E[i][j] = 1;
            }
            else
            {
                E[i][j] = 0;
            }
        }
    }

    // Ma tran 2E
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            ma_tran_2E[i][j] = 2 * E[i][j];
        }
    }

    double t1, t2;
    t1 = Chuan_1(size, A); // t1

    t2 = Chuan_vo_cung(size, A); // t2

    Chep_ma_tran(size, A, ma_tran1);
    Nhan_1_so_ma_tran(1 / (t1 * t2), size, ma_tran1); // tinh A/(t1*t2)=ma_tran1

    chuyen_vi(size, ma_tran1, Xo); // Tinh (A/(t1*t2))^T=Xo

    while (attempt < max_attempt)
    {
        Nhan_Ma_Tran(size, A, Xo, ma_tran2); // ma_tran2=A*Xo;

        Tru_Ma_Tran(size, ma_tran_2E, ma_tran2, ma_tran4); // 2*E-A*Xo=ma_tran4;

        Nhan_Ma_Tran(size, Xo, ma_tran4, Xo);
        /*printf("\n");
        In_Ma_Tran(size, Xo);
        printf("\n");*/

        // Tinh ma tran A*Xo;
        Nhan_Ma_Tran(size, A, Xo, ma_tran3); // ma_tran3=A*Xo

        // Tinh E-A*Xo;
        Tru_Ma_Tran(size, E, ma_tran3, ma_tran5);

        if (Chuan_2(size, ma_tran5) < 1)
        {
            attempt++;
        }
    }
    printf("\n");
    printf("Ma tran xap xi Xo : \n\n");
    In_Ma_Tran(size, Xo);
    printf("\n");
}
void Chep_ma_tran(int size, double A[MaxN][MaxN], double B[MaxN][MaxN])
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            B[i][j] = A[i][j];
        }
    }
}
void Tru_Ma_Tran(int size, double A[MaxN][MaxN], double B[MaxN][MaxN], double C[MaxN][MaxN]) // C la ma tran ket qua
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

void Ma_tran_nghich(int size, double A[MaxN][MaxN], double Xo[MaxN][MaxN], double ep)
{
    int k = 0;
    double q, t;
    double E[MaxN][MaxN], ma_tran1[MaxN][MaxN], ma_tran2[MaxN][MaxN], ma_tran3[MaxN][MaxN], ma_tran4[MaxN][MaxN], ma_tran6[MaxN][MaxN], ma_tran_2E[MaxN][MaxN]; // Ma tran don vi E
    double X[MaxN][MaxN];

    // Khoi tao ma tran E
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                E[i][j] = 1;
            }
            else
            {
                E[i][j] = 0;
            }
        }
    }

    // Tinh ma tran 2E
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            ma_tran_2E[i][j] = 2 * E[i][j];
        }
    }

    // Tinh A*Xo
    Nhan_Ma_Tran(size, A, Xo, ma_tran1); // ma_tran1=A*Xo;

    // Tinh E-A*Xo
    Tru_Ma_Tran(size, E, ma_tran1, ma_tran2); // ma_tran2=E-A*Xo;

    // Tinh q=||E-A*Xo||
    q = Chuan_2(size, ma_tran2);
    /*printf("\nq = %.15lf", q);
    printf("\n\n");*/

    // Tinh t=||Xo||
    t = Chuan_2(size, Xo);
    /*printf("t = %.15lf", t);
    printf("\n");*/

    // X=Xo
    Chep_ma_tran(size, Xo, X);

    // k=0;
    k = 0;
    do
    {
        Nhan_Ma_Tran(size, A, X, ma_tran3); // ma_tran3= A*X

        Tru_Ma_Tran(size, ma_tran_2E, ma_tran3, ma_tran4); // ma_tran4=2*E-A*X;

        Nhan_Ma_Tran(size, X, ma_tran4, X); // X=X*(2*E-A*X);

        k++;
        //printf("\n Lan lap thu %d \n", k);
        //In_Ma_Tran(size, X);
    } while (t * pow(q, pow(2, k)) / (1 - q) > ep);

    printf("\nMa tran nghich dao cua A : \n");
    printf("\n");
    In_Ma_Tran(size, X);
    printf("\n");
    printf(" So lan lap : %d", k);

    //Kiểm tra nhân ngược
    Nhan_Ma_Tran(size, A, X, ma_tran6);
    printf("\n");
    printf("\nKiem tra nhan nguoc : \n");
    printf("\n");
    In_Ma_Tran(size, ma_tran6);
    
}
double det_A(int size, double A[MaxN][MaxN])
{
    double A1[MaxN][MaxN], ratio, det = 1;
    int temp;
    Chep_ma_tran(size, A, A1);
    for (int i = 0; i < size; i++)
    {
        if (A1[i][i] == 0)
        {
            for (int j = i + 1; j < size; j++)
            {
                if (A1[j][j] != 0)
                {
                    for (int k = 0; k < size; k++)
                    {
                        temp = A1[i][k];
                        A1[i][k] = A1[j][k];
                        A1[j][k] = temp;
                    }
                }
            }
        }
        for (int j = i + 1; j < size; j++)
        {
            ratio = A1[j][i] / A1[i][i];
            for (int k = i; k < size; k++)
            {
                A1[j][k] = A1[j][k] - ratio * A1[i][k];
            }
        }
    }
    for (int i = 0; i < size; i++)
    {
        det = det * A1[i][i];
    }
    printf("\n Det A : %lf \n ", det);
    return det;
}
