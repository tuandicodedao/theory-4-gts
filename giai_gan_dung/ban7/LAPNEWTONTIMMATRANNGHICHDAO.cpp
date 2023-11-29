#include <stdio.h>
#include <math.h>

void NhapA(double a[100][100], int* m)
{
    FILE* fp;
    int i, j;
    fp = fopen("LAPNEWTON.txt", "r");

    fscanf(fp, "%d\n", m);
    fscanf(fp, "%d\n", m);

    printf("So hang cua ma tran: %d\n", *m);
    printf("So cot cua ma tran: %d\n", *m);

    for (i = 1; i <= *m; i++)
    {
        for (j = 1; j <= *m; j++)
        {
            fscanf(fp, "%lf\n", &a[i][j]);
        }
    }
}

void NhapX0(double a[100][100], int* m)
{
    FILE* fp;
    int i, j;
    fp = fopen("LAPNEWTON2.txt", "r");

    fscanf(fp, "%d\n", m);
    fscanf(fp, "%d\n", m);
    printf("Ma Tran Nghich Dao Tim Theo G-J:\n");

    for (i = 1; i <= *m; i++)
    {
        for (j = 1; j <= *m; j++)
        {
            fscanf(fp, "%lf\n", &a[i][j]);
        }
    }
}

void InMt(double a[100][100], int m)
{
    int i, j;
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= m; j++)
        {
            printf("%.10lf\t", a[i][j]);
        }
        printf("\n\n");
    }
}

double TimChuan(double a[100][100], int m)
{
    double chuan = 0;
    for (int i = 1; i <= m; i++)
    {
        double tong = 0;
        for (int j = 1; j <= m; j++)
        {
            tong += fabs(a[i][j]);
        }
        if (tong > chuan)
        {
            chuan = tong;
        }
    }
    return chuan;
}

int main()
{
    double A[100][100];
    double X0[100][100];
    int m;
    NhapA(A, &m);
    InMt(A, m);
    NhapX0(X0, &m);
    InMt(X0, m);
    printf("\n");
    double e = 1e-20;
    printf("epsilon: %.20lf\n", e);

    double AX0[100][100];
    int i, j, k;
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= m; j++)
        {
            AX0[i][j] = 0;
            for (k = 1; k <= m; k++)
            {
                AX0[i][j] += A[i][k] * X0[k][j];
            }
        }
    }
    printf("Tich cua AX0: \n");
    InMt(AX0, m);

    double E[100][100];
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= m; j++)
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
    printf("Ma Tran Don Vi: \n");
    InMt(E, m);

    double B[100][100];
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= m; j++)
        {
            B[i][j] = E[i][j] - AX0[i][j];
        }
    }
    printf("Hieu E-AX0: \n");
    InMt(B, m);

    double q = TimChuan(B, m);
    printf("Chuan cua E-AX0: %lf\n", q);

    double X[100][100];
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= m; j++)
        {
            X[i][j] = X0[i][j];
        }
    }
    double p = TimChuan(X, m);
    int l = 0;
    while (p * pow(q, pow(2, l)) / (1 - q) > e)
    {
        printf("Lan Lap %d: %.14lf\n", l + 1, p * pow(q, pow(2, l)) / (1 - q));

        double AX[100][100];
        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= m; j++)
            {
                AX[i][j] = 0;
                for (k = 1; k <= m; k++)
                {
                    AX[i][j] += A[i][k] * X[k][j];
                }
            }
        }

        double U[100][100];
        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= m; j++)
            {
                U[i][j] = E[i][j] * 2;
            }
        }

        double C[100][100];
        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= m; j++)
            {
                C[i][j] = U[i][j] - AX[i][j];
            }
        }

        double D[100][100];
        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= m; j++)
            {
                D[i][j] = 0;
                for (k = 1; k <= m; k++)
                {
                    D[i][j] += X[i][k] * C[k][j];
                }
            }
        }

        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= m; j++)
            {
                X[i][j] = D[i][j];
            }
        }

        l++;
        printf("Ma tran X:\n");
        InMt(X, m);
    }
    printf("Ma Tran Nghich Dao Cua A La: \n");
    InMt(X, m);

    printf("//////////////////////\n");
    printf("Kiem tra lai, neu tich cua A va X la ma tran don vi thi true!\n");
    double matrantich[100][100];
    for (i = 1; i <= m; i++)
    {
        for (j = 1; j <= m; j++)
        {
            matrantich[i][j] = 0;
            for (k = 1; k <= m; k++)
            {
                matrantich[i][j] += A[i][k] * X[k][j];
            }
        }
    }
    InMt(matrantich, m);
}
