#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define Na 1000 // block have 1000 bits

unsigned long long SEED = 312891;
unsigned long long RANV;
int RANI = 0;
double Ranq1();
void normal(double sig, double *n1, double *n2);
void table_symbol(int i, int a, int b, int c, int d);
void x_estimate_m(int i, double a, double b);
void table_receive(int i, double a, double b);

double **s;
int srow = 250;  // Na / 4 transmitt how many symbols (transmitter)
int scolumn = 2; //1 symbol vector have 2 index

double **m_estimate;
int mrow = 250;  // Na / 4 transmitt estimate symbols (receiver+estimate)
int mcolumn = 2; //1 symbol vector have 2 index

int **receive_sym;

int receivesymrow = 250; // Na/4 transmitt how many symbols
int receivesymcolumn = 4; //16QAM PRESENT 1symbol = 4bits

int main()
{
    //int Na = 1000;  // block have 1000 bits
    int Nf; // Nf receive 100;   // receive Nf error block
    int i, j;

    int *a; // denote as 1000 bits
    int alen = Na;
    double z;
    double x, y;
    a = (int *)malloc(alen * sizeof(int)); // generate a block of Na random bits

    int **sym;         // present symbol ( 4bits to 1symbols)
    int symrow = srow; // Na/4 transmitt how many symbols
    int symcolumn = 4; //16QAM PRESENT 1symbol = 4bits
    sym = (int **)malloc(symrow * sizeof(int *));
    for (i = 0; i < symrow; i++)
        sym[i] = (int *)malloc(symcolumn * sizeof(int));

    //double **s;
    //int srow = Na/4; // transmitt how many symbols
    //int scolumn = 2; //1 symbol vector have 2 index
    s = (double **)malloc(srow * sizeof(double *));
    for (i = 0; i < srow; i++)
        s[i] = (double *)malloc(scolumn * sizeof(double));
    
    double **x_receive;
    int xrow = 250;  // Na / 4 receive how many symbols (receiver)
    int xcolumn = 2; //1 symbol vector have 2 index
    x_receive = (double **)malloc(xrow * sizeof(double *));
    for (i = 0; i < xrow; i++)
        x_receive[i] = (double *)malloc(xcolumn * sizeof(double));

    m_estimate = (double **)malloc(mrow * sizeof(double *));
    for (i = 0; i < mrow; i++)
        m_estimate[i] = (double *)malloc(mcolumn * sizeof(double));

    receive_sym = (int **)malloc(receivesymrow * sizeof(int *));
    for (i = 0; i < receivesymrow; i++)
        receive_sym[i] = (int *)malloc(receivesymcolumn * sizeof(int));
    int *recieve_a;
    recieve_a = (int *)malloc(alen * sizeof(int));

    double ebn0s = 13;
    double sigma;
    //sigma = sqrt(5.0 / (pow(10, ebn0s / 10)));
    sigma = sqrt(1.25 / (pow(10, ebn0s / 10)));
    printf("%g,%g \n", pow(10, ebn0s / 10), sigma);

    double error_bit;
    double BER;
    int num = 0;

    while (Nf < 100) {
        num++;
        for (i = 0; i < Na; i++) { // using ranq1 to produces bits
            z = Ranq1();
            if (z >= 0.5)
                a[i] = 1;
            else
                a[i] = 0;
            // printf("a[%d] = %d;",i,a[i]);
        }
        for (i = 0; i < symrow; i++)
            for (j = 0; j < symcolumn; j++) {
                sym[i][j] = a[4 * i + j];
                if (i == 0 || i == 10)
                    printf("sym[%d][%d] = %d; ", i, j, sym[i][j]);
            }
        for (i = 0; i < srow; i++)
            table_symbol(i, sym[i][0], sym[i][1], sym[i][2], sym[i][3]);
        for (i = 0; i < srow; i++)
            for (j = 0; j < scolumn; j++)
                if (i == 0 || i == 10)
                    printf("s[%d][%d] = %g\n", i, j, s[i][j]);
        for(i = 0; i < srow; i++) {
            normal(sigma, &x, &y);
            x_receive[i][0] = s[i][0] + x;
            x_receive[i][1] = s[i][1] + y;
        }
        for (i = 0; i < xrow; i++)
            for (j = 0; j < xcolumn; j++)
                if (i == 0 || i == 10)
                    printf("x_receive[%d][%d] = %g\n", i, j, x_receive[i][j]);
        for (i = 0; i < xrow; i++) x_estimate_m(i,x_receive[i][0],x_receive[i][1]);
        for (i = 0; i < mrow; i++)
            for (j = 0; j < mcolumn; j++)
                if (i == 0 || i == 10)
                    printf("m_estimate[%d][%d] = %g\n", i, j,m_estimate[i][j]);
        for (i = 0; i < receivesymrow; i++) {
            table_receive(i, m_estimate[i][0],m_estimate[i][1]);
            if ( i == 0 || i == 1) printf("receive_sym[%d]= %d %d %d %d\n",i,receive_sym[i][0],receive_sym[i][1],receive_sym[i][2],receive_sym[i][3]);
        }
        for (i = 0; i < alen; i++) recieve_a[i] = receive_sym[i/4][i%4];
        printf("receive_a = %d %d %d %d\n",recieve_a[4],recieve_a[5],recieve_a[6],recieve_a[7]);
        for (i = 0; i < alen; i++)
            if (a[i] != recieve_a[i]) error_bit++;

        if (error_bit != 0)Nf++;
    }
    BER = error_bit / (num * 1000);
    printf("BER = %g \n", BER);
    FILE *outfp3; 
        outfp3 = fopen("QAM16_SNR_change_sigma.txt","a");
        fprintf(outfp3,"SNR = %g ; BER = %g \n", ebn0s, BER);
    fclose(outfp3);
    return 0;
}
double Ranq1()
{
    if (RANI == 0)
    {
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;

    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

void normal(double sig, double *n1, double *n2)
{
    double x1, x2;
    double s;
    //printf("sigma = %g\n", sig);
    do
    {
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);
    *n1 = sig * x1 * sqrt((-2.0 * log(s)) / s);
    *n2 = sig * x2 * sqrt((-2.0 * log(s)) / s);
}
void table_symbol(int i, int a, int b, int c, int d)
{
    if (a == 1 && b == 1 && c == 0 && d == 1)
    {
        s[i][0] = -3;
        s[i][1] = 3;
    }
    else if (a == 1 && b == 1 && c == 0 && d == 0)
    {
        s[i][0] = -1;
        s[i][1] = 3;
    }
    else if (a == 1 && b == 1 && c == 1 && d == 0)
    {
        s[i][0] = 1;
        s[i][1] = 3;
    }
    else if (a == 1 && b == 1 && c == 1 && d == 1)
    {
        s[i][0] = 3;
        s[i][1] = 3;
    }
    else if (a == 1 && b == 0 && c == 0 && d == 1)
    {
        s[i][0] = -3;
        s[i][1] = 1;
    }
    else if (a == 1 && b == 0 && c == 0 && d == 0)
    {
        s[i][0] = -1;
        s[i][1] = 1;
    }
    else if (a == 1 && b == 0 && c == 1 && d == 0)
    {
        s[i][0] = 1;
        s[i][1] = 1;
    }
    else if (a == 1 && b == 0 && c == 1 && d == 1)
    {
        s[i][0] = 3;
        s[i][1] = 1;
    }
    else if (a == 0 && b == 0 && c == 0 && d == 1)
    {
        s[i][0] = -3;
        s[i][1] = -1;
    }
    else if (a == 0 && b == 0 && c == 0 && d == 0)
    {
        s[i][0] = -1;
        s[i][1] = -1;
    }
    else if (a == 0 && b == 0 && c == 1 && d == 0)
    {
        s[i][0] = 1;
        s[i][1] = -1;
    }
    else if (a == 0 && b == 0 && c == 1 && d == 1)
    {
        s[i][0] = 3;
        s[i][1] = -1;
    }
    else if (a == 0 && b == 1 && c == 0 && d == 1)
    {
        s[i][0] = -3;
        s[i][1] = -3;
    }
    else if (a == 0 && b == 1 && c == 0 && d == 0)
    {
        s[i][0] = -1;
        s[i][1] = -3;
    }
    else if (a == 0 && b == 1 && c == 1 && d == 0)
    {
        s[i][0] = 1;
        s[i][1] = -3;
    }
    else if (a == 0 && b == 1 && c == 1 && d == 1)
    {
        s[i][0] = 3;
        s[i][1] = -3;
    }
    else
    {
        s[i][0] = 0;
        s[i][1] = -0;
    }
}
void x_estimate_m(int i, double a, double b){
    if (a >= 2) m_estimate[i][0] = 3;
    else if (a < 2 && a >= 0) m_estimate[i][0] = 1;
    else if (a < 0 && a >= -2) m_estimate[i][0] = -1;
    else if (a < -2) m_estimate[i][0] = -3;
    else m_estimate[i][0] = 0;
    if (b >= 2) m_estimate[i][1] = 3;
    else if (b < 2 && b >= 0) m_estimate[i][1] = 1;
    else if (b < 0 && b >= -2) m_estimate[i][1] = -1;
    else if (b < -2) m_estimate[i][1] = -3;
    else m_estimate[i][1] = 0;
}
void table_receive(int i, double a, double b) {
    if (a == -3 && b == 3) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 1;
    } else if (a == -1 && b == 3) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0;
    } else if (a == 1 && b == 3) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 0;
    } else if (a == 3 && b == 3) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 1;
    } else if (a == -3 && b == 1) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 1;
    } else if (a == -1 && b == 1) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0;
    } else if (a == 1 && b == 1) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 0;
    } else if (a == 3 && b == 1) {
        receive_sym[i][0] = 1;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 1;
    } else if (a == -3 && b == -1) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 1;
    } else if (a == -1 && b == -1) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0;
    } else if (a == 1 && b == -1) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 0;
    } else if (a == 3 && b == -1) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 1;
    } else if (a == -3 && b == -3) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 1;
    } else if (a == -1 && b == -3) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0;
    } else if (a == 1 && b == -3) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 0;
    } else if (a == 3 && b == -3) {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 1;
        receive_sym[i][2] = 1;
        receive_sym[i][3] = 1;
    } else {
        receive_sym[i][0] = 0;
        receive_sym[i][1] = 0;
        receive_sym[i][2] = 0;
        receive_sym[i][3] = 0; 
    }
}
