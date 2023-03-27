#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>

#define pi 3.14159265359
#define n 16//n=2^4
#define N 128//n=2^7

double cmod(complex double a)//取复数的模
{
    return sqrtl(powf(creal(a),2)+powf(cimag(a),2));
}

complex double f1(double t)
{
    return 0.7*sin(2*pi*2*t)+sin(2*pi*5*t);
}

complex double f2(double t)
{
    //srand((unsigned int)time(NULL));
    return 0.7*sin(2*pi*2*t)+sin(2*pi*5*t)+0.3*rand()/(double)(RAND_MAX);
}

complex double *FFT(int n_1,complex double f[n_1])//对函数FFT
{
    int i,j,k;
    complex double *g_0,*g_1;
    complex double w_n,w;
    if(n_1==1)
    {
        return f;
    }

    complex double f_0[n_1/2],f_1[n_1/2],*g;
    g=(complex double *)malloc(n_1*sizeof(complex double)); 
    w=1;
    w_n=cexpf(-_Complex_I*2*pi/n_1);
    for(i=0;i<n_1/2;i++)
    {
        f_0[i]=f[2*i];
        f_1[i]=f[2*i+1];
    }
    g_0=FFT(n_1/2,f_0);
    g_1=FFT(n_1/2,f_1);
    for(k=0;k<n_1/2;k++)
    {
        g[k]=(g_0[k]+w*g_1[k])/2;
        g[k+n_1/2]=(g_0[k]-w*g_1[k])/2;
        w=w*w_n;
    }
    return g;
}

complex double *IFFT(int n_1,complex double f[n_1])//对函数IFFT
{
    int i,j,k;
    complex double *g_0,*g_1;
    complex double w_n,w;
    if(n_1==1)
    {
        return f;
    }

    complex double f_0[n_1/2],f_1[n_1/2],*g;
    g=(complex double *)malloc(n_1*sizeof(complex double)); 
    w=1;
    w_n=cexpf(_Complex_I*2*pi/n_1);
    for(i=0;i<n_1/2;i++)
    {
        f_0[i]=f[2*i];
        f_1[i]=f[2*i+1];
    }
    g_0=IFFT(n_1/2,f_0);
    g_1=IFFT(n_1/2,f_1);
    for(k=0;k<n_1/2;k++)//反变换去掉'/2'
    {
        g[k]=(g_0[k]+w*g_1[k]);
        g[k+n_1/2]=(g_0[k]-w*g_1[k]);
        w=w*w_n;
    }
    return g;
}

int main()
{
    int i,j,k;
    double w;
    complex double a,b;
    complex double f_1_4[n],f_1_7[N],f_2[N],g_2_025[N];//变换前离散函数，g2的前25%
    complex double *g_1_4,*g_1_7,*g_2,*if_1_4,*if_1_7,*if_2,*if_2_025;//变换后与逆变换
    FILE *fp;

    // a=3.0+3.0*_Complex_I;
    // b=cexpf(pi/2*_Complex_I);
    // printf("b=%f+%fi\n",creal(b),cimag(b));
    // printf("%f\n",carg(b));
    // printf("%d\n",n);
    
    srand(n);
    // printf("%f\n",f2(0));
    // printf("%f\n",f2(0));
    // printf("%f\n",f2(0));

    for(k=0;k<n;k++)//初始化离散函数f1,f2
    {
        f_1_4[k]=f1((double)k/n);
    }
    for(k=0;k<N;k++)
    {
        f_1_7[k]=f1((double)k/N);
        f_2[k]=f2((double)k/N);
    }

    //f1 n=2^4变换
    g_1_4=FFT(n,f_1_4);
    if_1_4=IFFT(n,g_1_4);
    if((fp=fopen("E:\\C_Projects\\computation_methods\\HW6\\f_1_4.txt","w"))==NULL)
    {
        printf("can't open");
        exit(0);
    }
    printf("f1 n=2^4变换\nf1\t\t\tg1\t\t\t反变换f1\n");
    for(k=0;k<n;k++)//初始化离散函数f1,f2
    {
        printf("%f+%fi\t",creal(f_1_4[k]),cimag(f_1_4[k]));
        printf("%f+%fi\t",creal(g_1_4[k]),cimag(g_1_4[k]));
        printf("%f+%fi\n",creal(if_1_4[k]),cimag(if_1_4[k]));
        fprintf(fp,"%f\t%f\t%f\n",creal(f_1_4[k]),cmod(g_1_4[k]),creal(if_1_4[k]));
    }
    fclose(fp);
    printf("\n");

    //f1 n=2^7变换
    
    g_1_7=FFT(N,f_1_7);
    if_1_7=IFFT(N,g_1_7);
    if((fp=fopen("E:\\C_Projects\\computation_methods\\HW6\\f_1_7.txt","w"))==NULL)
    {
        printf("can't open");
        exit(0);
    }
    printf("f1 n=2^7变换\nf1\t\t\tg1\t\t\t反变换f1\n");
    for(k=0;k<N;k++)//初始化离散函数f1,f2
    {
        printf("%f+%fi\t",creal(f_1_7[k]),cimag(f_1_7[k]));
        printf("%f+%fi\t",creal(g_1_7[k]),cimag(g_1_7[k]));
        printf("%f+%fi\n",creal(if_1_7[k]),cimag(if_1_7[k]));
        fprintf(fp,"%f\t%f\t%f\n",creal(f_1_7[k]),cmod(g_1_7[k]),creal(if_1_7[k]));
    }
    fclose(fp);
    printf("\n"); 

    //f2 n=2^7变换
    g_2=FFT(N,f_2);
    if_2=IFFT(N,g_2);
    if((fp=fopen("E:\\C_Projects\\computation_methods\\HW6\\f_2.txt","w"))==NULL)
    {
        printf("can't open");
        exit(0);
    }
    printf("f2 n=2^7变换\nf2\t\t\tg2\t\t\t反变换f2\n");
    for(k=0;k<N;k++)//初始化离散函数f1,f2
    {
        printf("%f+%fi\t",creal(f_2[k]),cimag(f_2[k]));
        printf("%f+%fi\t",creal(g_2[k]),cimag(g_2[k]));
        printf("%f+%fi\n",creal(if_2[k]),cimag(if_2[k]));
        fprintf(fp,"%f\t%f\t%f\n",creal(f_2[k]),cmod(g_2[k]),creal(if_2[k]));
    }
    fclose(fp);
    printf("\n");

    for(k=0;k<N;k++)//将g2的前25%赋给g_2_025
    {
        if(k<N/8||k>7*N/8)
        {
            g_2_025[k]=g_2[k];
        }
        else
        {
            g_2_025[k]=0;
        }
    }
    if_2_025=IFFT(N,g_2_025);
    if((fp=fopen("E:\\C_Projects\\computation_methods\\HW6\\if_2_025.txt","w"))==NULL)
    {
        printf("can't open");
        exit(0);
    }
    for(k=0;k<N;k++)//g2的前25%逆变换
    {
        fprintf(fp,"%f\n",creal(if_2_025[k]));
    }
    fclose(fp);

    system("pause");
    return 0;
}