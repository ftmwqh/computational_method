#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define N 21

double *chase(int n,double a[n],double b[n-1],double c[n-1],double f[n])//追赶法求解三对角
{//其中a为对角，b为上一行，c为下一行，f为方程右边向量
    double *x,u[n],v[n-1],y[n];
    int i,j,k;
    x=(double *)malloc(n*sizeof(double)); 
    for(k=0;k<n;k++)//追
    {
        if(k==0)
        {
            u[k]=a[k];
            v[k]=b[k]/u[k];
            y[k]=f[k]/u[k];
        }
        else if (k==n-1)
        {
            u[k]=a[k]-c[k-1]*v[k-1];
            y[k]=(f[k]-c[k-1]*y[k-1])/u[k];
        }
        
        else
        {
            u[k]=a[k]-c[k-1]*v[k-1];
            v[k]=b[k]/u[k];
            y[k]=(f[k]-c[k-1]*y[k-1])/u[k];
        }     
    }
    for(k=n-1;k>=0;k--)//赶
    {
        if(k==n-1)
        {
            x[k]=y[k];
        }
        else
        {
            x[k]=y[k]-v[k]*x[k+1];
        }
        
    }
    return x;
}

int main()
{
    double *M;
    float y[N],x[N];
    double h[20],lambda[18],mu[18],A[19],d[19];
    double a[20],b[20],c[20],e[20];//S(x)的四个系数
    int i,j,k;

    FILE* fp = fopen("E:\\C_Projects\\computation_methods\\HW5\\point.txt", "r");
    if(fp==NULL)
    {
        printf("无文件");
        return -1;
    }
      
    for (i = 0; i < 21; i++)
    {//不能直接把int,float读入到double数组！！可以先到float，再到double
        fscanf(fp,"%f\t",&x[i]);
        fscanf(fp,"%f\n",&y[i]);
        if(i!=0)
        {
            h[i-1]=x[i]-x[i-1];
        }
    }
    fclose(fp);

    for(i=0;i<18;i++)
    {
        lambda[i]=h[i+1]/(h[i]+h[i+1]);
        mu[i]=1-h[i+2]/(h[i+1]+h[i+2]);
    }
    for(i=0;i<19;i++)
    {
        A[i]=2.0;
        d[i]=6/(h[i]+h[i+1])*((y[i+2]-y[i+1])/h[i+1]-(y[i+1]-y[i])/h[i]);
    }
    M=chase(19,A,lambda,mu,d);
    printf("追赶法求解得到M\n");
    for(i=0;i<19;i++)
    {
        printf("%f\t",M[i]);
    }
    printf("\n");
    //解方程得到各段系数
    printf("求解得到多项式如下\n");

    if((fp=fopen("E:\\C_Projects\\computation_methods\\HW5\\coefficient.txt","w"))==NULL)
    {
        printf("can't open");
        exit(0);
    }

    for(i=0;i<20;i++)
    {
        if(i==0)
        {
            a[i]=(M[i])/(6*(x[i+1]-x[i]));
            b[i]=-3*x[i]*a[i];
            c[i]=(y[i+1]-y[i]-(powf(x[i+1],3.0)-powf(x[i],3.0))*a[i]-(powf(x[i+1],2.0)-powf(x[i],2.0))*b[i])/(x[i+1]-x[i]);  
            e[i]=y[i]-x[i]*c[i]-powf(x[i],2.0)*b[i]-powf(x[i],3.0)*a[i];
        }
        if(i==19)
        {
            a[i]=(-M[i-1])/(6*(x[i+1]-x[i]));
            b[i]=0.5*M[i-1]-3*x[i]*a[i];
            c[i]=(y[i+1]-y[i]-(powf(x[i+1],3.0)-powf(x[i],3.0))*a[i]-(powf(x[i+1],2.0)-powf(x[i],2.0))*b[i])/(x[i+1]-x[i]);  
            e[i]=y[i]-x[i]*c[i]-powf(x[i],2.0)*b[i]-powf(x[i],3.0)*a[i];
        }
        else
        {
            a[i]=(M[i]-M[i-1])/(6*(x[i+1]-x[i]));
            b[i]=0.5*M[i-1]-3*x[i]*a[i];
            c[i]=(y[i+1]-y[i]-(powf(x[i+1],3.0)-powf(x[i],3.0))*a[i]-(powf(x[i+1],2.0)-powf(x[i],2.0))*b[i])/(x[i+1]-x[i]);  
            e[i]=y[i]-x[i]*c[i]-powf(x[i],2.0)*b[i]-powf(x[i],3.0)*a[i];
        }
        printf("S(x)=(%f)x^3+(%f)x^2+(%f)x+(%f),%f<x<%f\n",a[i],b[i],c[i],e[i],x[i],x[i+1]);
        fprintf(fp,"%f\t%f\t%f\t%f\n",a[i],b[i],c[i],e[i]);
    } 
    fclose(fp);

    //改变第十个压铁坐标
    y[9]=10.0;
    for(i=0;i<19;i++)
    {
        A[i]=2.0;
        d[i]=6/(h[i]+h[i+1])*((y[i+2]-y[i+1])/h[i+1]-(y[i+1]-y[i])/h[i]);
    }
    M=chase(19,A,lambda,mu,d);
    printf("追赶法求解得到M\n");
    for(i=0;i<19;i++)
    {
        printf("%f\t",M[i]);
    }
    printf("\n");
    //解方程得到各段系数
    printf("改变第10个压铁点坐标后的多项式如下\n");
    if((fp=fopen("E:\\C_Projects\\computation_methods\\HW5\\coefficient_10.txt","w"))==NULL)
    {
        printf("can't open");
        exit(0);
    }
    for(i=0;i<20;i++)
    {
        if(i==0)
        {
            a[i]=(M[i])/(6*(x[i+1]-x[i]));
            b[i]=-3*x[i]*a[i];
            c[i]=(y[i+1]-y[i]-(powf(x[i+1],3.0)-powf(x[i],3.0))*a[i]-(powf(x[i+1],2.0)-powf(x[i],2.0))*b[i])/(x[i+1]-x[i]);  
            e[i]=y[i]-x[i]*c[i]-powf(x[i],2.0)*b[i]-powf(x[i],3.0)*a[i];
        }
        if(i==19)
        {
            a[i]=(-M[i-1])/(6*(x[i+1]-x[i]));
            b[i]=0.5*M[i-1]-3*x[i]*a[i];
            c[i]=(y[i+1]-y[i]-(powf(x[i+1],3.0)-powf(x[i],3.0))*a[i]-(powf(x[i+1],2.0)-powf(x[i],2.0))*b[i])/(x[i+1]-x[i]);  
            e[i]=y[i]-x[i]*c[i]-powf(x[i],2.0)*b[i]-powf(x[i],3.0)*a[i];
        }
        else
        {
            a[i]=(M[i]-M[i-1])/(6*(x[i+1]-x[i]));
            b[i]=0.5*M[i-1]-3*x[i]*a[i];
            c[i]=(y[i+1]-y[i]-(powf(x[i+1],3.0)-powf(x[i],3.0))*a[i]-(powf(x[i+1],2.0)-powf(x[i],2.0))*b[i])/(x[i+1]-x[i]);  
            e[i]=y[i]-x[i]*c[i]-powf(x[i],2.0)*b[i]-powf(x[i],3.0)*a[i];
        }
        printf("S(x)=(%f)x^3+(%f)x^2+(%f)x+(%f),%f<x<%f\n",a[i],b[i],c[i],e[i],x[i],x[i+1]);
        fprintf(fp,"%f\t%f\t%f\t%f\n",a[i],b[i],c[i],e[i]);
    } 
    fclose(fp);

    system("pause");
    return 0;
}