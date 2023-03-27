#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>

#define pi 3.14159265359
#define e powf(10,-6)//控制精度
int reachx=0,reachy=0,Mx=0,My=0;//定义全局变量，用于记录每次计算是否达到精度,还有计算调用次数
long long n=0,N=0;

double ax(double t)//x方向加速度
{
    return sinl(t)/(sqrtl(t)+1);
}

double ay(double t)//y方向加速度
{
    return logl(t+1)/(t+1);
}

double vx(double t,int M)//计算给出Vx在t时刻的值
{
    int k,j,i,count;
    double R[2][M],h,sum;
    h=t;
    R[0][0]=(ax(0.0)+ax(t))*h/2;
    count=0;//使用count进行记录上次使用的R[0]还是R[1]
    for(k=2;k<=M;k++)
    {
        sum=0;
        h=h/2;//将h赋予下一次的区间长度
        for(i=1;i<=pow(2,k-2);i++)
        {
            sum+=ax(0.0+(2*i-1)*h);
        }
        R[(count+1)%2][0]=(R[count%2][0]+2*h*sum)/2.0;
 
        for(j=1;j<k;j++)
        {
            R[(count+1)%2][j]=R[(count+1)%2][j-1]+(R[(count+1)%2][j-1]-R[(count)%2][j-1])/(pow(4,j)-1);
        }
        if(fabsl(R[(count+1)%2][k-1]-R[count%2][k-2])<e) 
        {
            reachx+=1;
            n+=1;
            break;
        }
        count++;//计次，使下一次R使用的两列调换
    }
    //printf("%10.8f\n",fabsl(R[(count+1)%2][k-1]-R[count%2][k-2]));
    Mx+=1;
    N+=1;
    return R[(count+1)%2][k-1];
}

double vy(double t,int M)//计算给出Vy在t时刻的值
{
    int k,j,i,count;
    double R[2][M],h,sum;
    h=t;
    R[0][0]=(ay(0.0)+ay(t))*h/2;
    count=0;//使用count进行记录上次使用的R[0]还是R[1]
    for(k=2;k<=M;k++)
    {
        sum=0;
        h=h/2;//将h赋予下一次的区间长度
        for(i=1;i<=pow(2,k-2);i++)
        {
            sum+=ay(0.0+(2*i-1)*h);
        }
        R[(count+1)%2][0]=(R[count%2][0]+2*h*sum)/2.0;
 
        for(j=1;j<k;j++)
        {
            R[(count+1)%2][j]=R[(count+1)%2][j-1]+(R[(count+1)%2][j-1]-R[(count)%2][j-1])/(pow(4,j)-1);
        }
        if(fabsl(R[(count+1)%2][k-1]-R[count%2][k-2])<e) 
        {
            reachy+=1;
            n+=1;
            break;
        }
        count++;//计次，使下一次R使用的两列调换
    }
    My+=1;
    N+=1;
    return R[(count+1)%2][k-1];
}

double x(double t,int M)//计算x在t时刻的值
{
    int k,j,i,count;
    double R[2][M],h,sum;
    h=t;
    R[0][0]=(vx(0.0,M)+vx(t,M))*h/2;
    count=0;//使用count进行记录上次使用的R[0]还是R[1]
    for(k=2;k<=M;k++)
    {
        sum=0;
        h=h/2;//将h赋予下一次的区间长度
        for(i=1;i<=pow(2,k-2);i++)
        {
            sum+=vx(0.0+(2*i-1)*h,M);
        }
        R[(count+1)%2][0]=(R[count%2][0]+2*h*sum)/2.0;
 
        for(j=1;j<k;j++)
        {
            R[(count+1)%2][j]=R[(count+1)%2][j-1]+(R[(count+1)%2][j-1]-R[(count)%2][j-1])/(pow(4,j)-1);
        }
        if(fabsl(R[(count+1)%2][k-1]-R[count%2][k-2])<e) 
        {
            reachx+=1;//达到精度退出时，+1
            n+=1;
            break;
        }
    
        count++;//计次，使下一次R使用的两列调换
    }
    Mx+=1;
    N+=1;
    return R[(count+1)%2][k-1];
}

double y(double t,int M)//计算给出Vy在t时刻的值
{
    int k,j,i,count;
    double R[2][M],h,sum;
    h=t;
    R[0][0]=(vy(0.0,M)+vy(t,M))*h/2;
    count=0;//使用count进行记录上次使用的R[0]还是R[1]
    for(k=2;k<=M;k++)
    {
        sum=0;
        h=h/2;//将h赋予下一次的区间长度
        for(i=1;i<=pow(2,k-2);i++)
        {
            sum+=vy(0.0+(2*i-1)*h,M);
        }
        R[(count+1)%2][0]=(R[count%2][0]+2*h*sum)/2.0;
 
        for(j=1;j<k;j++)
        {
            R[(count+1)%2][j]=R[(count+1)%2][j-1]+(R[(count+1)%2][j-1]-R[(count)%2][j-1])/(pow(4,j)-1);
        }
        if(fabsl(R[(count+1)%2][k-1]-R[count%2][k-2])<e) 
        {
            reachy+=1;//达到精度退出时，+1
            n+=1;
            break;
        }
        count++;//计次，使下一次R使用的两列调换
    }
    My+=1;
    N+=1;
    return R[(count+1)%2][k-1];
}

int main()
{
    int i,j,k,M;
    double t,ratex[5]={0,0,0,0,0},ratey[5]={0,0,0,0,0};
    FILE *fp;
    M=8;
    if((fp=fopen("E:\\C_Projects\\computation_methods\\HW7\\xy.txt","w"))==NULL)
    {
        printf("can't open");
        exit(0);
    }  
    for(i=1;i<=100;i++)
    {
        fprintf(fp,"%f\t%f\n",x(0.1*i,M),y(0.1*i,M));
    }
    fclose(fp);
    reachx=0;
    reachy=0;
    Mx=0;
    My=0;
    n=0;
    N=0;

    printf("M\tx方向达到精度比例\ty方向达到精度比例\t总积分达到精度比例\n");
    for(i=1;i<=5;i++)
    {
        M=i*4;
        for(j=1;j<=100;j++)
        {
            x(0.1*j,M);
            y(0.1*j,M);
            
        }
        ratex[i-1]=reachx;
        ratey[i-1]=reachy;
        ratex[i-1]/=(double)Mx;//总共有10000次积分
        ratey[i-1]/=(double)My;
        printf("%d\t%f\t\t%f\t\t%f\n",M,ratex[i-1],ratey[i-1],(double)(n)/(double)(N));
        reachx=0;
        reachy=0;
        Mx=0;
        My=0;
        n=0;
        N=0;
    }


    printf("finish\n");

    system("pause");
    return 0;
}