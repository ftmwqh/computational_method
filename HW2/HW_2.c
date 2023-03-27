#include<stdio.h>
#include<math.h>
#include<stdlib.h>
double standard(double x,double a,double epsilon)//标准精确解
{
    return (1-a)/(1-exp(-1/epsilon))*(1-exp(-x/epsilon))+a*x;
}
double module_2(double a[],int n)//求解一个n维向量的p=2的模
{
    int i;
    double sum=0;
    for(i=0;i<n-1;i++)
    {   
        sum+=pow(a[i],2);

    }
    return sqrtl(sum);
}
double max(double a,double b)
{
	double z;
	z=a>b?a:b;
	return(z);
}
double mod_Ax_b(double y[],double a,int n,double epsilon)//计算Ax-b的模长
{
    int i,j;
    double h=1.0/(double)(n);
    double vector[n-1];//用于记录向量Ax-b
    for(i=0;i<n-1;i++)
    {
        if(i==0)
        {
            vector[i]=-(2*epsilon+h)*y[i]+(epsilon+h)*y[i+1]-a*pow(h,2);
        }
        else if (i==n-2)
        {
            vector[i]=-(2*epsilon+h)*y[i]+epsilon*y[i-1]+(epsilon+h)*1-a*pow(h,2);
        }
        else
        {
            vector[i]=-(2*epsilon+h)*y[i]+(epsilon+h)*y[i+1]+epsilon*y[i-1]-a*pow(h,2);
        }
        
    }
    return module_2(vector,n);


}
double mod_x_x(double x1[],double x2[],int n)//计算x(n)-x(n-1)的模
{
    int i;
    double x_x[n-1];
    for(i=0;i<n-1;i++)
    {
        x_x[i]=x1[i]-x2[i];
    }
    return module_2(x_x,n);
}
double Gauss_Seidel(double a,int n,double epsilon,double delta,char* filename)//G-S迭代法
{
    double h,y[n-1],y_former[n-1],now,var,inf_mod;//注意：此处由于不取端点，只有n-1维
    int i,j,count;
    count=0;
    h=1.0/(double)(n);
    for(i=0;i<n-1;i++)
    {
        y[i]=1;
    }
    while(1)
    {
        count++;
        for(j=0;j<n-1;j++)
        {
            y_former[j]=y[j];
        }
        for(i=0;i<n-1;i++)
        {
            if(i==0)
            {
                y[i]=((epsilon+h)*y[i+1]-a*pow(h,2))/(2*epsilon+h);
            }
            else if (i==n-2)
            {
                y[i]=(epsilon*y[i-1]-a*pow(h,2)+(epsilon+h)*1)/(2*epsilon+h);
            }
            else
            {
                y[i]=((epsilon+h)*y[i+1]+epsilon*y[i-1]-a*pow(h,2))/(2*epsilon+h);
            }        
        }
        if(max(mod_Ax_b(y,a,n,epsilon),mod_x_x(y,y_former,n))<delta)break;   
    }
    
    printf("times of iteration:%d\n",count);
    // for(i=0;i<n-1;i++)
    // {
    //     printf("mine:%.15f\t",y[i]);
    //     printf("standard:%.15f\n",standard((double)(i+1)/n,a,epsilon));
    // }
    var=0;
    for(i=0;i<n-1;i++)//计算方差
    {
        var+=pow(standard((double)(i+1)/n,a,epsilon)-y[i],2);
    }
    var=var/(n-1);
    printf("Variance:%.10f\t",var);
    
    now=fabsl(y[0]-standard((double)(1)/n,a,epsilon));//求取解与标准方程差的无穷范数
    for(i=0;i<n-1;i++)
    {
        if(fabsl(y[i]-standard((double)(i+1)/n,a,epsilon))>fabsl(now))
        {
            now=fabsl(y[i]-standard((double)(i+1)/n,a,epsilon));
        }
    }
    printf("Infinity_Module:%.10f\n",now);

    FILE*fp;
    if((fp=fopen(filename,"w"))==NULL)
    {
        printf("cant open");
        exit(0);
    }
    for(i=0;i<n-1;i++)
    {
        fprintf(fp,"%f\t",y[i]);
        fprintf(fp,"%f\n",standard((double)(i+1)/n,a,epsilon));
    }
    fclose(fp);

}

int Gause_Elimination(double a,int n,double epsilon,char* filename)//高斯列主元消元法
{
    double A[n-1][n-1],b[n-1],h,now,para,var,inf_mod;
    int i,j,k,max_row;
    h=1.0/(double)n;
    for(i=0;i<n-1;i++)
    {
        for(j=0;j<n-1;j++)
        {
            A[i][j]=0;
        }
    }
    for(i=0;i<n-1;i++)//初始化A,b矩阵
    {
        if(i==0)
        {
            A[i][i]=-(2*epsilon+h);
            A[i][i+1]=epsilon+h;
            b[i]=a*pow(h,2);
        }
        else if(i==n-2)
        {
            A[i][i]=-(2*epsilon+h);
            A[i][i-1]=epsilon;
            b[i]=a*pow(h,2)-(epsilon+h);
        }
        else
        {
            A[i][i+1]=epsilon+h;
            A[i][i]=-(2*epsilon+h);
            A[i][i-1]=epsilon;   
            b[i]=a*pow(h,2);         
        }
    }

    for(k=0;k<n-2;k++)
    {
        
        now=A[k][k];//求取最大列主元的行数max_row
        max_row=k;
        for(i=k;i<n-1;i++)
        {
            if(fabsl(A[i][k])>fabsl(now))
            {
                now=A[i][k];
                max_row=i;
            }
        }
        
        for(j=k;j<n-1;j++)//交换行得到最大列主元A,b矩阵
        {
            now=A[max_row][j];
            A[max_row][j]=A[k][j];
            A[k][j]=now;
            now=b[max_row];
            b[max_row]=b[k];
            b[k]=now;
        }
 
        for(i=k+1;i<n-1;i++)//消元过程
        {
            para=A[i][k]/A[k][k];
            A[i][k]=0;
            for(j=k+1;j<n-1;j++)
            {
                A[i][j]=A[i][j]-para*A[k][j];
            }
            b[i]=b[i]-para*b[k];
        }

    }

    b[n-1]=b[n-1]/A[n-1][n-1];//回代解放方程
    for(i=n-2;i>=0;i--)
    {
        for(j=i+1;j<n-1;j++)
        {
            b[i]=b[i]-A[i][j]*b[j];
        }
        b[i]=b[i]/A[i][i];
    }
    var=0;
    for(i=0;i<n-1;i++)//计算方差
    {
        var+=pow(standard((double)(i+1)/n,a,epsilon)-b[i],2);
    }
    var=var/(n-1);
    printf("Variance:%.10f\t",var);
    
    now=fabsl(b[0]-standard((double)(1)/n,a,epsilon));//求取解与标准方程差的无穷范数
    for(i=0;i<n-1;i++)
    {
        if(fabsl(b[i]-standard((double)(i+1)/n,a,epsilon))>fabsl(now))
        {
            now=fabsl(b[i]-standard((double)(i+1)/n,a,epsilon));
        }
    }
    printf("Infinity_Module:%.10f\n",now);
    // for(i=0;i<n-1;i++)          //输出对比结果
    // {
    //     printf("mine:%f\t",b[i]);
    //     printf("standard:%f\n",standard((double)(i+1)/n,a,epsilon));
    // }
    FILE*fp;
    if((fp=fopen(filename,"w"))==NULL)
    {
        printf("cant open");
        exit(0);
    }
    for(i=0;i<n-1;i++)
    {
        fprintf(fp,"%f\t",b[i]);
        fprintf(fp,"%f\n",standard((double)(i+1)/n,a,epsilon));
    }
    fclose(fp);


    return 0;
}
int main()
{
    double a,epsilon,n;
    double x[2]={1,2},y[2]={0,0001};
    int i=3;
    n=100;
    a=0.5;
    
    epsilon=1;
    printf("epsilon=%f\n",epsilon);
    printf("Gause Seidel:\n");
    Gauss_Seidel(a,n,epsilon,0.000000001,"G_S_1.txt");
    printf("Gauss Elimination:\n");
    Gause_Elimination(a,n,epsilon,"G_E_1.txt");
    printf("\n");
    
    epsilon=0.1;
    printf("epsilon=%f\n",epsilon);
    printf("Gause Seidel:\n");
    Gauss_Seidel(a,n,epsilon,0.0000001,"G_S_01.txt");
    printf("Gauss Elimination:\n");
    Gause_Elimination(a,n,epsilon,"G_E_01.txt");
    printf("\n");
    
    epsilon=0.01;
    printf("epsilon=%f\n",epsilon);
    printf("Gause Seidel:\n");
    Gauss_Seidel(a,n,epsilon,0.000001,"G_S_001.txt");
    printf("Gauss Elimination:\n");
    Gause_Elimination(a,n,epsilon,"G_E_001.txt");
    printf("\n");
    
    epsilon=0.0001;
    printf("epsilon=%f\n",epsilon);
    printf("Gause Seidel:\n");
    Gauss_Seidel(a,n,epsilon,0.000001,"G_S_00001.txt");
    printf("Gause Elimination:\n");
    Gause_Elimination(a,n,epsilon,"G_E_00001.txt");
    printf("\n");

    system("pause");
    return 0;
}