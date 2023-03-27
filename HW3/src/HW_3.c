#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double max(int row,double x[row])//求向量中的最大值，用于归一化
{
    int i;
    double max;
    max=0;
    for(i=0;i<row;i++)
    {
        if(fabsl(x[i])>max)
        {
            max=fabsl(x[i]);
        }
        
    }
    return max;
}

double **LU(int row,int col,double **A,char choice)//进行LU分解，结果返回L，U矩阵
{
    int i,j,k,r;
    double **L,**U;
    L=(double **)malloc(row*sizeof(double *)); //L,U矩阵分配内存
    for(i=0;i<row;i++)
    {
        L[i]=(double *)malloc(col*sizeof(double));
    }
    U=(double **)malloc(row*sizeof(double *));
    for(i=0;i<row;i++)
    {
        U[i]=(double *)malloc(col*sizeof(double));
    }
    for(i=0;i<row;i++)//初始化L对角元素为1
    {
        L[i][i]=1.0;
    }
    //进行分解
    for(k=0;k<row;k++)
    {
        for(j=k;j<row;j++)
        {
            U[k][j]=A[k][j];
            for(r=0;r<k;r++)
            {
                U[k][j]-=L[k][r]*U[r][j];
            }
        }
        for(i=k+1;i<col;i++)
        {
            L[i][k]=A[i][k];
            for(r=0;r<k;r++)
            {
                L[i][k]-=L[i][r]*U[r][k];
            }
            L[i][k]=L[i][k]/U[k][k];
        }
    }

    //根据输入参数判断返回L矩阵还是U矩阵
    if(choice=='L')
    {
        return L;
    }
    else if(choice=='U')
    {
        return U;
    }

}

double iteration(int row,int col,double **A,double y_0[],double epsilon)
{
    double *a[row],**L,**U;
    double x[row],y[row],z[row];
    double former_EV=0.0,EV=1.0,M;
    int i,j,k;
    for(i=0;i<row;i++)
        {
            a[i]=&A[i][0];
            y[i]=y_0[i];
        }   
    L=LU(row,col,a,'L');
    U=LU(row,col,a,'U');

    printf("Y");//制表
    for(i=0;i<row;i++)
    {
        printf("\t\t");
    }
    printf("X\n");
    //解LU方程
    while(fabsl(EV-former_EV)>epsilon)
    {
        former_EV=EV;
        for(i=0;i<row;i++)
        {
            z[i]=y[i];
            for(j=0;j<i;j++)
            {
                z[i]-=L[i][j]*z[j];
            }
        }
        for(i=row-1;i>=0;i--)
        {
            x[i]=z[i];
            for(j=i+1;j<row;j++)
            {
                x[i]-=U[i][j]*x[j];
            }
            x[i]=x[i]/U[i][i];
        }
        M=max(row,x);

        for(i=0;i<row;i++)//打印出表格
        {
            printf("%10.7f\t",y[i]);
        }
        for(i=0;i<row;i++)
        {
            printf("%14.7f\t",x[i]);
        }
        printf("\n");
        for(i=0;i<row;i++)
        {
            y[i]=x[i]/M;
        }
        EV=1.0/M;
        

    }
    printf("\nEigen Value:%10.9f\n",EV);
    printf("Eigen Vector:(");
    for(i=0;i<row;i++)
    {
        printf("%f",y[i]);
        if(i!=row-1) printf(",");
        else printf(")\n\n");
    }
}
int main()
{
    double A_1[5][5],y_01[5]={1,1,1,1,1},y_02[4]={1,1,1,1};
    double A_2[4][4]={{4.0,-1.0,1.0,3.0},{16.0,-2.0,-2.0,5.0},{16.0,-3.0,-1.0,7.0},{6.0,-4.0,2.0,9.0}};
    double *a_2[4],*a_1[5];
    int i,j;
    double epsilon=pow(10.0,-5);
    double b[1][3]={1,2,3};
    //初始化矩阵A_1,A_2
    for(i=0;i<5;i++)
    {
        for(j=0;j<5;j++)
        {
            A_1[i][j]=1.0/(9-i-j);
        }
    }
    //把A_1,A_2每行开头地址赋给a_1,2，以作为函数的参数输入
    for(i=0;i<4;i++)
    {
        a_2[i]=&A_2[i][0];
    }
    for(i=0;i<5;i++)
    {
        a_1[i]=&A_1[i][0];
    }   
    printf("A_1:\n");
    iteration(5,5,a_1,y_01,epsilon);
    printf("A_2:\n");
    iteration(4,4,a_2,y_02,epsilon);
    
    system("pause");
    return 0;
}