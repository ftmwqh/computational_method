#include<stdio.h>
#include<math.h>
#include<stdlib.h>
double two(int row,int col,double A[row][col])//先定义行列数，是一种二维参数方法
{
    int i,j;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            printf("%f\t",A[i][j]);
        }
        printf("\n");
    }
}
double max(int row,int col,double A[row][col])
{
    int i,j;
    double max;
    max=0;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            if(fabsl(A[i][j])>max)
            {
                max=A[i][j];
            }
        }
    }
    return max;
}
double **TD(int row,int col,double **A)//指针
{
    int i,j;
    double **L,**U;
    L=(double **)malloc(row*sizeof(double *));
    for(i=0;i<row;i++)
    {
        L[i]=(double *)malloc(col*sizeof(double));
    }
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            L[i][j]=A[i][j];
        }

    }
    return L;
}
double **LU(int row,int col,double **A,char choice)
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
        if(k==0)
        {
            for(j=k;j<row;j++)
            {
                U[k][j]=A[k][j];
            }
            for(i=k+1;i<col;i++)
            {
                L[i][k]=A[i][k]/U[k][k];
            }
        }
        else
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
    }


    if(choice=='L')
    {
        return L;
    }
    else if(choice=='U')
    {
        return U;
    }

}
int main()
{
    double A_1[5][5];
    double A_2[4][4]={{4.0,-1.0,1.0,3.0},{16.0,-2.0,-2.0,5.0},{16.0,-3.0,-1.0,7.0},{6.0,-4.0,2.0,9.0}};
    double *a_2[4],*a_1[5],**a;
    int i,j;
    double epsilon;
    double b[1][3]={1,2,3};
    for(i=0;i<5;i++)
    {
        for(j=0;j<5;j++)
        {
            A_1[i][j]=1.0/(9-i-j);
        }
    }
    // two(5,5,A_1);
    // two(4,4,A_2);
    // two(1,3,b);
    for(i=0;i<4;i++)
    {
        a_2[i]=&A_2[i][0];
    }
    for(i=0;i<5;i++)
    {
        a_1[i]=&A_1[i][0];
    }   
    

    a=LU(4,4,a_2,'U');
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            printf("%f\t",a[i][j]);
        }
        printf("\n");
    }


    // printf("%f",max(5,5,A_2));
    system("pause");
    return 0;
}