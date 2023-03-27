#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

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

double max_index(int row,int col,double **A,char choice)//求最大非对角元的行、列指标
{
    int i,j,p,q;
    double max;
    max=0;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            if(i!=j)
            {
                if(fabsl(A[i][j])>max)
                {
                    max=fabsl(A[i][j]);
                    p=i;
                    q=j;
                }
            }
        }
    }
    if(choice=='p')//返回行指标
    {
        return p;
    }
    if(choice=='q')//返回列指标
    {
        return q;
    }
}

double module_2(double a[],int row)//求解一个n维向量的p=2的模
{
    int i;
    double sum=0;
    for(i=0;i<row-1;i++)
    {   
        sum+=pow(a[i],2);
    }
    return sqrtl(sum);
}

double square_sum(int row,int col,double **A)//求非对角元平方和
{
    double sum;
    int i,j;
    sum=0;
    for(i=0;i<row;i++)
    {
        for(j=0;j<col;j++)
        {
            if(i!=j)
            {
                sum+=A[i][j]*A[i][j];
            }
        }
    }
    return sum;
}

double **Transpose(int row,int col,double **A)//矩阵转置
{
    int i,j;
    double **A_T;
    A_T=(double **)malloc(col*sizeof(double *)); //为转置后的矩阵A_T分配内存
    for(i=0;i<col;i++)
    {
        A_T[i]=(double *)malloc(row*sizeof(double));
    }
    
    for(i=0;i<row;i++)//转置矩阵赋值
    {
        for(j=0;j<col;j++)
        {
            A_T[j][i]=A[i][j];
        }
    }
    return A_T;

}

double **multiply(int row_1,int col_1,double **A_1,int row_2,int col_2,double **A_2)
{
    int i,j,k;
    double **result;
    result=(double **)malloc(row_1*sizeof(double *)); //为乘积分配内存row_1行col_2列
    for(i=0;i<row_1;i++)
    {
        result[i]=(double *)malloc(col_2*sizeof(double));
    }
    for(i=0;i<row_1;i++)
    {
        for(j=0;j<col_2;j++)
        {
            result[i][j]=0;
            for(k=0;k<col_1;k++)
            {
                result[i][j]+=A_1[i][k]*A_2[k][j];
            }
        }
    }
    return result;

}//函数multiply、Transpose会改变所赋给的地址，之后操作需要重新进行赋地址的操作

double det_4(int row,int col,double **A)
{
    int i,j,m,n,k=0;
    double l=0,b[24];
    for(i=0;i<4;i++) //先不考虑（-1）的n+j次方，计算出每一个值，四阶方程式共有24个
    {
        for(j=0;j<4;j++)
        {if(j==i)
        continue;
            for(m=0;m<4;m++)
            {if((m==i)||(m==j))
            continue;
                for(n=0;n<4;n++)
                {if((n==m)||(n==j)||(n==i))
                continue;
                    b[k]=A[0][i]*A[1][j]*A[2][m]*A[3][n];
                    k=k+1;
                }
            }
        }
    }
        j=0;
    for(i=1;i<24;i++)    //考虑（-1）的n+j次方，那24个值得正负规律为+--++--++--++--++--++--+
    {
    if(j<2)
    {
        l=l-b[i];
    }
    if(j>=2)
    {
        l=l+b[i];
    }
    j++;
    if(j==4)
    {
        j=0;
    }
    
    }
    return (l+b[0]);
}

double **Jacobi(int row,int col,double **A,double epsilon,char choice)
{
    double **EV,**Q;//由于使用指针**A,直接操作会改变外部A的值，所以采用EV进行对角化操作
    double *ev[row],*e_vec[row],*q_i[row];
    double sub[row][col],Q_i[row][col];//sub用于记录EV前一个的值,Q_i为当前所乘的矩阵
    double s,t,t1,t2,c,d;
    int i,j,p,q;
    EV=(double **)malloc(row*sizeof(double *)); //为EV分配内存
    for(i=0;i<row;i++)
    {
        EV[i]=(double *)malloc(col*sizeof(double));
    }

    Q=(double **)malloc(row*sizeof(double *)); //Q即为特征向量矩阵，是每一次Q_i的乘积
    for(i=0;i<row;i++)
    {
        Q[i]=(double *)malloc(col*sizeof(double));
    }

    for(i=0;i<row;i++)//令EV=A
    {
        for(j=0;j<col;j++)
        {
            EV[i][j]=A[i][j];
        }
    }

    for(i=0;i<row;i++)//ev指向EV每行地址,e_vec-Q,q_i-Q_i
    {
        ev[i]=&EV[i][0];
        e_vec[i]=&Q[i][0];
        q_i[i]=&Q_i[i][0];
    }

    for(i=0;i<row;i++)//初始化Q,Q_i为单位矩阵。
    {
        for(j=0;j<col;j++)
        {
            if(i==j)
            {
                Q[i][j]=1.0;
                Q_i[i][j]=1.0;
            }
            else
            {
                Q[i][j]=0;
                Q_i[i][j]=0;
            }
        }
    }
    
    while (square_sum(row,col,ev)>epsilon)
    {
        for(i=0;i<row;i++)
        {
            for(j=0;j<col;j++)
            {
                sub[i][j]=EV[i][j];
            }
        }

        p=max_index(row,col,ev,'p');
        q=max_index(row,col,ev,'q');
        s=(EV[q][q]-EV[p][p])/(2*EV[p][q]);
        if(s==0)
        {
            t=1.0;
        }
        else
        {
            t1=-s-sqrtl(pow(s,2)+1.0);
            t2=-s+sqrtl(pow(s,2)+1.0);
            if(fabsl(t1)>fabsl(t2))
            {
                t=t2;
            }
            else
            {
                t=t1;
            }
        }
        c=1.0/sqrtl(1+pow(t,2));
        d=t/sqrtl(1+pow(t,2));

        for(i=0;i<row;i++)
        {
            if((i!=p)&&(i!=q))
            {
                EV[i][p]=c*sub[p][i]-d*sub[q][i];
                EV[p][i]=EV[i][p];
                EV[i][q]=c*sub[q][i]+d*sub[p][i];
                EV[q][i]=EV[i][q];
            }
        }
        EV[p][p]=sub[p][p]-t*sub[p][q];
        EV[q][q]=sub[q][q]+t*sub[p][q];
        EV[p][q]=0.0;
        EV[q][p]=0.0;
        if(choice=='E')
        {
            printf("非对角元平方和:%9.8f\n",square_sum(row,col,ev));
        }
        for(i=0;i<row;i++)//初始化Q_i为单位矩阵。
        {
            for(j=0;j<col;j++)
            {
                if(i==j)
                {
                    Q_i[i][j]=1.0;
                }
                else
                {
                    Q_i[i][j]=0;
                }
            }
        }
        Q_i[p][p]=c;
        Q_i[p][q]=d;
        Q_i[q][p]=-d;
        Q_i[q][q]=c;//令Q_i为当前的变换矩阵
        Q=multiply(row,col,e_vec,col,row,q_i);//Q为把所有Q_i相乘
        for(i=0;i<row;i++)//ev指向EV每行地址,e_vec-Q,q_i-Q_i
        {
            e_vec[i]=&Q[i][0];//函数multiply会改变Q的地址，需要重新将Q地址传给e_vec
        }
    }
    // printf("eigen vector\n");
    // for(i=0;i<row;i++)//初始化Q_i为单位矩阵。
    //     {
    //         for(j=0;j<col;j++)
    //         {
    //             printf("%f\t",Q[i][j]);
    //         }
    //         printf("\n");
    //     }
    if(choice=='E')//E则返回特征值对角阵
    {
        return EV;
    }
    if(choice=='Q')//Q则返回特征向量矩阵
    {
        return Q;//Q为特征向量列组成的矩阵
    }
}

int main()
{
    double A_1[4][3]={{0,1,0},{1,1,0},{1,0,0},{0,0,1}},SIGMA[4][4];
    double *a_1[4],*a_T[4],*aat[4],*a_1t[3],*ata[3],*u[4];
    double **AAT,**ATA,**A_1T,**EV_1,**U,**EV_2,**V;
    int i,j,k;
    double epsilon=pow(10.0,-6),bigger;
    //初始化矩阵A_1,A_2
    printf("初始随机矩阵:\n");
    for(i=0;i<4;i++)
    {
        for(j=0;j<3;j++)
        {
            A_1[i][j]=rand()/(double)(RAND_MAX);//赋予[0,1]的随机数
            printf("%f\t",A_1[i][j]);
        }
        printf("\n");
    }

    //把A_1,A_2每行开头地址赋给a_1,2，以作为函数的参数输入
    for(i=0;i<4;i++)
    {
        a_1[i]=&A_1[i][0];
    }   
    A_1T=Transpose(4,3,a_1);//如果没有地址就运行会导致停止；将A_1转置
    for(i=0;i<3;i++)
    {
        a_1t[i]=&A_1T[i][0];//地址赋给a_1t
    }
    AAT=multiply(4,3,a_1,3,4,a_1t);//A*AT得到AAT矩阵
    printf("A*A^T:\n");
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            printf("%f\t",AAT[i][j]);
        }
        printf("\n");
    }
    for(i=0;i<4;i++)
    {
        aat[i]=&AAT[i][0];//地址赋给aat
    }
    
    EV_1=Jacobi(4,4,aat,epsilon,'E');//AAT的特征值矩阵
    printf("A*A^T的Jacobi特征值矩阵:\n");
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            printf("%f\t",EV_1[i][j]);
        }
        printf("\n");
    }
    U=Jacobi(4,4,aat,epsilon,'Q');//AAT的特征列向量矩阵U
    

    ATA=multiply(3,4,a_1t,4,3,a_1);//A*AT得到AAT矩阵
    printf("A^T*A:\n");
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            printf("%f\t",ATA[i][j]);
        }
        printf("\n");
    }

    for(i=0;i<3;i++)
    {
        ata[i]=&ATA[i][0];//地址赋给ata
    }
    EV_2=Jacobi(3,3,ata,epsilon,'E');//ATA的特征值矩阵
    printf("A^T*A的Jacobi特征值矩阵:\n");
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            printf("%f\t",EV_2[i][j]);
        }
        printf("\n");
    }
    
    

    for(i=0;i<4;i++)
    {
        for (j=0;j < 4; j++)
        {
            if(i==j)
            {
                SIGMA[i][j]=EV_1[i][j];
            }
            else
            {
                SIGMA[i][j]=0.0;
            }
        }
        
    }
    
    for(i=0;i<4;i++)//冒泡法对SIGMA排序，并同时交换U的列顺序
    {
        for(j=0;j<3-i;j++)
        {
            if(SIGMA[j+1][j+1]>SIGMA[j][j])
            {
                bigger=SIGMA[j+1][j+1];
                SIGMA[j+1][j+1]=SIGMA[j][j];
                SIGMA[j][j]=bigger;
                for(k=0;k<4;k++)//根据特征值大小排序后对相应特征向量U排序
                {
                    bigger=U[k][j+1];
                    U[k][j+1]=U[k][j];
                    U[k][j]=bigger;
                }
            }
        }
    }

    for(i=0;i<4;i++)
    {
        u[i]=&U[i][0];
    }
    V=multiply(3,4,a_1t,4,4,u);
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            V[j][i]/=sqrt(SIGMA[i][i]);
        }
    }
   
    // for(i=0;i<3;i++)//冒泡法对SIGMA排序，并同时交换V的列顺序
    // {
    //     for(j=0;j<2-i;j++)
    //     {
    //         if(EV_2[j+1][j+1]>EV_2[j][j])
    //         {
    //             bigger=EV_2[j+1][j+1];
    //             EV_2[j+1][j+1]=EV_2[j][j];
    //             EV_2[j][j]=bigger;
    //             for(k=0;k<3;k++)//根据特征值大小排序后对相应特征向量V排序
    //             {
    //                 bigger=V[k][j+1];
    //                 V[k][j+1]=V[k][j];
    //                 V[k][j]=bigger;
    //             }
    //         }
    //     }
    // }
   
    printf("U:\n");
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            printf("%f\t",U[i][j]);
        }
        printf("\n");
    }

    printf("SIGMA:\n");
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            printf("%f\t",sqrt(SIGMA[i][j]));
        }
        printf("\n");
    }
    printf("%f\t%f\t%f\n",0,0,0);
    printf("V^T:\n");
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            printf("%f\t",V[j][i]);
        }
        printf("\n");
    }
    for(j=0;j<4;j++)
    {
        for(i=0;i<4;i++)
        {
            AAT[i][i]-=EV_1[j][j];
        }
        printf("det(AAT-lambda_i*I)=%14.12f\n",det_4(4,4,aat));//求AAT行列式
        for(i=0;i<4;i++)
        {
            AAT[i][i]+=EV_1[j][j];//矩阵恢复
        }
    }

    //第二问
    printf("以下为第二题iris\n");
    float iris[150][5];//文件中为float,不能读到double中
    double iris_0[50][4],iris_1[50][4],iris_2[50][4],iris_4[150][4];
    double **iris_4T,**Cov,**EV_cov,**Vec_cov;
    double *i0[50],*i1[50],*i2[50],*i4[150],*i4t[4],*cov[4],x[150],y[150];
    double ave=0,vx,vy;
    FILE* fp = fopen("E:\\C_Projects\\computation_methods\\HW4\\iris.txt", "r");
    if(fp==NULL)
    {
        printf("无文件");
        return -1;
    }
      
    for (i = 0; i < 150; i++)
    {
        for (j = 0; j < 4; j++)
        {
            fscanf(fp, "%f,", &iris[i][j]);//注意读取格式","隔开，结尾是"\n"
        }
        fscanf(fp,"%f\n",&iris[i][4]);
    }
     

    for(i=0;i<4;i++)//去中心化后的行向量矩阵
    {
        ave=0;
        for(j=0;j<150;j++)
        {
            ave+=iris[j][i]/150.0;
        }
        for(j=0;j<150;j++)
        {
            iris_4[j][i]=iris[j][i]-ave;
        }
    }
    
    for(i=0;i<150;i++)
    {
        i4[i]=&iris_4[i][0];
    }   
    iris_4T=Transpose(150,4,i4);
    for(i=0;i<4;i++)
    {
        i4t[i]=&iris_4T[i][0];
    }   

    Cov=multiply(4,150,i4t,150,4,i4);
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            Cov[i][j]/=150.0;
        }
    }
    printf("协方差矩阵：\n");
    for(i=0;i<4;i++)
    {
        cov[i]=&Cov[i][0];
        for(j=0;j<4;j++)
        {
            printf("%f\t",Cov[i][j]);
        }
        printf("\n");
    }   
    EV_cov=Jacobi(4,4,cov,epsilon,'E');
    Vec_cov=Jacobi(4,4,cov,epsilon,'Q');
    for(i=0;i<4;i++)//冒泡法对EV_cov，并同时交换Vec_cov的列顺序,得到最大2个特征值对应列向量
    {
        for(j=0;j<3-i;j++)
        {
            if(EV_cov[j+1][j+1]>EV_cov[j][j])
            {
                bigger=EV_cov[j+1][j+1];
                EV_cov[j+1][j+1]=EV_cov[j][j];
                EV_cov[j][j]=bigger;
                for(k=0;k<4;k++)//根据特征值大小排序后对相应特征向量V排序
                {
                    bigger=Vec_cov[k][j+1];
                    Vec_cov[k][j+1]=Vec_cov[k][j];
                    Vec_cov[k][j]=bigger;
                }
            }
        }
    }
    printf("协方差特征值矩阵：\n");
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            printf("%f\t",EV_cov[i][j]);
        }
        printf("\n");
    }
    printf("协方差特征值列向量矩阵：\n");
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            printf("%f\t",Vec_cov[i][j]);
        }
        printf("\n");
    }

    if((fp=fopen("E:\\C_Projects\\computation_methods\\HW4/coordinate.txt","w"))==NULL)
    {
        printf("can't open");
        exit(0);
    }
    for(i=0;i<150;i++)//投影操作，用x[150],y[150]记录两个方向的坐标
    {
        vx=0;
        vy=0;
        for(j=0;j<4;j++)
        {
            vx+=Vec_cov[j][0]*iris_4[i][j];
            vy+=Vec_cov[j][1]*iris_4[i][j];
        }
        x[i]=vx;
        y[i]=vy;
        fprintf(fp,"%f\t",x[i]);
        fprintf(fp,"%f\n",y[i]);
    }
    fclose(fp);
    system("pause");
    return 0;
}