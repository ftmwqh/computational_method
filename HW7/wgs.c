#include<math.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>



int M = 8;             //设定最大迭代次数M和精度控制值e
double e = 0.000001;
long long N = 0;        //统计达到精度的比例
long long n = 0;

double ax(double t)  //x方向加速度函数
{
	return sin(t) / (sqrt(t) + 1);
}

double ay(double t)  //y方向加速度函数
{
	return log(t + 1) / (t + 1);
}

double Romberg(double begin,double end,double (*f)(double x))//Romberg积分，参数为积分区间起始点和被积函数
{
	N = N + 1;

	double R[21][21] = { 0 };
	double hk = 0, h = end - begin;

	R[1][1] = (f(begin) + f(end)) * h / 2;

	for (int k = 2; k <= M; k++)
	{
		hk = h / pow(2, k - 1);

		R[k][1] = R[k - 1][1];
		for (int i = 1; i <= pow(2, k - 2); i++)
			R[k][1] = R[k][1] + 2 * hk * f(begin + (2 * i - 1) * hk);
		R[k][1] = R[k][1] / 2;

		for (int j = 2; j <= k; j++)
			R[k][j] = R[k][j - 1] + (R[k][j - 1] - R[k - 1][j - 1]) / (pow(4, j - 1) - 1);

		if (fabsl(R[k][k] - R[k - 1][k - 1]) < e)
		{
			n = n + 1;
			return R[k][k];
		}
	}
	return R[M][M];
}

double vx(double t)   //x方向上速度函数，给定时间t，调用Romberg函数求得
{
	
		return Romberg(0, t, ax);
}

double vy(double t)   //y方向上速度函数
{
	
		return Romberg(0, t, ay);
}


int main()
{
	double x[100] = { 0 }, y[100] = {0};//记录x方向和y方向上的位移的数组

	for (int i = 0; i < 100; i++)
	{
		x[i] = Romberg(i*0.1, (i + 1) * 0.1, vx);  //求0-10每个小区间(间隔为0.1)的位移如t=0.1-0.2,9.8-9.9的位移
		y[i] = Romberg(i*0.1, (i + 1) * 0.1, vy);
	}

	for (int i = 1; i < 100; i++)
	{
		x[i] = x[i] + x[i - 1];                    //求0-10每个点处的位移,如t=0.1处的位移
		y[i] = y[i] + y[i - 1];
	}
	


	printf("%f",(double)n/N);                  //输出满足精度的比例
}