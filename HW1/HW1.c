#include<stdio.h>
#include<math.h>
#include<stdlib.h>
double max(double a,double b)
{
	double z;
	z=a>b?a:b;
	return(z);
}
/*
double equation(double t,double xp,double xq,double yq)
{
	double t_1,t_2;
	t_1=t/(-1-xp*sqrt(pow(t,2)+1));
	t_2=(yq*sqrt(pow(t,2)+1)-t)/(-1-xq*sqrt(pow(t,2)+1));
	return(2*t/(1-pow(t,2))-(t_2-t_1)/(1+t_1*t_2));
}
*/

double equation_c(double c,double xp,double xq,double yq)
{
	double c_1,c_2,s,s_1,s_2;
	s=sqrt(1-pow(c,2));
	c_1=(-c-xp)/(sqrt(pow(s,2)+pow((c+xp),2)));
	c_2=(-c-xq)/(sqrt(pow(yq-s,2)+pow(c+xq,2)));
	s_1=sqrt(1-pow(c_1,2));
	s_2=sqrt(1-pow(c_2,2));
	return (1-2*pow(c,2)+c_2*c_1+s_1*s_2);
}
double equation_s(double s,double xp,double xq,double yq)
{
	double c_1,c_2,c,s_1,s_2;
	c=sqrt(1-pow(s,2));
	c_1=(-c-xp)/(sqrt(pow(s,2)+pow((c+xp),2)));
	c_2=(-c-xq)/(sqrt(pow(yq-s,2)+pow(c+xq,2)));
	s_1=sqrt(1-pow(c_1,2));
	s_2=sqrt(1-pow(c_2,2));
	return (2*s*c-s_2*c_1+s_1*c_2);
}
double solve_eq_s(double xp,double xq,double yq,double epsilon)//epsilon为所需精度 
{
	double a,b,y;
	b=yq/sqrt(pow(xq,2)+pow(yq,2));
	a=0.000000000000000;//
	if((equation_s(/*sqrt(1-pow(a,2))*/a,xp,xq,yq)*equation_s(/*sqrt(1-pow(b,2))*/b,xp,xq,yq))<0)//
	{
		while((b-a)>epsilon)
		{
			y=(a+b)/2;
			if(equation_s(/*sqrt(1-pow(y,2))*/y,xp,xq,yq)>0)//
			{
				b=y;
				
				if(max((b-a),(sqrt(1-pow(a,2))-sqrt(1-pow(b,2))))<epsilon)
				{
					y=(a+b)/2;
				
					return y;//返回反射点y坐标 
					break;
				}
			}
			else
			{
				a=y;
				if(max((b-a),(sqrt(1-pow(a,2))-sqrt(1-pow(b,2))))<epsilon)
				{
					y=(a+b)/2;
					return y;//返回反射点y坐标 
					break;
				}
			}
		}
	}
    return 0;
}
double solve_eq_c(double xp,double xq,double yq,double epsilon)//epsilon为所需精度 
{
	double a,b,y;
	b=yq/sqrt(pow(xq,2)+pow(yq,2));
	a=0.000000000000000;//
	if((equation_c(sqrt(1-pow(a,2)),xp,xq,yq)*equation_c(sqrt(1-pow(b,2)),xp,xq,yq))<0)//
	{
		while((b-a)>epsilon)
		{
			y=(a+b)/2;
			if(equation_c(sqrt(1-pow(y,2)),xp,xq,yq)>0)//
			{
				b=y;
				
				if(max((b-a),(sqrt(1-pow(a,2))-sqrt(1-pow(b,2))))<epsilon)
				{
					y=(a+b)/2;
				
					return y;//返回反射点y坐标 
					break;
				}
			}
			else
			{
				a=y;
				if(max((b-a),(sqrt(1-pow(a,2))-sqrt(1-pow(b,2))))<epsilon)
				{
					y=(a+b)/2;
					return y;//返回反射点y坐标 
					break;
				}
			}
		}
	}
    return 0;
}
struct Line
{
	double a;
	double b;
	double c;//ax+by=c格式 
};
double solve_linear_y(double xp,double xq,double yq,double t)
{
	double k_1,k_2,Ry;
	k_1=t/(-1-xp*sqrt(pow(t,2)+1));
	k_2=-t;
	struct Line line[2]={{k_1,-1,k_1*xp},{k_2,-1,k_2*xq-yq}};
	Ry=(line[1].c-line[1].a/line[0].a*line[0].c)/(line[1].b-line[1].a/line[0].a*line[0].b);
	return Ry;
}
double solve_linear_x(double xp,double xq,double yq,double t)
{
	double k_1,k_2,Rx;
	k_1=t/(-1-xp*sqrt(pow(t,2)+1));
	k_2=-t;
	struct Line line[2]={{k_1,-1,k_1*xp},{k_2,-1,k_2*xq-yq}};
	Rx=(line[1].c-line[1].b/line[0].b*line[0].c)/(line[1].a-line[1].b/line[0].b*line[0].a);
	return Rx;
}
int main()
{
	int i,j;
	double a[2][3];
	FILE* f1;
	if ((f1 = fopen("test.txt", "r")) == NULL)
	{
		printf("cannot open file\n");
		exit(0);
	}
	for(i=0;!feof(f1);i++)//a[i][0],[1],[2]分别为xP,xQ,yQ 
	{
		for(j=0;j<3;j++)
		{
			fscanf(f1,"%lf",&a[i][j]);
		}
		
	}
	fclose(f1);
	double T[2][2],R[2][2];
	for(i=0;i<2;i++)
	{
		T[i][1]=solve_eq_c(a[i][0],a[i][1],a[i][2],0.000000000001);
		T[i][0]=-sqrt(1-pow(T[i][1],2));
		R[i][0]=solve_linear_x(a[i][0],a[i][1],a[i][2],-T[i][1]/T[i][0]);
		R[i][1]=solve_linear_y(a[i][0],a[i][1],a[i][2],-T[i][1]/T[i][0]);
		printf("T(%f,%f)\t",T[i][0],T[i][1]);
		printf("R(%f,%f)\n",R[i][0],R[i][1]);
	}
	printf("\n");


	double b[8][3];
	FILE* f2;
	if ((f2 = fopen("data.txt", "r")) == NULL)
	{
		printf("cannot open file\n");
		exit(0);
	}
	for(i=0;!feof(f2);i++)//a[i][0],[1],[2]分别为xP,xQ,yQ 
	{
		for(j=0;j<3;j++)
		{
			fscanf(f2,"%lf",&b[i][j]);
		}
	}
	fclose(f2);
	double T_data[8][2],R_data[8][2];
	for(i=0;i<8;i++)
	{
		if(b[i][0]>-1.01||b[i][2]<0.01)
		{
			T_data[i][1]=solve_eq_s(b[i][0],b[i][1],b[i][2],0.00000000000001);
		}
		else
		{
			T_data[i][1]=solve_eq_c(b[i][0],b[i][1],b[i][2],0.00000000000001);
		}
		T_data[i][0]=-sqrt(1-pow(T_data[i][1],2));
		R_data[i][0]=solve_linear_x(b[i][0],b[i][1],b[i][2],-T_data[i][1]/T_data[i][0]);
		R_data[i][1]=solve_linear_y(b[i][0],b[i][1],b[i][2],-T_data[i][1]/T_data[i][0]);
		printf("T(%20.15f,%20.15f)\t",T_data[i][0],T_data[i][1]);
		printf("R(%20.15f,%20.15f)\n",R_data[i][0],R_data[i][1]);
	}
	//if(!)
    system("pause");
    return 0;


}
