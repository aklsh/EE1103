/*
***********************************************************************************************************************************************
EE1103-Numerical Methods
Assignment on implementation of cubic spline using LU decomposition
This program will:
	-use LU decomposition instead of Gauss Elimination to find the coefficients of the cubic polynomials
Done by:  Ch Sarvani (EE18B120), Abhijeet Ajithkumar (EE18B121), Akilesh Kannan (EE18B122)		
***********************************************************************************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define num 502   
#define n 101
#define N n-1

//To decompose given tridiagonal matrix into upper and lower triangular matrices
void LUdecompose(double tri[N][N],double l[N][N],double u[N][N])
{
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
        {
        	u[i][j]=tri[i][j];
        	if(i<j)
            	l[i][j]=0;
        	else if(i==j)
            	l[i][j]=1;

        }
    for(int m=0;m<N;m++)
        for(int i=m+1;i<N;i++)
    {
        double factor=u[i][m]/u[m][m];
        for(int j=0;j<N;j++)
           {
           		u[i][j]-=(factor*u[m][j]);
           		l[i][m]=factor;
           }
    }

}

//To solve for the second derivatives from L and U
void solver(double l[N][N],double u[N][N],double sigTemp[N],double B[N])
{
    double Y[N];
    int i,j;
    for(i=0; i<N; i++)
    {
        Y[i]=B[i];
        for(j=0; j<i; j++)
        {
            Y[i]-=l[i][j]*Y[j];
        }
    }
    for(i=N-1; i>=0; i--)
    {
        sigTemp[i]= Y[i];
        for(j=i+1; j<n; j++)
        {
            sigTemp[i]-=u[i][j]*sigTemp[j];
        }
        sigTemp[i]/=u[i][i];
    }

}

//To calculate the coefficients of each of the cubic polynomials using the second derivatives
void cSCoeffCalc(double h[N+1], double sig[n], double y[n], double a[n], double b[n], double c[n], double d[n])
{
    int i;

    for(i=0;i<n;i++)
    {
        d[i]=y[i];
        b[i]=sig[i]/2.0;
        a[i]=(sig[i+1]-sig[i])/(h[i]*6.0);
        c[i]=(y[i+1]-y[i])/h[i]-h[i]*(2*sig[i]+sig[i+1])/6.0;
    }
}

//To generate the simplified tridiagonal matrix
void tridiagonalCubicSplineGen(double h[N+1], double tri[N][N])
{
    int i;
    for(i=0;i<n-1;i++)
    {
        tri[i][i]=2*(h[i]+h[i+1]);
    }

    for(i=0;i<n-2;i++)
    {
        tri[i][i+1]=h[i+1];
        tri[i+1][i]=h[i+1];
    }

}

//To calculate the interpolated value
double calcintpol(int x, double a[],double b[], double c[], double d[])
{
    int m = x/5;
    double calcval;

    calcval = a[m]*pow((x-(5*m + 1)),3) + b[m]*pow((x-(5*m + 1)),2) + c[m]*(x-(5*m + 1)) + d[m];
    return calcval;
}


void cspline(double x[n],double y[n],double arr[num][2])
{
	double h[N+1];
	
	for(int i=0;i<=N;i++)
		h[i]=x[i+1]-x[i];
		
    double l[N][N],u[N][N];
    double a[n]; //array to store the ai's
    double b[n]; //array to store the bi's
    double c[n]; //array to store the ci's
    double d[n]; //array to store the di's
    double sig[n],sigTemp[N]; //array to store Si's
  
    sig[0]=0;	//since we have assumed the second derivative to be 0 at end points.
    sig[N]=0;
  
    double tri[N][N];
    tridiagonalCubicSplineGen(h,tri);
	
    double B[N];
    for(int i=1;i<N+1;i++)
        B[i-1]=(y[i+1]-y[i])*6/(double)h[i]-(y[i]-y[i-1])*6/(double)h[i-1];
		
    LUdecompose(tri,l,u);
	
    solver(l,u,sigTemp,B);
	
    for(int i=1;i<n;i++)
    {
        sig[i]=sigTemp[i-1];
    }
	
    cSCoeffCalc(h,sig,y,a,b,c,d);
    
	//To calculate the interpolated value, and hence compute root mean square error
	double cspintval[num];
	double rmse,sum2;

    for(int r=0;r<502;r++)
    {
		cspintval[r]=calcintpol(r+1,a,b,c,d);
		sum2+=pow((cspintval[r]-arr[r][1]),2);
    }
    rmse=sqrt(sum2/num);
    printf("%lf is the error value.",rmse);

}

int main()
{
	FILE *fptr = fopen("out1_test0.csv","r"); //opening the data file

	if(fptr==NULL)
		printf("Error");

	double arr[num][2],x[n],y[n]; 	//arrays x and y store the downsampled data
	int count=0;
	for(int z=0;z<num;z++)
		arr[z][0] = z+1;

	for(int i=0;i<num;i++)
	{
		fscanf(fptr,"%lf",&(arr[i][1]));
		count++;
		if(count%5==1)
			{
				x[i/5]=arr[i][0];
				y[i/5]=arr[i][1];
			}

		if(feof(fptr))
			break;
	}
	
	cspline(x,y,arr);
}
