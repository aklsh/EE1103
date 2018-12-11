/*
***********************************************************************************************************************************************
EE1103-Numerical Methods
Assignment on LU Decomposition
This program will:
	-decompose the given coefficient matrix into its corresponding lower and upper matrices (by Gaussian elimination)
	-check if the two matrices give the original matrix on multiplication
	-find roots by means of forward and back substitution
Done by:  Ch Sarvani (EE18B120), Abhijeet Ajithkumar (EE18B121), Akilesh Kannan (EE18B122)		
***********************************************************************************************************************************************
*/

#include <stdio.h>
#include <string.h>
#define N 3

// Matrix initialisation
double a[N][N],u[N][N],l[N][N],b[N],y[N],x[N];

// To decompose given matrix A into L and U 
void LUdecompose()
{
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
        {
        u[i][j]=a[i][j];
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

//To multiply L and U, and to check if A is obtained
void multiplier()
{
    printf("Printing matrix A, as obtained by multiplication:\n\n");
    for(int i=0;i<N;i++)
	{
	   for(int j=0;j<N;j++)	
	   {
		double sum=0;
		for(int k=0;k<N;k++)
	    	   sum+=l[i][k]*u[k][j];
		printf("%lf	",sum);
	   }
	   printf("\n"); 
	}
}

//To find the roots of the equation by means of forward and back substitutions
void solution()
{
    y[0]=b[0];
    y[1]=b[1]-y[0]*a[1][0];
    y[2]=b[2]-y[0]*a[2][0]-y[1]*a[2][1];
    x[2]=y[2]/u[2][2];
    x[1]=(y[1]-u[1][2]*x[2])/u[1][1];
    x[0]=(y[0]-u[0][1]*x[1]-u[0][2]*x[2])/u[0][0];
}

void main()
{
    //initialization
    
    a[0][0]=3;
    a[0][1]=-0.1;
    a[0][2]=-0.2;
    a[1][0]=0.1;
    a[1][1]=7;
    a[1][2]=-0.3;
    a[2][0]=0.3;
    a[2][1]=-0.2;
    a[2][2]=10;
    b[0]=7.85;
    b[1]=-19.3;
    b[2]=71.4;
    LUdecompose();
    multiplier();
    solution();

//printing solutions
    printf("\nTherefore, x1=%lf,x2=%lf and x3=%lf.\n",x[0],x[1],x[2]);
}











































































