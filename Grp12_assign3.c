//EE1103 Assignment-3-Downsampling and Interpolation
//Done by Group 12 - Ch Sarvani (ee18b120), Abhijeet Ajithkumar (ee18b121), Akilesh Kannan (ee18b122)
//This program will implement the 3 interpolation techniques - Lagrange's method, Newton's Divided Differences, Cubic spline interpolation

//Commands to use in terminal: (These have to be typed when the working directory is the directory in which this file is in)

/*

gcc Grp12_Assignment3.c -lm -o Grp12
./Grp12

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double lagrange(double x,double arr[][2], int r) //function for calculating the interpolated value using lagrange polynomial - will return the interpolated value
{
	double intpolvallag=0;

	for(int i=0;i<r;i++)
	{
		double numer=1, denom=1;

		for(int j=0;j<r;j++)
			if(j!=i)
			{
				numer=numer*(x-arr[j][0]);
				denom=denom*(arr[i][0]-arr[j][0]);
			}
		intpolvallag+=((numer/denom)*arr[i][1]);
	}

	return intpolvallag;
}

double newton(double subarr[][2], int num, int k,double arr[][2]) //function for calculating the root mean squared error of interpolated data using Newton's method
{
	int m=k;
	double A[101][101],err=0;
//To calculate differences
	for(int l=0;l<k;l++)
    	for(int p=0;p<m;p++)
   		{
        	if(l==0)
            	A[l][p]=subarr[p][1];
        	else
            	A[l][p]=(A[l-1][p+1]-A[l-1][p])/(50*l);
        	m--;
    	}
//Error calculation
	for(int q=0;q<num;q++)
	{
        double prod,sum=0;
        double sum2=0;
    
        for(int z=0;z<k;z++)
        {
          	prod=A[z][0];
       
            for(int w=0;w<z;w++)
                prod*=(q+1-(50*(w+1)));
                sum2+=prod;
        }
        sum+=sum2;
		err+=pow((arr[q][1]-sum),2);
	}

	err/=num;

	return sqrt(err);
} 

void GauElim(int m, int n, double a[m][n], double x[n-1])//performing Gauss Elimination
{
    int i,j,k;
    for(i=0;i<m-1;i++)
    {
        for(k=i+1;k<m;k++)
	   {
            double  term=a[k][i]/ a[i][i];
            for(j=0;j<n;j++)
	        {
                a[k][j]=a[k][j]-term*a[i][j];
            }
        }
         
    }
    
    for(i=m-1;i>=0;i--)
    {
        x[i]=a[i][n-1];
        
        for(j=i+1;j<n-1;j++)
        {
            x[i]=x[i]-a[i][j]*x[j];
        }
        
        x[i]=x[i]/a[i][i];
    }           
}

void cSCoeffCalc(int n, double h[n], double sig[n+1], double y[n+1], double a[n], double b[n], double c[n], double d[n])//Cubic Spline coefficients calculator
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

void tridiagonalCubicSplineGen(int n, double h[n], double a[n-1][n], double y[n+1])//Function to generate the tridiagonal matrix 
{
    int i;
    for(i=0;i<n-1;i++)
    {
        a[i][i]=2*(h[i]+h[i+1]);
    }
   
    for(i=0;i<n-2;i++)
    {
        a[i][i+1]=h[i+1];
        a[i+1][i]=h[i+1];
    }
    
    for(i=1;i<n;i++)
    {
        a[i-1][n-1]=(y[i+1]-y[i])*6/(double)h[i]-(y[i]-y[i-1])*6/(double)h[i-1];
    }
} 

double calcintpol(int x, double a[],double b[], double c[], double d[])//calculates the interpolated value using Cubic spline method
{
    int m = x/5;
    double calcval;

    calcval = a[m]*pow((x-(5*m + 1)),3) + b[m]*pow((x-(5*m + 1)),2) + c[m]*(x-(5*m + 1)) + d[m];

    return calcval;
}


void cspline(double arr[][2], int num)
{
    double x[101]; //array to store the x-axis points
    double y[101]; //array to store the y-axis points
    double h[100];   //array to store the successive interval widths

    int count=0;
    
    for(int i=0;i<num;i++)
    {    
        count++;
        
        if(count%5==1)
            {
                y[i/5]=arr[i][1];
                x[i/5]=arr[i][0];
                h[i/5]=5;
            }
    }
    
    int n = num/5 + 1;
    
    double a[n]; //array to store the ai's
    double b[n]; //array to store the bi's
    double c[n]; //array to store the ci's
    double d[n]; //array to store the di's
    double sig[n+1]; //array to store Si's
    double sigTemp[n-1]; //array to store the Si's except S0 and Sn
    sig[0]=0;
    sig[n]=0;
    double tri[n-1][n]; //matrix to store the tridiagonal system of equations that will solve for Si's
    tridiagonalCubicSplineGen(n,h,tri,y); //to initialize tri[n-1][n]
    
    GauElim(n-1,n,tri,sigTemp);
    for(int i=1;i<n;i++)
    {
        sig[i]=sigTemp[i-1];
    }
    //calculate the values of ai's, bi's, ci's, and di's
    cSCoeffCalc(n,h,sig,y,a,b,c,d);
   

	double cspintval[502];
	double rmse,sum2;

    for(int r=0;r<502;r++)
    {
	cspintval[r]=calcintpol(r+1,a,b,c,d);
	sum2+=pow((cspintval[r]-arr[r][1]),2);
    }

    rmse=sqrt(sum2/num);

    printf("%lf is the root mean square error by using cubic spline interpolated data\n", rmse);
 
}

int main()
{
	FILE *fptr = fopen("out1_test0.csv","r"); //opening the data file
	
	if(fptr==NULL)
		printf("Error");
	
	int num=502;

	int r;

	r = num/50 + 1;

	double arr[502][2];  // creating array to store all the data points
	
	double subarr[11][2]; // creating array to store downsampled points
	
	int count=0;

	for(int z=0;z<num;z++)
		arr[z][0] = z+1;
	
	for(int i=0;i<num;i++)
	{
		fscanf(fptr,"%lf",&(arr[i][1]));		

		count++;
		
		if(count%50==1)
			{
				subarr[i/50][1]=arr[i][1];
				subarr[i/50][0]=arr[i][0];
			}
		
		if(feof(fptr))  
			break;
	}

	double vallag[502];

	double error2lag[502];
	
	for(double x=1;x<=num;x++)
		{
			vallag[(int)x-1] = lagrange(x,subarr,r);
		}

	double sumerrorsqlag=0;

	for(int t=0;t<num;t++)
	{
		double diff = arr[t][1] - vallag[t];
		error2lag[t] = pow(diff,2);
		sumerrorsqlag+=error2lag[t];
	}

	double rmselag,rmsenew;

	rmselag = sqrt(sumerrorsqlag/num);
	rmsenew = newton(subarr,num,r,arr);
	cspline(arr,num);

	//for(int i=0;i<num;i++)
	//printf("%lf \n",vallag[i]);

	printf("%lf is the root mean square error of the lagrange's polynomial interpolated data.\n", rmselag);
	printf("%lf is the root mean square error of the newton's polynomial interpolated data.\n", rmsenew);


	return 0;
}
