//EE1103 - Assignment 4 - Integration
//Group 12 - Ch Sarvani (ee18b120), Abhijit Ajithkumar (ee18b121), Akilesh Kannan (ee18b122)
//Commands to run in terminal:

/*

gcc Grp12Assgn4.c -lm -o grp12a4

./grp12a4 7

// The number 7 in the above command can be replaced with any number (max. 25). The number is the number of iterations for romberg's method that the user requires.

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int countlinefile()
{
	FILE *ptrf = fopen("out1_test0.csv","r");

	int count =0;
	char s;

	while(!feof(ptrf))
	{
		fscanf(ptrf,"%c",&s);
		if(s=='\n')
			count++;
	}

	return count;

}

double x[101]; //array to store the x-axis points
double y[101]; //array to store the y-axis points
double h[100];   //array to store the successive interval widths

double traparea=0; //stores area of graph calculated using trapezium method

void trapz(double y[], int num) //trapezium rule
{
	for(int a=0;a<num-1;a++)
		traparea+= (2.5)*(y[a+1]+y[a]);
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


double cspline(double xreq) //preforming cubic spline interpolation 
{
    FILE *fptr = fopen("out1_test0.csv","r"); //opening the data file
	
	if(fptr==NULL)
		printf("Error");
	
	int num=countlinefile();
	
	int n = num/5 + 1;

	double arr[num][2];  // creating array to store all the data points
	
	for(int z=0;z<num;z++)
		arr[z][0] = z+1;
	
	for(int i=0;i<num;i++)
	{
		fscanf(fptr,"%lf",&(arr[i][1]));

		if(feof(fptr))  
			break;
	}

	fclose(fptr);

   
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
    
    if(xreq==0)
    	return 0;
    	
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

	double cspintval;

	cspintval=calcintpol(xreq,a,b,c,d);
	
	return cspintval;
}

void romberg(double f(double), double a, double b, double **R, int n)  //romber'g method
{
	double h = b-a, sum; 
	
	trapz(y,b/5 +1);

	R[0][0] = h*(f(a) + f(b))*0.5;
   
    for (int i=1; i<=n; i++)
    { 
     	h = h*0.5;
     	sum = 0;
     	
     	for (int k=0; k<pow(2,i); k++)
       		sum+=f(a + k*h);
     	
     	R[i][0] = 0.5*(R[i-1][0] + sum*h);  
     
     	for (int j=1; j<=i; j++)
       		R[i][j] = (pow(4,j)*R[i][j-1] - R[i-1][j-1])/(pow(4,j) - 1); 
    }

	printf("The romberg method areas after each iteration are given below in triangle form:\nOur best guess will be A[%d][%d]. Try executing the code again with a different number of iterations to see more accurate values for romberg's method\n",n-1,n-1);
    for(int c=0;c<n;c++)
 	{	
 		for(int d=0;d<=c;d++)
   			printf("%lf \t", R[c][d]);

   		printf("\n");
   	}

   	printf("\n%lf is trapezium rule area \n", traparea);
}


int main(int argc, char** argv)
{
	if(argc!=2)
		{
			fprintf(stderr,"example usage: %s 5\n", argv[0]);
			exit(0);
		}

	int n=25; //max no.of iterations
	int numline = countlinefile()-1;//number of data points in file
	
	double y=cspline(0.0); //dummy variable to call cspline function to read the file and assign x and y values
	double cspline(double);
	double **R;// 2-d array to store the areas calculated after each romberg iteration

	R = calloc((n), sizeof(double *));
  	
	for (int i=0; i<=n; i++)
    		R[i] = calloc((n), sizeof(double));

	int numiter=atoi(argv[1]); //number of iterations to which accuracy is needed.
	
	romberg(cspline,1,numline,R,numiter);

	return 0;

}
