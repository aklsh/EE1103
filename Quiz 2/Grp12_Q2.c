/***********************************************************************************************************************************************************************************************
Quiz 2 - EE1103 - Numerical Methods

The program's objective is to solve the differential equation governing the dynamics of a magnetic moment in presence of an applied magnetic field.
Upon execution, the program asks user for input of a question number as given in Quiz-2 (https://courses.iitm.ac.in/pluginfile.php/167556/mod_assign/introattachment/0/Q2.pdf?forcedownload=1),
after which the program will output the required plot(s) in GNUPLOT window(s).

Name of Program - Grp12_Q2.c
Author(s) - Ch Sarvani (ee18b120), Abhijeet Ajithkumar (ee18b121), Akilesh Kannan (ee18b122)
Date created - 19-October-2018

Input(s): question number (1 or 2 or 3 or 4)

How to execute:
				gcc Grp12_Q2.c -lm -o grp12q2
				./grp12q2

Outputs: GNUPLOT windows with the plot(s) of trajectories followed by a magnetic moment in presence of an applied magnetic field, obtained by solving the DE(s) in 2 different methods
*************************************************************************************************************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


//dummy arrays used in functions - these are declared as global variables because compiler will not be able to assign the array returned from the function -
//- to the required variable if the array returned is declared locally, inside the function.

double cp[3]; 
double add[3];
double cm[3]; 

double Mg[2000000][3];                                                    

void delay(int number_of_seconds) 
{   
    int milli_seconds = 1000 * number_of_seconds;                    // Converting time into milli_seconds 

    clock_t start_time = clock();                                    
//Storing start time using clock(), a standard library function returns the number of seconds ticked from the beginning of the program
   
    while (clock() < start_time + milli_seconds);                    // looping till required time is not acheived
} 

//to generate random number using polar form of box-muller transform
double rndgen() 
{
	static double U1,U2;
	static int p =0;

	double rnd;
	if(p==0)
	{
		U1=(rand()+1.0)/(RAND_MAX+2.0);
		U2=rand()/(RAND_MAX+1.0);
		rnd=sqrt(-2.0*log(U1))*cos(2*3.1415*U2)/2000.0;				//the division by 2000.0 is done to make the random number small
	}

	else
		rnd=sqrt(-2.0*log(U1))*sin(2*3.1415*U2)/2000.0;				//the division by 2000.0 is done to make the random number small
	
	p = 1-p;

	return rnd;
}

//to evaluate the cross product of 2 vectors stored in form of arrays and return the resultant vector as an array
double* crossprod(double a[], double b[])                           
{
	cp[0] = a[1]*b[2] - b[1]*a[2]; 
	cp[1] = a[2]*b[0] - b[2]*a[0];
	cp[2] = a[0]*b[1] - b[0]*a[1];

	return cp;
}

// to evaluate and return the dot product of 2 vectors (for correlation calculation)
double dotprod(double a[], double b[])
{
	double dp;

	dp+=(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);

	return dp;
}

//to multiply a vector by a scalar constant and return the result as an array 
double* mulconst(double a[], double h)                             
{
	for(int i=0;i<3;i++)
		cm[i] = h*a[i];

	return cm;
}

//to add two vectors and return the resultant vector as an array
double* addvectors(double a[], double b[])                         
{
	for(int i=0;i<3;i++)
		add[i] = a[i]+b[i];

	return add;
}

//to solve the DE using RK4 method and to plot the obtained trajectory in gnuplot
void rksolver(double Mo[], double H[], double a, double step)     
{
	
	FILE* fileptr = fopen("rkdata.txt", "w");

	memcpy(Mg[0],Mo, 3*sizeof(double));
//memcpy() is a standard library function used to copy a block of memory from one location to another. 
//memcpy() is an easier way of copying arrays.

	fprintf(fileptr,"%lf %lf %lf\n", Mg[0][0], Mg[0][1], Mg[0][2]);

	for(unsigned long int j=1;j<200000*0.001/step;j++)
	{
		double *mxh = malloc(3 * sizeof(double));
		double *mxmxh = malloc(3 * sizeof(double));

		memcpy(mxh,crossprod(Mg[j-1],H), 3*sizeof(double));                          
		memcpy(mxmxh,crossprod(Mg[j-1],mxh), 3*sizeof(double));

		double* k1 = malloc(3 * sizeof(double));
		double* k2 = malloc(3 * sizeof(double));
		double* k3 = malloc(3 * sizeof(double));
		double* k4 = malloc(3 * sizeof(double));

		//calculating the values of k1,k2,k3,k4 required to solve the DE using RK4 method

		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));

		double* temp = malloc(3*sizeof(double));
		memcpy(temp,Mg[j-1], 3*sizeof(double));
		
		memcpy(k1,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(Mg[j-1],addvectors(temp, mulconst(k1,0.5)), 3*sizeof(double));
		
		memcpy(mxh,crossprod(Mg[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(Mg[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
		
		memcpy(k2,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(Mg[j-1],addvectors(temp, mulconst(k2,0.5)), 3*sizeof(double));
		
		memcpy(mxh,crossprod(Mg[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(Mg[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
		
		memcpy(k3,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(Mg[j-1],addvectors(temp, mulconst(k3,1.0)), 3*sizeof(double));
		
		memcpy(mxh, crossprod(Mg[j-1],H), 3*sizeof(double));
		memcpy(mxmxh, crossprod(Mg[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
	
		memcpy(k4,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(Mg[j-1],temp, 3*sizeof(double));

		Mg[j][0] = (1.0/6.0)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])+Mg[j-1][0];
		Mg[j][1] = (1.0/6.0)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])+Mg[j-1][1];
		Mg[j][2] = (1.0/6.0)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])+Mg[j-1][2];

		double mag;
		mag = sqrt(Mg[j][0]*Mg[j][0] + Mg[j][1]*Mg[j][1] + Mg[j][2]*Mg[j][2]);	            //finding magnitude for normalisation

		Mg[j][0]/=mag;
		Mg[j][1]/=mag;
		Mg[j][2]/=mag;							//normalisation of values to stay on the unit sphere

		fprintf(fileptr,"%lf %lf %lf\n", Mg[j][0],Mg[j][1],Mg[j][2]);	
	}

	if(step!=0.0001)	//to not generate the plot for std=0.0001 as it is not necessary. We only collect the data for std=0.0001 because it is used in another question.
	{
		//plotting the data in gnuplot
		FILE* gnuplotptr = popen("gnuplot -persistent", "w");
		fprintf(gnuplotptr,"set title 'Trajectory of M by RK45 solver with step size of %lf '\nset xlabel 'Mx'\nset ylabel 'My'\nset zlabel 'Mz'\n", step);
		fprintf(gnuplotptr,"set view equal xyz\n splot 'rkdata.txt' u 1:2:3 w l ls 3 title 'trajectory in 3d space'\n");

		fflush(gnuplotptr);
		pclose(gnuplotptr);
	}

	fclose(fileptr);
	
}

//to solve the DE using RK4 method and to plot the graph of alpha vs switchtime in gnuplot
void rksolverforalpha(double Mo[], double H[], double a, double step)                   
{
	FILE* fileptr = fopen("alphatime.txt", "a");
	
	double** M = malloc((200000*0.001/step) * sizeof(double**));

	for(unsigned long int i=0;i<200/step;i++)
		M[i] = (double*) malloc(3 * sizeof(double));

	memcpy(M[0],Mo, 3*sizeof(double));
	double switchtime;
	
	for(unsigned long int j=1;j<200000*0.001/step;j++)
	{
		double *mxh = malloc(3 * sizeof(double));
		double *mxmxh = malloc(3 * sizeof(double));

		memcpy(mxh,crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(M[j-1],mxh), 3*sizeof(double));

		double* k1 = malloc(3 * sizeof(double));
		double* k2 = malloc(3 * sizeof(double));
		double* k3 = malloc(3 * sizeof(double));
		double* k4 = malloc(3 * sizeof(double));

		//calculating the values of k1,k2,k3,k4 required to solve the DE using RK4 method

		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));		
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));	

		double* temp = malloc(3*sizeof(double));
		memcpy(temp,M[j-1], 3*sizeof(double));
		
		memcpy(k1,addvectors(mxh,mxmxh), 3*sizeof(double));	

		memcpy(M[j-1],addvectors(temp, mulconst(k1,0.5)), 3*sizeof(double));
		
		memcpy(mxh,crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(M[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
		
		memcpy(k2,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(M[j-1],addvectors(temp, mulconst(k2,0.5)), 3*sizeof(double));
		
		memcpy(mxh,crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(M[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
		
		memcpy(k3,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(M[j-1],addvectors(temp, mulconst(k3,1.0)), 3*sizeof(double));
		
		memcpy(mxh, crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh, crossprod(M[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
	
		memcpy(k4,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(M[j-1],temp, 3*sizeof(double));

		M[j][0] = (1.0/6.0)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])+M[j-1][0];
		M[j][1] = (1.0/6.0)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])+M[j-1][1];
		M[j][2] = (1.0/6.0)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])+M[j-1][2];

		double mag;
		mag = sqrt(M[j][0]*M[j][0] + M[j][1]*M[j][1] + M[j][2]*M[j][2]);	         //finding magnitude for normalisation

		M[j][0]/=mag;
		M[j][1]/=mag;
		M[j][2]/=mag;							//normalisation of values to stay n the unit sphere

		if(M[j][0]>=0)
		{
			switchtime = j/1000.0;	//finding switchtime - the division by 1000.0 is done as we have our step-size as 0.001
			break;
		}

	}

	if(switchtime>0.001)
		fprintf(fileptr,"%lf %.10lf\n",-1*a,switchtime);	
// for the Landau-Lifshitz equation given in the quiz paper, the value of alpha used has to be negative

	fclose(fileptr);
}

//to solve the DE using euler's method and to plot the obtained trajectory in gnuplot
void euler(double Mo[], double H[], double a, double step)                              
{
	FILE* gnuplotptr = popen("gnuplot -persistent", "w");
	FILE* fileptr = fopen("eulerdata.txt", "w");
	
	fprintf(gnuplotptr,"set title 'Trajectory of M by euler method with step size %lf'\nset xlabel 'Mx'\nset ylabel 'My'\nset zlabel 'Mz'\n", step);

	double** M = malloc((200000*0.001/step) * sizeof(double**));

	for(unsigned long int i=0;i<200/step;i++)
		M[i] = (double*) malloc(3 * sizeof(double));

	memcpy(M[0],Mo, 3*sizeof(double));

	fprintf(fileptr,"%lf %lf %lf\n", M[0][0], M[0][1], M[0][2]);

	for(unsigned long int j=1;j<200000*0.001/step;j++)
	{
		double *mxh = malloc(3 * sizeof(double));
		double *mxmxh = malloc(3 * sizeof(double));

		memcpy(mxh,crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(M[j-1],mxh), 3*sizeof(double));

		M[j][0] = (-1 * (mxh[0]) - a * mxmxh[0])*step + M[j-1][0];
		M[j][1] = (-1 * (mxh[1]) - a * mxmxh[1])*step + M[j-1][1];
		M[j][2] = (-1 * (mxh[2]) - a * mxmxh[2])*step + M[j-1][2];

		double mag;
		mag = sqrt(M[j][0]*M[j][0] + M[j][1]*M[j][1] + M[j][2]*M[j][2]);		//finding magnitude of vector for normalisation

		M[j][0]/=mag;
		M[j][1]/=mag;
		M[j][2]/=mag;							//normalisation of values to stay on the unit sphere


		fprintf(fileptr,"%lf %lf %lf\n", M[j][0],M[j][1],M[j][2]);	

		free(mxh);
		free(mxmxh);	
	}

	fprintf(gnuplotptr,"set view equal xyz\n splot 'eulerdata.txt' u 1:2:3 w l ls 3 title 'trajectory in 3d space'\n");

	fflush(gnuplotptr);

	fclose(fileptr);
	pclose(gnuplotptr);
}
              
//to solve the DE using the RK4 method for a noisy system
//to plot the obtained trajectory in gnuplot
//to determine the correlation between the original and the noisy simulations
void noise(double Mo[], double H[], double a, double step,double std)
{
	srand(0);
	
	FILE* fileptrnoise = fopen("noise.txt", "w");
	FILE* filecor = fopen("correlation.txt", "a");
	
	double** M = malloc((200000*0.001/step) * sizeof(double**));

	for(unsigned long int i=0;i<200.0/step;i++)
		M[i] = (double*) malloc(3 * sizeof(double));

	memcpy(M[0],Mo, 3*sizeof(double));

	double cor;
	
	for(unsigned long int j=1;j<200000*0.001/step;j++)
	{
		double *mxh = malloc(3 * sizeof(double));
		double *mxmxh = malloc(3 * sizeof(double));

		memcpy(mxh,crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(M[j-1],mxh), 3*sizeof(double));

		double* k1 = malloc(3 * sizeof(double));
		double* k2 = malloc(3 * sizeof(double));
		double* k3 = malloc(3 * sizeof(double));
		double* k4 = malloc(3 * sizeof(double));

		//calculating the values of k1,k2,k3,k4 required to solve the DE using RK4 method

		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));

		double* temp = malloc(3*sizeof(double));
		memcpy(temp,M[j-1], 3*sizeof(double));
		
		memcpy(k1,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(M[j-1],addvectors(temp, mulconst(k1,0.5)), 3*sizeof(double));
		
		memcpy(mxh,crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(M[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
		
		memcpy(k2,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(M[j-1],addvectors(temp, mulconst(k2,0.5)), 3*sizeof(double));
		
		memcpy(mxh,crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh,crossprod(M[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
		
		memcpy(k3,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(M[j-1],addvectors(temp, mulconst(k3,1.0)), 3*sizeof(double));
		
		memcpy(mxh, crossprod(M[j-1],H), 3*sizeof(double));
		memcpy(mxmxh, crossprod(M[j-1],mxh), 3*sizeof(double));
		
		memcpy(mxh,mulconst(mxh,-1*step), 3*sizeof(double));
		memcpy(mxmxh,mulconst(mxmxh,-1*a*step), 3*sizeof(double));
	
		memcpy(k4,addvectors(mxh,mxmxh), 3*sizeof(double));

		memcpy(M[j-1],temp, 3*sizeof(double));

		//adding noise to values obtained through rksolver

		M[j][0] = (1.0/6.0)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])+M[j-1][0] + std*rndgen();
		M[j][1] = (1.0/6.0)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])+M[j-1][1] + std*rndgen();
		M[j][2] = (1.0/6.0)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])+M[j-1][2] + std*rndgen();			

		double mag;
		mag = sqrt(M[j][0]*M[j][0] + M[j][1]*M[j][1] + M[j][2]*M[j][2]);			//finding magnitude for normalisation

		M[j][0]/=mag;
		M[j][1]/=mag;
		M[j][2]/=mag;								//normalisation of values to stay n the unit sphere

		fprintf(fileptrnoise,"%lf %lf %lf\n",M[j][0],M[j][1],M[j][2]);

		cor+=dotprod(M[j],Mg[j]);
	}

	cor/=(200.0/step);

	fprintf(filecor,"%lf %lf\n", std/2000.0,cor);		//we are printing std/2000.0 onto the file because we had already multiplied by 1/2000.0 inside the rndgen() function. So, the effective standard deviation will be st/2000.0

if(std==0.1||std==0.9||std==1.3) //plotting the trajectory only for these standard deviation values - only to show the user how the trajectory looks like with noise
	{
		FILE* gnuplotptr = popen("gnuplot -persistent", "w");
		fprintf(gnuplotptr,"set title 'Trajectory of M with noise - standard deviation %lf'\nset xlabel 'Mx'\nset ylabel 'My'\nset zlabel 'Mz'\n", std/2000.0);
		fprintf(gnuplotptr,"set view equal xyz\n splot 'noise.txt' u 1:2:3 w l ls 3 title 'trajectory in 3d space'\n");
		fflush(gnuplotptr);
		pclose(gnuplotptr);
	}

	fclose(fileptrnoise);
	fclose(filecor);
	
}

int main()
{
	double Mo[] = {-0.99,0.0,0.0};
	double H[] = {1.0,0.01,0.0};

	double a = 0.05;
	int ch=1;
	double std;

	do
	{
		int choice;

		printf("\nEnter the question number as given in the question paper to get the required plots for that question:\n");
		scanf("%d", &choice);

		switch(choice)								//switch statement for easy access to each part of the program
		{
			case 1:								//euler method
				
				printf("The following 8 gnuplot windows show the trajectory obtained for various step sizes using euler's method, which are mentioned in the title of the plot\n");
				
				for(int i=0;i<32768;i++)		//to generate a pause
					for(int j=0;j<32768;j++);
				
				//Calling the euler function for various step sizes
				euler(Mo,H,a,0.001);
				delay(30);
				euler(Mo,H,a,0.005);
				delay(30);
				euler(Mo,H,a,0.01);
				delay(30);
				euler(Mo,H,a,0.05);
				delay(30);
				euler(Mo,H,a,0.1);
				delay(30);
				euler(Mo,H,a,0.12);
				delay(30);
				euler(Mo,H,a,0.15);
				delay(30);
				euler(Mo,H,a,0.2);
				delay(30);

				printf("\nFrom these plots we can see that the maximum step size for a stable solution is somewhere between 0.1 - 0.2.\n");
				break;

			case 2:								//runge-kutta (RK4) method

				printf("The following 9 gnuplot windows show the trajectory obtained for various step sizes using Runge-Kutta (RK4) method, which are mentioned in the title of the plot\n");
				
				for(int i=0;i<32768;i++)		//to generate a pause
					for(int j=0;j<32768;j++);

				//Calling the rksolver function for various step sizes

				rksolver(Mo,H,a,0.001);
				delay(30);
				rksolver(Mo,H,a,0.02);
				delay(30);
				rksolver(Mo,H,a,0.1);
				delay(30);
				rksolver(Mo,H,a,0.2);
				delay(30);
				rksolver(Mo,H,a,0.4);
				delay(30);
				rksolver(Mo,H,a,0.45);
				delay(30);
				rksolver(Mo,H,a,0.5);
				delay(30);
				rksolver(Mo,H,a,0.6);
				delay(30);
				rksolver(Mo,H,a,0.0001);  				//calling this(stepsize=0.0001) at last because we need this data for q4. Plot will not be generated.
				delay(30);

				printf("\nFrom these plots we can see that the maximum step size for a stable solution is somewhere around 0.4.\n");
				break;

			case 3:								//time for magnetisation to switch

				printf("\nPlease wait for some time as the process of finding switchtime for 100 alpha values is highly time consuming.\n");

				for(double i=0;i<100;i++)
				{
					double a1 = pow(10, -1*(4-(0.04*i)));
					rksolverforalpha(Mo,H,a1,0.001);
				}

				//plotting the data in gnuplot

				FILE* gnuplotptr = popen("gnuplot -persistent", "w");
				fprintf(gnuplotptr,"set title 'plot of switching time vs. alpha'\n set xlabel 'alpha (a)'\n set ylabel 'switchtime (t)'\n");
				fprintf(gnuplotptr,"plot 'alphatime.txt' u 1:2 w l title 't = f(a)'\n");

				fflush(gnuplotptr);
				pclose(gnuplotptr);

				break;

			case 4:								//noise to simulate a real physical system

				printf("\nPlease wait for some time. Each plot takes time to generate.\n");

				for(int i=1;i<=50;i+=2)
				{
					std = i/10.0;  
					noise(Mo,H,0.05,0.0001,std);
				}

				//plotting the data in gnuplot

				FILE* gnuplotcor = popen("gnuplot -persistent","w");
				fprintf(gnuplotcor,"set title 'plot of correlation vs. standard deviation'\n set xlabel 'standard deviation (std)'\n set ylabel 'correlation (cor)'\n");
				fprintf(gnuplotcor,"plot 'correlation.txt' u 1:2 w linespoints title 'cor = f(std)' ls 1\n");

				fflush(gnuplotcor);
				pclose(gnuplotcor);
				break;

			default:
				printf("\nEnter a valid question number (1-4 only allowed)");
		}

		printf("\nDo you want to continue to look at other questions? Enter 1 to do so. Else, press any other number on the keyboard.");
		scanf("%d", &ch);


		if(ch!=1)
		{
			//deleting the files created during runtime of program. This section can be commented if user wants to look at the files after the termination of program
			remove("alphatime.txt");
			remove("rkdata.txt");
			remove("eulerdata.txt");	
			remove("noise.txt");
			remove("correlation.txt");
		}

	}while(ch==1);

	return 0;
}
