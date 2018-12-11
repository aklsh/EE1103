#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 110315245
#define C 12345						//assigning values for A, C, M
#define M 4294967296							

/* Commands used in terminal to execute the programs:-
	
	student@pclab-desktop:~/Desktop$ gcc GrpC3_Quiz1.c -lm -o GrpC3_Quiz1
	student@pclab-desktop:~/Desktop$ ./GrpC3_Quiz1 13456 45690

Disclaimer: The numbers used in the above command maybe replaced by other numbers as per the requirement of the user. 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Answers to the question 1(b) and 1(c)
	1(b):
		To generate numbers between 0 and 1, we can use the function |sin(x)| of the randomly generated values.
		
	The Program is given below:
			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc,char **argv)				//values accepted from user for number of random numbers and seed value
{
	unsigned long int n = atol(argv[1]);		//n contains value for the number of random numbers to be generated.
	
	long int *ptr = malloc((n+1) * sizeof(long int));	//ptr points to the starting location of the dynamically allocated array
	
	long double *ptr1 = malloc((n+1) * sizeof(long double));	//new array allocated dynamically to store values between 0 and 1 based on the originally generated random numbers.
	
	*(ptr) = atol(argv[2]);		//assigning seed value to first element of array

	for (unsigned long int i=0 ; i<n ; i++)		//loop to generate random numbers and printing them
		 *(ptr+i+1) = (A * (*(ptr+i)) + C) % M;		//generating random numbers using given expression
 
	for (unsigned long int i=0 ; i<n+1 ; i++)
		{
		 	*(ptr1+i)=fabs(sin(*(ptr+i)));	//using |sin(x)| function to generate numbers from 0 to 1

			printf("%Le " , *(ptr1+i));	//printing random numbers generated in the iteration
		}
	
	long double sum1 = 0 , sum12 = 0 , avg1 , stdev1;		//declaring variables for the calculation of mean and standard deviation
   	
	for(unsigned long int i=0 ; i < (n+1) ; i++)		//for loop to calculate sum and sum of squares of the random numbers generated 
		{
			sum1 += *(ptr1 + i);
			sum12 += (*(ptr1 + i)) * (*(ptr1 + i));
		}
	
	avg1 = sum1/(n+1);		//calculating average of the random numbers generated
	
	stdev1 = sqrt((sum12/(n+1)) - (avg1*avg1));		//calculating standard deviation of the random numbers
	
	printf("\n The average of the random numbers generated is %Le \n The standard deviation of the random numbers generated is %Le \n" , avg1 , stdev1);					//printing mean and standard deviation calculated

	
	free(ptr1);
	free(ptr);
		
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
	1(c):
		The main disadvantage is that there is a upper limit on the values that can be generated - As we are performing modulo 			m, the max value that can be generated is m-1, and the minimum is 0. Also, the random numbers generated using the given 		expression lack slightly in randomness that they are alternately odd and even.
		
		If we decrease the value of m, the upper limit of random number generated is reduced, thus reducing the randomness of the 			numbers generated if there is a large number of values that have to be generated.
		
	
*/ 	  
	  


int main(int argc,char **argv)				//values accepted from user for number of random numbers and seed value

{
	unsigned long int n = atol(argv[1]);		//n contains value for the number of random numbers to be generated.
	unsigned long int *ptr = malloc((n+1) * sizeof(unsigned long int));	/*ptr points to the starting location of the dynamically 											   allocated array*/
	
	*(ptr) = atol(argv[2]);		//assigning seed value to first element of array

	for (unsigned long int i=0 ; i<n ; i++)		//loop to generate random numbers and printing them
		{
			 printf("%ld \n" , *(ptr+i));	//printing random numbers generated in the previous iteration
			 *(ptr+i+1) = (A * (*(ptr+i)) + C) % M;		//generating random numbers using given expression
   		}
   		
   	unsigned long int sum = 0 , sum2 = 0 , avg , stdev;	//declaring variables for the calculation of mean and standard deviation
   	
	for(unsigned long int i=0 ; i < (n+1) ; i++)		//for loop to calculate sum and sum of squares of the random numbers generated 
		{
			sum += *(ptr + i);
			sum2 += (*(ptr + i)) * (*(ptr + i));
		}
	
	avg = sum/(n+1);		//calculating average of the random numbers generated
	
	stdev = sqrt((sum2/(n+1)) - (avg*avg));		//calculating standard deviation of the random numbers
	
	printf("\n The average of the random numbers generated is %ld \n The standard deviation of the random numbers generated is %ld \n" , avg , stdev);					//printing mean and standard deviation calculated

	free(ptr);			//freeing the allocated memory
	
}
   		
   			  

   

