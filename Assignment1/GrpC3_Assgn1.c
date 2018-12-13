		// EE1103 - Assignment 1 -  Min, Max, Mean, StdDev

//Done by Anirudh R(EE18B103), Akilesh Kannan(EE18B122), Aravint Annamalai(EE18B125)

/* Commands used in Linux terminal : 

Assuming that file formatted data, i.e., only times taken :

more ping.txt | awk 'BEGIN{max=0;min=10000;sum=0;sum2=0 {if($1>max)max=$1;if($1<min)min=$1;sum+=$1;sum2+=$1*$1} END {print max,min,sum/NR,sqrt(sum2/NR-(sum/NR)*(sum/NR))}'

To compile and run the C file in terminal:

gcc GrpC3-Assgn1-14Aug.c -lm; ./a.out
*/


#include<stdio.h>
#include<math.h>

int main()
{
	float sum2=0, sum=0, max=0, min=10000, x, stddev, avg;
	int nr=0;
	/*
	  sum2 for storing value of sum of squares of time(s)
	  sum for storing value of sum of time(s)
	  max for storing value of maximum time
	  min for storing value of minimum time
	  stddev for storing value of standard deviation
	  avg for storing value of mean (average)
	  nr for storing the value of number of rows of data (no. of data values)
	*/
	FILE *datfilptr = fopen("ping.txt","r");        //datfilptr is the data file pointer and this statement is aimed at opening the file ping.txt

	
	printf("\tDISCLAIMER: All times mentioned here are measured in milliseconds\n\n");
		if(datfilptr == NULL)
		{
			printf("Error opening the file\n");           //error message in case file does not open
		}

		else
		{
			 for( ; fscanf(datfilptr,"%f",&x) && !feof(datfilptr) ; )     //reading from file and entering into the loop
			 {
				sum += x;              //updating sum
				sum2 += x*x;           //updating sum of squares

				nr = nr+1;             //updating number of rows

				if(max<x)
					max=x;         //finding maximum time to send and receive a packet

				if(min>x)
					min=x;         //finding minimum time to send and receive a packet
			 }

			fclose(datfilptr);             //closing file

			avg = sum/nr;                  //computing average
			stddev = sqrt((sum2/nr)-(avg*avg));    //computing standard deviation

                        //round trip here means a packet being sent and received

			printf("\nthe average time taken for 1 round trip is: %f", avg);      
			printf("\nthe maximum time taken for 1 round trip is: %f", max);      
			printf("\nthe minimum time taken for 1 round trip is: %f", min);      
			printf("\nthe standard deviation for time taken for 1 round trip is: %f", stddev);   
		}

	return 0;
}
