//EE1103 Assignment-3-Downsampling and Interpolation
//Done by Group 12 - Ch Sarvani (ee18b120), Abhijeet Ajithkumar (ee18b121), Akilesh Kannan (ee18b122)
//This program is for generating a file with the downsampled data (taking every 5th point of original data points) 

#include <stdio.h>

int main()
{
	FILE *fptr = fopen("out1_test0.csv","r");
	
	if(fptr==NULL)
		printf("Error");
	
	int num=502;

	int r;

	r = num/5 + 1;

	double arr[502][2];
	
	double subarr[101][2];
	
	int count=0;

	for(int z=0;z<num;z++)
		arr[z][0] = z+1;
	
	for(int i=0;i<num;i++)
	{
		fscanf(fptr,"%lf",&(arr[i][1]));		

		count++;
		
		if(count%5==1)
			{
				subarr[i/5][1]=arr[i][1];
				subarr[i/5][0]=arr[i][0];
			}
		
		if(feof(fptr))  
			break;
	}
	for(int k=0;k<r;k++)
		printf("%lf %lf\n",subarr[k][0],subarr[k][1]);

	return 0;
}

/* How to plot the downsampled data in gnuplot:

Commands to be used in terminal:

gcc downsampleplot.c -o downsamplefile
./downsamplefile > downsample.txt
gnuplot
gnuplot> plot "downsample.txt" w lines 		// command 1

To check if the downsampled data plotted is close enought to actual plot, type the following command instead of command 1:
gnuplot> plot "downsample.txt" w lines, "out1_test0.csv" w lines

*/