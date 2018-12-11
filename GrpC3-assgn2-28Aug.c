						//EE1103 - Assignment 2 - Random Numbers
		//Assignment done by Group C3 - Akilesh Kannan(EE18B122), Aravint Annamalai(EE18B125), Anirudh R(EE18B103)

//Commands to be used in Terminal to run this program:
/*

<all these commands are to be executed after changing to the directory where this source code file is stored>

gcc GrpC3-assgn2-28aug.c -o GrpC3-assgn2-28aug
./GrpC3-assgn2-28aug 100 8 6

** the numbers used in the above statement can be changed and are used only as an example **

*/
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

int rndgen(int l, int u) 			//function to generate random numbers in range [l,u]
{
	return (rand()%(u-l+1) + l);
}

int main(int argc, char *argv[])
{
	srand(time(0)); 		//to seed the initial value for random number generator
	
	int N,M,P,*ptrorig,*ptrsep,*newarr;;

	N = atoi(argv[1]); 		//User input for total number of bits
	M = atoi(argv[2]); 		//User input for offset value (from 1st bit)
	P = atoi(argv[3]); 		//User input for number of separations, i.e. number of extra bits apart from offset bit to form the new substring

	ptrorig = malloc(N * sizeof(int)); 		//Dynamic memory allocation of N bits

	for(int i=0;i<N;i++)
		*(ptrorig+i) = rndgen(0,1); 		//assigning value to each of the N bits using a random number generator

	ptrsep = malloc((P+1) * sizeof(int)); 		//Dynamic memory allocation for the separation values

	int sumsep = 0;

	for(int i=0;i<P+1;i++) 		//Generating the separation values randomly
	{
		if(i==0)
			*(ptrsep) = 0;
		else	
			*(ptrsep + i) = rndgen(1,(N-P-M-sumsep+i-1)/10 + 1);
		
		//the 2nd argument in the above function call is made in such a way that there will never be a case where the separations  
		//between any 2 adjacent bits will overshoot the total number of bits available
		// /10 + 1 is done so that there won't be an accumulation of 1s in the end - can be edited out and it will not have any consequence on the output
		sumsep += *(ptrsep + i);
	}

	newarr = malloc((P+1) * sizeof(int)); 		//Dynamic memory allocation of substring

	int sumsep1 = 0;

	for(int r=0;r<P+1;r++) 		//storing the substring to be matched for in a separate memory allocation
		{
			sumsep1 += *(ptrsep + r);

			*(newarr + r) = *(ptrorig + M + sumsep1);
		}

	printf("Summary: total number of bits is %d, original offset value is %d, number of spacings is %d\n", N,M,P);

 	printf("\nThe offset(s) where hamming distance is 0 is: \n");     //Printing offsets for the substrings where the hamming distance is minimum(i.e.)0

	for(int h=0; h < N-sumsep; h++) 		//Finding hamming distance
	{
		int hamdist = 0;
		int sumsep2 = 0;

		for(int t=0; t<P+1 ;t++)
			{
				sumsep2 += *(ptrsep + t);

				if(*(newarr + t) != *(ptrorig + h + sumsep2))
					hamdist++;
			}

                  if(hamdist==0)
			printf("%d\n", h);
                       

	}
		
	//freeing of the pointers to mallocs to avoid memory leakage

	free(ptrorig); 		
	free(ptrsep);
	free(newarr);

	printf("\n");
	

	return 0;
}
