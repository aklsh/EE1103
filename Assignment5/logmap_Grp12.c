//EE1103 - Assignment 5(a) - Logistic Map generation
//Done by Group 12 - Ch Sarvani (ee18b120), Abhijeet Ajithkumar (ee18b121), Akilesh Kannan (ee18b122)

//Commands to used in terminal:
/*
gcc logmap_Grp12.c -o lgmapg12
./lgmapg12
*/

#include <stdio.h>
#include <stdlib.h>
  
double width=1500, height=1500; //declaring size and width of the .pgm file
int pixel[3000][3000];  //pixel array

double start=2.8; //starting value of r in the iterative equation: x(n+1) = r * x(n) * (1 - x(n)). Can be changed.

void logistic(double x, double r)  //function for generating the x(i)'s and accordingly incrementing corresponding pixel values
{
    double xval = x;

    for(int j=1;j<10000;j++)
    {
        xval=r*(xval)*(1-xval);

        //greater the intensity (black being 0 and white being 255) of the pixel in the .pgm file, greater number of times a particular number occurs in the list of x(i)'s.
       
        if(pixel[(int)(height*xval)][(int)(width*(r-start))]<255) //since 255 (white) is the max value that a pixel can take
           pixel[(int)(height*xval)][(int)(width*(r-start))]++;

    }    

}

int main()
{
    double xstart=0.5;  //the initial value, x(0) for the iterations. Can be changed.
    
    for (int n=0;n<height;++n) 
        for (int m=0;m<width;++m) 
            pixel[n][m] = 0;  //initialising all pixel values to 0 (making all of them black)

    for(double i=0;i<=(4-start)*width;i++)
            logistic(xstart, start+(i/width));  //calling function to generate the list of x(i)'s for a particular r v

	printf("The range of r used to generate the logistic map is [2.8,4].\nA file named logistic_map_Grp12.pgm is created in the pwd.\n");
    FILE *f = fopen("logistic_map_Grp12.pgm", "w"); //declaring an object for the graph file.
    
    fprintf(f, "P2\n# comment\n%d %d\n255\n",(int)width,(int)height); //declaring the height and width of the file
           
    for (int n=0 ; n<height ; ++n)
        for (int m=0 ; m<width ; ++m) 
            fprintf(f, "%d ", pixel[n][m]); //writing onto .pgm file
    
    fclose(f);

    return 0;
}
