#include <stdio.h>
#include <string.h>
#include <stdlib.h>
int main()
{
    char buffer[1024] ;
    char *record,*line;
    int i=0,j=0;
    double temp = 0.0;
    double mat[100][3];
    FILE *fstream = fopen("data.txt","r");
    for(int i = 0; i < 100; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            fscanf(fstream,"%f",&temp);
            mat[i][j] = temp;
        }
    }
    return 0 ;
}