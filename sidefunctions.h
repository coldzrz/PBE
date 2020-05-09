# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
#include <omp.h>
#define __USE_C99_MATH

#include <stdbool.h>

void readparameters ( void);

bool check(double a[],double b[],int n);

void regioncount(double radial[],double height[],double umat[],double x[], double y[]);

void add(double res[],double a[],double b[],double coeff,int n);

void copy(double a[],double b[],int n);

double mtm(double a[],double b[],int n);

void regioncountmaple(int region[],double radial[],double height[],double x[], double y[]);

double *r8vec_linspace_new ( int n, double a, double b );

void print2d(int n,double a[],char *data_filename,double fa1,double fa2);

void timestamp ( void );

void print2dint(int n,int a[],char *data_filename);

void test();

double potential(double r);

void printmatrix(double val[],int col[],int row[],int count);
