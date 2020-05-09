#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
# include <stdbool.h>


	extern double w;
	extern double pos;
	extern double neg;
	extern double mid;
	extern int ilimit;
	extern int imin;
	extern int dt;

void memheight(double radial[],double height[]);
double f1(double s,double psi,double r,double h,double ps,double pr,double sigma,double a);
double f2(double s,double psi,double r,double h,double ps,double pr,double sigma,double a);
double f3(double s,double psi,double r,double h,double ps,double pr,double sigma,double a);
double f4(double s,double psi,double r,double h,double ps,double pr,double sigma,double a);
double f5(double s,double psi,double r,double h,double ps,double pr,double sigma,double a);
void linspace( int dt,double x[], double start,double end);
int rk4(double s[],double psi[],double r[],double h[],double ps[],double pr[],double k[],double l[],double m[],double n[],double o[] ,double angle,double start,double end,int dt,double sigma,double a,double dx);
double trapz(double l[],double psi[],double r[],double h[],double ps[],double pr[],double dx,int dt,double a,double sigma);
double trapzarea(double lagrange[],double s[],double r[],double psi[],double dx,int dt,double angle,double a);
bool checkrk(double psi[],double dx,int dt);
bool findmin(double psi[]);
int numofmin(double psi[]);
double findw(double s[],double psi[],double r[],double h[],double ps[],double pr[],double k[],double l[],double m[],double n[],double o[] ,double angle,double sigma,double a,double dx,double start,double end);
