#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include "functions.h"
#include "sidefunctions.h"
#include "rk4.h"


	double w = -3;
	double pos = +1.0;
	double neg = -7.0;
	double mid = -3.0;
	int ilimit=0;
	int imin=0;
	int dt = 100000;

void memheight(double radial[], double height[])
{
	double alpha = angle; //wrapping degree 
	double a = radius+thickness/2.0 ;	//radius of particle+ thickness of membrane
	double sig = sigma;	//surface tension
	dt =nm;	//number of discreatization points
	double dx=0.0;
	double start = startx;	//start of the arc length, start of r-axis
	double end = endx;	//end of the arc length, end of r-axis
	double area=0;


	FILE *file;
	dx=(end-start)/dt;
	double *k,*l,*m,*n,*o;
	double *s,*psi,*r,*h,*ps,*pr,*lagrange,integral=0.0,integralwo=0.0;
	int i=0.0;
	k = ( double * ) malloc ( 5 * sizeof ( double ) );
	l = ( double * ) malloc ( 5 * sizeof ( double ) );
	m = ( double * ) malloc ( 5 * sizeof ( double ) );
	n = ( double * ) malloc ( 5 * sizeof ( double ) );
	o = ( double * ) malloc ( 5 * sizeof ( double ) );
	s   = ( double * ) malloc ( dt * sizeof ( double ) );
	psi = ( double * ) malloc ( dt * sizeof ( double ) );
	r   = ( double * ) malloc ( dt * sizeof ( double ) );
	h   = ( double * ) malloc ( dt * sizeof ( double ) );
	ps  = ( double * ) malloc ( dt * sizeof ( double ) );
	pr  = ( double * ) malloc ( dt * sizeof ( double ) );
	lagrange  = ( double * ) malloc ( dt * sizeof ( double ) );
	for(i=0;i<5;i++)
	{k[i]=0.0;l[i]=0.0;m[i]=0.0;n[i]=0.0;o[i]=0.0;}
	for(i=0;i<dt;i++)
	{s[i]=0.0;psi[i]=0.0;r[i]=0.0;h[i]=0.0;ps[i]=0.0;pr[i]=0.0;lagrange[i]=0.0;}
	//loop to stretch the membrane to have a radius same as box size by changing the arc length
	i=0;
	while(1)
	{
	i++;
	if(i>1000)
	{
	fprintf(errorfile,"System exited because cant solve shape equations\n");
	break;
	}
	w=findw(s,psi,r,h,ps,pr,k,l,m,n,o,alpha,sig,a,dx,start,end);
	if(r[nm-1]<endx)
	end+=endx-r[nm-1];
	else if(r[nm-1]>endx)
	end-=r[nm-1]-endx;
	if(fabs(r[nm-1]-endx)*1e9<5e-5)
	break;
	dx=(end-start)/dt;
	//printf("end\t%g\t%g\t%g\tW\t%g\n",end*1e9,r[dt-1]*1e9,(r[nm-1]-endx)*1e9,w);
	}
	//new set of family curves, by a factor
	double temp;
	temp=h[0];
	for(i=0;i<nm;i++)
	{
		h[i] -= temp;
	}
	//printf("start of height %g %g\n",h[500],height[500]);
	for(i=0;i<nm;i++)
	{
		h[i] *= memfactor;
	}
	//printf("start of height %g %g\n",h[500],height[500]);
	for(i=0;i<nm;i++)
	{
		h[i] += temp;
	}

	
	//bending energy with sigma=sig
	//printf("Sig %g\n",sig);
	integral=trapz(lagrange,psi,r,h,ps,pr,dx,dt,a,sig);
	//printf("%0.2g\t%0.9g\t%0.9g\t%0.9g\n",alpha,integral,w,r[dt-1]*1e9);
	
	//bending energy with sigma=0
	sig=0.0;
	//printf("Sig %g\n",sig);
	integralwo=trapz(lagrange,psi,r,h,ps,pr,dx,dt,a,sig);
	//printf("%0.2g\t%0.9g\t%0.9g\t%0.9g\n",alpha,integralwo,w,r[dt-1]*1e9);
	area=trapzarea(lagrange,s,r,psi,dx,dt,angle,a);
	memarea=area;
	EnergyBen=integral*kappa*acos(-1);	//original energy is divided by kappa*pi
	EnergyBenwosigma=integralwo*kappa*acos(-1);	
	
	

	
	
	
	//Writting the membrane shape to file
	file = fopen ( "ShapeEquation.dat", "wt" );
	fprintf(file,"#1:s\t2:psi(s)\t3:r(s)\t4:h(s)\t5:P_psi(s)\t6:P_r(s)\t\t\t1,2,3 are in nm\n");
	fprintf(file,"#Shpae equations solution for\n");
	fprintf(file,"#Sigma=%g\tKappa=%g\tArc length=%g and Radius=%g, area=%.14g\n",sigma,kappa,end*1e9,r[nm-1]*1e9,area);
	fprintf(file,"#Wrapping degree=%g\tParticle size=%g\tAdhesion strength=%g\n",angle,radius*1e9,w);
	fprintf(file,"#Energy with sigma=%g\tEnergy w/o sigma=%g\n",integral,integralwo);
	for ( i = 0; i < dt ; i++ )//output
	fprintf(file,"%g\t%g\t%g\t%g\t%g\t%g\n",s[i]*1e9,psi[i],r[i]*1e9,h[i]*1e9,ps[i],pr[i]);
	fclose ( file );
	printf ("	Shape equations solved to file ShapeEquation.dat\n");
	for(i=0;i<nm;i++)
	{
		radial[i] = r[i];
		if(isNP == 1 || isNP ==2)
			height[i] = h[i]+parshift  	  +memshift;
		else
			height[i] = h[i]+parshift + a +memshift;
	}
	if (isNP!=0 && isMem!=0 && angle <1)	
	if(((height[0]-thickness/2.0)<(parshift+radius)&&(height[0]-thickness/2.0)>(parshift-radius))||((height[0]+thickness/2.0)<(parshift+radius)&&(height[0]+thickness/2.0)>(parshift-radius)) )
	{
		fprintf(errorfile,"System exited,Membrane and particle overlap because of incorrect membrane and particle shift\n");exit(0);
	}

	//test();
	
	free(k);  
	free(l);  
	free(m);  
	free(n); 
	free(o); 
	free(s);
	free(psi);  
	free(r);    
	free(h);    
	free(ps); 
	free(pr);  
	free(lagrange);  
	return;
}
double findw(double s[],double psi[],double r[],double h[],double ps[],double pr[],double k[],double l[],double m[],double n[],double o[] ,double alpha,double sig, double a,double dx,double start,double end)
{
	//iterative function to find the value of adhesion strength W which make the membrane flat at end. w=1-sqrt(W)
	int flag,i,j,count=0;
	double temp,radian;
	if(alpha<90)
		{neg=-7.0;pos=1.0;mid=w;}
	else
		{neg=1.0;pos=-7.0;mid=w;}
	radian = alpha*acos(-1)/180.0;
	w=mid;
	temp=DBL_MAX;
	j=0;
	while(1)
	{
	j++;
	count++;
	flag=rk4(s,psi,r,h,ps,pr,k,l,m,n,o,radian,start,end,dt,sig,a,dx);
	//printf("%g\t%g\t%g\t%g\t%d\t%d\n",pos,mid,neg,temp,imin,ilimit);
	findmin(psi);
	if(flag==1)
	{
	if(psi[imin]<0)
	{neg=mid;mid= (pos+neg)/2.0;w=mid;}
	else
	{pos=mid;mid= (pos+neg)/2.0;w=mid;}
	}
	else if(flag==0)
	{
	if(psi[imin]>0)
	{pos=mid;mid= (pos+neg)/2.0;w=mid;}
	else
	{neg=mid;mid= (pos+neg)/2.0;w=mid;}	
	}
	else if(flag==2)
	{
		if(psi[imin]>0)
		{
		pos=mid;mid= (pos+neg)/2.0;w=mid;
		}
		else
		{
		neg=mid;mid= (pos+neg)/2.0;w=mid;
		}

	}
	i=numofmin(psi);
	if(fabs(temp-mid)<1e-7 && i==1)
	break;
	else
	temp=mid;
	if(j>50)
	if(alpha<90)
		{neg=-7.0;pos=1.0;j=0;temp=DBL_MAX;}
	else
		{neg=1.0;pos=-7.0;j=0;temp=DBL_MAX;}
	if(count>1000)
	{fprintf(errorfile,"System exited, couldnt solve shape equations because it couldnt find correct value of adhesion strength\n");exit(0);}
	}
	return mid;
}
bool findmin(double psi[])
{
	//find minima location in the array
	double max=DBL_MAX;
	int i,j;//imin;
	for(i=1;i<ilimit;i++)
		if(psi[i]<max)
		{max=psi[i];imin=i;}
	if(imin<dt-1)
	return false;
	else 
	return true;
}
int numofmin(double psi[])
{
	//find number of minima in the array
	int i,j=0;
	for(i=1;i<ilimit-1;i++)
	if(psi[i]<psi[i-1] &&psi[i]<psi[i+1])
	j++;
	if(psi[dt-1]<psi[dt-2])
	j++;
	if(psi[0]<psi[1])
	j++;
	return j;
}
double trapz(double lagrange[],double psi[],double r[],double h[],double ps[],double pr[],double dx,int dt,double a,double sig)
{
	//use Trapezoid rule integral to calculate integral of L
	double integral=0.0;
	int i;
	for ( i = 0; i < dt ; i++ )
	lagrange[i]=ps[i]*ps[i]*0.25/r[i] +2*sig*r[i]*(1-cos(psi[i]))/a/a;
	integral=0.0;
	for( i = 1; i<dt ; i++)
	integral+= (lagrange[i-1]+lagrange[i])*dx/2;
	return integral;
}
double trapzarea(double lagrange[],double s[],double r[],double psi[],double dx,int dt,double angle,double a)
{
	//use Trapezoid rule integral to calculate membrane area
	double integral=0.0;
	int i;
	for ( i = 0; i < dt ; i++ )
	lagrange[i]=r[i];
	integral=0.0;
	for( i = 1; i<dt ; i++)
	integral+= (lagrange[i-1]*sqrt(1+tan(psi[i-1])*tan(psi[i-1]))+lagrange[i]*sqrt(1+tan(psi[i])*tan(psi[i])))*dx/2;
	return integral*2.0*acos(-1);
}
int rk4(double s[],double psi[],double r[],double h[],double ps[],double pr[],double k[],double l[],double m[],double n[],double o[] ,double radian,double start,double end,int dt,double sig,double a,double dx)
{
	//using Runge-Kutta method to solve the shape equations 
	int i,j;
	linspace(dt,s,start,end);
	psi[0]=radian;
	 r[0] =a*sin(radian);
	 h[0] =-a*cos(radian);
	ps[0] =2.0*sin(radian)*(w+1);
	pr[0] =tan(radian)*(1+2.0*sig*(1-cos(radian))-w*w)/a;
	k[0]=psi[0];
	l[0]=  r[0];
	m[0]=  h[0];
	n[0]= ps[0];
	o[0]= pr[0];
	for(i=0;i<dt-1;i++)
	{
	k[1]=dx*f1(s[i],psi[i],r[i],h[i],ps[i],pr[i],sig,a);
	l[1]=dx*f2(s[i],psi[i],r[i],h[i],ps[i],pr[i],sig,a);
	m[1]=dx*f3(s[i],psi[i],r[i],h[i],ps[i],pr[i],sig,a);
	n[1]=dx*f4(s[i],psi[i],r[i],h[i],ps[i],pr[i],sig,a);
	o[1]=dx*f5(s[i],psi[i],r[i],h[i],ps[i],pr[i],sig,a);
	
	k[2]=dx*f1(s[i]+dx/2.0,psi[i]+k[1]/2.0,r[i]+l[1]/2.0,h[i]+m[1]/2.0,ps[i]+n[1]/2.0,pr[i]+o[1]/2.0,sig,a);
	l[2]=dx*f2(s[i]+dx/2.0,psi[i]+k[1]/2.0,r[i]+l[1]/2.0,h[i]+m[1]/2.0,ps[i]+n[1]/2.0,pr[i]+o[1]/2.0,sig,a);
	m[2]=dx*f3(s[i]+dx/2.0,psi[i]+k[1]/2.0,r[i]+l[1]/2.0,h[i]+m[1]/2.0,ps[i]+n[1]/2.0,pr[i]+o[1]/2.0,sig,a);
	n[2]=dx*f4(s[i]+dx/2.0,psi[i]+k[1]/2.0,r[i]+l[1]/2.0,h[i]+m[1]/2.0,ps[i]+n[1]/2.0,pr[i]+o[1]/2.0,sig,a);
	o[2]=dx*f5(s[i]+dx/2.0,psi[i]+k[1]/2.0,r[i]+l[1]/2.0,h[i]+m[1]/2.0,ps[i]+n[1]/2.0,pr[i]+o[1]/2.0,sig,a);
	
	k[3]=dx*f1(s[i]+dx/2.0,psi[i]+k[2]/2.0,r[i]+l[2]/2.0,h[i]+m[2]/2.0,ps[i]+n[2]/2.0,pr[i]+o[2]/2.0,sig,a);
	l[3]=dx*f2(s[i]+dx/2.0,psi[i]+k[2]/2.0,r[i]+l[2]/2.0,h[i]+m[2]/2.0,ps[i]+n[2]/2.0,pr[i]+o[2]/2.0,sig,a);
	m[3]=dx*f3(s[i]+dx/2.0,psi[i]+k[2]/2.0,r[i]+l[2]/2.0,h[i]+m[2]/2.0,ps[i]+n[2]/2.0,pr[i]+o[2]/2.0,sig,a);
	n[3]=dx*f4(s[i]+dx/2.0,psi[i]+k[2]/2.0,r[i]+l[2]/2.0,h[i]+m[2]/2.0,ps[i]+n[2]/2.0,pr[i]+o[2]/2.0,sig,a);
	o[3]=dx*f5(s[i]+dx/2.0,psi[i]+k[2]/2.0,r[i]+l[2]/2.0,h[i]+m[2]/2.0,ps[i]+n[2]/2.0,pr[i]+o[2]/2.0,sig,a);
	
	k[4]=dx*f1(s[i]+dx,psi[i]+k[3],r[i]+l[3],h[i]+m[3],ps[i]+n[3],pr[i]+o[3],sig,a);
	l[4]=dx*f2(s[i]+dx,psi[i]+k[3],r[i]+l[3],h[i]+m[3],ps[i]+n[3],pr[i]+o[3],sig,a);
	m[4]=dx*f3(s[i]+dx,psi[i]+k[3],r[i]+l[3],h[i]+m[3],ps[i]+n[3],pr[i]+o[3],sig,a);
	n[4]=dx*f4(s[i]+dx,psi[i]+k[3],r[i]+l[3],h[i]+m[3],ps[i]+n[3],pr[i]+o[3],sig,a);
	o[4]=dx*f5(s[i]+dx,psi[i]+k[3],r[i]+l[3],h[i]+m[3],ps[i]+n[3],pr[i]+o[3],sig,a);
	
	psi[i+1]=psi[i]+(k[1]+2.0*k[2]+2.0*k[3]+k[4])/6.0;
	r[i+1]  = r[i] +(l[1]+2.0*l[2]+2.0*l[3]+l[4])/6.0;
	h[i+1]  = h[i] +(m[1]+2.0*m[2]+2.0*m[3]+m[4])/6.0;
	ps[i+1] =ps[i] +(n[1]+2.0*n[2]+2.0*n[3]+n[4])/6.0;
	pr[i+1] =pr[i] +(o[1]+2.0*o[2]+2.0*o[3]+o[4])/6.0;
	if(psi[i+1]>acos(-1)||psi[i+1]<-acos(-1))

	{//printf("Reach a psi[%d]=%g\n",i,psi[i]);psi[i+1]=0.0;ilimit=i-2;
	if(psi[ilimit-1]>0)
	return 0;
	else
	return 1;
	}
	ilimit=i;
	}
	return 2;
}
bool checkrk(double psi[],double dx,int dt)
{
	//check that last point is very close to zero "flat membrane"
	int i,j,k;
	double err=1e-6;
	double nerr=-1e-6;
	double temp=0.0;
	temp = (psi[dt-1]-psi[dt-2])/dx;
	if(fabs(temp)<err)
		if( nerr<psi[dt-1] && psi[dt-1]<err )
			return false;
	return true;
}
/******************************************************************************/
double f1(double s,double psi,double r,double h,double ps,double pr,double sig,double a)	//1st equation: psi.dot =P_psi/2r -sin(psi)/r
{
	return ps/2/r -sin(psi)/r;
}
double f2(double s,double psi,double r,double h,double ps,double pr,double sig,double a)	//2nd equation: r.dot =cos psi
{
	return cos(psi);
}
double f3(double s,double psi,double r,double h,double ps,double pr,double sig,double a)	//3rd equation: h.dot =sin psi
{
	return sin(psi);
}
double f4(double s,double psi,double r,double h,double ps,double pr,double sig,double a)	//4th equation: P_psi.dot = (P_psi/r -P_h)cos psi +(2*sigma*r/a^2 +P_r)sin psi ,  P_h=0 
{
	return (ps/r)*cos(psi)+(2.0*sig*r/a/a +pr)*sin(psi);
}
double f5(double s,double psi,double r,double h,double ps,double pr,double sig,double a)	//5th equation: P_r.dot =P_psi/2r -sin(psi)/r
{
	return ps*(ps/4.0/r -sin(psi)/r)/r +2.0*sig*(1-cos(psi))/a/a;
}
void linspace( int dt,double x [],double start ,double end)	//discreatize s over the region between start-->end and number of points dt
{
	int i,n;
	n=dt;
	if ( n == 1 )
	{
		x[0] = ( start + end ) / 2.0;
	}
	else
	{
		for ( i = 0; i < n; i++ )
		{
			x[i] = ( ( double ) ( n -1- i ) * start + ( double ) ( i ) * end ) / ( double ) ( n-1  );
		}
	}
	return;
}
