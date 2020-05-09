#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "functions.h"
#include "sidefunctions.h"
# include "rk4.h"

	double angle = 80;	//wrapping degree
	double sigma = 0.5;	//surface tension
	double kappa = 50.0;	//bending rigidity
	double startx = 0.0e-9;	//start point in the r-axis
	double endx = 10.0e-9;	//end point in the r-axis
	double starty = 0.0e-9;	//start point in the z-axis
	double endy = 20.0e-9;	//end point in the z-axis
	double thickness = 2.0e-9;	//thickness of the membrane
	double hgthickness = 0.5e-9;//membrane headgroup thickness
	double pthickness = 0.5e-9;//particle outer layer thickness
	double radius = 5.0e-9; //radius of the nanoparticle
	int nx = 1000; //number of points in the r direction
	int ny = 1000; //number of points in the z direction
	int nm = 1000; //number of points for membrane
	int ncount[10]={0,0,0,0,0,0,0,0,0,0}; //point count for each region
	double qp = 1;//*;	//1.6/pi =0.50929581789407//18.735362998
	double epi = 1.0;//dielectric constant of particle inner area
	double epo = 10.0;//dielectric constant of particle outer shell 
	double ewu = 80.0;//dielectric constant of bulk above membrane
	double ewb = 80.0;//dielectric constant of bulk below membrane
	double emt = 2.0;//dielectric constant of membrane tailgroup
	double emhu = 5.0;//dielectric constant of membrane headgroup
	double emhb = 5.0;
	double F = 96354.33212; //Faraday's' constant 6.022e23*1.6e-19
	double T = 300; // Temperature
	double Kb = 1.38064852e-23; //KbT 1.38064852*10e-23 boltzmann constant* temperature
	double KbT = 4.11640356238e-21; //KbT 1.38064852*10e-23 boltzmann constant* temperature
	double RT = 2477.750204;//8.314459*298
	double ec = 1.60217662e-19;
	double h = 0.0;	//discreatization
	double hx = 0.0;	//discreatization
	double hy = 0.0;	//discreatization
	double slu = 0.0; //debye length of upper solution
	double slb = 0.0; //debye length of lower solution
	int drepeat=10;	//smoothing repeation on dielectric constant
	int qrepeat=2;	//smoothing repeation on charge density
	int prepeat=2;	//smoothing repeation on potential pot0 
	int nbkup = 2;	//ions number on upper solution
	int nbkbt = 2;	//ions number on bottom solution
	int nmemu = 1;	//ions number on membrane headgroup upper
	int nmemb = 1;	//ions number on membrane headgroup lower
	double zbkup[6]={1,-1,0,0,0,0};	//valence ions on upper solution
	double cbkup[6]={0.15,0.15,0,0,0,0};	//concentration of ions in upper solution
	double cbkupM[6]={22.9897,35.453,0,0,0,0};	//concentration of ions in upper solution
	double zbkbt[6]={1,-1,0,0,0,0};		//valence ions on bottom solution
	double cbkbt[6]={0.15,0.15,0,0,0,0};	//concentration of ions in bottom solution
	double cbkbtM[6]={22.9897,35.453,0,0,0,0};	//concentration of ions in upper solution
	double zmemu[6]={1,-1,0,0,0,0};	//valence ions on membrane headgroup
	double zmemb[6]={1,-1,0,0,0,0};	//valence ions on membrane headgroup
	double cmemu[6]={1,0.15,0,0,0,0};	//concentration of ions in upper membrane headgroup
	double cmemb[6]={1,0.15,0,0,0,0};	//concentration of ions in lower membrane headgroup
	double cmemuM[6]={22.9897,35.453,35.453,0,0,0};	//mass of ions in upper headgroup
	double cmembM[6]={35.453,22.9897,35.453,0,0,0};	//mass of ions in lower headgroup
	double pbfactor = 96485332.12;	//farady constant *( 1000 to change Liter to m^3) 
	double EnergyTot=0.0;	//total energy of the system
	double EnergyBen=0.0;	//bending energy of the membrane
	double EnergyBenwosigma=0.0;	//bending energy of the membrane
	double EnergyEle=0.0;	//electrostatic energy of the system
	double EnergyEleflat= 0.0; //electrostatic energy of flat membrane
	double memarea =0.0;	//membrane area
	FILE *errorfile;
	int id0 = 0;	//Inside particle 
	int id1 = 1;	//Particle outter layer
	int id2 = 2;	//Solution 1, upper solution
	int id3 = 3;	//Solution 2, lower solution
	int id4 = 4;	//Headgroup
	int id5 = 5;	//Tailgroup
	int id6 = 6;	//Headgroup
	int id7 = 7;	//Tailgroup
	int id8 = 8;	//Tailgroup
	int id9 = 9;	//Tailgroup
	int isMem = 1;	//Membrane is in system or not
	int isNP = 1;	//is Particle in system or not
	int isFd2d = 1;	//start FD2D functions
	int isPlot = 1;	//start gnuplot functions
	int isclean = 1; //remove unwanted files
	int isPrint = 1; //print area,epsilon,charge,efield,pot0
	int pot2d = 1; //change 3d pot to 2d 
	double memfactor = 1.0;	//factor to multiply membrane height with for shape family h=h*memfactor
	double memshift = 0.0e-9; // a factor to shift membrane up or down from the solved position
	double parshift = 0.0e-9; // a factor to shift paritcle up or down from the solved position
	clock_t times = 0;
	int status = 0;
	double jacobfact = 0.0;

/******************************************************************************/
double electricfield(double efieldr[],double efieldz[],double pot[],double temp[])
/******************************************************************************/
{
	//calculate the electrostatic energy of the system
	int i,j,k,n;
	FILE *data_unit;
	double res=0.0, iu=0,ib=0,maxnorm=-1e-20;
	for(i=1;i<ny-1;i++)
		for(j=1;j<nx-1;j++)
		{
			efieldr[j+i*nx]=-(pot[j+1+  i*nx]-pot[j-1+  i*nx])/2.0/hx;
			efieldz[j+i*nx]=-(pot[j+(i+1)*nx]-pot[j+(i-1)*nx])/2.0/hy;
		}
	for(j=0;j<nx*ny;j++)
	temp[j]=0.0;
	smoothq(5,efieldr,temp);
	smoothq(5,efieldz,temp);
	for(i=0;i<nx*ny;i++)
		if(sqrt(efieldr[i]*efieldr[i]+efieldz[i]*efieldz[i])>maxnorm)
		maxnorm=sqrt(efieldr[i]*efieldr[i]+efieldz[i]*efieldz[i]);
	data_unit = fopen ( "Efield.dat", "wt" );
	fprintf ( data_unit, "#Maxnorm=%g\n",maxnorm);
	for ( i = 0; i < ny; i++ )
	{
		for ( j = 0; j < nx; j++ )
		{
			//if(sqrt(efieldr[j+i*nx]*efieldr[j+i*nx]+efieldz[j+i*nx]*efieldz[j+i*nx])>0.1*maxnorm)
				fprintf ( data_unit, "%g\t%g\t%g\t%g\n",(double)(j*endx*1e9/(nx-1)),(double)(i*endy*1e9/(ny-1)),efieldr[j+i*nx],efieldz[j+i*nx]);
			//else 
			//	fprintf ( data_unit, "%g\t%g\t%g\t%g\n",(double)(j*endx*1e9/(nx-1)),(double)(i*endy*1e9/(ny-1)),0.0,0.0);
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );
	printf ( "\n" );
	printf ( "	Matrix (2D) printed to file '%s'\n", "Efield.dat" );
}
/******************************************************************************/
double debye()
/******************************************************************************/
{
//calculate the electrostatic energy of the system
	int i,j,k,n;
	double res=0.0, iu=0,ib=0;
	for(k=0;k<nbkup;k++)
	{
		iu+=zbkup[k]*zbkup[k]*cbkup[k];
	}
	for(k=0;k<nbkbt;k++)
	{
		ib+=zbkbt[k]*zbkbt[k]*cbkbt[k];
	}
	if(nbkup!=0)
	slu=sqrt(8.8541878128e-12*ewu*KbT/iu/F/ec/1000);
	if(nbkbt!=0)
	slb=sqrt(8.8541878128e-12*ewb*KbT/ib/F/ec/1000);
	if(nbkup!=0)
	printf("Screening length %g\t%g\n",slu,slb);
}
/******************************************************************************/
double lambda(int id,int i)
/******************************************************************************/
{
	double temp,lcube;
	double blankc=6.62607004e-34; //J.s
	double av=6.0221366516752e26; //avocadro number *1000 to get mass in KG
	temp = blankc/sqrt(2*acos(-1)*cbkbtM[i]*KbT/av);
	if(id==id2)
		temp = blankc/sqrt(2*acos(-1)*cbkupM[i]*KbT/av);
	else if(id==id3) 
		temp = blankc/sqrt(2*acos(-1)*cbkbtM[i]*KbT/av);
	else if(id==id4)
		temp = blankc/sqrt(2*acos(-1)*cmemuM[i]*KbT/av);
	else if(id==id6)
		temp = blankc/sqrt(2*acos(-1)*cmembM[i]*KbT/av);			
	lcube = temp*temp*temp;
	return lcube; //unit in m^3
}
/******************************************************************************/
double freeenergy(double pot[],double qregion[],int region[],double xvec[],double pot0[])
/******************************************************************************/
{	//double c1=0.0,c2=0.0,c3=0.0,c4=0.0,c5=0.0,c6=0.0,c7=0.0,c8=0.0;
	//double a1=0.0,a2=0.0,a3=0.0,a4=0.0,a5=0.0,a6=0.0,a7=0.0,a8=0.0;
	//double temp1=0.0,temp2=0.0;
	double term1=0.0, term2=0.0, term3=0.0, term4=0.0,sum=0.0;
	int i,j,k;
	for(i=0;i<ny;i++)
		for(j=0;j<nx;j++)
		{	
			//term1
			term1+=(f(region[j+i*nx],pot[j+i*nx],pot0[j+i*nx])+qregion[j+i*nx])*pot[j+i*nx]*xvec[j]/2;
			
			//term2
			if(region[j+i*nx]==id2)
				for(k=0;k<nbkup;k++)
					term2 += 1000*cbkup[k]*6.02214076e23*xvec[j]*KbT;
			else if(region[j+i*nx]==id3) 
				for(k=0;k<nbkbt;k++)
					term2 += 1000*cbkbt[k]*6.02214076e23*xvec[j]*KbT;
			else if(region[j+i*nx]==id4) 
				for(k=0;k<nmemu;k++)
					term2 += 1000*cmemu[k]*6.02214076e23*xvec[j]*KbT;
			else if(region[j+i*nx]==id6) 
				for(k=0;k<nmemb;k++)
					term2 += 1000*cmemb[k]*6.02214076e23*xvec[j]*KbT;
			//term3
			
			//if(region[j+i*nx]==id2 && cbkup[0] !=0.0)
			//	{
			//		c1 += (1000*cbkup[0]*exp(-zbkup[0]*F*pot[j+i*nx]/RT))*(log(6.02214076e23*lambda(region[j+i*nx],0)*1000*cbkup[0]*exp(-zbkup[0]*F*pot[j+i*nx]/RT))-1)*xvec[j];
			//		c2 += (1000*cbkup[1]*exp(-zbkup[1]*F*pot[j+i*nx]/RT))*(log(6.02214076e23*lambda(region[j+i*nx],1)*1000*cbkup[1]*exp(-zbkup[1]*F*pot[j+i*nx]/RT))-1)*xvec[j];
			//	}
			//else if(region[j+i*nx]==id3 && cbkbt[0] !=0.0) 
			//	{
			//		c3 += (1000*cbkbt[0]*exp(-zbkbt[0]*F*pot[j+i*nx]/RT))*(log(6.02214076e23*lambda(region[j+i*nx],0)*1000*cbkbt[0]*exp(-zbkbt[0]*F*pot[j+i*nx]/RT))-1)*xvec[j];
			//		c4 += (1000*cbkbt[1]*exp(-zbkbt[1]*F*pot[j+i*nx]/RT))*(log(6.02214076e23*lambda(region[j+i*nx],1)*1000*cbkbt[1]*exp(-zbkbt[1]*F*pot[j+i*nx]/RT))-1)*xvec[j];
			//	}
			//else if(region[j+i*nx]==id4 && cmemu[0] !=0.0) 
			//	{
			//		c5 += (1000*cmemu[0]*exp(-zmemu[0]*F*(pot[j+i*nx]-v0)/RT))*(log(6.02214076e23*lambda(region[j+i*nx],0)*1000*cmemu[0]*exp(-zmemu[0]*F*(pot[j+i*nx]-v0)/RT))-1)*xvec[j];
			//		c6 += (1000*cmemu[1]*exp(-zmemu[1]*F*(pot[j+i*nx]-v0)/RT))*(log(6.02214076e23*lambda(region[j+i*nx],1)*1000*cmemu[1]*exp(-zmemu[1]*F*(pot[j+i*nx]-v0)/RT))-1)*xvec[j];
			//	}
			//else if(region[j+i*nx]==id6 && cmemb[0] !=0.0) 
			//	{
			//		c7 += (1000*cmemb[0]*exp(-zmemb[0]*F*(pot[j+i*nx]-v1)/RT))*(log(6.02214076e23*lambda(region[j+i*nx],0)*1000*cmemb[0]*exp(-zmemb[0]*F*(pot[j+i*nx]-v1)/RT))-1)*xvec[j];
			//		c8 += (1000*cmemb[1]*exp(-zmemb[1]*F*(pot[j+i*nx]-v1)/RT))*(log(6.02214076e23*lambda(region[j+i*nx],1)*1000*cmemb[1]*exp(-zmemb[1]*F*(pot[j+i*nx]-v1)/RT))-1)*xvec[j];
			//	}
//


			if(region[j+i*nx]==id2 && cbkup[0] !=0.0)
				for(k=0;k<nbkup;k++)
					term3 += KbT*(6.02214076e23*1000*cbkup[k]*exp(-zbkup[k]*F*(pot[j+i*nx])/RT))*(log(6.02214076e23*lambda(region[j+i*nx],k)*1000*cbkup[k]*exp(-zbkup[k]*F*(pot[j+i*nx])/RT))-1)*xvec[j];
			else if(region[j+i*nx]==id3 && cbkbt[0] !=0.0) 
				for(k=0;k<nbkbt;k++)
					term3 += KbT*(6.02214076e23*1000*cbkbt[k]*exp(-zbkbt[k]*F*(pot[j+i*nx])/RT))*(log(6.02214076e23*lambda(region[j+i*nx],k)*1000*cbkbt[k]*exp(-zbkbt[k]*F*(pot[j+i*nx])/RT))-1)*xvec[j];
			else if(region[j+i*nx]==id4 && cmemu[0] !=0.0) 
				for(k=0;k<nmemu;k++)
					term3 += KbT*(6.02214076e23*1000*cmemu[k]*exp(-zmemu[k]*F*(pot[j+i*nx]-pot0[j+i*nx])/RT))*(log(6.02214076e23*lambda(region[j+i*nx],k)*1000*cmemu[k]*exp(-zmemu[k]*F*(pot[j+i*nx]-pot0[j+i*nx])/RT))-1)*xvec[j];
			else if(region[j+i*nx]==id6 && cmemb[0] !=0.0) 
				for(k=0;k<nmemb;k++)
					term3 += KbT*(6.02214076e23*1000*cmemb[k]*exp(-zmemb[k]*F*(pot[j+i*nx]-pot0[j+i*nx])/RT))*(log(6.02214076e23*lambda(region[j+i*nx],k)*1000*cmemb[k]*exp(-zmemb[k]*F*(pot[j+i*nx]-pot0[j+i*nx])/RT))-1)*xvec[j];
			//term4

			//if(region[j+i*nx]==id2 && cbkup[0] !=0.0)
			//		{
			//			a1 += (KbT*log(lambda(region[j+i*nx],0)*1000*cbkup[0]*6.02214076e23) + KbT * log(1000*cbkup[0]*6.02214076e23*exp(-zbkup[0]*F*pot[j+i*nx]/RT)))*1000*cbkup[0]*exp(-zbkup[0]*F*pot[j+i*nx]/RT)*xvec[j];
			//			a2 += (KbT*log(lambda(region[j+i*nx],1)*1000*cbkup[1]*6.02214076e23) + KbT * log(1000*cbkup[1]*6.02214076e23*exp(-zbkup[1]*F*pot[j+i*nx]/RT)))*1000*cbkup[1]*exp(-zbkup[1]*F*pot[j+i*nx]/RT)*xvec[j];
			//		}
			//else if(region[j+i*nx]==id3 && cbkbt[0] !=0.0)
			//		{
			//			a3 += (KbT*log(lambda(region[j+i*nx],0)*1000*cbkbt[0]*6.02214076e23) + KbT * log(1000*cbkbt[0]*6.02214076e23*exp(-zbkbt[0]*F*pot[j+i*nx]/RT)))*1000*cbkbt[0]*exp(-zbkbt[0]*F*pot[j+i*nx]/RT)*xvec[j];
			//			a4 += (KbT*log(lambda(region[j+i*nx],1)*1000*cbkbt[1]*6.02214076e23) + KbT * log(1000*cbkbt[1]*6.02214076e23*exp(-zbkbt[1]*F*pot[j+i*nx]/RT)))*1000*cbkbt[1]*exp(-zbkbt[1]*F*pot[j+i*nx]/RT)*xvec[j];
			//		}
			//else if(region[j+i*nx]==id4 && cmemu[0] !=0.0) 
			//		{
			//			a5 += (KbT*log(lambda(region[j+i*nx],0)*1000*cmemu[0]*6.02214076e23) + KbT * log(1000*cmemu[0]*6.02214076e23*exp(-zmemu[0]*F*(pot[j+i*nx]-v0)/RT)))*1000*cmemu[0]*exp(-zmemu[0]*F*(pot[j+i*nx]-v0)/RT)*xvec[j];
			//			a6 += (KbT*log(lambda(region[j+i*nx],1)*1000*cmemu[1]*6.02214076e23) + KbT * log(1000*cmemu[1]*6.02214076e23*exp(-zmemu[1]*F*(pot[j+i*nx]-v0)/RT)))*1000*cmemu[1]*exp(-zmemu[1]*F*(pot[j+i*nx]-v0)/RT)*xvec[j];
			//		}
			//else if(region[j+i*nx]==id6 && cmemb[0] !=0.0) 
			//		{
			//			a7 += (KbT*log(lambda(region[j+i*nx],0)*1000*cmemb[0]*6.02214076e23) + KbT * log(1000*cmemb[0]*6.02214076e23*exp(-zmemb[0]*F*(pot[j+i*nx]-v1)/RT)))*1000*cmemb[0]*exp(-zmemb[0]*F*(pot[j+i*nx]-v1)/RT)*xvec[j];
			//			a8 += (KbT*log(lambda(region[j+i*nx],1)*1000*cmemb[1]*6.02214076e23) + KbT * log(1000*cmemb[1]*6.02214076e23*exp(-zmemb[1]*F*(pot[j+i*nx]-v1)/RT)))*1000*cmemb[1]*exp(-zmemb[1]*F*(pot[j+i*nx]-v1)/RT)*xvec[j];
			//		}
//
//


			if(region[j+i*nx]==id2 && cbkup[0] !=0.0)
				for(k=0;k<nbkup;k++)
					term4 +=   (log(lambda(region[j+i*nx],k)*1000*cbkup[k]*6.02214076e23)+log(exp(-zbkup[k]*F*(pot[j+i*nx])/RT)))*RT*1000*cbkup[k]*exp(-zbkup[k]*F*(pot[j+i*nx])/RT)*xvec[j];
			else if(region[j+i*nx]==id3 && cbkbt[0] !=0.0) 
				for(k=0;k<nbkbt;k++)
					term4 += (log(lambda(region[j+i*nx],k)*1000*cbkbt[k]*6.02214076e23) + log(exp(-zbkbt[k]*F*(pot[j+i*nx])/RT)))*RT*1000*cbkbt[k]*exp(-zbkbt[k]*F*(pot[j+i*nx])/RT)*xvec[j];
			else if(region[j+i*nx]==id4 && cmemu[0] !=0.0) 
				for(k=0;k<nmemu;k++)
					term4 += (log(lambda(region[j+i*nx],k)*1000*cmemu[k]*6.02214076e23) + log(exp(-zmemu[k]*F*(pot[j+i*nx]-pot0[j+i*nx])/RT)))*RT*1000*cmemu[k]*exp(-zmemu[k]*F*(pot[j+i*nx]-pot0[j+i*nx])/RT)*xvec[j];
			else if(region[j+i*nx]==id6 && cmemb[0] !=0.0) 
				for(k=0;k<nmemb;k++)
					term4 += (log(lambda(region[j+i*nx],k)*1000*cmemb[k]*6.02214076e23) + log(exp(-zmemb[k]*F*(pot[j+i*nx]-pot0[j+i*nx])/RT)))*RT*1000*cmemb[k]*exp(-zmemb[k]*F*(pot[j+i*nx]-pot0[j+i*nx])/RT)*xvec[j];
		}	

		//AV #= 6.02214076e23
		term1*=2*acos(-1)*hx*hy/KbT;
		term2*=2*acos(-1)*hx*hy/KbT;
		term3*=2*acos(-1)*hx*hy/KbT;
		term4*=2*acos(-1)*hx*hy/KbT;
		//printf("term4=%g\n",term4);
		//temp1*=2*acos(-1)*hx*hy/KbT;
		//temp2*=2*acos(-1)*hx*hy/KbT;
		//term4=temp1;
		//printf("c1=%g,c2=%g,c3=%g,c4=%g,c5=%g,c6=%g,c7=%g,c8=%g,sum=%g,final=%g\n",c1,c2,c3,c4,c5,c6,c7,c8,c1+c2+c3+c4+c5+c6+c7+c8,0.0);
		//printf("a1=%g,a2=%g,a3=%g,a4=%g,a5=%g,a6=%g,a7=%g,a8=%g,sum=%g,final=%g\n",a1,a2,a3,a4,a5,a6,a7,a8,a1+a2+a3+a4+a5+a6+a7+a8,0.0);
		//c1*=2*acos(-1)*hx*hy*KbT*6.02214076e23/KbT;
		//c2*=2*acos(-1)*hx*hy*KbT*6.02214076e23/KbT;
		//c3*=2*acos(-1)*hx*hy*KbT*6.02214076e23/KbT;
		//c4*=2*acos(-1)*hx*hy*KbT*6.02214076e23/KbT;
		//c5*=2*acos(-1)*hx*hy*KbT*6.02214076e23/KbT;
		//c6*=2*acos(-1)*hx*hy*KbT*6.02214076e23/KbT;
		//c7*=2*acos(-1)*hx*hy*KbT*6.02214076e23/KbT;
		//c8*=2*acos(-1)*hx*hy*KbT*6.02214076e23/KbT;
//
		//a1*=2*acos(-1)*hx*hy*6.02214076e23/KbT;
		//a2*=2*acos(-1)*hx*hy*6.02214076e23/KbT;
		//a3*=2*acos(-1)*hx*hy*6.02214076e23/KbT;
		//a4*=2*acos(-1)*hx*hy*6.02214076e23/KbT;
		//a5*=2*acos(-1)*hx*hy*6.02214076e23/KbT;
		//a6*=2*acos(-1)*hx*hy*6.02214076e23/KbT;
		//a7*=2*acos(-1)*hx*hy*6.02214076e23/KbT;
		//a8*=2*acos(-1)*hx*hy*6.02214076e23/KbT;
		sum=term1+term2+term3-term4;
		printf("Term1=%g\nTerm2=%g\nTerm3=%g\nTerm4=%g\nsum3terms=%g\nsumofterms=%g\n",term1,term2,term3,term4,term2+term3-term4,term1+term2+term3-term4);
		//printf("c1=%g,c2=%g,c3=%g,c4=%g,c5=%g,c6=%g,c7=%g,c8=%g,sum=%g,final=%g\n",c1,c2,c3,c4,c5,c6,c7,c8,c1+c2+c3+c4+c5+c6+c7+c8,0.0);
		//printf("a1=%g,a2=%g,a3=%g,a4=%g,a5=%g,a6=%g,a7=%g,a8=%g,sum=%g,final=%g\n",a1,a2,a3,a4,a5,a6,a7,a8,a1+a2+a3+a4+a5+a6+a7+a8,0.0);
		//printf("temp1=%g\t temp2=%g \t sum=%g\n",temp1,temp2,temp1+temp2);

		return sum;
	
}
/******************************************************************************/

void addrhs(double acc[],int icc[],int ccc[],double rhs[],double iniguess[],int count,int n,int region[],double qregion[],double pot0[])

/******************************************************************************/
{
	//create right hand side of PBE , b=-{A*V+h*h*F(v)}
	int i,j,k;
	double sum=0.0;
	//cc_v ( n, count, icc, ccc, acc, iniguess, rhs); //compute b=A*V
	//for(i=0;i<nx*ny;i++)
	//	sum+=rhs[i]*iniguess[i];
	for(i=0;i<ny*nx;i++)	//compute b=h*h*F(v)
			{
			rhs[i]+=f(region[i],iniguess[i],pot0[i]);
			rhs[i]+=qregion[i];
			}
	for(i=0;i<ny*nx;i++)	//b=-b
			rhs[i]=-rhs[i];
	
	
	//printf("Definite %g\n",sum);
	return;
}
/******************************************************************************/

void jacobian ( int n, double val[],int row[],int col[] ,int count,double iniguess[],int region[],double pot0[])

/******************************************************************************/
{
	//Compute the Jacobian matrix by adding +h*h*F`(v) to A
			int j;
			for(j=0;j<count;j++)
			if(row[j]==col[j])
			{val[j]+=ff(region[row[j]],iniguess[row[j]],pot0[row[j]]) * (1.0+jacobfact);}

			jacobfact=jacobfact*jacobfact;

			//{val[j]+=hx*hy*ff(region[row[j]],iniguess[row[j]]);}
}
/******************************************************************************/

void smooth ( int rep,double dregion [],double temp[])

/******************************************************************************/
{
	//Use geometrical averaging to smooth the dielectric constant on the system with rep repeations
	int i,j,k,n;
	n=nx*ny;
	for(i=0;i<ny;i++)
		for(j=0;j<nx;j++)
		temp[j+i]=0.0;

	for(k=0;k<rep;k++)
	{
		for(i=1;i<ny-1;i++)
			for(j=1;j<nx-1;j++)
				temp[j+i*nx]=pow((dregion[j+i*nx-nx]*dregion[j+i*nx-1]*dregion[j+i*nx]*dregion[j+i*nx+1]*dregion[j+i*nx+nx]),0.2);
		j=0;
		for(i=1;i<ny-1;i++)
				temp[j+i*nx]=pow((dregion[j+i*nx-nx]*dregion[j+i*nx]*dregion[j+i*nx+1]*dregion[j+i*nx+nx]),0.25);
		j=nx-1;
		for(i=1;i<ny-1;i++)
				temp[j+i*nx]=pow((dregion[j+i*nx-nx]*dregion[j+i*nx]*dregion[j+i*nx-1]*dregion[j+i*nx+nx]),0.25);
		for(i=1;i<ny-1;i++)
			for(j=0;j<nx;j++)
				dregion[j+i*nx]=temp[j+i*nx];
	}
}
/******************************************************************************/

void smoothq ( int rep,double dregion [],double temp[])

/******************************************************************************/
{
	//Use mathematical averaging to smooth the dielectric constant on the system with rep repeations
	int i,j,k,n;
	n=nx*ny;
	for(i=0;i<ny;i++)
		for(j=0;j<nx;j++)
		temp[j+i]=0.0;
	for(k=0;k<rep;k++)
	{
		for(i=1;i<ny-1;i++)
			for(j=1;j<nx-1;j++)
				temp[j+i*nx]=(dregion[j+i*nx-nx]+dregion[j+i*nx-1]+dregion[j+i*nx]+dregion[j+i*nx+1]+dregion[j+i*nx+nx])/5.0;
		j=0;
		for(i=1;i<ny-1;i++)
				temp[j+i*nx]=(dregion[j+i*nx-nx]+dregion[j+i*nx]+dregion[j+i*nx+1]+dregion[j+i*nx+nx])/4.0;
		for(i=1;i<ny-1;i++)
			for(j=0;j<nx-1;j++)
				dregion[j+i*nx]=temp[j+i*nx];
	}
}
/******************************************************************************/
void init ( double x[], double y[], int n, 
	double a[],int b[],int c[], double rhs[] ,double height [] ,double dregion[],int region[])

/******************************************************************************/
{
	//Triple format a 'value',b 'row index',c 'column index'
	//Initialize the system in sparse triplet matrix format
	double dc0,dcc,dce,dcn,dcs,dcw;
	double dr;
	double dz;
	int ic,jc;
	int kc,ke,kn,ks,kw;
	int k=0;
	int i,j;
	dr = hx;
	dz = hy;
	for ( i = 0; i < ny ; i++ )
	{
		for ( j = 0; j < nx ; j++ )
		{
//Bottom boundary.
			if(i == 0)
			{
				kc = i * nx + j;
				a[k] = 1.0;
				b[k] = kc;
				c[k] = kc;
				rhs[kc] = 0.0;
				k++;
			}
//Left boundary
			else if((j == 0 && i!=0 && i!=ny-1) )
			{	
				kc = i * nx + j;
				kn = kc + nx;
				ks = kc - nx;
				ke = kc + 1;
				dcc = dregion[kc];
				dce = dregion[ke];
				dcn = dregion[kn];
				dcs = dregion[ks];
				a[k] = (-0.25*(dcn-dcs)+ dcc)/(hy*hy);
				b[k] = kc;
				c[k] = ks;
				k++;
				a[k] = -4.0*dcc/(hx*hx)-2.0*dcc/(hy*hy);
				b[k] = kc;
				c[k] = kc;
				k++;
				a[k] = 4.0*dcc/(hx*hx);
				b[k] = kc;
				c[k] = ke;
				k++;
				a[k] = (0.25*(dcn-dcs)+dcc)/(hy*hy);
				b[k] = kc;
				c[k] = kn;
				k++;
				rhs[kc] = 0.0;
			}
//Right boundary.
			else if(j == nx - 1 && i!=0 && i!=ny-1)
			{
				/*kc = i * nx + j;
				a[k] = 1.0;
				b[k] = kc;
				c[k] = kc;
				k++;
				if(region[kc]==4 ||region[kc]==5)
				a[k] = -1.0;
				else
				a[k] = 0.0;
				b[k] = kc;
				c[k] = kc-1;
				k++;
				rhs[kc] = 0.0;*/
				kc = i * nx + j;
				kn = kc + nx;
				ks = kc - nx;
				kw = kc - 1;
				dcc = dregion[kc];
				dcw = dregion[kw];
				dcn = dregion[kn];
				dcs = dregion[ks];
				a[k] = (-0.25*(dcn-dcs)+ dcc)/(hy*hy);
				b[k] = kc;
				c[k] = ks;
				k++;
				a[k] = 2.0*dcc/(hx*hx);
				b[k] = kc;
				c[k] = kw;
				k++;
				a[k] = -2.0*dcc/(hx*hx)-2.0*dcc/(hy*hy);
				b[k] = kc;
				c[k] = kc;
				k++;
				a[k] = (0.25*(dcn-dcs)+dcc)/(hy*hy);
				b[k] = kc;
				c[k] = kn;
				k++;
				rhs[kc] = 0.0;	
			}
//Upper boundary.

			else if(i == ny - 1)
			{
				kc = i * nx + j;
				a[k] = 1.0;
				b[k] = kc;
				c[k] = kc;
				k++;
				rhs[kc] = 0.0;
			}
//Inner points
			else {
			kc = i * nx + j;
			ke = kc + 1;
			kw = kc - 1;
			kn = kc + nx;
			ks = kc - nx;
			dcc = dregion[kc];
			dce = dregion[ke];
			dcw = dregion[kw];
			dcn = dregion[kn];
			dcs = dregion[ks];
			a[k] = (-0.25*(dcn-dcs) +dcc)/(hy*hy);	//south
			b[k] = kc;
			c[k] = ks;
			k++;
			a[k] = (-0.25*(dce-dcw) +dcc*(-0.5*dr/x[j+1] +1))/(hx*hx);	//west
			b[k] = kc;
			c[k] = kw;
			k++;
			a[k] = -2.0*dcc/(hx*hx)  -2.0*dcc/(hy*hy);	//center
			b[k] = kc;
			c[k] = kc;
			k++;
			a[k] = ( 0.25*(dce-dcw) +dcc*( 0.5*dr/x[j+1] +1))/(hx*hx);	//east
			b[k] = kc;
			c[k] = ke;
			k++;
			a[k] = ( 0.25*(dcn-dcs) +dcc)/(hy*hy);	//north
			b[k] = kc;
			c[k] = kn;
			k++;
			rhs[kc] = 0.0;//f ( x[j], y[i] );
			}
		}	
	}
	return;
}
/******************************************************************************/
double d ( int nr)
/******************************************************************************/
//D return the dielectric coefficient.
{
	double value;
	if(nr==-1)
	{printf("wrong location\n");
	exit(0);}
	if(nr==id0)//Inside the particle
		value=epi;
	if(nr==id1)//particle outer layer
		value=epo;
	if(nr==id2)//Solution above membrane
		value=ewu;
	if(nr==id3)//Solution under the membrane
		value=ewb;
	if(nr==id4)//upper Head group
		value=emhu;
	if(nr==id5)//Tail group
		value=emt;
	if(nr==id6)//bottom headgroup
		value=emhb;
	return value;
}
/******************************************************************************/

double f ( int nr,double v,double v0)

/******************************************************************************/
//F evaluates the electric source term.
{
	double value=0.0;
	int i;
	if(nr==-1)
	{printf("wrong location\n");
	exit(0);}
	if(nr==id0)//Inside the particle 
		value=0.0;
	if(nr==id1)//particle outer layer, can be found in qregion matrix
		value=0.0;
	if(nr==id2)//Solution above membrane
		{
		for(i=0;i<nbkup;i++)
			value+=pbfactor*zbkup[i]*cbkup[i]*exp(-zbkup[i]*F*(v)/RT);
		}
	if(nr==id3)//Solution under the membrane
		{
		for(i=0;i<nbkbt;i++)
			value+=pbfactor*zbkbt[i]*cbkbt[i]*exp(-zbkbt[i]*F*(v)/RT);
		}
	if(nr==id4)//Head group
		{
		for(i=0;i<nmemu;i++)
			value+=pbfactor*zmemu[i]*cmemu[i]*exp(-zmemu[i]*F*(v-v0)/RT);
		}
	if(nr==id5)//Tail group
		value=0.0;
	if(nr==id6)//Head group
		{
		for(i=0;i<nmemb;i++)
			value+=pbfactor*zmemb[i]*cmemb[i]*exp(-zmemb[i]*F*(v-v0)/RT);
		}
	return value;
}
/******************************************************************************/

double ff ( int nr,double v,double v0)

/******************************************************************************/

//FF evaluates the derivative of electric source term.

{
	double value;
	int i;
	if(nr==-1)
	{printf("wrong location\n");
	exit(0);}
	value=0.0;
	if(nr==id0)//Inside the particle
		value=0.0;
	if(nr==id1)//particle outer layer
		value=0.0;
	if(nr==id2)//Solution above membrane
		{
		for(i=0;i<nbkup;i++)
			value+=-pbfactor*zbkup[i]*zbkup[i]*cbkup[i]*F*exp(-zbkup[i]*F*(v)/RT)/RT;
		}
	if(nr==id3)//Solution under the membrane
		{
		for(i=0;i<nbkbt;i++)
			value+=-pbfactor*zbkbt[i]*zbkbt[i]*cbkbt[i]*F*exp(-zbkbt[i]*F*(v)/RT)/RT;
		}
	if(nr==id4)//Head group
		{
		for(i=0;i<nmemu;i++)
			value+=-pbfactor*zmemu[i]*zmemu[i]*cmemu[i]*F*exp(-zmemu[i]*F*(v-v0)/RT)/RT;
		}
	if(nr==id5)//Tail group
		value=0.0;
	if(nr==id6)//Head group
		{
		for(i=0;i<nmemb;i++)
			value+=-pbfactor*zmemb[i]*zmemb[i]*cmemb[i]*F*exp(-zmemb[i]*F*(v-v0)/RT)/RT;
		}
	return value;
}
/******************************************************************************/

void r8sp_mv ( int m, int n, int nz_num, int row[], int col[], 
  double a[], double x[] , double b[])

/******************************************************************************/
/*
  Purpose:

    R8SP_MV multiplies a R8SP matrix times a vector.

  Discussion:

    The R8SP storage format stores the row, column and value of each nonzero

    It is possible that a pair of indices (I,J) may occur more than
    once.  Presumably, in this case, the intent is that the actual value
    of A(I,J) is the sum of all such entries.  This is not a good thing
    to do, but I seem to have come across this in MATLAB.

    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 February 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.

    Input, int NZ_NUM, the number of nonzero elements in the matrix.

    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
    of the nonzero elements.

    Input, double A[NZ_NUM], the nonzero elements of the matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double R8SP_MV[M], the product vector A*X.
*/
{
  int i;
  int j;
  int k;

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( k = 0; k < nz_num; k++ )
  {
    i = row[k];
    j = col[k];
    b[i] = b[i] + a[k] * x[j];
  }
}
