#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include "fd2d.h"
#include "functions.h"
#include "sidefunctions.h"
#include "rk4.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/IterativeLinearSolvers>
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::ConjugateGradient;
 
using std::cout;
using std::endl;

void fd2d (double radial[],double height[])

/******************************************************************************/

{
	FILE *data_unit;
	int i=0;
	int j=0;
	int n=0;
	int count=0;
	double *val;	//triple format matrix store value
	int *row;	//triple format matrix store row index
	int *col;	//triple format matrix store column index
	double *rhs;
	double *xvec;
	double *yvec;
	double *iniguess, *efieldr,*efieldz, *pot0;
	int *region;
	double *dregion,*temp;
	double *qregion;
	double norm=0.0;
	double tenergy=0.0;
	//count =nx*nx*5-12*nx+8;
	//count =nx*nx*5-14*nx+12;
	//count =nx*nx*5-13*nx+10;
	//count = 5*nx*ny-5*ny-8*nx+10+ny-2;
	count = 5*(nx-2)*(ny-2) +4*(ny-2)+4*(ny-2)+1*(nx)+1*(nx);
	printf("Number of NNZ = %d\n",count);
	n = nx * ny;
	xvec = r8vec_linspace_new ( nx, startx, endx );
	yvec = r8vec_linspace_new ( ny, starty, endy );
	val  = ( double * ) malloc ( count * sizeof ( double ) );
	row = ( int * ) malloc ( count * sizeof ( int ) );
	col = ( int * ) malloc ( count * sizeof ( int ) );
	rhs = ( double * ) malloc ( n * sizeof ( double ) );
	iniguess = ( double * ) malloc ( n * sizeof ( double ) );
	efieldr = ( double * ) malloc ( n * sizeof ( double ) );
	efieldz = ( double * ) malloc ( n * sizeof ( double ) );
	pot0 = ( double * ) malloc ( n * sizeof ( double ) );
	region  = (int *) malloc (n *sizeof (int));
	dregion = (double *) malloc (n *sizeof (double));
	temp = (double *) malloc (n *sizeof (double));
	qregion = (double *) malloc (n *sizeof (double));
	
	
	//Initialize initial guess with value
	for(j=0;j<n;j++)
	{efieldr[j]=0.0;efieldz[j]=0.0;pot0[j]=0.0;iniguess[j]=0.0;rhs[i]=0.0;region[i]=0;dregion[i]=0.0;temp[i]=0.0;qregion[j]=0.0;}
	
	
	//create region matrix which contain the ID of each grid point
	regioncountmaple( region,radial,height, xvec, yvec);
	if(isPrint == 1)
	print2dint(nx,region,"area.dat");
	if(isPlot!=0)
	status = system("gnuplot area.gp");
	if(isclean == 1)
	    status = system("rm area.gp area.dat");
	//create dregion matrix which contain the dielectric value of each grid point
	for(j=0;j<n;j++)
	     if(region[j]==id0) dregion[j]=ewu *8.8541878128e-12; 
	else if(region[j]==id1) dregion[j]=ewu *8.8541878128e-12;	
	else if(region[j]==id2) dregion[j]=ewu *8.8541878128e-12; 
	else if(region[j]==id3) dregion[j]=ewb *8.8541878128e-12; 
	else if(region[j]==id4) dregion[j]=emhu*8.8541878128e-12; 
	else if(region[j]==id5) dregion[j]=emt *8.8541878128e-12;
	else if(region[j]==id6) dregion[j]=emhb*8.8541878128e-12;
	//smoothing dielectric constant with neighboring points
	smooth(drepeat,dregion,temp);
	for(j=0;j<n;j++)
	     if(region[j]==id0) dregion[j]=epi *8.8541878128e-12; 
	else if(region[j]==id1) dregion[j]=epo *8.8541878128e-12;
	smooth(10,dregion,temp);
	if(isPrint == 1)
	print2d(nx,dregion,"Epsilon.dat",1e9,1.0/8.854e-12);
	if(isPlot!=0)
	status = system("gnuplot epsilon.gp");
	if(isclean == 1)
	    status = system("rm epsilon.gp Epsilon.dat");
	//create qregion matrix which contain the charge value of each grid point
	qp=qp*ec;	//total charge * elementary charge
	qp=qp*3.0/(4.0*acos(-1)*((radius*radius*radius)-(radius-pthickness)*(radius-pthickness)*(radius-pthickness)));//*1.106183685
	for(j=0;j<n;j++)
		if(region[j]==1) qregion[j]=qp;
		else  qregion[j]=0.0;
	smoothq(qrepeat,qregion,temp);
	//calculate debye-screening length
	debye();
		
	//print total number of points on each region
	printf("%d \t %d \t %d \t %d \t %d \t %d \t%d\t%d\n",ncount[0],ncount[1],ncount[2],ncount[3],ncount[4],ncount[5],ncount[6],ncount[0]+ncount[1]+ncount[2]+ncount[3]+ncount[4]+ncount[5]+ncount[6]);
	//v0=0.0;
	//v1=0.0;
	//calculate v0 of the system
	for(j=0;j<n;j++)
	{
		 if(region[j]==id0) pot0[j]=0.0; 
	else if(region[j]==id1) pot0[j]=0.0;	
	else if(region[j]==id2) pot0[j]=0.0; 
	else if(region[j]==id3) pot0[j]=0.0; 
	else if(region[j]==id4)
	{
		for(i=0;i<nmemu;i++)
		pot0[j]+=1000*zmemu[i]*cmemu[i]*6.02214076e23*hgthickness*ec*slu*0.5/dregion[j];
	}
	else if(region[j]==id5) 
	{
		for(i=0;i<nmemu;i++)
		pot0[j]+=1000*zmemu[i]*cmemu[i]*6.02214076e23*hgthickness*ec*slu*0.5/dregion[j];
	}
	else if(region[j]==id6)
	{
		for(i=0;i<nmemb;i++)
		pot0[j]+=1000*zmemb[i]*cmemb[i]*6.02214076e23*hgthickness*ec*slu*0.5/dregion[j];
	}
	}
	//smooth of v0
	smoothq(prepeat,pot0,temp);
	if(isPrint == 1)
	print2d(nx,pot0,"pot0.dat",1e9,1.0);
	if(isclean == 1)
	    status = system("rm epsilon.gp Epsilon.dat");

	typedef Eigen::Triplet<double> Tp;
	typedef SparseMatrix<double> SpMat;
	Eigen::BiCGSTAB<SpMat> solver;

	SpMat A(n, n);
	VectorXd x(n), beigen(n),x0(n);
	//initial guess of the membrane and inside the particle
	for(j=0;j<n;j++)
	iniguess[j]=pot0[j];
	//if(region[j]==0)
	//{
	//	iniguess[j]=qp*(radius*radius-(radius-pthickness)*(radius-pthickness))/2.0/dregion[j] - qp*4.0*acos(-1)*((radius*radius*radius)-(radius-pthickness)*(radius-pthickness)*//(radius-pthickness))/slu/4.0/acos(-1)/dregion[j]/3.0/(1.0+radius/slu);
	//}
	//else
	//start of newton raphson method
	j=0;
	if(isFd2d == 1)
	while(1)
	{
	//initate of PBE on each grid point
	init ( xvec , yvec , n , val , row , col , rhs , height,dregion,region );
	r8sp_mv ( n, n, count, row, col, val, iniguess,rhs );
	addrhs(val,row,col,rhs,iniguess,count,n,region,qregion,pot0);
	jacobian ( n, val ,row ,col,count, iniguess,region,pot0);
	//transform triplet format to sparsematrix format of eigen
	std::vector< Tp > triplets;
	triplets.reserve(count);
	
	for(i=0;i<count;i++)
	{
		triplets.push_back(Tp(row[i], col[i], val[i]));
	}
	//transform rhs to eigen format
	for(i=0;i<n;i++)
		beigen[i] = rhs[i];
	A.setFromTriplets(triplets.begin(), triplets.end());
	solver.compute(A);
	//solve Ax=b
	x = solver.solve(beigen);
	if(solver.info()!=Eigen::Success) {
	// decomposition failed
	std::cout << "Failure decomposing\n";}
	//r8sp_cg ( n,count,row,col,val,rhs,iniguess1 );
	//calculate the norm of change in the newton raphson method
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error()      << std::endl;
	//add x1+=x0
	for(i=0;i<n;i++)
		{
		iniguess[i]+=x[i];
		}

	tenergy=EnergyEle;
	EnergyEle=freeenergy(iniguess,qregion,region,xvec,pot0);
	if(fabs(EnergyEle-tenergy)<1e-5)break;
	if(j>1500)break;
	j++;
	printf("%d\n",j);
	}
	electricfield(efieldr,efieldz,iniguess,temp);
	//EnergyTot=EnergyBen+EnergyEle;
	if (isNP!=1)
	EnergyEleflat=EnergyEle;
	for(j=0;j<n;j++)
	{qregion[j]+=f(j,iniguess[j],pot0[j]);}
	if (isPrint == 1)
	print2d(nx,qregion,"Rhocharge.dat",1e9,1.0);
	if(isPlot!=0)
	status = system("gnuplot charge.gp");
	if(isclean == 1)
	    status = system("rm charge.gp Rhocharge.dat");
	data_unit=fopen("erg.dat","a");
	//fprintf(data_unit,"#1:Sim_time\t#2:Angle\t#3:E_Bending\t#4:E_Bending_NoSigma\t#5:E_Electrostatic\t#6:E_flatMemEle\t#7:nx\t#8:ny\t#9:r_end\t#10:z_end\t#11:A_freeMem\t#12:A_memTot\n");
	fprintf(data_unit,"%g\t%g\t%.10g\t%.10g\t%.10g\t%.10g\t%d\t%d\t%g\t%g\t%.14g\t%.14g\n",((double) (clock() - times)) / CLOCKS_PER_SEC,angle,EnergyBen,EnergyBenwosigma,EnergyEle,EnergyEleflat,nx,ny,endx*1e9,endy*1e9,memarea,memarea+2*acos(-1)*(radius+thickness/2.0)*(radius+thickness/2.0)*(1.0-cos(angle*acos(-1)/180.0)));
	fclose(data_unit);
	data_unit = fopen ( "Pot.dat", "wt" );
	fprintf ( data_unit, "#Angle=%g\n",angle);
	fprintf ( data_unit, "#\n");
	fprintf ( data_unit, "#\n");
	fprintf ( data_unit, "#Particle radius =%g\n",radius);
	fprintf ( data_unit, "#Box size R_end*Z_end=%5g*%5g\n",endx,endy);
	fprintf ( data_unit, "#h%g=\n",h);
	for ( i = 0; i < ny; i++ )
	{
		for ( j = 0; j < nx; j++ )
		{
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, iniguess[j+i*nx]);
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );
	status = system("gnuplot Pot.gp");
	printf ( "\n" );
	printf ( "	Created data file '%s'\n", "Pot.dat" );



	//printf ( "	Creating C files '%s'\n");
/*
	char command[100];
	int status;
	//Outside Solution

	data_unit = fopen ( "cbkuppos.dat", "wt" );
	for ( i = 0; i < ny; i++ )    //       REVERT THIS LINE 
	//i=ny/2;
	{
		for ( j = 0; j < nx; j++ )
		{	
			if(region[j+i*nx]==id2)
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, cbkup[0]*exp(-zbkup[0]*F*iniguess[j+i*nx]/RT));
			else
			{
				fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, 0.0);
			}
			
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );
	
    sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' cbkuppos.dat >  cbkuppos2dr.dat", endy*1e9/2);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' cbkuppos.dat >  cbkuppos2dz.dat");
    status = system(command);
    status = system("awk  NF  cbkuppos2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp >cbkuppos2dz.dat;rm temp");

	data_unit = fopen ( "cbkupneg.dat", "wt" );
	for ( i = 0; i < ny; i++ )    //       REVERT THIS LINE 
	//i=ny/2;
	{
		for ( j = 0; j < nx; j++ )
		{
			if(region[j+i*nx]==id2)
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, cbkup[1]*exp(-zbkup[1]*F*iniguess[j+i*nx]/RT));
			else
			{
				fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, 0.0);
			}
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );

	sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' cbkupneg.dat >  cbkupneg2dr.dat", endy*1e9/2);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' cbkupneg.dat >  cbkupneg2dz.dat");
    status = system(command);
    status = system("awk  NF  cbkupneg2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp > cbkupneg2dz.dat;rm temp");

	//Inside Solution

	data_unit = fopen ( "cbkbtpos.dat", "wt" );
	for ( i = 0; i < ny; i++ )    //       REVERT THIS LINE 
	//i=ny/2;
	{
		for ( j = 0; j < nx; j++ )
		{	
			if(region[j+i*nx]==id3)
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, cbkbt[0]*exp(-zbkbt[0]*F*iniguess[j+i*nx]/RT));
			else
			{
				fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, 0.0);
			}
			
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );

		sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' cbkbtpos.dat >  cbkbtpos2dr.dat", endy*1e9/2);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' cbkbtpos.dat >  cbkbtpos2dz.dat");
    status = system(command);
    status = system("awk  NF  cbkbtpos2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp > cbkbtpos2dz.dat;rm temp");

	data_unit = fopen ( "cbkbtneg.dat", "wt" );
	for ( i = 0; i < ny; i++ )    //       REVERT THIS LINE 
	//i=ny/2;
	{
		for ( j = 0; j < nx; j++ )
		{
			if(region[j+i*nx]==id3)
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, cbkbt[1]*exp(-zbkbt[1]*F*iniguess[j+i*nx]/RT));
			else
			{
				fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, 0.0);
			}
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );


		sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' cbkbtneg.dat >  cbkbtneg2dr.dat", endy*1e9/2);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' cbkbtneg.dat >  cbkbtneg2dz.dat");
    status = system(command);
    status = system("awk  NF  cbkbtneg2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp > cbkbtneg2dz.dat;rm temp");

	//Outside Membrane

	data_unit = fopen ( "cmemupos.dat", "wt" );
	for ( i = 0; i < ny; i++ )    //       REVERT THIS LINE 
	//i=ny/2;
	{
		for ( j = 0; j < nx; j++ )
		{	
			if(region[j+i*nx]==id4)
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, cmemu[0]*exp(-zmemu[0]*F*(iniguess[j+i*nx]-v0)/RT));
			else
			{
				fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, 0.0);
			}
			
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );

		sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' cmemupos.dat >  cmemupos2dr.dat", endy*1e9/2);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' cmemupos.dat >  cmemupos2dz.dat");
    status = system(command);
    status = system("awk  NF  cmemupos2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp > cmemupos2dz.dat;rm temp");


	data_unit = fopen ( "cmemuneg.dat", "wt" );
	for ( i = 0; i < ny; i++ )    //       REVERT THIS LINE 
	//i=ny/2;
	{
		for ( j = 0; j < nx; j++ )
		{
			if(region[j+i*nx]==id4)
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, cmemu[1]*exp(-zmemu[1]*F*(iniguess[j+i*nx]-v0)/RT));
			else
			{
				fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, 0.0);
			}
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );



	sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' cmemuneg.dat >  cmemuneg2dr.dat", endy*1e9/2);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' cmemuneg.dat >  cmemuneg2dz.dat");
    status = system(command);
    status = system("awk  NF  cmemuneg2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp > cmemuneg2dz.dat;rm temp");



	//Inside Membrane

	data_unit = fopen ( "cmembpos.dat", "wt" );
	for ( i = 0; i < ny; i++ )    //       REVERT THIS LINE 
	//i=ny/2;
	{
		for ( j = 0; j < nx; j++ )
		{	
			if(region[j+i*nx]==id6)
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, cmemb[0]*exp(-zmemb[0]*F*(iniguess[j+i*nx]-v1)/RT));
			else
			{
				fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, 0.0);
			}
			
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );


	sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' cmembpos.dat >  cmembpos2dr.dat", endy*1e9/2);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' cmembpos.dat >  cmembpos2dz.dat");
    status = system(command);
    status = system("awk  NF  cmembpos2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp > cmembpos2dz.dat;rm temp");


	data_unit = fopen ( "cmembneg.dat", "wt" );
	for ( i = 0; i < ny; i++ )    //       REVERT THIS LINE 
	//i=ny/2;
	{
		for ( j = 0; j < nx; j++ )
		{
			if(region[j+i*nx]==id6)
			fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, cmemb[1]*exp(-zmemb[1]*F*(iniguess[j+i*nx]-v1)/RT));
			else
			{
				fprintf ( data_unit, "%14g	%14g	%14g\n",xvec[j]*1e9, yvec[i]*1e9, 0.0);
			}
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );

	sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' cmembneg.dat >  cmembneg2dr.dat", endy*1e9/2);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' cmembneg.dat >  cmembneg2dz.dat");
    status = system(command);
    status = system("awk  NF  cmembneg2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp > cmembneg2dz.dat;rm temp");
	
	
	
	printf ( "\n" );
	printf ( "	Created data file '%s'\n", "Pot.dat" );
*/




//Free memory
	free ( efieldr );
	free ( efieldz );
	free ( iniguess );
	free ( pot0 );
	free ( rhs );
	free ( xvec );
	free ( yvec );
	free ( val );
	free ( row );
	free ( col );
	free (  region );
	free ( dregion );
	free ( qregion );
	free ( temp );
}
