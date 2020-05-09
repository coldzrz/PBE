# include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
# include <stdio.h>
# include <math.h>
# include <time.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "fd2d.h"
#include "functions.h"
#include "sidefunctions.h"
#include "rk4.h"
void gnuploteps(void);
void gnuplotpdf(void);
void gnuplotpng(void);
/******************************************************************************/
using namespace Eigen;
using namespace std;
int main ( )

/******************************************************************************/
 // MAIN is the main program for FD2D_PBE.
{	
  times = clock();
  timestamp ( );
  readparameters();
	int i,j;
	//Change RT value to the correct temperature
  RT=T*8.314459727525675;//calculated from kB*Av =1.38064852e-23*t6.02214076e23
	KbT=T*Kb;
	errorfile=fopen("Errors.txt","a");
  if(isNP==0 && angle <1.0)
  parshift=2.0;
  parshift=endy/parshift;
	h = (endx-startx)/nx;
	hx = (endx-startx)/nx;
	hy = (endy-starty)/ny;
	printf("At start hx=%g\thy=%g\n",(endx-startx)/nx,(endy-starty)/ny);
	double  *radial,*height;
	radial = ( double * ) malloc ( nm* sizeof ( double ) );
	height = ( double * ) malloc ( nm* sizeof ( double ) );
	for(i=0;i<nm;i++)
	{
	radial[i]=0.0;
	radial[i]=0.0;
	}
  memheight(radial,height); //solve shape equations
  if(isPlot == 2)
  {
    gnuploteps();
    //status = system("gnuplot *gp");
  }
  else if(isPlot == 1)
  {
    gnuplotpng();
    //status = system("gnuplot *gp");
  }
  //if(isNP == 0)
  //{qp=0;angle=0.0000001;}
  fd2d (radial,height);
  if(isclean == 1)
  {
    status = system("rm epsilon.gp charge.gp area.gp area.dat Epsilon.dat Rhocharge.dat");
  }
  if(pot2d != 0)
  {
    char command[100];
    sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' Pot.dat >  pot2dr.dat", endy*1e9/parshift);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' Pot.dat >  pot2dz.dat");
    status = system(command);
    status = system("awk  NF  pot2dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp >pot2dz.dat;rm temp");
    sprintf(command, "awk '{if (sqrt(($2-%lf)^2) < \"0.2\") print $1 \"\t\" $3}' Pot0.dat >  pot02dr.dat", endy*1e9/parshift);
    status = system(command);
    sprintf(command, " awk '{if ($1 < \"1.0\")  print $2 \"\t\" $3}' Pot0.dat >  pot02dz.dat");
    status = system(command);
    status = system("awk  NF  pot02dz.dat|awk 'NR>3 ' >temp");
    status = system("cat temp >pot02dz.dat;rm temp");
    if(pot2d == 2)
    status = system("rm Pot.dat");
    // awk '{if (sqrt(($2-50)^2) < \"0.2\") print $1 \"\t\" $3}' Pot.dat >  pot2d.dat
  }
  if(isNP==2)
  {
    isNP=0;
    printf("Before\nnx=%d, ny=%d, endx=%g, endy=%g, hx=%g, hy=%g\n",nx,ny,endx,endy,hx,hy);
    endx=sqrt((memarea+2*acos(-1)*(radius+thickness/2.0)*(radius+thickness/2.0)*(1.0-cos(angle*acos(-1)/180.0)))/acos(-1));
    endy=sqrt((memarea+2*acos(-1)*(radius+thickness/2.0)*(radius+thickness/2.0)*(1.0-cos(angle*acos(-1)/180.0)))/acos(-1));
    nx=endx/hx;
    ny=endy/hy;
    hx = (endx-startx)/nx;
    hy = (endy-starty)/ny;
    printf("After\nnx=%d, ny=%d, endx=%g, endy=%g, hx=%g, hy=%g\n",nx,ny,endx,endy,hx,hy);

    qp=0;angle=0.00001;
    for(i=0;i<nm;i++)
    {
    radial[i]=0.0;
    radial[i]=0.0;
    }
  }
  if(j==0 && isNP!=0)
  {
    system("mv Pot.dat fullPot.dat 2>/dev/null");
    system("mv pot2dr.dat fullpot2dr.dat 2>/dev/null");
    system("mv pot2dz.dat fullpot2dz.dat 2>/dev/null");
    system("mv Pot.png fullPot.png 2>/dev/null");
    system("mv Epsilon.dat fullEpsilon.dat 2>/dev/null");
    system("mv area.dat fullarea.dat 2>/dev/null");
    system("mv Rhocharge.dat fullRhocharge.dat 2>/dev/null");
    system("mv epsilon.png fullEpsilon.png 2>/dev/null");
    system("mv area.png fullarea.png 2>/dev/null");
    system("mv Rhocharge.png fullRhocharge.png 2>/dev/null");
  }
  else
  {
    system("mv Pot.dat flatPot.dat 2>/dev/null");
    system("mv pot2dr.dat flatpot2dr.dat 2>/dev/null");
    system("mv pot2dz.dat flatpot2dz.dat 2>/dev/null");
    system("mv Pot.png flatPot.png 2>/dev/null");
    system("mv Epsilon.dat flatEpsilon.dat 2>/dev/null");
    system("mv area.dat flatarea.dat 2>/dev/null");
    system("mv Rhocharge.dat flatRhocharge.dat 2>/dev/null");
    system("mv epsilon.png flatEpsilon.png 2>/dev/null");
    system("mv area.png flatarea.png 2>/dev/null");
    system("mv Rhocharge.png flatRhocharge.png 2>/dev/null");
  }
  system("mkdir ./datafiles");
  system("mv *.dat ./datafiles");
  system("mv datafiles/PBEparameters.dat .");
  system("mv datafiles/erg.dat .");
  system("mv *.png datafiles");
  system("mv *.gp datafiles");
  fclose(errorfile);
  errorfile=fopen("Energy.dat","a");
	fprintf(errorfile,"#1:Sim_time\t#2:Angle\t#3:E_Bending\t#4:E_Bending_NoSigma\t#5:E_Electrostatic\t#6:E_flatMemEle\t#7:nx\t#8:ny\t#9:r_end\t#10:z_end\t#11:A_freeMem\t#12:A_memTot\n");
  fclose(errorfile);
  int statu;
  statu = system("echo `sed -n -e 1p erg.dat` \" \" ` sed -n -e 2p erg.dat`|awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$18\"\t\"$7\"\t\"$8\"\t\"$9\"\t\"$10\"\t\"$11\"\t\"$12}' >>Energy.dat");

	//int status = system("gnuplot *gp");
	printf("At end hx=%g\thy=%g\n",(endx-startx)/nx,(endy-starty)/ny);
	timestamp ( );

  //if(isNP!=2)
  //{
  //  isNP=0;
  //  endx=sqrt((memarea+2*acos(-1)*(radius+thickness/2.0)*(radius+thickness/2.0)*(1.0-cos(angle*acos(-1)/180.0)))/acos(-1));
  //  endy=sqrt((memarea+2*acos(-1)*(radius+thickness/2.0)*(radius+thickness/2.0)*(1.0-cos(angle*acos(-1)/180.0)))/acos(-1));
  //  nx=endx/hx;
  //  ny=endy/hy;
  //	hx = (endx-startx)/nx;
  //	hy = (endy-starty)/ny;
  //  qp=0;angle=0.00001;
  //  for(i=0;i<nm;i++)
  //  {
  //  radial[i]=0.0;
  //  radial[i]=0.0;
  //  }
  //  memheight(radial,height); //solve shape equations
  //  fd2d (radial,height);
  //}
  
  
  free(radial);
	free(height);
	return 0;
}

/******************************************************************************/

void gnuploteps(void)	//create gnuplot scripts to generate pdf files of potential, charge density, dielectric regions, and points ID

/******************************************************************************/
{
  char command_filename[] = "Pot.gp";
  FILE *command_unit;
  char data_filename[] = "Pot.eps";
  FILE *data_unit;
  //Potential file
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set term post eps enhanced color font 'Arial,30' \n" );
  fprintf ( command_unit, "set output '%s'	\n",data_filename );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set border lw 3\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set mxtics 3\n" );
  fprintf ( command_unit, "set mytics 3\n" );
  fprintf ( command_unit, "set xlabel '{/:Bold r (nm)}'\n" );
  fprintf ( command_unit, "set ylabel '{/:Bold z (nm)}'\n" );
  fprintf ( command_unit, "set cblabel '{/:Bold V (k_{B}T/e)}'\n" );
  fprintf ( command_unit, "set xtics %f,%f,%f\n",startx*1e9+5,(endx*1e9-5-(startx*1e9+5))/3.0,endx*1e9-5 );
  fprintf ( command_unit, "set ytics %f,%f,%f\n",starty*1e9,(endy*1e9-starty*1e9)/4.0,endy*1e9-1 );
  fprintf ( command_unit, "#set cbtics 0.2,0.4,1.4\n" );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot 'Pot.dat' u 1:2:($3*%f) notitle \n",ec/KbT);
  fprintf ( command_unit, "system('epstopdf Pot.eps') \n");
  fprintf ( command_unit, "#system('epstopdf Pot.eps')\n");
  fprintf ( command_unit, "system('rm Pot.eps')\n");
  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );


  //Dielectric file
  strcpy(command_filename,"epsilon.gp");
  //strcpy(data_filename,"epsilon.eps");
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set term post eps enhanced color font 'Arial,30' \n" );
  fprintf ( command_unit, "set output 'epsilon.eps'	\n" );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set border lw 3\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set mxtics 3\n" );
  fprintf ( command_unit, "set mytics 3\n" );
  fprintf ( command_unit, "set xlabel '{/:Bold r (nm)}'\n" );
  fprintf ( command_unit, "set ylabel '{/:Bold z (nm)}'\n" );
  fprintf ( command_unit, "set cblabel '{/:Bold {/Symbol e}_r}'\n" );
  fprintf ( command_unit, "set xtics %f,%f,%f\n",startx*1e9+5,(endx*1e9-5-(startx*1e9+5))/3.0,endx*1e9-5 );
  fprintf ( command_unit, "set ytics %f,%f,%f\n",starty*1e9,(endy*1e9-starty*1e9)/4.0,endy*1e9-1 );
  fprintf ( command_unit, "set cbtics %f,%f,%f\n",0,20,80 );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot 'Epsilon.dat' u 1:2:3 notitle \n");
  fprintf ( command_unit, "system('epstopdf epsilon.eps') \n");
  fprintf ( command_unit, "#system('epstopdf epsilon.eps')\n");
  fprintf ( command_unit, "system('rm epsilon.eps')\n");
  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );


  //Charge file
  strcpy(command_filename,"Rhocharge.gp");
  command_unit = fopen ( command_filename, "wt" );
  //strcpy(data_filename,"Rhocharge.eps");
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set term post eps enhanced color font 'Arial,30' \n" );
  fprintf ( command_unit, "set output 'charge.eps'	\n");
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set border lw 3\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set mxtics 3\n" );
  fprintf ( command_unit, "set mytics 3\n" );
  fprintf ( command_unit, "set xlabel '{/:Bold r (nm)}'\n" );
  fprintf ( command_unit, "set ylabel '{/:Bold z (nm)}'\n" );
  fprintf ( command_unit, "set cblabel '{/:Bold {/Symbol r} (C/nm^{3})}'\n" );
  fprintf ( command_unit, "set xtics %f,%f,%f\n",startx*1e9+5,(endx*1e9-5-(startx*1e9+5))/3.0,endx*1e9-5 );
  fprintf ( command_unit, "set ytics %f,%f,%f\n",starty*1e9,(endy*1e9-starty*1e9)/4.0,endy*1e9-1 );
  fprintf ( command_unit, "#set cbtics 0.2,0.4,1.4\n" );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot 'Rhocharge.dat' u 1:2:3 notitle \n");
  fprintf ( command_unit, "system('epstopdf charge.eps') \n");
  fprintf ( command_unit, "#system('epstopdf charge.eps')\n");
  fprintf ( command_unit, "system('rm charge.eps')\n");
  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );





  //System file
  strcpy(command_filename,"area.gp");
  //strcpy(data_filename,"area.eps");
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set term post eps enhanced color font 'Arial,30' \n" );
  fprintf ( command_unit, "set output 'area.eps'	\n");
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set border lw 3\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set mxtics 3\n" );
  fprintf ( command_unit, "set mytics 3\n" );
  fprintf ( command_unit, "set xlabel '{/:Bold r (nm)}'\n" );
  fprintf ( command_unit, "set ylabel '{/:Bold z (nm)}'\n" );
  fprintf ( command_unit, "set cblabel '{/:Bold Grid point ID}'\n" );
  fprintf ( command_unit, "set xtics %f,%f,%f\n",startx*1e9+5,(endx*1e9-5-(startx*1e9+5))/3.0,endx*1e9-5 );
  fprintf ( command_unit, "set ytics %f,%f,%f\n",starty*1e9,(endy*1e9-starty*1e9)/4.0,endy*1e9-1 );
  fprintf ( command_unit, "#set cbtics 0.2,0.4,1.4\n" );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot 'area.dat' u 1:2:3 notitle \n");
  fprintf ( command_unit, "system('epstopdf area.eps') \n");
  fprintf ( command_unit, "#system('epstopdf area.eps')\n");
  fprintf ( command_unit, "system('rm area.eps')\n");
  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );



}
/******************************************************************************/

void gnuplotpdf(void)	//create gnuplot scripts to generate pdf files of potential, charge density, dielectric regions, and points ID

/******************************************************************************/
{
  char command_filename[] = "Pot.gp";
  FILE *command_unit;
  char data_filename[] = "test01_data.txt";
  FILE *data_unit;
  //Potential file
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set samples 100000\n" );
  fprintf ( command_unit, "set term epslatex standalone solid color  colortext \n" );
  fprintf ( command_unit, "filename = 'potential'	\n" );
  fprintf ( command_unit, "set output filename.'.tex'\n" );
  fprintf ( command_unit, "#set yrange [0:10]\n" );
  fprintf ( command_unit, "#set xrange [0:20]\n" );
  fprintf ( command_unit, "#set cbrange [-0.1:0]\n" );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size 1,1\n" );
  fprintf ( command_unit, "set xlabel 'z (nm)'\n" );
  fprintf ( command_unit, "set ylabel 'r (nm)'\n" );
  fprintf ( command_unit, "set zlabel '<---U(X,Y)--->'\n" );
  fprintf ( command_unit, "set title 'Potential'\n" );
  fprintf ( command_unit, "set lmargin screen 0.12\n" );
  fprintf ( command_unit, "set rmargin screen 0.83\n" );
  fprintf ( command_unit, "#set bmargin screen 0.05\n" );
  fprintf ( command_unit, "#set tmargin screen 1.0\n" );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot 'Pot.dat' u 1:2:($3*%f) notitle \n",ec/KbT);
  fprintf ( command_unit, "set output \n");
  fprintf ( command_unit, "system 'xelatex -interaction=nonstopmode '.filename.'.tex' \n");
  fprintf ( command_unit, "system 'rm '.filename.'{.log,.aux,.tex,-inc.eps}'\n");
  fprintf ( command_unit, "set t pop \n");
  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );

  //Dielectric file
  strcpy(command_filename,"epsilon.gp");
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set samples 100000\n" );
  fprintf ( command_unit, "set term epslatex standalone solid color  colortext \n" );
  fprintf ( command_unit, "filename = 'dielectric'	\n" );
  fprintf ( command_unit, "set output filename.'.tex'\n" );
  fprintf ( command_unit, "#set yrange [0:10]\n" );
  fprintf ( command_unit, "#set xrange [0:20]\n" );
  fprintf ( command_unit, "#set cbrange [-0.1:0]\n" );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size 1,1\n" );
  fprintf ( command_unit, "set xlabel 'z (nm)'\n" );
  fprintf ( command_unit, "set ylabel 'r (nm)'\n" );
  fprintf ( command_unit, "set zlabel '<---U(X,Y)--->'\n" );
  fprintf ( command_unit, "set title 'Dielectric'\n" );
  fprintf ( command_unit, "set lmargin screen 0.12\n" );
  fprintf ( command_unit, "set rmargin screen 0.83\n" );
  fprintf ( command_unit, "#set bmargin screen 0.05\n" );
  fprintf ( command_unit, "#set tmargin screen 1.0\n" );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot 'Epsilon.dat' u 1:2:3 notitle \n");
  fprintf ( command_unit, "set output \n");
  fprintf ( command_unit, "system 'xelatex -interaction=nonstopmode '.filename.'.tex' \n");
  fprintf ( command_unit, "system 'rm '.filename.'{.log,.aux,.tex,-inc.eps}'\n");
  fprintf ( command_unit, "set t pop \n");

  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );

  //Charge file
  strcpy(command_filename,"charge.gp");
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set samples 100000\n" );
  fprintf ( command_unit, "set term epslatex standalone solid color  colortext \n" );
  fprintf ( command_unit, "filename = 'charge'	\n" );
  fprintf ( command_unit, "set output filename.'.tex'\n" );
  fprintf ( command_unit, "#set yrange [0:10]\n" );
  fprintf ( command_unit, "#set xrange [0:20]\n" );
  fprintf ( command_unit, "#set cbrange [-0.1:0]\n" );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size 1,1\n" );
  fprintf ( command_unit, "set xlabel 'z (nm)'\n" );
  fprintf ( command_unit, "set ylabel 'r (nm)'\n" );
  fprintf ( command_unit, "set zlabel '<---U(X,Y)--->'\n" );
  fprintf ( command_unit, "set title 'Charge Density'\n" );
  fprintf ( command_unit, "set lmargin screen 0.12\n" );
  fprintf ( command_unit, "set rmargin screen 0.83\n" );
  fprintf ( command_unit, "#set bmargin screen 0.05\n" );
  fprintf ( command_unit, "#set tmargin screen 1.0\n" );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot 'Rhocharge.dat' u 1:2:3 notitle \n");
  fprintf ( command_unit, "set output \n");
  fprintf ( command_unit, "system 'xelatex -interaction=nonstopmode '.filename.'.tex' \n");
  fprintf ( command_unit, "system 'rm '.filename.'{.log,.aux,.tex,-inc.eps}'\n");
  fprintf ( command_unit, "set t pop \n");

  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );
  
  //SystemID file
  strcpy(command_filename,"area.gp");
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set samples 100000\n" );
  fprintf ( command_unit, "set term epslatex standalone solid color  colortext \n" );
  fprintf ( command_unit, "filename = 'System'	\n" );
  fprintf ( command_unit, "set output filename.'.tex'\n" );
  fprintf ( command_unit, "#set yrange [0:10]\n" );
  fprintf ( command_unit, "#set xrange [0:20]\n" );
  fprintf ( command_unit, "#set cbrange [-0.1:0]\n" );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size 1,1\n" );
  fprintf ( command_unit, "set xlabel 'z (nm)'\n" );
  fprintf ( command_unit, "set ylabel 'r (nm)'\n" );
  fprintf ( command_unit, "set zlabel '<---U(X,Y)--->'\n" );
  fprintf ( command_unit, "set title 'Charge Density'\n" );
  fprintf ( command_unit, "set lmargin screen 0.12\n" );
  fprintf ( command_unit, "set rmargin screen 0.83\n" );
  fprintf ( command_unit, "#set bmargin screen 0.05\n" );
  fprintf ( command_unit, "#set tmargin screen 1.0\n" );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot 'area.dat' u 2:1:3 notitle \n");
  fprintf ( command_unit, "set output \n");
  fprintf ( command_unit, "system 'xelatex -interaction=nonstopmode '.filename.'.tex' \n");
  fprintf ( command_unit, "system 'rm '.filename.'{.log,.aux,.tex,-inc.eps}'\n");
  fprintf ( command_unit, "set t pop \n");

  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );
}

/******************************************************************************/

void gnuplotpng(void)	//create gnuplot scripts to generate png files of potential, charge density, dielectric regions, and points ID

/******************************************************************************/
{
  char command_filename[] = "Pot.gp";
  FILE *command_unit;
  char data_filename[] = "Pot.png";
  FILE *data_unit;
  //Potential file
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set term pngcairo size 4000,2200 enhanced color font 'Arial,90' \n" );
  fprintf ( command_unit, "set output '%s'	\n",data_filename );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set border lw 2\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set mxtics 3\n" );
  fprintf ( command_unit, "set mytics 3\n" );
  fprintf ( command_unit, "set tics scale 0.6\n" );
  fprintf ( command_unit, "set border linecolor rgb 'grey75'\n" );
  fprintf ( command_unit, "set tic textcolor rgb 'black'\n" );
  fprintf ( command_unit, "set xlabel '{/:Bold r (nm)}'\n" );
  fprintf ( command_unit, "set ylabel '{/:Bold z (nm)}'\n" );
  fprintf ( command_unit, "set cblabel '{/:Bold V (k_{B}T/e)}'\n" );
  fprintf ( command_unit, "set xtics %f,%f,%f\n",startx*1e9,(endx*1e9-startx*1e9)/4.0,endx*1e9-1 );
  fprintf ( command_unit, "set ytics %f,%f,%f\n",starty*1e9,(endy*1e9-starty*1e9)/4.0,endy*1e9-1 );
  fprintf ( command_unit, "#set cbtics 0.2,0.4,1.4\n" );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot [%f:%f][%f:%f]'Pot.dat' u 1:2:($3*%f) notitle \n",startx*1e9+1,endx*1e9-1,starty*1e9+1,endy*1e9-1,ec/KbT);
  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );

  //Dielectric file
  strcpy(command_filename,"epsilon.gp");
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set term pngcairo size 4000,2200 enhanced color font 'Arial,90' \n" );
  fprintf ( command_unit, "set output 'Epsilon.png'	\n" );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set border lw 2\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set mxtics 3\n" );
  fprintf ( command_unit, "set mytics 3\n" );
  fprintf ( command_unit, "set tics scale 0.6\n" );
  fprintf ( command_unit, "set border linecolor rgb 'grey75'\n" );
  fprintf ( command_unit, "set tic textcolor rgb 'black'\n" );
  fprintf ( command_unit, "set xlabel '{/:Bold r (nm)}'\n" );
  fprintf ( command_unit, "set ylabel '{/:Bold z (nm)}'\n" );
  fprintf ( command_unit, "set cblabel '{/:Bold {/Symbol e}_r}'\n" );
  fprintf ( command_unit, "set xtics %f,%f,%f\n",startx*1e9,(endx*1e9-startx*1e9)/4.0,endx*1e9-1 );
  fprintf ( command_unit, "set ytics %f,%f,%f\n",starty*1e9,(endy*1e9-starty*1e9)/4.0,endy*1e9-1 );
  fprintf ( command_unit, "set cbtics %d,%d,%d\n",0,20,79 );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot [%f:%f][%f:%f]'Epsilon.dat' u 1:2:3 notitle \n",startx*1e9+1,endx*1e9-1,starty*1e9+1,endy*1e9-1);

  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );

  //Charge file
  strcpy(command_filename,"charge.gp");
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set term pngcairo size 4000,2200 enhanced color font 'Arial,90' \n" );
  fprintf ( command_unit, "set output 'Rhocharge.png'	\n" );
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set border lw 2\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set mxtics 3\n" );
  fprintf ( command_unit, "set mytics 3\n" );
  fprintf ( command_unit, "set tics scale 0.6\n" );
  fprintf ( command_unit, "set border linecolor rgb 'grey75'\n" );
  fprintf ( command_unit, "set tic textcolor rgb 'black'\n" );
  fprintf ( command_unit, "set xlabel '{/:Bold r (nm)}'\n" );
  fprintf ( command_unit, "set ylabel '{/:Bold z (nm)}'\n" );
  fprintf ( command_unit, "set cblabel '{/:Bold {/Symbol r} (C/nm^{3})}'\n" );
  fprintf ( command_unit, "set xtics %f,%f,%f\n",startx*1e9,(endx*1e9-startx*1e9)/4.0,endx*1e9-1 );
  fprintf ( command_unit, "set ytics %f,%f,%f\n",starty*1e9,(endy*1e9-starty*1e9)/4.0,endy*1e9-1 );
  fprintf ( command_unit, "#set cbtics %f,%f,%f\n",0,qp/4.0,qp );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot [%f:%f][%f:%f]'Rhocharge.dat' u 1:2:3 notitle \n",startx*1e9+1,endx*1e9-1,starty*1e9+1,endy*1e9-1);

  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );
  
  //SystemID file
  strcpy(command_filename,"area.gp");
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "reset\n" );
  fprintf ( command_unit, "set term pngcairo size 2000,1100 enhanced color font 'Arial,50' \n" );
  fprintf ( command_unit, "set output 'Area.png'\n");
  fprintf ( command_unit, "set pm3d map interpolate 2,2\n" );
  fprintf ( command_unit, "set size ratio -1\n" );
  fprintf ( command_unit, "set border lw 2\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set mxtics 3\n" );
  fprintf ( command_unit, "set mytics 3\n" );
  fprintf ( command_unit, "set tics scale 0.6\n" );
  fprintf ( command_unit, "set border linecolor rgb 'grey75'\n" );
  fprintf ( command_unit, "set tic textcolor rgb 'black'\n" );
  fprintf ( command_unit, "set xlabel '{/:Bold r (nm)}'\n" );
  fprintf ( command_unit, "set ylabel '{/:Bold z (nm)}'\n" );
  fprintf ( command_unit, "set cblabel '{/:Bold Grid point ID}'\n" );
  fprintf ( command_unit, "set xtics %f,%f,%f\n",startx*1e9,(endx*1e9-startx*1e9)/4.0,endx*1e9-1 );
  fprintf ( command_unit, "set ytics %f,%f,%f\n",starty*1e9,(endy*1e9-starty*1e9)/4.0,endy*1e9-1 );
  fprintf ( command_unit, "set cbtics %d,%d,%d\n",0,1,id9 );
  fprintf ( command_unit, "set palette defined (0  0.0 0.0 0.5, \\\n" );
  fprintf ( command_unit, "                     1  0.0 0.0 1.0, \\\n" );
  fprintf ( command_unit, "                     2  0.0 0.5 1.0, \\\n" );
  fprintf ( command_unit, "                     3  0.0 1.0 1.0, \\\n" );
  fprintf ( command_unit, "                     4  0.5 1.0 0.5, \\\n" );
  fprintf ( command_unit, "                     5  1.0 1.0 0.0, \\\n" );
  fprintf ( command_unit, "                     6  1.0 0.5 0.0, \\\n" );
  fprintf ( command_unit, "                     7  1.0 0.0 0.0, \\\n" );
  fprintf ( command_unit, "                     8  0.5 0.0 0.0 ) \n" );
  fprintf ( command_unit, "splot [%f:%f][%f:%f]'area.dat' u 1:2:3 notitle \n",startx*1e9+1,endx*1e9-1,starty*1e9+1,endy*1e9-1);
  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );
}
