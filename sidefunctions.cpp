# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>
#include <string.h>
# include "sidefunctions.h"
# include "functions.h"
# include <stdbool.h>
# include "rk4.h"
#include <float.h>
double potential(double r)
{

if(r<radius)
return -0.287753;
else
return -1.438768184/r;
}


/******************************************************************************/
void regioncountmaple(int region[],double radial[],double height[],double x[], double y[])

/******************************************************************************/
{
	//This function goes over all grid points and give them id:{0,1,2,3,4,5} for inside particle, particle outter layer, upper solution, lower solution, headgroup and tailgroup respectivily
	int i,j,k,n,total,a,b,c,d,e,f,g;
	double xj,yi,R,Ro,Ri,r,hgt,pt,memr,memd;
	hgt=hgthickness;
	pt=pthickness;
	memr=0.0;
	memd=0.0;
	n=nx*ny;
	a = id0;
	b = id1;
	c = id2;
	d = id3;
	e = id4;
	f = id5;
	g = id6;
	if (isNP == 0)
		{a =id2;b=id2;radius=0.0;pt=0.0;}
	if (isMem == 0)
		{e=id2;f=id2;d=id2;}
	R=(radius)*(radius);
	Ri=(radius-pt)*(radius-pt);
	//R=(radius-h*2)*(radius-h*2);
	//Ri=(radius-pt-h*2)*(radius-pt-h*2);
	for(i=0;i<ny;i++)
		for(j=0;j<nx;j++)
		region[j+i*nx]=-1;
	ncount[0]=0,ncount[1]=0,ncount[2]=0,ncount[3]=0,ncount[4]=0,ncount[5]=0,ncount[6]=0,ncount[7]=0;
	for(i=0;i<ny;i++)
	{
		for(j=0;j<nx;j++)
		{
			xj=x[j];
			yi=y[i];
			r=(xj)*(xj)+(yi-parshift)*(yi-parshift);
			//find ID0 and ID1 of particle
			if(r<=Ri)
			{region[j+i*nx]=id0;}
			else if(r<=R)
			{region[j+i*nx]=id1;}
			else
			//Find ID4 and ID6 of the membrane
			for(k=0;k<nm;k++)
			{
				memd=sqrt((height[k]-yi)*(height[k]-yi)+(radial[k]-xj)*(radial[k]-xj));
				if(memd<=thickness/2.0)
				region[j+i*nx]=id4;
				if(memd<=(thickness/2.0-hgt))
				{region[j+i*nx]=id5;break;}
			}
			//Add membrnae wrapped around particle
			for(k=0;k<=angle+0.01*angle;k++) //
			{
				if (isNP == 0 && angle <1.0)
					break;
				if(region[j+i*nx]==id5)
				continue;
				memd=sqrt((parshift-(radius+thickness/2.0)*cos(k*acos(-1)/180)-yi)*(parshift-(radius+thickness/2.0)*cos(k*acos(-1)/180)-yi)
				+((radius+thickness/2.0)*sin(k*acos(-1)/180)-xj)*((radius+thickness/2.0)*sin(k*acos(-1)/180)-xj));
				if(memd<=thickness/2.0)//(thickness/2.0-hgt)<memd &&
				region[j+i*nx]=id4;
				if(memd<=(thickness/2.0-hgt))
				{region[j+i*nx]=id5;break;}
				if(region[j+i*nx]==-1)
				{
					memd=sqrt((parshift-yi)*(parshift-yi)+(-xj)*(-xj));
					if(memd<=(radius+thickness/2.0) && yi <(parshift-(radius+thickness/2.0)*cos((angle-0.1*angle)*acos(-1)/180) ))
						region[j+i*nx]=id4;	
				}
			}
			if(region[j+i*nx]!=id0 && region[j+i*nx]!=id1&& region[j+i*nx]!=e&& region[j+i*nx]!=id5)
			region[j+i*nx]=id2;
			
			ncount[region[j+i*nx]]++;
		}
	}
	ncount[region[0]]--;
	region[0]=id3;
	ncount[region[0]]++;
	j=0;
	//find ID3	
	for(i=0;i<ny-1;i++)
	{
		for(j=0;j<nx-1;j++)
		{
			if(region[j+i*nx]==id3 && j+(i)*nx > -1 && j+(i)*nx < nx*ny+1 )
				{
					if (j-1+(i-1)*nx > -1 && j-1+(i-1)*nx < nx*ny+1 && region[j-1+(i-1)*nx]==id2 )
					{
						ncount[region[j-1+(i-1)*nx]]--;region[j-1+(i-1)*nx]=id3;ncount[region[j-1+(i-1)*nx]]++;
					}
					if ( j-1+(i)*nx > -1 && j-1+(i)*nx < nx*ny+1 && region[j-1+(i)*nx]==id2 )
					{
						ncount[region[j-1+(i)*nx]]--;region[j-1+(i)*nx]=id3;ncount[region[j-1+(i)*nx]]++;
					}
					if (j-1+(i+1)*nx > -1 && j-1+(i+1)*nx < nx*ny+1 && region[j-1+(i+1)*nx]==id2  )
					{
						ncount[region[j-1+(i+1)*nx]]--;region[j-1+(i+1)*nx]=id3;ncount[region[j-1+(i+1)*nx]]++;
					}
					if (j+(i-1)*nx > -1 && j+(i-1)*nx < nx*ny+1 && region[j+(i-1)*nx]==id2  )
					{
						ncount[region[j+(i-1)*nx]]--;region[j+(i-1)*nx]=id3;ncount[region[j+(i-1)*nx]]++;
					}
					if ( j+(i+1)*nx > -1 && j+(i+1)*nx < nx*ny+1 && region[j+(i+1)*nx]==id2 )
					{
						ncount[region[j+(i+1)*nx]]--;region[j+(i+1)*nx]=id3;ncount[region[j+(i+1)*nx]]++;
					}
					if (j+1+(i-1)*nx > -1 && j+1+(i-1)*nx < nx*ny+1 && region[j+1+(i-1)*nx]==id2 )
					{
						ncount[region[j+1+(i-1)*nx]]--;region[j+1+(i-1)*nx]=id3;ncount[region[j+1+(i-1)*nx]]++;
					}
					if ( j+1+(i)*nx > -1 && j+1+(i)*nx < nx*ny+1 && region[j+1+(i)*nx]==id2 )
					{
						ncount[region[j+1+(i)*nx]]--;region[j+1+(i)*nx]=id3;ncount[region[j+1+(i)*nx]]++;
					}
					if (j+1+(i+1)*nx > -1 && j+1+(i+1)*nx < nx*ny+1 && region[j+1+(i+1)*nx]==id2  )
					{
						ncount[region[j+1+(i+1)*nx]]--;region[j+1+(i+1)*nx]=id3;ncount[region[j+1+(i+1)*nx]]++;
					}
					
				}	
		}
	}
	//find ID3
	for(i=ny-1;i>0;i--)
	{
		for(j=nx-1;j>0;j--)
		{
			if(region[j+i*nx]==id3 && j+(i)*nx > -1 && j+(i)*nx < nx*ny+1 )
				{
					if (j-1+(i-1)*nx > -1 && j-1+(i-1)*nx < nx*ny+1 && region[j-1+(i-1)*nx]==id2 )
					{
						ncount[region[j-1+(i-1)*nx]]--;region[j-1+(i-1)*nx]=id3;ncount[region[j-1+(i-1)*nx]]++;
					}
					if ( j-1+(i)*nx > -1 && j-1+(i)*nx < nx*ny+1 && region[j-1+(i)*nx]==id2 )
					{
						ncount[region[j-1+(i)*nx]]--;region[j-1+(i)*nx]=id3;ncount[region[j-1+(i)*nx]]++;
					}
					if (j-1+(i+1)*nx > -1 && j-1+(i+1)*nx < nx*ny+1 && region[j-1+(i+1)*nx]==id2  )
					{
						ncount[region[j-1+(i+1)*nx]]--;region[j-1+(i+1)*nx]=id3;ncount[region[j-1+(i+1)*nx]]++;
					}
					if (j+(i-1)*nx > -1 && j+(i-1)*nx < nx*ny+1 && region[j+(i-1)*nx]==id2  )
					{
						ncount[region[j+(i-1)*nx]]--;region[j+(i-1)*nx]=id3;ncount[region[j+(i-1)*nx]]++;
					}
					if ( j+(i+1)*nx > -1 && j+(i+1)*nx < nx*ny+1 && region[j+(i+1)*nx]==id2 )
					{
						ncount[region[j+(i+1)*nx]]--;region[j+(i+1)*nx]=id3;ncount[region[j+(i+1)*nx]]++;
					}
					if (j+1+(i-1)*nx > -1 && j+1+(i-1)*nx < nx*ny+1 && region[j+1+(i-1)*nx]==id2 )
					{
						ncount[region[j+1+(i-1)*nx]]--;region[j+1+(i-1)*nx]=id3;ncount[region[j+1+(i-1)*nx]]++;
					}
					if ( j+1+(i)*nx > -1 && j+1+(i)*nx < nx*ny+1 && region[j+1+(i)*nx]==id2 )
					{
						ncount[region[j+1+(i)*nx]]--;region[j+1+(i)*nx]=id3;ncount[region[j+1+(i)*nx]]++;
					}
					if (j+1+(i+1)*nx > -1 && j+1+(i+1)*nx < nx*ny+1 && region[j+1+(i+1)*nx]==id2  )
					{
						ncount[region[j+1+(i+1)*nx]]--;region[j+1+(i+1)*nx]=id3;ncount[region[j+1+(i+1)*nx]]++;
					}
					
				}	
		}
	}
	//find ID6
	for(i=1;i<ny-1;i++)
	{
		for(j=1;j<nx-1;j++)
		{
			if(region[j+i*nx]==id3 ||region[j+i*nx]==id6)
				{
					if (region[j-1+(i-1)*nx]==id4)
					{
						ncount[region[j-1+(i-1)*nx]]--;region[j-1+(i-1)*nx]=id6;ncount[region[j-1+(i-1)*nx]]++;
					}
					if (region[j-1+(i)*nx]==id4)
					{
						ncount[region[j-1+(i)*nx]]--;region[j-1+(i)*nx]=id6;ncount[region[j-1+(i)*nx]]++;
					}
					if (region[j-1+(i+1)*nx]==id4)
					{
						ncount[region[j-1+(i+1)*nx]]--;region[j-1+(i+1)*nx]=id6;ncount[region[j-1+(i+1)*nx]]++;
					}
					if (region[j+(i-1)*nx]==id4)
					{
						ncount[region[j+(i-1)*nx]]--;region[j+(i-1)*nx]=id6;ncount[region[j+(i-1)*nx]]++;
					}
					if (region[j+(i+1)*nx]==id4)
					{
						ncount[region[j+(i+1)*nx]]--;region[j+(i+1)*nx]=id6;ncount[region[j+(i+1)*nx]]++;
					}
					if (region[j+1+(i-1)*nx]==id4)
					{
						ncount[region[j+1+(i-1)*nx]]--;region[j+1+(i-1)*nx]=id6;ncount[region[j+1+(i-1)*nx]]++;
					}
					if (region[j+1+(i)*nx]==id4)
					{
						ncount[region[j+1+(i)*nx]]--;region[j+1+(i)*nx]=id6;ncount[region[j+1+(i)*nx]]++;
					}
					if (region[j+1+(i+1)*nx]==id4)
					{
						ncount[region[j+1+(i+1)*nx]]--;region[j+1+(i+1)*nx]=id6;ncount[region[j+1+(i+1)*nx]]++;
					}
					
				}	
		}
	}
	//find ID6
	for(i=ny-2;i>1;i--)
	{
		for(j=nx-2;j>1;j--)
		{
			if(region[j+i*nx]==id3 ||region[j+i*nx]==id6)
				{
					if (region[j-1+(i-1)*nx]==id4)
					{
						ncount[region[j-1+(i-1)*nx]]--;region[j-1+(i-1)*nx]=id6;ncount[region[j-1+(i-1)*nx]]++;
					}
					if (region[j-1+(i)*nx]==id4)
					{
						ncount[region[j-1+(i)*nx]]--;region[j-1+(i)*nx]=id6;ncount[region[j-1+(i)*nx]]++;
					}
					if (region[j-1+(i+1)*nx]==id4)
					{
						ncount[region[j-1+(i+1)*nx]]--;region[j-1+(i+1)*nx]=id6;ncount[region[j-1+(i+1)*nx]]++;
					}
					if (region[j+(i-1)*nx]==id4)
					{
						ncount[region[j+(i-1)*nx]]--;region[j+(i-1)*nx]=id6;ncount[region[j+(i-1)*nx]]++;
					}
					if (region[j+(i+1)*nx]==id4)
					{
						ncount[region[j+(i+1)*nx]]--;region[j+(i+1)*nx]=id6;ncount[region[j+(i+1)*nx]]++;
					}
					if (region[j+1+(i-1)*nx]==id4)
					{
						ncount[region[j+1+(i-1)*nx]]--;region[j+1+(i-1)*nx]=id6;ncount[region[j+1+(i-1)*nx]]++;
					}
					if (region[j+1+(i)*nx]==id4)
					{
						ncount[region[j+1+(i)*nx]]--;region[j+1+(i)*nx]=id6;ncount[region[j+1+(i)*nx]]++;
					}
					if (region[j+1+(i+1)*nx]==id4)
					{
						ncount[region[j+1+(i+1)*nx]]--;region[j+1+(i+1)*nx]=id6;ncount[region[j+1+(i+1)*nx]]++;
					}
					
				}	
		}
	}
	
	//To remove particle only and keep membrane shape fixed
	if(isNP==0)
	for(i=0;i<ny-1;i++)
	{
		for(j=0;j<nx-1;j++)
		{
			if (region[j+(i)*nx]==id0 ||  region[j+(i)*nx]==id1 )
					{
						ncount[region[j+(i)*nx]]--;region[j+(i)*nx]=id2;ncount[region[j+(i)*nx]]++;
					}
		}
	}
	if(isMem==0)
	for(i=0;i<ny;i++)
	{
		for(j=0;j<nx;j++)
		{
			if (region[j+(i)*nx]!=id0 &&  region[j+(i)*nx]!=id1 )
					{
						ncount[region[j+(i)*nx]]--;region[j+(i)*nx]=id2;ncount[region[j+(i)*nx]]++;
					}
		}
	}
	
	//To have qrepeat of particle outer layer thickness changed to inner particle
	if(isNP!=0)
	for(i=0;i<ny;i++)
	{
		for(j=0;j<nx;j++)
		{
			xj=x[j];
			yi=y[i];
			r=(xj)*(xj)+(yi-parshift)*(yi-parshift);
			//find ID0 and ID1 of particle
			if(r<=R && r >= (R-hx*qrepeat*hx*qrepeat))
			{ncount[region[j+(i)*nx]]--;region[j+(i)*nx]=id0;ncount[region[j+(i)*nx]]++;}
		}
	}
		printf("Finished classifying regions\n");


}/******************************************************************************/
double mtm(double a[],double b[],int n)
{
	//compute vTv
	int i,j,k;
	double sum=0.0,total=0.0;
	for(i=0;i<n;i++)
	{
		sum+=a[i]*b[i];
	}
	return sum;
}
void add(double res[],double a[],double b[],double coeff,int n)
{
	//add matrix with another multiplied by a coefficient
	int i,j,k;
	double sum=0.0,total=0.0;
	for(i=0;i<n;i++)
	{
		res[i]=a[i]+coeff*b[i];
	}
}
void copy(double a[],double b[],int n)
{
	//copy double matrix
	int i;
	for(i=0;i<n;i++)
	{
		b[i]=a[i];
	}
}
bool check(double a[],double b[],int n)
{
	//check that the change in each point is less than tolerance
	int counter=0,i;
	double err=1e-3;
	for(i=0;i<n;i++)
			if(fabs(b[i]-a[i])<err)
				counter++;
	//printf("Error %d\n",counter);
	if(counter>=n)
		return false;
		else 
		return true;
}
/******************************************************************************/

void printrhs(double rhs[],int n,char *data_filename)

/******************************************************************************/
{
	//print vector to file
	int i,j;
	//char data_filename[] = "rhs.txt";
	FILE *data_unit;
	data_unit = fopen ( data_filename, "wt" );
	for ( i = 0; i < n; i++ )
	{
		fprintf ( data_unit, "%g\n",rhs[i]);
	}
	fclose ( data_unit );
	printf ( "\n" );
	printf ( "	Created graphics data file '%s'\n", data_filename );
}
/******************************************************************************/

double *r8vec_linspace_new ( int n, double a, double b )

/******************************************************************************/
{
	int i;
	double *x;
	x = ( double * ) malloc ( n * sizeof ( double ) );

	if ( n == 1 )
	{
		x[0] = ( a + b ) / 2.0;
	}
	else
	{
		for ( i = 0; i < n; i++ )
		{
			x[i] = ( ( double ) ( n - 1.0 - i ) * a  + ( double ) ( i ) * b )  / ( double ) ( n - 1.0);
		}
	}
	return x;
}
/******************************************************************************/

void print2dint(int n,int a[],char *data_filename)

/******************************************************************************/
{
//print region id for PB (particle, membrane, solution)
	int i,j;
	FILE *data_unit;
	data_unit = fopen ( data_filename, "wt" );
	for ( i = 0; i < ny; i++ )
	{
		for ( j = 0; j < nx; j++ )
		{
			fprintf ( data_unit, "%g\t%g\t%d\n",(double)j*endx*1e9/(nx-1),(double)i*endy*1e9/(ny-1),a[j+i*n]);
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );
	printf ( "\n" );
	printf ( "	Matrix (2D) printed to file '%s'\n", data_filename );
}
/******************************************************************************/

void print2d(int n,double a[],char *data_filename,double fa1,double fa2)

/******************************************************************************/
{
//Normal 2D matrix
	int i,j;
	FILE *data_unit;
	data_unit = fopen ( data_filename, "wt" );
	for ( i = 0; i < ny; i++ )
	{
		for ( j = 0; j < nx; j++ )
		{
			fprintf ( data_unit, "%g\t%g\t%g\n",(double)(j*endx*fa1/(nx-1)),(double)(i*endy*fa1/(ny-1)),a[j+i*n]*fa2);
		}
		fprintf ( data_unit, "\n" );
	}
	fclose ( data_unit );
	printf ( "\n" );
	printf ( "	Matrix (2D) printed to file '%s'\n", data_filename );
}
/******************************************************************************/

void printmatrix(double val[],int col[],int row[],int count)

/******************************************************************************/
{
	char filename[20];
	sprintf(filename, "laplacematrix.dat");
	FILE *file = fopen(filename, "w");
	int i,j;
	for(i=0;i<count;i++)
	{
			fprintf(file,"%0.5f\t%d\t%d\n",val[i],col[i],row[i]);
	}
	fclose(file);
} 
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;
	double cpu_time_used;


	cpu_time_used = ((double) (clock() - times)) / CLOCKS_PER_SEC;
	printf("Time used %f seconds\n",cpu_time_used);

	now = time ( NULL );
	tm = localtime ( &now );

	len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

	fprintf ( stdout, "%s\n", time_buffer );

	return;
# undef TIME_SIZE
}
/******************************************************************************/

/******************************************************************************/
void test()

/******************************************************************************/
{
printf("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n");
}

/******************************************************************************/

void readparameters ( void)

/******************************************************************************/
{
FILE * fr = fopen("PBEparameters.dat", "rt");


if(fr == NULL){printf("file %s not found", "PBEparameters.txt");}

char tmpstr1[16];
char tmpstr2[16];
    char tempbuff[1000];

    while(!feof(fr)) 
    {
         if (fgets(tempbuff,1000,fr)) {

            sscanf(tempbuff, "%15s : %15[^;];", tmpstr1, tmpstr2); //will read a line till it hit a semicolumn
            //printf("%s = \t",  tmpstr1);
            //printf("%s\n",  tmpstr2);

            if (strcmp(tmpstr1,"nx")==0) {
                 nx = atoi(tmpstr2);
				 printf("nx = %d\n",nx);
            } 
            else if (strcmp(tmpstr1,"ny")==0) {
                 ny = atoi(tmpstr2);
				 printf("ny = %d\n",ny);
            }
            else if (strcmp(tmpstr1,"nm")==0) {
                 nm = atoi(tmpstr2);
				 printf("nm = %d\n",nm);
            }
            else if (strcmp(tmpstr1,"startx")==0) {
                 startx = atof(tmpstr2);
				 printf("startx = %5g\n",startx);
            }
            else if (strcmp(tmpstr1,"endx")==0) {
                 endx = atof(tmpstr2);
				 printf("endx = %5g\n",endx);
            }
            else if (strcmp(tmpstr1,"starty")==0) {
                 starty = atof(tmpstr2);
				 printf("starty = %5g\n",starty);
            }
            else if (strcmp(tmpstr1,"endy")==0) {
                 endy = atof(tmpstr2);
				 printf("endy = %5g\n",endy);
            }
            else if (strcmp(tmpstr1,"angle")==0) {
                 angle = atof(tmpstr2);
				 printf("angle = %5g\n",angle);
            }
            else if (strcmp(tmpstr1,"sig")==0) {
                 sigma = atof(tmpstr2);
				 printf("sigma = %5g\n",sigma);
            }
            else if (strcmp(tmpstr1,"kappa")==0) {
                 kappa = atof(tmpstr2);
				 printf("kappa = %5g\n",kappa);
            }
            else if (strcmp(tmpstr1,"thickness")==0) {
                 thickness = atof(tmpstr2);
				 printf("thickness = %5g\n",thickness);
            }
            else if (strcmp(tmpstr1,"hgthickness")==0) {
                 hgthickness = atof(tmpstr2);
				 printf("hgthickness = %5g\n",hgthickness);
            }
            else if (strcmp(tmpstr1,"pthickness")==0) {
                 pthickness = atof(tmpstr2);
				 printf("pthickness = %5g\n",pthickness);
            }
            else if (strcmp(tmpstr1,"radius")==0) {
                 radius = atof(tmpstr2);
				 printf("radius = %5g\n",radius);
            }
            else if (strcmp(tmpstr1,"qp")==0) {
                 qp = atof(tmpstr2);
				 printf("qp = %5g\n",qp);
            }
            else if (strcmp(tmpstr1,"epi")==0) {
                 epi = atof(tmpstr2);
				 printf("epi = %5g\n",epi);
            }
            else if (strcmp(tmpstr1,"epo")==0) {
                 epo = atof(tmpstr2);
				 printf("epo = %5g\n",epo);
            }
            else if (strcmp(tmpstr1,"ewu")==0) {
                 ewu = atof(tmpstr2);
				 printf("ewu = %5g\n",ewu);
            }
            else if (strcmp(tmpstr1,"ewb")==0) {
                 ewb = atof(tmpstr2);
				 printf("ewb = %5g\n",ewb);
            }
            else if (strcmp(tmpstr1,"emt")==0) {
                 emt = atof(tmpstr2);
				 printf("emt = %5g\n",emt);
            }
            else if (strcmp(tmpstr1,"emhb")==0) {
                 emhb = atof(tmpstr2);
				 printf("emhb = %5g\n",emhb);
            }
			else if (strcmp(tmpstr1,"emhu")==0) {
                 emhu = atof(tmpstr2);
				 printf("emhu = %5g\n",emhu);
            }
            else if (strcmp(tmpstr1,"drepeat")==0) {
                 drepeat = atoi(tmpstr2);
				 printf("drepeat = %d\n",drepeat);
            }
            else if (strcmp(tmpstr1,"qrepeat")==0) {
                 qrepeat = atoi(tmpstr2);
				 printf("qrepeat = %d\n",qrepeat);
            }
			else if (strcmp(tmpstr1,"prepeat")==0) {
                 qrepeat = atoi(tmpstr2);
				 printf("qrepeat = %d\n",qrepeat);
            }
            else if (strcmp(tmpstr1,"nbkup")==0) {
                 nbkup = atoi(tmpstr2);
				 printf("nbkup = %d\n",nbkup);
            }
            else if (strcmp(tmpstr1,"zbkup1")==0) {
                 zbkup[0] = atof(tmpstr2);
				 printf("zbkup1 = %5g\n",zbkup[0]);
            }
            else if (strcmp(tmpstr1,"zbkup2")==0) {
                 zbkup[1] = atof(tmpstr2);
				 printf("zbkup2 = %5g\n",zbkup[1]);
            }
			else if (strcmp(tmpstr1,"zbkup3")==0) {
                 zbkup[2] = atof(tmpstr2);
				 printf("zbkup3 = %5g\n",zbkup[2]);
            }
			else if (strcmp(tmpstr1,"zbkup4")==0) {
                 zbkup[3] = atof(tmpstr2);
				 printf("zbkup4 = %5g\n",zbkup[3]);
            }
            else if (strcmp(tmpstr1,"cbkup1")==0) {
                 cbkup[0] = atof(tmpstr2);
				 printf("cbkup1 = %5g\n",cbkup[0]);
            }
            else if (strcmp(tmpstr1,"cbkup2")==0) {
                 cbkup[1] = atof(tmpstr2);
				 printf("cbkup2 = %5g\n",cbkup[1]);
            }
			else if (strcmp(tmpstr1,"cbkup3")==0) {
                 cbkup[2] = atof(tmpstr2);
				 printf("cbkup3 = %5g\n",cbkup[2]);
            }
			else if (strcmp(tmpstr1,"cbkup4")==0) {
                 cbkup[3] = atof(tmpstr2);
				 printf("cbkup4 = %5g\n",cbkup[3]);
            }
            else if (strcmp(tmpstr1,"nbkbt")==0) {
                 nbkbt = atoi(tmpstr2);
				 printf("nbkbt = %d\n",nbkbt);
            }
            else if (strcmp(tmpstr1,"zbkbt1")==0) {
                 zbkbt[0] = atof(tmpstr2);
				 printf("zbkbt1 = %5g\n",zbkbt[0]);
            }
            else if (strcmp(tmpstr1,"zbkbt2")==0) {
                 zbkbt[1] = atof(tmpstr2);
				 printf("zbkbt2 = %5g\n",zbkbt[1]);
            }
			else if (strcmp(tmpstr1,"zbkbt3")==0) {
                 zbkbt[2] = atof(tmpstr2);
				 printf("zbkbt3 = %5g\n",zbkbt[2]);
            }
			else if (strcmp(tmpstr1,"zbkbt4")==0) {
                 zbkbt[3] = atof(tmpstr2);
				 printf("zbkbt4 = %5g\n",zbkbt[3]);
            }
            else if (strcmp(tmpstr1,"cbkbt1")==0) {
                 cbkbt[0] = atof(tmpstr2);
				 printf("cbkbt1 = %5g\n",cbkbt[0]);
            }
            else if (strcmp(tmpstr1,"cbkbt2")==0) {
                 cbkbt[1] = atof(tmpstr2);
				 printf("cbkbt2 = %5g\n",cbkbt[1]);
            }
			else if (strcmp(tmpstr1,"cbkbt3")==0) {
                 cbkbt[2] = atof(tmpstr2);
				 printf("cbkbt3 = %5g\n",cbkbt[2]);
            }
			else if (strcmp(tmpstr1,"cbkbt4")==0) {
                 cbkbt[3] = atof(tmpstr2);
				 printf("cbkbt4 = %5g\n",cbkbt[3]);
            }
            else if (strcmp(tmpstr1,"nmemb")==0) {
                 nmemb = atoi(tmpstr2);
				 printf("nmemb = %d\n",nmemb);
            }
            else if (strcmp(tmpstr1,"zmemb1")==0) {
                 zmemb[0] = atof(tmpstr2);
				 printf("zmemb1 = %5g\n",zmemb[0]);
            }
            else if (strcmp(tmpstr1,"zmemb2")==0) {
                 zmemb[1] = atof(tmpstr2);
				 printf("zmemb2 = %5g\n",zmemb[1]);
            }
			else if (strcmp(tmpstr1,"zmemb3")==0) {
                 zmemb[2] = atof(tmpstr2);
				 printf("zmemb3 = %5g\n",zmemb[2]);
            }
			else if (strcmp(tmpstr1,"zmemb4")==0) {
                 zmemb[3] = atof(tmpstr2);
				 printf("zmemb4 = %5g\n",zmemb[3]);
            }
            else if (strcmp(tmpstr1,"cmemb1")==0) {
                 cmemb[0] = atof(tmpstr2);
				 printf("cmemb1 = %5g\n",cmemb[0]);
            }
            else if (strcmp(tmpstr1,"cmemb2")==0) {
                 cmemb[1] = atof(tmpstr2);
				 printf("cmemb2 = %5g\n",cmemb[1]);
            }
			else if (strcmp(tmpstr1,"cmemb3")==0) {
                 cmemb[2] = atof(tmpstr2);
				 printf("cmemb3 = %5g\n",cmemb[2]);
            }
			else if (strcmp(tmpstr1,"cmemb4")==0) {
                 cmemb[3] = atof(tmpstr2);
				 printf("cmemb4 = %5g\n",cmemb[3]);
            }
			else if (strcmp(tmpstr1,"cmembM1")==0) {
                 cmembM[0] = atof(tmpstr2);
				 printf("cmembM1 = %5g\n",cmembM[0]);
            }
			else if (strcmp(tmpstr1,"cmembM2")==0) {
                 cmembM[1] = atof(tmpstr2);
				 printf("cmembM2 = %5g\n",cmembM[1]);
            }
			else if (strcmp(tmpstr1,"cmembM3")==0) {
                 cmembM[2] = atof(tmpstr2);
				 printf("cmembM3 = %5g\n",cmembM[2]);
            }
			else if (strcmp(tmpstr1,"cmembM4")==0) {
                 cmembM[3] = atof(tmpstr2);
				 printf("cmembM4 = %5g\n",cmembM[3]);
            }
			else if (strcmp(tmpstr1,"nmemu")==0) {
                 nmemu = atoi(tmpstr2);
				 printf("nmemu = %d\n",nmemu);
            }
            else if (strcmp(tmpstr1,"zmemu1")==0) {
                 zmemu[0] = atof(tmpstr2);
				 printf("zmemu1 = %5g\n",zmemu[0]);
            }
            else if (strcmp(tmpstr1,"zmemu2")==0) {
                 zmemu[1] = atof(tmpstr2);
				 printf("zmemu2 = %5g\n",zmemu[1]);
            }
			else if (strcmp(tmpstr1,"zmemu3")==0) {
                 zmemu[2] = atof(tmpstr2);
				 printf("zmemu3 = %5g\n",zmemu[2]);
            }
			else if (strcmp(tmpstr1,"zmemu4")==0) {
                 zmemu[3] = atof(tmpstr2);
				 printf("zmemu4 = %5g\n",zmemu[3]);
            }
            else if (strcmp(tmpstr1,"cmemu1")==0) {
                 cmemu[0] = atof(tmpstr2);
				 printf("cmemu1 = %5g\n",cmemu[0]);
            }
            else if (strcmp(tmpstr1,"cmemu2")==0) {
                 cmemu[1] = atof(tmpstr2);
				 printf("cmemu2 = %5g\n",cmemu[1]);
            }
			else if (strcmp(tmpstr1,"cmemu3")==0) {
                 cmemu[2] = atof(tmpstr2);
				 printf("cmemu3 = %5g\n",cmemu[2]);
            }
			else if (strcmp(tmpstr1,"cmemu4")==0) {
                 cmemu[3] = atof(tmpstr2);
				 printf("cmemu4 = %5g\n",cmemu[3]);
            }
			else if (strcmp(tmpstr1,"cmemuM1")==0) {
                 cmemuM[0] = atof(tmpstr2);
				 printf("cmemuM1 = %5g\n",cmemuM[0]);
            }
			else if (strcmp(tmpstr1,"cmemuM2")==0) {
                 cmemuM[1] = atof(tmpstr2);
				 printf("cmemuM2 = %5g\n",cmemuM[1]);
            }
			else if (strcmp(tmpstr1,"cmemuM3")==0) {
                 cmemuM[2] = atof(tmpstr2);
				 printf("cmemuM3 = %5g\n",cmemuM[2]);
            }
			else if (strcmp(tmpstr1,"cmemuM4")==0) {
                 cmemuM[3] = atof(tmpstr2);
				 printf("cmemuM4 = %5g\n",cmemuM[3]);
            }
			else if (strcmp(tmpstr1,"T")==0) {
                 T = atof(tmpstr2);
				 printf("T = %5g\n",T);
            }
            else if (strcmp(tmpstr1,"id0")==0) {
                 id0 = atoi(tmpstr2);
				 printf("id0 = %d\n",id0);
            }
            else if (strcmp(tmpstr1,"id1")==0) {
                 id1 = atoi(tmpstr2);
				 printf("id1 = %d\n",id1);
            }
            else if (strcmp(tmpstr1,"id2")==0) {
                 id2 = atoi(tmpstr2);
				 printf("id2 = %d\n",id2);
            }
            else if (strcmp(tmpstr1,"id3")==0) {
                 id3 = atoi(tmpstr2);
				 printf("id3 = %d\n",id3);
            }
			else if (strcmp(tmpstr1,"id4")==0) {
                 id4 = atoi(tmpstr2);
				 printf("id4 = %d\n",id4);
            }
			else if (strcmp(tmpstr1,"id5")==0) {
                 id5 = atoi(tmpstr2);
				 printf("id5 = %d\n",id5);
            }
			else if (strcmp(tmpstr1,"id6")==0) {
                 id6 = atoi(tmpstr2);
				 printf("id6 = %d\n",id6);
            }
			else if (strcmp(tmpstr1,"id7")==0) {
                 id7 = atoi(tmpstr2);
				 printf("id7 = %d\n",id7);
            }
			else if (strcmp(tmpstr1,"isMem")==0) {
                 isMem = atoi(tmpstr2);
				 printf("isMem = %d\n",isMem);
            }
			else if (strcmp(tmpstr1,"isNP")==0) {
                 isNP = atoi(tmpstr2);
				 printf("isNP = %d\n",isNP);
            }
			else if (strcmp(tmpstr1,"isPlot")==0) {
                 isPlot = atoi(tmpstr2);
				 printf("isPlot = %d\n",isPlot);
            }
			else if (strcmp(tmpstr1,"isFd2d")==0) {
                 isFd2d = atoi(tmpstr2);
				 printf("isFd2d = %d\n",isFd2d);
            }
			else if (strcmp(tmpstr1,"isclean")==0) {
                 isclean = atoi(tmpstr2);
				 printf("isclean = %d\n",isclean);
            }
			else if (strcmp(tmpstr1,"isPrint")==0) {
                 isPrint = atoi(tmpstr2);
				 printf("isPrint = %d\n",isPrint);
            }
			else if (strcmp(tmpstr1,"pot2d")==0) {
                 pot2d = atoi(tmpstr2);
				 printf("pot2d = %d\n",pot2d);
            }
			else if (strcmp(tmpstr1,"memfactor")==0) {
                 memfactor = atof(tmpstr2);
				 printf("memfactor = %5g\n",memfactor);
            }
			else if (strcmp(tmpstr1,"memshift")==0) {
                 memshift = atof(tmpstr2);
				 printf("memshift = %5g\n",memshift);
            }
			else if (strcmp(tmpstr1,"parshift")==0) {
                 parshift = atof(tmpstr2);
				 printf("parshift = %5g\n",parshift);
            }
            else{
                printf("Unrecongized parameter : \"%s\"\n", tmpstr1);
                free (tmpstr1);
				free (tmpstr2);
				free (tempbuff);
                fclose(fr);
         }
		}
	}
//free (tmpstr1);
//free (tmpstr2);
//free (tempbuff);

}

