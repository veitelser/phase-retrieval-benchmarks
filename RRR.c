/*
-----------------------------------
RRR: Simple Phase Retrieval
-----------------------------------

This program is deliberately minimalist so as not to obscure
the structure of the algorithm. It needs the FFTW3 library:

	http://www.fftw.org

The program is set up for solving the benchmarks described in:

	"Benchmark problems for phase retrieval", V. Elser, T.-Y. Lan & T. Bendory

To compile:

	gcc -O2 RRR.c -lm -lfftw3 -o RRR

To run:

	./RRR [datafile] [supp] [powgoal] [beta] [iterlimit] [trials] [resultfile] &

datafile:	one of the benchmark datafiles (data100E, data140E, ...)
supp:		support size = 8*N, N = 100, 140, ... is the number of atoms
powgoal:	fractional power in support
beta:		RRR parameter
iterlimit:	RRR iteration limit (long int)
trials:		number of random starts
resultfile:	ASCII file of the iteration counts for each trial

Example:

	./RRR data/data100E 800 .95 .5 1000 5 results100E &

Solutions are written to a file named sol (M x M table of floats).

*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

// 2D grid: M x M
#define M 128

double x[M*M],x1[M*M],x2[M*M],v[M*M],fmag[M*(M/2+1)];
double complex f2[M*(M/2+1)];
double beta,solf0,solpow;
int supp,datapow;

fftw_plan forward_plan,backward_plan;


// FFTW setup for 2D transforms
void setup()
{
forward_plan=fftw_plan_dft_r2c_2d(M,M,x2,f2,FFTW_MEASURE);
backward_plan=fftw_plan_dft_c2r_2d(M,M,f2,x2,FFTW_MEASURE);
}

// read data, set Fourier magnitudes, compute Fourier power
// first column of data should already be symmetrized
// normalization of magnitudes compensates for FFTW's missing normalization
void getdata(char *datafile)
{
FILE *fp;
int i,j,k,data;

fp=fopen(datafile,"r");

datapow=0;
k=0;

for(i=0;i<M;++i)
	{
	fscanf(fp,"%d",&data);
	datapow+=data;

	fmag[k++]=sqrt((double)data)/M;

	for(j=1;j<M/2;++j)
		{
		fscanf(fp,"%d",&data);
		datapow+=2*data;

		fmag[k++]=sqrt((double)data)/M;
		}

	fmag[k++]=.0;
	}

fclose(fp);
}

// quick-select algorithm
// recursively finds value in array 'val' below element at a given 'rank'
double qselect(double *val, int len, int rank)
{
int p,q;
double tmp;

# define SWAP(p,q){tmp=val[p]; val[p]=val[q]; val[q]=tmp;}

for(p=q=0;p<len-1;++p)
	{
	if(val[p]<val[len-1])
		continue;

	SWAP(p,q);
	++q;
	}

SWAP(q,len-1);

return q==rank ? val[q] : q>rank ? qselect(val,q,rank) : qselect(val+q,len-q,rank-q);
}

// support-size projection: xin -> x1
// uses quick-select with rank = supp
void proj1(double *xin)
{
int k;
double low;

for(k=0;k<M*M;++k)
	v[k]=xin[k];

low=qselect(v,M*M,supp);

// also set to zero any negative values in support
if(low<.0)
	low=.0;

for(k=0;k<M*M;++k)
	{
	if(xin[k]<=low)
		x1[k]=.0;
	else
		x1[k]=xin[k];
	}
}

// Fourier magnitude projection: x2 -> f2 -> x2
void proj2()
{
int k;
double f2mag;

// x2 -> f2
fftw_execute(forward_plan);

// project to f2[0]>=0 (x2 is positive) and normalize
if(creal(f2[0])<.0)
	f2[0]=.0;
else
	f2[0]/=M*M;

// rescale Fourier magnitudes and normalize
for(k=1;k<M*(M/2+1);++k)
	{
	f2mag=cabs(f2[k]);
	if(f2mag==.0)
		f2[k]=fmag[k]/M;
	else
		f2[k]*=fmag[k]/f2mag;
	}

// f2 -> x2
fftw_execute(backward_plan);
}

// relaxed-reflect-reflect (RRR) iteration
void RRR()
{
int k;

proj1(x);

for(k=0;k<M*M;++k)
	x2[k]=2.*x1[k]-x[k];

proj2();

// solf0 is the [0] Fourier coefficient of the candidate solution
// solpow is the Fourier power of the candidate solution
solf0=.0;
solpow=.0;

for(k=0;k<M*M;++k)
	{
	x[k]+=beta*(x2[k]-x1[k]);

	// support of candidate solution x2 is defined by positive values of x1
	// x2 is synthesized from data magnitudes
	if(x1[k]>.0)
		{
		solf0+=x2[k];
		solpow+=x2[k]*x2[k];
		}
	}

solf0/=M;
}

// initialize x
void init()
{
int k;

// start with positive random numbers
for(k=0;k<M*M;++k)
	x2[k]=((double) rand())/RAND_MAX;

// scale with data magnitudes
proj2();

for(k=0;k<M*M;++k)
	x[k]=x2[k];
}

// write x2 to the file 'sol'
void printsol()
{
int i,j;
FILE *fp;

fp=fopen("sol","w");

for(i=0;i<M;++i)
	{
	for(j=0;j<M;++j)
		fprintf(fp,"%12.6f",x2[i*M+j]);
	fprintf(fp,"\n");
	}

fclose(fp);
}


int main(int argc,char* argv[])
{
char *datafile,*resultfile;
int try,trials,succ;
long long iter,iterlimit,itercount;
double pow,powgoal;
FILE *fp;
clock_t start;
double elapsed;

if(argc==8)
	{
	datafile=argv[1];
	supp=atoi(argv[2]);
	powgoal=atof(argv[3]);
	beta=atof(argv[4]);
	iterlimit=atol(argv[5]);
	trials=atoi(argv[6]);
	resultfile=argv[7];
	}
else
	{
	fprintf(stderr,"expected seven arguments: datafile, supp, powgoal, beta, iterlimit, trials, resultfile\n");
	return 1;
	}

getdata(datafile);
setup();

fp=fopen(resultfile,"w");
fprintf(fp,"datafile: %s\n",datafile);
fprintf(fp,"support: %d  power goal: %f  beta: %f\n\n",supp,powgoal,beta);
fprintf(fp,"trial    iterations\n");
fclose(fp);

succ=0;
itercount=0;

srand(time(0));
start=clock();

for(try=1;try<=trials;++try)
	{
	// randomly initialize for each trial
	init();

	for(iter=0;iter<iterlimit;++iter)
		{
		RRR();

		// solf0 and solpow are computed in RRR()
		pow=solpow/(solf0*solf0+datapow);

		// terminate iterations and record results when power criterion is satisfied
		if(pow>powgoal)
			{
			printsol();

			fp=fopen(resultfile,"a");
			fprintf(fp,"%5d%14lld\n",try,iter);
			fclose(fp);

			++succ;
			itercount+=iter;

			goto next;
			}
		}

	fp=fopen(resultfile,"a");
	fprintf(fp,"%5d%14d\n",try,-1);
	fclose(fp);

	itercount+=iterlimit;

	next:;
	}

elapsed=((double)(clock()-start))/CLOCKS_PER_SEC;

fp=fopen(resultfile,"a");

if(succ)
	fprintf(fp,"\niterations per solution: %f\n",((double)itercount)/succ);
else
	fprintf(fp,"\nno solutions found\n");

fprintf(fp,"\niterations per sec: %f\n",((double)itercount)/elapsed);

fclose(fp);

return 0;
}
