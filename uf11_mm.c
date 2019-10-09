/****************************************************************************/
// UF11 test probem from cec2009 competition. 
// this function is a rotated and extended version of DTLZ2
// For detailed formulation see page 13 of: 
// Zhang, Q., Zhou, A., Zhao, S., Suganthan, P. N., Liu, W., & Tiwari, S.
// (2008). Multiobjective optimization test instances for the CEC 2009 
// special session and competition. University of Essex, Colchester, UK 
// and Nanyang technological University, Singapore, special session on
// performance assessment of multi-objective optimization algorithms, 
// technical report, 264.
/****************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "borgmm.h"

#define PI 3.14159265358979323846
#define NVARS 30
#define NOBJS 5

int nvars = NVARS;
int nobjs = NOBJS;

void uf11(double* vars, double* objs, double* consts)
{	
	
	int i=0,j=0;    
	int k = nvars - nobjs + 1;
	double g = 0;
	double sum = 0;
	//double z[nvars],zz[nvars],p[nvars],psum[nobjs],M[nvars][nvars],lamda_l[nvars];
	double *z, *zz, *p, *psum, **M, *lamda_l;
	double M_10D[10][10]={{0.0346, -0.7523, 0.3561, -0.2958, 0.4675,0,0,0,0,0},{0.8159, -0.0423, 0.4063, 0.3455, -0.2192,0,0,0,0,0},{-0.3499, 0.3421, 0.8227, -0.2190, -0.1889,0,0,0,0,0},{-0.0963, -0.4747, -0.0998, -0.2429, -0.8345,0,0,0,0,0},{-0.4487, -0.2998, 0.1460, 0.8283, -0.0363,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1}};
	double lamda_l_10D[10]={0.313,0.312,0.321,0.316,0.456,1,1,1,1,1};
	
	double M_30D[30][30]={{0.0128,0.2165,0.4374,-0.0800,0.0886,-0.2015,0.1071,0.2886,0.2354,0.2785,-0.1748,0.2147,0.1649,-0.3043,0.5316,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{0.4813,   0.2420,    -0.3663,   -0.0420,   -0.0088,   -0.4945,   -0.3073, 0.1990, 0.0441, -0.0627, 0.0191, 0.3880, -0.0618, -0.0319, -0.1833,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{0.4816,   -0.2254,    0.0663,    0.4801,    0.2009,   -0.0008,   -0.1501,    0.0269,   -0.2037,0.4334,   -0.2157,   -0.3175,   -0.0923,    0.1451,    0.1118,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{-0.0876,   -0.2667,   -0.0063,    0.2114,    0.4506,    0.0823,   -0.0125,    0.2313,    0.0840,-0.2376,    0.1938,   -0.0030,    0.3391,    0.0863,    0.1231,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{-0.1025,    0.4011,   -0.0117,    0.2076,    0.2585,    0.1124,   -0.0288,    0.3095,   -0.6146,-0.2376,    0.1938,   -0.0030,    0.3391,    0.0863,    0.1231,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{0.4543,   -0.2761,   -0.2985,   -0.2837,    0.0634,    0.1070,    0.2996,   -0.2690,   -0.1634,-0.1452,    0.1799,   -0.0014,    0.2394,   -0.2745,    0.3969,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{-0.1422,   -0.4364,    0.0751,   -0.2235,    0.3966,   -0.0252,    0.0908,    0.0477,   -0.2254,0.1801,   -0.0552,    0.5770,   -0.0396,    0.3765,   -0.0522,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{0.3542,   -0.2245,    0.3497,   -0.1609,   -0.1107,    0.0079,    0.2241,    0.4517,    0.1309,-0.3355,   -0.1123,   -0.1831,    0.3000,    0.2045,   -0.3191,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{0.0005,    0.0377,   -0.2808,   -0.0641,    0.1316,    0.2191,    0.0207,    0.3308,    0.4117,0.3839,    0.5775,   -0.1219,    0.1192,    0.2435,    0.0414,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{-0.1177,   -0.0001,   -0.1992,   -0.4533,    0.4234,   -0.0191,   -0.3740,    0.1325,    0.0972,-0.2042,   -0.3493,   -0.4018,   -0.1087,    0.0918,    0.2217,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{0.1818,    0.3022,   -0.1388,   -0.2380,   -0.0773,    0.6463,    0.0450,    0.1030,   -0.0958,0.2837,   -0.3969,    0.1779,   -0.0251,   -0.1543,   -0.2452,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000,0.000},
	{-0.1889,   -0.4397,   -0.2206,    0.0981,   -0.5203,    0.1325,   -0.3427,    0.4242,   -0.1271,-0.0291,   -0.0795,    0.1213,    0.0565,   -0.1092,    0.2720,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{-0.1808,   -0.0624,   -0.2689,    0.2289,    0.1128,   -0.0844,   -0.0549,   -0.2202,    0.2450,0.0825,   -0.3319,    0.0513,    0.7523,    0.0043,   -0.1472,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{-0.0983,    0.0611,   -0.4145,    0.3017,    0.0410,   -0.0703,    0.6250,    0.2449,    0.1307,-0.1714,   -0.3045,    0.0218,   -0.2837,    0.1408,    0.1633,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0.2026,    0.0324,    0.1496,    0.3129,    0.1437,    0.4331,   -0.2629,   -0.1498,    0.3746,-0.4366,    0.0163,    0.3316,   -0.0697,    0.1833,    0.2412,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}};
	
    double lamda_l_30D[30]={0.113,0.105,0.117,0.119,0.108,0.110,0.101,0.107,0.111,0.109,0.120,0.108,0.101,0.105,0.116,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000};
	
	z       = new double[nvars];
	zz      = new double[nvars];
	p       = new double[nvars];
	psum    = new double[nobjs];
	M       = new double*[nvars]; for(i=0; i<nvars; i++) M[i] = new double[nvars];
	lamda_l = new double[nvars];
	
	if (nvars==10)
	{
		for (i=0;i<nvars;i++)
		{
			for (j=0;j<nvars;j++)
			{
				M[i][j]=M_10D[i][j];	
			}
			lamda_l[i]=lamda_l_10D[i];
		}
	}
	else 
	{
		for (i=0;i<nvars;i++)
		{
			for (j=0;j<nvars;j++)
			{
				M[i][j]=M_30D[i][j];	
			}
			lamda_l[i]=lamda_l_30D[i];
		}
	}
	
	for (i=0;i<nvars;i++)
	{
		z[i]=0;
		for (j=0;j<nvars;j++)
		{
			z[i]+=M[i][j]*vars[j];
		}		
		if (z[i]>=0 && z[i]<=1)
		{
			zz[i]=z[i];
			p[i]=0;
		}
		else if (z[i]<0)
		{			
			zz[i]=-lamda_l[i] * z[i];	
			p[i]=-z[i];
		}
		else
		{
			zz[i]=1-lamda_l[i]*(z[i]-1);
			p[i]=z[i]-1;		
		}
	}	
	for(j=0;j<nobjs;j++)
	{
		psum[j] = 0;
	}
	
	for (i = nvars - k + 1; i <= nvars; i++)
	{
		g += pow(zz[i-1]-0.5,2);
		for(j=0;j<nobjs;j++)
		{
			psum[j]= sqrt( pow(psum[j],2) + pow(p[i-1],2) );
		}
	}
	
	for (i = 1; i <= nobjs; i++)
	{
		double ff = (1 + g);
		for (j = nobjs - i; j >= 1; j--)
		{
			ff *= cos(zz[j-1] * PI / 2.0);
			psum[i-1] = sqrt( pow(psum[i-1],2) + pow(p[j-1],2) );
		}
		if (i > 1)
		{
			ff *= sin(zz[(nobjs - i + 1) - 1] * PI / 2.0);
			psum[i-1] = sqrt( pow(psum[i-1],2) + pow(p[(nobjs - i + 1) - 1],2) );
		}
		
		objs[i-1] = 2.0/(1+exp(-psum[i-1])) * (ff+1);
	}

	delete []z; delete []zz; delete []p;delete []psum; delete []lamda_l;
	for(i=0; i<nvars; i++) delete []M[i]; delete []M;

	// sleep for .001 seconds (units of usleep are microseconds)
	usleep(1000);
}


int main(int argc, char* argv[]) {
	 // setting random seed
  	unsigned int seed = atoi(argv[1]);
	int j;
	int rank;
	int NFE = 1000000;
	char outputFilename[256];
	FILE* outputFile = NULL;
	char runtime[256];
	char timing[256];

	// All multi-master runs need to call startup, specify the number of
	// islands (one master per island), and set the runtime limits.
	BORG_Algorithm_ms_startup(&argc, &argv);
	BORG_Algorithm_ms_islands(1);
	BORG_Algorithm_ms_max_time(.125);
	BORG_Algorithm_ms_max_evaluations(NFE); 
	// Enable global Latin hypercube initialization to ensure each island
	// gets a well sampled distribution of solutions.
	BORG_Algorithm_ms_initialization(INITIALIZATION_LATIN_GLOBAL);

	// Define the problem.  Problems are defined the same way as the
	// serial example (see dtlz2_serial.c).
	BORG_Problem problem = BORG_Problem_create(nvars, nobjs, 0, uf11);

	double lower[NVARS] = {1.773,-1.846,-1.053,-2.370,-1.603,-1.878,-1.677,-0.935,-1.891,-0.964,-0.885,-1.690,-2.235,-1.541,-0.720,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000};

	double upper[NVARS] = {1.403,1.562,2.009,0.976,1.490,1.334,1.074,2.354,1.462,2.372,2.267,1.309,0.842,1.665,2.476,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000};

	for (j=0; j<nvars; j++) {
			BORG_Problem_set_bounds(problem, j, lower[j], upper[j]);
	}

	for (j=0; j<nobjs; j++) {
			BORG_Problem_set_epsilon(problem, j, 0.01);
	}

	
	// Save runtime dynamics to a file.  Each master node will
	// write its local runtime dynamics to this file.  The %%d
	// gets replaced by the index of the master.
	BORG_Algorithm_output_frequency((int)NFE/20);
	sprintf(outputFilename, "./sets/uf11_S%d.set", seed); 
	sprintf(runtime, "./runtime/uf11/uf11_S%d_M%%d.runtime", seed);
	BORG_Algorithm_output_runtime(runtime);

	// Get the rank of this process.  The rank is used to ensure each
	// parallel process uses a different random seed.
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// When running experiments, we want to run the algorithm multiple
	// times and average the results.
    BORG_Random_seed(37*seed*(rank+1));


	// Run the multi-master Borg MOEA on the problem.
	BORG_Archive result = BORG_Algorithm_ms_run(problem);

	// Only the controller process will return a non-NULL result.
	// The controller aggregates all of the Pareto optimal
	// solutions generated by each master.  Then print the Pareto
	// optimal solutions to the screen.
	if (result != NULL) {
		outputFile = fopen(outputFilename, "w");
		if (!outputFile) {
			BORG_Debug("Unable to open final output file\n");
		}
		BORG_Archive_print(result, outputFile);
		BORG_Archive_destroy(result);
		fclose(outputFile);
	}
	// Shutdown the parallel processes and exit.
	BORG_Algorithm_ms_shutdown();
	BORG_Problem_destroy(problem);
	return EXIT_SUCCESS;
}
