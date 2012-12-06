/**********************************************************************
 * gillespie_ssim_sparse.cpp
 *
 *  a MEX file implementing the Gillespie stochastic solver method in 
 *  MATLAB optimized for sparse reaction matricies
 * 
 *  Gillespie D (1977) Exact Stochastic Simulation of Coupled Chemical 
 *  Reactions.The Journal of Physical Chemistry 81(25):2340–2361.
 *
 *  Usage: 
 *  result = gillespie_ssim_sparse(Reactants,Products,Rates,InitCon,dt,tMax)
 *
 *  Reactants & Products are M x N matricies with stoichiometries for 
 *  all M reactions and N species. Species not used in a given reaction
 *  must be set to zero.
 *
 *  Rates is a column vector containing M reaction propensities
 *
 *  InitCon is a column vector containing initial numbers of all 
 *  N particle species
 *
 *  tMax is the maximum simulated time 
 *
 *  dt is the time step for recording simulation state  
 * 
 *  result is a tMax/dt x (N+1) matrix of the species counts 
 *  at every time dt. The first column gives the specific times 
 *  at which each row was recorded.
 *
 *  
 *  Currently requires MersenneTwister.h from 
 *  http://www-personal.umich.edu/~wagnerr/MersenneTwister.html
 *  to compile.
 *
 *  
 *  Brian Kolterman
 *  11/2010
 **********************************************************************/

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "MersenneTwister.h"

// input arguments

#define REAC        prhs[0]
#define PROD        prhs[1]
#define RATE        prhs[2]
#define INI_CON     prhs[3]
#define DT          prhs[4]
#define TMAX        prhs[5]

// output matrix
 
#define RESULT      plhs[0]

void NchooseK(double &n, double &k, double &result);
double factorial(double n);

       
static MTRand mtrand; 


void gillespie_ssim(double rea[], double pro[], double rat[], double ini[], double dt, double tmax, double out[], mwSize s, mwSize r, mwSize tsteps)

{
    
    // init sim
    
    int i, j, mu, count, countP, iStep, nStep, last_time, iRe, iSp;
    double tau, randN, test, temp[10], *a_cumsum, *a_stepsum, *conc;   
    double *a, *h;      // reaction probabilities
    double t = 0.0;       // simulation time
    
    conc = (double*)mxCalloc(s,sizeof(double));
    a = (double*)mxCalloc(r,sizeof(double));
    a_cumsum = (double*)mxCalloc(r,sizeof(double));
    
    nStep = (int)ceil(sqrt((double)r));                  // used to reduce mu 
    a_stepsum = (double*)mxCalloc(nStep,sizeof(double)); // search space
    
    h = (double*)mxCalloc(r,sizeof(double));
  
    // Create sparse index for reactions
    
    int (*indSpar)[10] = new int[r][10];

    int (*indSparPro)[10] = new int[r][10];
    
    for (iRe = 0; iRe < r; iRe++)   {
        
        count = 0;
        countP = 0;
        
        for (iSp = 0; iSp < s; iSp++)   {
            
            j = iSp*r + iRe;
            
            if (rea[j] != 0)    {
                
                count++;
                indSpar[iRe][count] = iSp;
            }
            
            if (pro[j] != 0)    {
                
                countP++;
                indSparPro[iRe][countP] = iSp;
            }
        }
        
        indSpar[iRe][0] = count;
        indSparPro[iRe][0] = countP;
    }
            
              
    out[0] = t;
    last_time = 0;
     
    for (i = 0; i < s; i++) {
        
        j = (i+1)*tsteps;
        out[j] = ini[i];
        conc[i] = ini[i];
    }

    
    // start simulation loop 
    
    while (t < tmax) {
        
        iStep = 0; 
        // check if all species have been used up and halt if so
        
        temp[0] = 0;
        
        for (iSp = 0; iSp < s; iSp++)   {
            
            temp[0] += conc[iSp];
            
        
        if (temp[0] == 0)   {
        
            printf("All species extinguished\n");
            
            break;
        
        }
        
        // find reaction probabilities
        
        for (iRe = 0; iRe < r; iRe++) {
            
            count = 0;
            temp[count] = 0;
            
            for (iSp = 1; iSp <= indSpar[iRe][0]; iSp++) {
                
                // 2D MATLAB matrix index in column format
                j = indSpar[iRe][iSp]*r + iRe;
                if (rea[j] == 0) continue;
                
                if (conc[indSpar[iRe][iSp]] == 0) {
                    temp[0] = 0;    // make sure this reaction can't go
                    break;
                }
                
                if (rea[j] == 1) { 
                    
                    temp[count] = conc[indSpar[iRe][iSp]];
                    
                }   else    {
                    
                    NchooseK(conc[indSpar[iRe][iSp]],rea[j],temp[count]);
                    
                }
                
                count++;
                
            }
            
            
            //printf("temp = %g\n",temp[0]);
            h[iRe] = temp[0];
            
            for (i = 1; i < count; i++) {
                
                h[iRe] *= temp[i];
                 
            }
        
           
            a[iRe] = rat[iRe]*h[iRe];
           
            a_cumsum[iRe] = a[iRe];
            
            if (iRe != 0)    {
                
                a_cumsum[iRe] += a_cumsum[iRe-1];
           
           }
        
        // Get reaction time (tau) and number (mu)
        
        randN = mtrand.randDblExc();
        
        tau = (1/a_cumsum[r-1])*log((1/randN));
        
        randN = mtrand.randDblExc();
        test = randN*a_cumsum[r-1];
        
        for (mu = 0; mu < r ; mu++) {
            
           if (test <= a_cumsum[mu]) break;
          
        }
        
        // increase time and adjust species population levels
        
        t += tau;
        
        for (iSp = 1; iSp <= indSpar[mu][0]; iSp++)    {
            
            j = indSpar[mu][iSp]*r + mu;
            conc[indSpar[mu][iSp]] = conc[indSpar[mu][iSp]] - rea[j] + pro[j];

        }
        
        for (iSp = 1; iSp <= indSparPro[mu][0]; iSp++)    {
            
            j = indSparPro[mu][iSp]*r + mu;
            conc[indSparPro[mu][iSp]] = conc[indSparPro[mu][iSp]] - rea[j] + pro[j];
        }
        
        // record sim state when dt has elapsed
        
        if ((t - out[last_time]) > dt)  {
            
            last_time++;
            
            for (i = 1; i <= s; i++) {
                
                j = i*tsteps + last_time;
                out[j] = conc[i-1];
                out[last_time] = t;
            }
        
        }
        
    }   // end sim loop
    
    // record final state if not done so
    
    if ((t - out[last_time]) != 0)  {
        
        last_time++;
            
        for (i = 1; i <= s; i++) {
                
            j = i*tsteps + last_time;
            out[j] = conc[i-1];
            out[last_time] = t;
            
        }    
    }
    
    
    delete [] indSpar;
    delete [] indSparPro;
    mxFree(a);
    mxFree(a_cumsum);
    mxFree(h);
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  
    int i, j, ind, ind2;
    mwSize nSpecies, nReactions, nTsteps, iT;
    double tMax, dt, temp;
    double *rea, *pro, *rat, *ini, *out, *res;
    mxArray *OUT;
    
    // seed the rng with N = 624 random uint32s

    MTRand mtrands;
    MTRand::uint32 seed[ MTRand::N ];
    for( int n = 0; n < MTRand::N; ++n ) seed[n] = 23 * mtrands.randInt(); 

    mtrand.seed(seed);
    
    dt = mxGetScalar(DT);
    tMax = mxGetScalar(TMAX);
    
    nSpecies = (mwSize)mxGetN(REAC);
    nReactions = (mwSize)mxGetM(REAC);
	nTsteps = (mwSize)tMax/dt;
    

    // create output matrix for use in simulation routine
    
    OUT = mxCreateDoubleMatrix(nTsteps, (nSpecies+1), mxREAL); 
    
    
    
    // Assign pointers  
    
    rea = mxGetPr(REAC);
    pro = mxGetPr(PROD);
    rat = mxGetPr(RATE);
    ini = mxGetPr(INI_CON);
    out = mxGetPr(OUT);
    
    
    
    // run sim
    
    gillespie_ssim(rea,pro,rat,ini,dt,tMax,out,nSpecies,nReactions,nTsteps);
	
    
    
    
    
    // reformat output to remove unused rows
    
    for (iT = 1; iT < nTsteps; iT++)   {
        
        if (out[iT] == 0) break;
        
    }
   
    RESULT = mxCreateDoubleMatrix(iT, (nSpecies+1), mxREAL);
    res = mxGetPr(RESULT);
    
    for (j = 0; j <= nSpecies; j++) {
    
        for (i = 0; i < iT; i++)    {
            
            ind = j*iT + i;
            ind2 = j*nTsteps + i;
            res[ind] = out[ind2];
            
        }
    }
    
    mxDestroyArray(OUT);
    
    return;

}


void NchooseK(double &n, double &k, double &result) {
        
    result = factorial(n)/(factorial(n-k)*factorial(k));
    
}

double factorial(double n)  {
    
    if (n == 0) return 1;
    
    return factorial(n-1)*n;
        
        
}