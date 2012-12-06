/**********************************************************************
 * gillespie_ssim_ORexpModel.cpp
 *
 *  a MEX file simulating the OR model in arxiv.org/abs/1201.2933
 *  in MATLAB based on the Gillespie stochastic solver method
 *  
 *
 *  Usage: 
 *  result = gillespie_ssim_ORmod(nORs,nPORs,rTFprod,rTFkill,rTFbasal, ...
 *              rDnaKill,rTFass,rTFdiss,alpha,dt,tMax,nExp,nMax)
 *
 *  rDiss(n) = rTFdiss*exp(-alpha(n-1))
 *
 *  nORs = # functional OR genes
 *  nPORs = # pseudogenes
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

#define NORS        prhs[0]
#define NPOR        prhs[1]
#define TFPROD      prhs[2]
#define TFKILL      prhs[3]
#define TFBASAL     prhs[4]
#define DNAKILL     prhs[5]
#define TFASS       prhs[6]
#define TFDISS      prhs[7]
#define ALPHA       prhs[8]
#define DT          prhs[9]
#define TMAX        prhs[10]
#define NEXP        prhs[11]
#define NMAX        prhs[12]

// output matrix
 
#define RESULT      plhs[0]

void NchooseK(double &n, double &k, double &result);
double factorial(double n);

       
static MTRand mtrand; 


void gillespie_ssim(double nexp, double nmax, double tfp, double tfk, double tfb, double dnak, double ass, double diss, double alpha, double dt, double tmax, double out[], mwSize nors, mwSize npor, mwSize tsteps)

{
    
    // init sim
    
    int i, j, mu, count, last_time, iRe, iSp, iTF;
    double tau, randN, test, rbas, *a_cumsum, *conc, *ratDiss;   
    double *a;      // reaction probabilities
    double t = 0.0;       // simulation time
    mwSize r;
    
    r = (2*(nors+npor)) + 2;   // number of possible reactions 
    iTF = nors + npor;
    
    conc = (double*)mxCalloc((nors+npor+1),sizeof(double));
    a = (double*)mxCalloc(r,sizeof(double));
    a_cumsum = (double*)mxCalloc(r,sizeof(double));
    ratDiss = (double*)mxCalloc((mwSize)(nmax+1),sizeof(double));
    
    for (i = 0; i <= nmax; i++)  {
           
        ratDiss[i] = diss*exp((-1.0)*alpha*(double)i);
    }
  
    
    rbas = tfp/tfb;
            
               
    out[0] = t;
    last_time = 0;
     
    for (i = 0; i <= npor; i++) {
        
        j = (i+1)*tsteps;
        out[j] = 0;
        conc[i] = 0;
    }

    
    // start simulation loop 
    
     while (t < tmax) {

        
        // find reaction probabilities for OR gain and loss
        
        iRe = 0;
        count = 0;
        
        for (iSp = 0; iSp < nors; iSp++) {
            
            if (conc[iSp] < nmax) {
                
                a[iRe] = ass*conc[iTF];
          
            } else {
                
                a[iRe] = 0;
            }
            
            iRe++;
            
            if (conc[iSp] !=0)  {
                
                a[iRe] = ratDiss[(int)conc[iSp]];
                
            }   else    {
                
                a[iRe] = 0;
            }
            
            iRe++;
            
            if (conc[iSp] >= nexp)  {
                
                count++;
            }
            
        }
        
        // pseudogenes (don't contribute to number of expressed genes)
        
        for (iSp = nors; iSp < (npor+nors); iSp++) {
            
            if (conc[iSp] < nmax) {
                
                a[iRe] = ass*conc[iTF];
          
            } else {
                
                a[iRe] = 0;
            }
            
            iRe++;
            
            if (conc[iSp] !=0)  {
                
                a[iRe] = ratDiss[(int)conc[iSp]];
                
            }   else    {
                
                a[iRe] = 0;
            }
            
            iRe++;
            
            
            
        }
        
        a[iRe] = tfp*(1.0/(pow(rbas,count))); // TF production 
        
        iRe++;
        
        a[iRe] = conc[iTF]*tfk; //TF killin
        
        
        
        a_cumsum[0] = a[0];
        
        for (iRe = 1; iRe < r; iRe++)   {
            
            a_cumsum[iRe] = a[iRe];
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
        
        
        iSp = (int)mu/2;
        
       // printf("mu= %i  tau= %g  cumsum= %g  iSp = %i\n",mu,tau,a_cumsum[r-1],iSp); 
        
        if (iSp < iTF)  {
            
            if ((mu % 2) == 0)     {
                
                conc[iSp]++;
                conc[iTF]--;
                
            }   else {
                
                conc[iSp]--;
                conc[iTF]++;
            }
        }   else    {
            
            if ((mu % 2) == 0) {
                
                conc[iTF]++;
            
            }   else    {
                
                conc[iTF]--;
            }
            
        }
        
        
        
        // record sim state when dt has elapsed
        
        if ((t - out[last_time]) > dt)  {
            
            last_time++;
            
            for (i = 1; i <= (nors+npor); i++) {
                
                j = i*tsteps + last_time;
                out[j] = conc[i-1];
                out[last_time] = t;
            }
        
        }
        
    }   // end sim loop
    
    // record final state if not done so
    
    if ((t - out[last_time]) != 0 && (t != t))  {
        
        last_time++;
            
        for (i = 1; i <= (nors+npor); i++) {
                
            j = i*tsteps + last_time;
            out[j] = conc[i-1];
            out[last_time] = t;
            
        }    
    }
    
    
    
    mxFree(a);
    mxFree(a_cumsum);
    mxFree(ratDiss);
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  
    int i, j, ind, ind2;
    mwSize nORs, nPORs, nTsteps, iT;
    double tfp, tfk, tfb, dnak, ass, diss, tMax, dt, temp;
    double *res, *out, alpha, nmax, nexp;
    mxArray *OUT;
    
    // seed the rng with N = 624 random uint32s

    MTRand mtrands;
    MTRand::uint32 seed[ MTRand::N ];
    for( int n = 0; n < MTRand::N; ++n ) seed[n] = 23 * mtrands.randInt(); 

    mtrand.seed(seed);
    
    nORs = (mwSize)mxGetScalar(NORS);
    nPORs = (mwSize)mxGetScalar(NPOR);
    tfp = mxGetScalar(TFPROD);
    tfk = mxGetScalar(TFKILL);
    tfb = mxGetScalar(TFBASAL);
    dnak = mxGetScalar(DNAKILL);
    ass = mxGetScalar(TFASS);
    diss = mxGetScalar(TFDISS);
    alpha = mxGetScalar(ALPHA);
    dt = mxGetScalar(DT);
    tMax = mxGetScalar(TMAX);
    nmax = mxGetScalar(NMAX);
    nexp = mxGetScalar(NEXP);
    
	nTsteps = (mwSize)tMax/dt;
    

    // create output matrix for use in simulation routine
    
    OUT = mxCreateDoubleMatrix(nTsteps, (nORs+nPORs+2), mxREAL); 
    
    
    
    // Assign pointer for output
    
    out = mxGetPr(OUT);
    
    
    
    // run sim
    
    gillespie_ssim(nexp,nmax,tfp,tfk,tfb,dnak,ass,diss,alpha,dt,tMax,out,nORs,nPORs,nTsteps);
	
    
    
    
    
    // reformat output to remove unused rows
    
    for (iT = 1; iT < nTsteps; iT++)   {
        
        if ((out[iT] == 0) || (out[iT] != out[iT])) break;
        
    }
   
    RESULT = mxCreateDoubleMatrix(iT, (nORs+nPORs+2), mxREAL);
    res = mxGetPr(RESULT);
    
    for (j = 0; j <= (nORs+nPORs+1); j++) {
    
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