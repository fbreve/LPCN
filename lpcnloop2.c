#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	maxiter_IN      prhs[0]
#define	npart_IN        prhs[1]
#define	nclass_IN       prhs[2]
#define	omega_IN        prhs[3]
#define partnode_IN     prhs[4]
#define slabel_IN       prhs[5]
#define nsize_IN        prhs[6]
#define nlist_IN        prhs[7]
#define ndist_IN        prhs[8]
#define pot_IN          prhs[9]

/* Output Arguments */

#define ph2_ttiter_OUT  plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )    
{    
    int maxiter, npart, nclass; // escalares int
    unsigned int *partnode;
    unsigned char *nsize;    
    unsigned short int *slabel;
    unsigned int *nlist; // matrizes de int
    double *ndist, *pot;  // matrizes de double
    int qtnode, neibmax;
    double omega;
    
    /* Check for proper number of arguments */
    
    
    if (nrhs != 10) { 
	    mexErrMsgTxt("10 argumentos de entrada requeridos."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Only 1 output argument allowed."); 
    }
    
    maxiter = (int) mxGetScalar(maxiter_IN);
    npart = (int) mxGetScalar(npart_IN);
    nclass = (int) mxGetScalar(nclass_IN);
    omega = mxGetScalar(omega_IN);
    partnode = (unsigned int *) mxGetData(partnode_IN);
    slabel = (unsigned short int *) mxGetData(slabel_IN);
    nsize = (unsigned char *) mxGetData(nsize_IN);    
    nlist = (unsigned int *) mxGetData(nlist_IN);    
    ndist = mxGetPr(ndist_IN);
    pot = mxGetPr(pot_IN);
    
    qtnode = (int) mxGetM(slabel_IN);
    neibmax = (int) mxGetN(nlist_IN);  // quantidade máxima de vizinhos que um nó tem   
            
    // non-Windows users should probably use /dev/random or /dev/urandom instead of rand_s
    //unsigned int seed;
    //errno_t err;
    //err = rand_s(&seed);
    //if (err != 0) printf_s("The rand_s function failed!\n");
    //srand(seed);
    double maxmmpot = 0;
    double *prob = malloc(sizeof(double)*neibmax); // vetor de probabilidades de visitar vizinho    
    double *nc = malloc(sizeof(double)*nclass);
    double *newpot = malloc(sizeof(double)*qtnode*nclass);
    for(int i=0; i<qtnode*nclass; i++) newpot[i]=pot[i];
    bool *labeled = malloc(sizeof(bool)*npart);
    for(int i=0; i<npart; i++) labeled[i]=false;
    int labeledc = 0;
    int i;
    for(i=0; i<maxiter; i++)
    {
        for(int j=0; j<npart; j++)
        {
            for(int i2=0; i2<nclass; i2++) nc[i2]=0;
            
            double sumweight=0;
            for(int i2=0; i2<nsize[partnode[j]-1]; i2++)
            {
                for(int i3=0; i3<nclass; i3++)
                    nc[i3] = nc[i3] + pot[(i3*qtnode + nlist[(i2*qtnode + partnode[j]-1)]-1)] * ndist[(i2*qtnode + partnode[j]-1)];
                sumweight += ndist[(i2*qtnode + partnode[j]-1)];
            }
            
            for(int i2=0; i2<nclass; i2++)
                newpot[(i2*qtnode + partnode[j]-1)] = nc[i2] / sumweight;
        }
            
        for(int j=0; j<npart; j++) 
            for(int i2=0; i2<nclass; i2++)
                pot[(i2*qtnode + partnode[j]-1)] = newpot[(i2*qtnode + partnode[j]-1)];                    
        
        if (i % 10 == 0)
        {
            double mmpot = 0;
            for(int i2=0; i2<qtnode; i2++)
            {
                double mpot=0;
                for(int i3=0; i3<nclass; i3++)
                    if(pot[i3*qtnode + i2]>mpot) mpot = pot[i3*qtnode + i2];
                mmpot += mpot;
            }
            mmpot /= qtnode;
            
//             if (i % 10000 == 0)
//             {
//                 printf("Iter: %i  Meanpot: %0.4f\n",i,mmpot);
//                 mexEvalString("drawnow");
//             }
                      
            if (mmpot - maxmmpot > omega) maxmmpot = mmpot;
            else break;
        }
    }
    free(prob);
    free(nc);
    free(newpot);
    free(labeled);
    
    ph2_ttiter_OUT = mxCreateDoubleScalar(i);
    
    return;    
}
