#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	maxiter_IN       prhs[0]
#define	nnonlabeled_IN   prhs[1]
#define indnonlabeled_IN prhs[2]
#define	omega_IN         prhs[3]
#define	potval_IN        prhs[4]
#define	k_IN             prhs[5]
#define	nlist_IN         prhs[6]
#define ndist_IN         prhs[7]

/* Output Arguments */

#define ph1_ttiter_OUT   plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{    
    int maxiter, nnonlabeled; // escalares int
    unsigned short int k; // escalar de uint16       
    unsigned int *indnonlabeled, *nlist; // vetores de uint32
    double *potval, *ndist;  // matrizes de double
    int qtnode, nclass;
    double omega;
    
    /* Check for proper number of arguments */
    
    
    if (nrhs != 8) { 
	    mexErrMsgTxt("8 argumentos de entrada requeridos."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Only 1 output argument allowed."); 
    }
    
    maxiter = (int) mxGetScalar(maxiter_IN);
    nnonlabeled = (int) mxGetScalar(nnonlabeled_IN);
    indnonlabeled = (unsigned int *) mxGetData(indnonlabeled_IN);
    omega = mxGetScalar(omega_IN);   
    k = (unsigned short int) mxGetScalar(k_IN);                
    nlist = (unsigned int *) mxGetData(nlist_IN);    
    potval = mxGetPr(potval_IN);
    ndist = mxGetPr(ndist_IN);
    
    qtnode = (int) mxGetM(potval_IN);    
    nclass = (int) mxGetN(potval_IN);      
       
    double maxmmpot = 0;
    double *newpot = malloc(sizeof(double) * nnonlabeled * nclass);
    int i;
    for(i=0; i<maxiter; i++)
    {
        
        for(int j=0; j<nnonlabeled; j++)
        {            
            // inicialmente potenciais novos s�o zerados
            //newpot[j] = 0;
            //newpot[j + nnonlabeled] = 0;        
            for(int j2=0; j2<nclass; j2++) newpot[j + nnonlabeled*j2] = 0;            
            // peso acumulado de todos os vizinhos
            double accweight=0;
            // para cada vizinho do n� n�o rotulado
            
            for(int ki=0; ki<k; ki++)
            {           
                // vamos pegar um vizinho
                int neib = nlist[j + ki*nnonlabeled]-1;
                // vamos somar os potenciais dos vizinhos no novo potencial do nosso n� n�o rotulado
                // lembrar que em C acessa-se matriz por [LINHA + COLUNA * QTDE DE LINHAS]
                //newpot[j]               += potval[neib]          * ndist[j + ki*nnonlabeled];
                //newpot[j + nnonlabeled] += potval[neib + qtnode] * ndist[j + ki*nnonlabeled];                
                for(int j2=0; j2<nclass; j2++) newpot[j + nnonlabeled*j2] += potval[neib + qtnode*j2] * ndist[j + ki*nnonlabeled];                                
                accweight += ndist[j + ki*nnonlabeled];
            }
            // dividindo os potenciais acumulados pela quantidade de vizinhos, para obter a m�dia
            //newpot[j] /= accweight;
            //newpot[j + nnonlabeled] /= accweight;            
            for(int j2=0; j2 <nclass; j2++) newpot[j + nnonlabeled*j2] /= accweight;                                    
        }
        
        // colocar os novos potenciais na lista de potenciais
        
        for(int j=0; j<nnonlabeled; j++) 
        {
            int ppj = indnonlabeled[j]-1;
            //potval[ppj] = newpot[j];
            //potval[ppj + qtnode] = newpot[j + nnonlabeled];
            for(int j2=0; j2 <nclass; j2++) potval[ppj + qtnode*j2] = newpot[j + nnonlabeled*j2];
        }

        // vamos testar converg�ncia                 
        if (i % 10 == 0)
        {
            // vari�vel para guardar a m�dia de maior potencial                      
            double mmpot = 0;
            for(int j=0; j<nnonlabeled; j++)
            {
                double mpot=0;
                // vamos pegar o n�
                int ppj = indnonlabeled[j]-1;
                // vamos achar o maior potencial dentro do n�
                for(int i3=0; i3<nclass; i3++) if(potval[ppj + qtnode*i3]>mpot) mpot = potval[ppj + qtnode*i3];
                // e ent�o som�-lo
                mmpot += mpot;
            }
            // e por fim dividir pela quantidade de n�s para obter a m�dia
            mmpot /= nnonlabeled;            
            
            //printf("Iter: %i  Meanpot: %0.4f\n",i,mmpot);            
            // se da �ltima maior m�dia para a atual aumentou mais que 0.001
            if (mmpot - maxmmpot > omega) maxmmpot = mmpot;   
            else break;
        }             
    }      
    
    free(newpot);
    
    ph1_ttiter_OUT = mxCreateDoubleScalar(i);
    return;
    
}
