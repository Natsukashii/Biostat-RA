#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"




static const double one=1.0;
static const double zero=0.0;
static const int inc_one = 1;
static const double EPS = 1e-4;
// static const double SCALE_FACTOR = 0.9;
static const double MAX_DEV_EXPLAINED = 0.99;

int sign(double x){
	if (x>0){ 
		return 1;
	} else if (x==0.0){
		return 0;
	} else{
		return -1;
	}
}

/* Computes difference between vectors x and y, both of length n
   and stores it in vector z */
void vectorDifference(int *n, double *x, double *y, double *z) {
  for(int i=0; i<*n; i++) {
    z[i] = x[i] - y[i];
  }
}

double implicitFunction(double c, double m, double *D, double *Vr, double lambda) {
  double out = -1;
  for(int i=0; i<m; i++) {
    out += pow(Vr[i]/(D[i]*c + lambda),2);
  }
  return out;
}

double lineSearch(double m, double *D, double *Vr, double lambda) {
  double tol = 0.0001;
  double start = 1;  
  double end = 2;
  //double scale = 2;
  double start_val = implicitFunction(start, m, D, Vr, lambda);
  double end_val = implicitFunction(end, m, D, Vr, lambda);
  while(sign(start_val) == sign(end_val)) {
    if(start_val > end_val && sign(start_val) > 0) {
      end*=2.0;
      end_val = implicitFunction(end, m, D, Vr, lambda);
    } else {
      start*=0.5;
      start_val = implicitFunction(start, m, D, Vr, lambda);
    }
  }
  double midpt = 0.5*(start + end);
  double midpt_val = implicitFunction(midpt, m, D, Vr, lambda);
  while(end - start > tol) {
    if(sign(midpt_val) == sign(start_val)) {
      start = midpt;
      start_val = midpt_val;
    } else {
      end = midpt;
      end_val = midpt_val;
    }
    midpt = 0.5*(start + end);
    midpt_val = implicitFunction(midpt, m, D, Vr, lambda);
  }
  return midpt; 
}



// Beta vector update
int updateBeta(int j, int *n, double *y, double *U, double *fit, double *lambdas, double *lambdas_alpha, double *psis, double *D, int *degrees, int *cum_degrees, double *betas, double *alphas, double *z, int *family, double *W) {
  double *resid = calloc(*n, sizeof(double));
  double *UjBj = calloc(*n, sizeof(double));
  double *DinvUjr = calloc(degrees[j], sizeof(double));
  double *old_beta = calloc(degrees[j], sizeof(double));
  int nonzero = 1;
  
  // Store old value of beta
  for(int i=0; i<degrees[j]; i++) {
    old_beta[i] = betas[cum_degrees[j] + i];
  }
  if(*family == 0) {  // Linear regression updates
    // Calculate residual
    vectorDifference(n,y,fit,resid);
    // Adjust resdiual for contribution of beta_j
    dgemv_("N", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, betas+cum_degrees[j], &inc_one, &zero, UjBj, &inc_one);
    for(int i=0; i<*n; i++) {
      resid[i] += UjBj[i];
    }
    // Compute D^{-1/2} * U_j^T * resid
    dgemv_("T", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, resid, 
      &inc_one, &zero, DinvUjr, &inc_one);
    for(int i=0; i<degrees[j]; i++) {
      DinvUjr[i] = DinvUjr[i] / sqrt(D[cum_degrees[j]+i]);
    }
    double norm_DinvUjr = sqrt(ddot_(degrees+j, DinvUjr, &inc_one, DinvUjr, &inc_one));
    // Check if beta_j is 0
    if(norm_DinvUjr > lambdas[j]) {
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j] + i] = 0.0;
      }
      nonzero = 0;
    } else {
      double *Dstar = calloc(degrees[j], sizeof(double));
      // Form D*
      for(int i=0; i<degrees[j]; i++) {
        Dstar[i] = psis[j] + 1/D[cum_degrees[j]+i];
      }
      Dstar[0] -= psis[j];  // Adjust first entry
      // Line search to find c
      double c = lineSearch(degrees[j], Dstar, DinvUjr, lambdas[j]);
      // Update beta using equation (4.3)
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j]+i] = (DinvUjr[i]/sqrt(D[cum_degrees[j]+i]))/(Dstar[i] + lambdas[j]/c);
      }
      free(Dstar);
    }
  }  else if(*family == 1) {  // Logistic regression updates without approximation
  	double *VT = calloc((*n)*(degrees[j]), sizeof(double));
  	double *VTr = calloc(degrees[j], sizeof(double));
  	double *QTVTr = calloc(degrees[j], sizeof(double));
  	
    // Calculate residual
    vectorDifference(n,z,fit,resid);
    
    // Adjust resdiual (W^{1/2} * resid) for contribution of beta_j
    dgemv_("N", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, betas+cum_degrees[j], &inc_one, &zero, UjBj, &inc_one); 
    for(int i=0; i<*n; i++) {
      resid[i] += UjBj[i];
      resid[i] = resid[i] * sqrt(W[i]);
    }
    
    // Compute V^T = D^{-1/2} * U_j^T * W^{1/2}
    for(int k=0; k<degrees[j]; k++) {
    	for(int i=0; i<*n; i++) 
	  VT[k*(*n)+i] = U[(*n)*(cum_degrees[j])+k*(*n)+i] * sqrt(W[i]) / sqrt(D[cum_degrees[j]+k]);
    }
    
    // Compute V^T * r = D^{-1/2} * U_j^T * W * resid
    dgemv_("N", degrees+j, n, &one, VT, degrees+j, resid, 
      &inc_one, &zero, VTr, &inc_one);
      
    double norm_VTr = sqrt(ddot_(degrees+j, VTr, &inc_one, VTr, &inc_one));
    // Check if beta_j is 0
    if(norm_VTr == lambdas[j]) {
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j] + i] = 0.0;
      }
      nonzero = 0;
    } else { 
      double *Dtilde = calloc(pow((degrees[j]),2), sizeof(double));
      // Caculate V^T * V 
      dsyrk_("U", "N", degrees+j, n, &one, VT, n, &zero, Dtilde, degrees+j);
      // Caculate D_tilde = V^T * V + psi * I_{(-1)}
      for(int i=0; i<degrees[j]; i++) {
        Dtilde[i*(degrees[j])+i] += psis[j];
      }
      Dtilde[0] -= psis[j];  // Adjust first entry
      
      double *Q = Dtilde;
      double *Lam = calloc(degrees[j], sizeof(double));
      int lwork = 3*degrees[j]-1;
      double *work = calloc(lwork, sizeof(double));
      int info = 100;
      // Eigen-decomposition 
      dsyev_("V", "U", degrees+j, Q, degrees+j, Lam, work, &lwork, &info);
      
      // Caculate Q^T * V^T * r
      dgemv_("T", degrees+j, degrees+j, &one, Q, degrees+j, VTr, &inc_one, &zero, QTVTr, &inc_one);
      
      // Line search to find c
      double c = lineSearch(degrees[j], Lam, QTVTr, lambdas[j]);
      
  	  // Caculate (Dtilde + lambda_tilde/c * I)^{(-1)}
      for(int i=0; i<degrees[j]; i++) {
        Dtilde[i*(degrees[j])+i] += lambdas[j]/c;
      }
      int *ipiv = calloc(degrees[j], sizeof(int));
      dgetrf_(degrees+j, degrees+j, Dtilde, degrees+j, ipiv, &info);
      dgetri_(degrees+j, Dtilde, degrees+j, ipiv, work, degrees+j, &info);
      
      // Caculate new beta
	  double *new_beta = malloc(degrees[j]*sizeof(sizeof(double)));
	  dgemv_("N", degrees+j, degrees+j, &one, Dtilde, degrees+j, VTr, &inc_one, &zero, new_beta, &inc_one); 
      for(int i=0; i<degrees[j]; i++) {
        new_beta[i] = new_beta[i]/sqrt(D[cum_degrees[j]+i]);
      }
      
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j]+i] = new_beta[i];
      }
      free(Dtilde);
      free(Q);
      free(Lam);
      free(work);
      free(ipiv);
      free(new_beta);
	}
    free(VT);
    free(VTr);
    free(QTVTr);
  }
  /** Fitted value update step */
  double *beta_diff = calloc(degrees[j], sizeof(double));
  int equal_Flag = 0;
  for(int i=0; i<degrees[j]; i++) {
    beta_diff[i] = betas[cum_degrees[j]+i] - old_beta[i];
    if(old_beta[i] != betas[cum_degrees[j]+i]) {
      equal_Flag = 1;
    }
  }
  if(equal_Flag != 0) {
    double *Ujbeta_diff = calloc(*n, sizeof(double));
    dgemv_("N", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, beta_diff, 
      &inc_one, &zero, Ujbeta_diff, &inc_one);
    for(int i=0; i<*n; i++) {
      fit[i] += Ujbeta_diff[i];
    }
    free(Ujbeta_diff);
  }
  free(resid);
  free(UjBj);
  free(old_beta);
  free(DinvUjr);
  free(beta_diff);
  return nonzero;
}
