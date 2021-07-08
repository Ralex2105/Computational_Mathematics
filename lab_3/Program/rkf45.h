#ifndef DIFFERENCIAL_RKF45
#define DIFFERENCIAL_RKF45


void RKF45(void (*F)(double t1,double * y1,double * dy1),
		   int NEQN,double *Y,double *T,double *TOUT,
		   double *RELERR,double *ABSERR,
		   int *IFLAG,double *WORK,int *IWORK);

#endif

