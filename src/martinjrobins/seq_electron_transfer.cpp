#include <stdio.h>
#include <iostream>
#include <math.h>
#include "sinusoidal_voltammetry.hpp"
#include "utilities.hpp"

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


#include <boost/math/constants/constants.hpp>

/* Problem Constants */

#define NEQ   3
#define NOUT  12

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Prototypes of functions called by IDA */

int resrob(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data);

static int grob(realtype t, N_Vector yy, N_Vector yp,
                realtype *gout, void *user_data);

int jacrob(long int Neq, realtype tt,  realtype cj, 
           N_Vector yy, N_Vector yp, N_Vector resvec,
           DlsMat JJ, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Prototypes of private functions */
static void PrintOutput(void *mem, realtype t, N_Vector y);
static void PrintFinalStats(void *mem);
static void check_flag(void *flagvalue, const char *funcname, int opt);

struct param_struct {
    double k01,k02;
    double alpha1,alpha2;
    double E01,E02;
    double Ru;
    double Cdl,CdlE,CdlE2,CdlE3;
    double Estart,Ereverse,dE,omega,phase;
    int reverse;
    Efun E;
};

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

void seq_electron_transfer(map& params, vector& Itot, vector& t) {
    param_struct p;
    p.k01 = get(params,std::string("k01"),35.0);
    p.k02 = get(params,std::string("k02"),65.0);
    p.alpha1 = get(params,std::string("alpha1"),0.5);
    p.alpha2 = get(params,std::string("alpha2"),0.5);
    p.E01 = get(params,std::string("E01"),0.25);
    p.E02 = get(params,std::string("E02"),-0.25);
    p.Ru = get(params,std::string("Ru"),2.74);
    p.Cdl = get(params,std::string("Cdl"),0.0037);
    p.CdlE = get(params,std::string("CdlE"),0.0);
    p.CdlE2 = get(params,std::string("CdlE2"),0.0);
    p.CdlE3 = get(params,std::string("CdlE3"),0.0);
    p.Estart = get(params,std::string("Estart"),-10.0);
    p.Ereverse = get(params,std::string("Ereverse"),10.0);

    const double pi = boost::math::constants::pi<double>();
    p.omega = get(params,std::string("omega"),2*pi);
    p.phase = get(params,std::string("phase"),0.0);
    p.dE = get(params,std::string("dE"),0.1);
    p.reverse = 0;
    
    p.E = Efun(p.Estart,p.Ereverse,p.dE,p.omega,p.phase,1.0);

#ifndef NDEBUG
    std::cout << "Running seq_electron_transfer with parameters:"<<std::endl;
    std::cout << "\tk01 = "<<p.k01<<std::endl;
    std::cout << "\tk02 = "<<p.k02<<std::endl;
    std::cout << "\talpha1 = "<<p.alpha1<<std::endl;
    std::cout << "\talpha2 = "<<p.alpha2<<std::endl;
    std::cout << "\tE01 = "<<p.E01<<std::endl;
    std::cout << "\tE02 = "<<p.E02<<std::endl;
    std::cout << "\tRu = "<<p.Ru<<std::endl;
    std::cout << "\tCdl = "<<p.Cdl<<std::endl;
    std::cout << "\tCdlE = "<<p.CdlE<<std::endl;
    std::cout << "\tCdlE2 = "<<p.CdlE2<<std::endl;
    std::cout << "\tCdlE3 = "<<p.CdlE3<<std::endl;
    std::cout << "\tEstart = "<<p.Estart<<std::endl;
    std::cout << "\tEreverse = "<<p.Ereverse<<std::endl;
    std::cout << "\tomega = "<<p.omega<<std::endl;
    std::cout << "\tdE= "<<p.dE<<std::endl;
#endif

    //IDA variables
    void *mem = NULL;
    N_Vector yy,yp,avtol;
    yy = yp = avtol = NULL;
    realtype *yval,*ypval,*atval;
    yval = ypval = atval = NULL;

    /* Allocate N-vectors. */
    yy = N_VNew_Serial(NEQ);
    check_flag((void *)yy, "N_VNew_Serial", 0);
    yp = N_VNew_Serial(NEQ);
    check_flag((void *)yp, "N_VNew_Serial", 0);
    avtol = N_VNew_Serial(NEQ);
    check_flag((void *)avtol, "N_VNew_Serial", 0);

    /* Create and initialize  y, y', and absolute tolerance vectors. */
    yval  = NV_DATA_S(yy);
    ypval = NV_DATA_S(yp);
    yval[0] = ONE;
    yval[1] = ZERO;
    yval[2] = ZERO;
    const realtype E = p.E(0);
    const realtype dEdt = p.E.ddt(0);
    yval[2] = p.Cdl*(ONE + p.CdlE*E + p.CdlE2*pow(E,2)+ p.CdlE3*pow(E,3))*dEdt;

#ifndef NDEBUG
    std::cout << "starting with initial conditions:"<<std::endl;
    std::cout << "\tyval = ("<<yval[0]<<","<<yval[1]<<","<<yval[2]<<")"<<std::endl;
    std::cout << "\typval = ("<<ypval[0]<<","<<ypval[1]<<","<<ypval[2]<<")"<<std::endl;
#endif

    realtype rtol = RCONST(1.0e-8);

    atval = NV_DATA_S(avtol);
    atval[0] = RCONST(1.0e-5);
    atval[1] = RCONST(1.0e-5);
    atval[2] = RCONST(1.0e-1);

    /* Integration limits */
    realtype tret;
    int Nt = t.size();
    if (Nt == 0) {
        Nt = 10000;
        const double Tmax = 2*(p.Ereverse-p.Estart);
        const realtype dt = Tmax/Nt;
        t.resize(Nt);
        for (int i=0; i<Nt; i++) {
            t[i] = i*dt;
        }
    }
    
    Itot.resize(Nt,0);
    Itot[0] = NV_DATA_S(yy)[2];

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate();
    check_flag((void *)mem, "IDACreate", 0);
    int retval = IDAInit(mem, resrob, t[0], yy, yp);
    check_flag(&retval, "IDAInit", 1);
    /* Call IDASVtolerances to set tolerances */
    retval = IDASVtolerances(mem, rtol, avtol);
    check_flag(&retval, "IDASVtolerances", 1);

    /* Free avtol */
    N_VDestroy_Serial(avtol);

    /* Call IDADense and set up the linear solver. */
    retval = IDADense(mem, NEQ);
    check_flag(&retval, "IDADense", 1);
    retval = IDADlsSetDenseJacFn(mem, jacrob);
    check_flag(&retval, "IDADlsSetDenseJacFn", 1);

    /* Set user_data */
    retval = IDASetUserData(mem, (void *)(&p)); 
    check_flag(&retval, "IDASetUserData", 1);

    /* In loop, call IDASolve, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached. */
    const realtype dE = p.Ereverse-p.Estart;
    for (int iout=1; iout < Nt ; iout++) {

        if (!p.reverse && t[iout] > dE) {
            retval = IDASolve(mem, dE, &tret, yy, yp, IDA_NORMAL);
            check_flag(&retval, "IDASolve", 1);
            p.reverse = 1;
            retval = IDASolve(mem, t[iout], &tret, yy, yp, IDA_NORMAL);
            check_flag(&retval, "IDASolve", 1);
        } else {
            retval = IDASolve(mem, t[iout], &tret, yy, yp, IDA_NORMAL);
            check_flag(&retval, "IDASolve", 1);
        }

        //PrintOutput(mem,tret,yy);

        if (retval == IDA_SUCCESS) {
            t[iout] = tret;
            Itot[iout] = NV_DATA_S(yy)[2];
            //std::cout << "t = "<<tret<<std::endl;
            //std::cout << "I = "<<NV_DATA_S(yy)[2]<<std::endl;
            //std::cout << "theta1 = "<<NV_DATA_S(yy)[0]<<" theta2 = "<<NV_DATA_S(yy)[1]<<std::endl;
        } else {
            char buffer[100];
            sprintf(buffer,"retval from IDASolve != IDA_SUCCESS. retval = %d",retval); 
            throw std::runtime_error(buffer);
        } 
    }

#ifndef NDEBUG
    PrintFinalStats(mem);
#endif
    /* Free memory */
    IDAFree(&mem);
    N_VDestroy_Serial(yy);
    N_VDestroy_Serial(yp);
}

/*
 *--------------------------------------------------------------------
 * Functions called by IDA
 *--------------------------------------------------------------------
 */

/*
 * Define the system residual function. 
 */

int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
    realtype *yval = NV_DATA_S(yy); 
    realtype *ypval = NV_DATA_S(yp); 
    realtype *rval = NV_DATA_S(rr);

    const param_struct &p = *((param_struct *)user_data);

    realtype E,dEdt;
    E = p.E(tres);
    dEdt = p.E.ddt(tres);
    
    //std::cout << "E = "<<E<<" dEdt = "<<dEdt<<std::endl;
    const realtype E2 = E*E;
    const realtype E3 = E2*E;
    const realtype expval1 = E - p.E01 - p.Ru*yval[2];
    const realtype expval2 = E - p.E02 - p.Ru*yval[2];
    const realtype dtheta1dt = p.k01*((ONE-yval[0]-yval[1])*std::exp((ONE-p.alpha1)*expval1) - yval[0]*std::exp(-p.alpha1*expval1));
    const realtype dtheta2dt = p.k02*((ONE-yval[0]-yval[1])*std::exp(-p.alpha2*expval2) - yval[1]*std::exp((ONE-p.alpha2)*expval2));

    rval[0] = dtheta1dt - ypval[0];
    rval[1] = dtheta2dt - ypval[1];
    rval[2] = p.Cdl*(ONE + p.CdlE*E + p.CdlE2*E2 + p.CdlE3*E3)*(dEdt-p.Ru*ypval[2]) + (ypval[0]- ypval[1]) - yval[2];

    //std::cout <<"res = ("<<rval[0]<<","<<rval[1]<<","<<rval[2]<<")"<<std::endl;

    return(0);
}

/*
 * Define the Jacobian function. 
 */

int jacrob(long int Neq, realtype tt,  realtype cj, 
           N_Vector yy, N_Vector yp, N_Vector resvec,
           DlsMat JJ, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
    realtype *yval;
  
    yval = NV_DATA_S(yy);

    const param_struct &p = *((param_struct *)user_data);

    realtype E = p.E(tt);

    const realtype E2 = E*E;
    const realtype E3 = E2*E;
    const double expval1 = E - p.E01 - p.Ru*yval[2];
    const double expval2 = E - p.E02 - p.Ru*yval[2];
    const double exp11 = std::exp((1-p.alpha1)*expval1);
    const double exp12 = std::exp(-p.alpha1*expval1);
    const double exp21 = std::exp(-p.alpha2*expval2);
    const double exp22 = std::exp((1-p.alpha2)*expval2);

    DENSE_ELEM(JJ,0,0) = p.k01*(-exp11-exp12) - cj;
    DENSE_ELEM(JJ,0,1) = p.k01*(-exp11);
    DENSE_ELEM(JJ,0,2) = p.k01*((1-yval[0]-yval[1])*(p.Ru*(p.alpha1-1))*exp11 - yval[0]*p.alpha1*p.Ru*exp12);

    DENSE_ELEM(JJ,1,0) = p.k02*(-exp21);
    DENSE_ELEM(JJ,1,1) = p.k02*(-exp21-exp22) - cj;
    DENSE_ELEM(JJ,1,2) = p.k02*((1-yval[0]-yval[1])*(p.Ru*p.alpha2)*exp21 + yval[1]*(1-p.alpha2)*p.Ru*exp22);

    DENSE_ELEM(JJ,2,0) = cj;
    DENSE_ELEM(JJ,2,1) = -cj;
    DENSE_ELEM(JJ,2,2) = -ONE - cj*p.Cdl*(ONE + p.CdlE*E + p.CdlE2*E2 + p.CdlE3*E3)*p.Ru;

    return(0);
}

/*
 *--------------------------------------------------------------------
 * Private functions
 *--------------------------------------------------------------------
 */

static void PrintOutput(void *mem, realtype t, N_Vector y)
{
  realtype *yval;
  int retval, kused;
  long int nst;
  realtype hused;

  yval  = NV_DATA_S(y);

  retval = IDAGetLastOrder(mem, &kused);
  check_flag(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_flag(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.4Le %12.4Le %12.4Le %12.4Le | %3ld  %1d %12.4Le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
  printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#endif
}

/*
 * Print final integrator statistics
 */

static void PrintFinalStats(void *mem)
{
  int retval;
  long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;

  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_flag(&retval, "IDAGetNumResEvals", 1);
  retval = IDADlsGetNumJacEvals(mem, &nje);
  check_flag(&retval, "IDADlsGetNumJacEvals", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDADlsGetNumResEvals(mem, &nreLS);
  check_flag(&retval, "IDADlsGetNumResEvals", 1);
  retval = IDAGetNumGEvals(mem, &nge);
  check_flag(&retval, "IDAGetNumGEvals", 1);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);
  printf("Number of root fn. evaluations     = %ld\n", nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static void check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        char buffer[100];
        sprintf(buffer,"SUNDIALS_ERROR: %s() failed - returned NULL pointer\n",funcname);
        throw std::runtime_error(buffer);
    } else if (opt == 1) {
        /* Check if flag < 0 */
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            char buffer[100];
            sprintf(buffer,"\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
                  funcname, *errflag);
            throw std::runtime_error(buffer);
        }
    } else if (opt == 2 && flagvalue == NULL) {
        /* Check if function returned NULL pointer - no memory allocated */
        char buffer[100];
        sprintf(buffer,"\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
                funcname);
        throw std::runtime_error(buffer);
    }
}
