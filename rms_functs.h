#ifndef RMS_FUNCTS_H
#define RMS_FUNCTS_H

#include <vector>
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"

void covar_evectors(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    TMatrixDSym& cov,
    TVectorD& eigenValues,
    TMatrixD& eigenVectors
);

double CoreMeanAlpha(double a, double b, double c);

double CoreMeanAlpha2(double a, double b, double c);

double DeltaAlpha2_Eigen(double a, double b, double c,
                         const TVectorD& eigenVec,
                         double eigenVal);

double CoreAlpha2Error(double a, double b, double c,
                       const TMatrixDSym& cov);


                       

double x2_expect(double a, double b, double c);

double rms2(double a, double b, double c);

double var_par(double a, double b, double c, 
                    const TVectorD& eigenVec,
                    double eigenVal);

double rms_err_f(double a, double b, double c,
                        const TMatrixDSym& cov);
#endif


/* Compile like this:
root -l
.L rms_functs.cxx+
.L full_acop_plot.C+
full_acop_plot();
*/