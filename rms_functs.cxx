#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"   // For Form()
#include <vector>
#include <stdexcept>

//################# this gives the eigenvectors #################
void covar_evectors(
    const vector<double>& x,
    const vector<double>& y,
    const vector<double>& z,
    TMatrixDSym& cov,
    TVectorD& eigenValues,
    TMatrixD& eigenVectors
    ){
    // Sanity checks
    if (x.size() != y.size() || x.size() != z.size())
        throw std::runtime_error("Input vectors must have the same size");

    const int N = x.size();
    if (N < 2)
        throw std::runtime_error("Need at least two entries to compute covariance");

    // Compute means
    double mx = 0.0, my = 0.0, mz = 0.0;
    for (int i = 0; i < N; ++i) {
        mx += x[i];
        my += y[i];
        mz += z[i];
    }
    mx /= N;
    my /= N;
    mz /= N;

    // Resize outputs
    cov.ResizeTo(3,3);
    cov.Zero();

    eigenValues.ResizeTo(3);
    eigenVectors.ResizeTo(3,3);

    // Fill covariance matrix
    for (int i = 0; i < N; ++i) {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        const double dz = z[i] - mz;

        cov(0,0) += dx * dx;
        cov(0,1) += dx * dy;
        cov(0,2) += dx * dz;

        cov(1,1) += dy * dy;
        cov(1,2) += dy * dz;

        cov(2,2) += dz * dz;
    }

    // Sample covariance normalization
    cov *= 1.0 / (N - 1);

    // Enforce symmetry explicitly
    cov(1,0) = cov(0,1);
    cov(2,0) = cov(0,2);
    cov(2,1) = cov(1,2);

    // Diagonalize
    TMatrixDSymEigen eigenSolver(cov);
    eigenValues  = eigenSolver.GetEigenValues();
    eigenVectors = eigenSolver.GetEigenVectors();

    for (int i = 0; i < 3; ++i) {
        double norm = 0.0;
        for (int j = 0; j < 3; ++j)
            norm += eigenVectors(j,i) * eigenVectors(j,i);

        norm = std::sqrt(norm);

        if (norm > 0) {
            for (int j = 0; j < 3; ++j)
                eigenVectors(j,i) /= norm;
        }
    }
}

// ################# Calculate <x^2> #################
double x2_expect(double a, double b, double c){
    TF1 variance(Form("variance_%p", &variance), "[0]*exp(-x/[1] + [2]*pow(x,0.25)) * x*x", 0.05*10e-5, 1);
    variance.SetParameters(a, b, c);

    return variance.Integral(0.5*10e-5, 1);
}


// ################# Calculate <x^2>/<x> = rms^2 #################
double rms2(double a, double b, double c){
    TF1 f(Form("f_%p", &f), "[0]*exp(-x/[1] + [2]*pow(x,0.25))", 0.05*10e-5, 1);
    f.SetParameters(a, b, c);

    double f_int = f.Integral(0.5*10e-5, 1);
    return x2_expect(a, b, c) / f_int;
}


// ################# Vary a single parameter #################
double var_par(double a, double b, double c, const TVectorD& eigenVec, double eigenVal){
    const double del = sqrt(eigenVal);

    double a_var = a + del * eigenVec(0);
    double b_var = b + del * eigenVec(1);
    double c_var = c + del * eigenVec(2);

    return rms2(a_var, b_var, c_var) - rms2(a, b, c);
}


// ################# Combine variations for error #################
double rms_err_f(double a, double b, double c, const TMatrixDSym& cov){
    TMatrixDSymEigen eig(cov);
    TVectorD eigenValues  = eig.GetEigenValues();
    TMatrixD eigenVectors = eig.GetEigenVectors();

    vector<double> varied(3); //a, b, c
    
    for (int i = 0; i < 3; ++i) {
        TVectorD v(3);
        v(0) = eigenVectors(0,i);
        v(1) = eigenVectors(1,i);
        v(2) = eigenVectors(2,i);

        double d = var_par(a, b, c, v, eigenValues(i));
        varied[i] = d * d; //squaring them here for convience
    }
    return sqrt(varied[0] + varied[1] + varied[2]);
}








//################# this gives <x^2>/<x> = rms #################
double CoreMeanAlpha2(double a, double b, double c){
    TF1 core_tmp(Form("core_tmp_%p", &core_tmp), "[0]*exp(-x/[1] + [2]*pow(x,0.25))", 0.5*10e-5, 1);

    core_tmp.SetParameters(a, b, c);

    double norm = core_tmp.Integral(0.5*10e-5, 1);
    if (norm <= 0.0) return 0.0;

    TF1 core_x2(Form("core_x2_%p", &core_x2), "[0]*exp(-x/[1] + [2]*pow(x,0.25)) * x*x", 0.5*10e-5, 1);
    core_x2.SetParameters(a, b, c);

    return core_x2.Integral(0.5*10e-5, 1) / norm;
}


////################# this gives <x^2> //#################
double x2_expt(double a, double b, double c){
    TF1 expect2(Form("expect2_%p", &expect2), "[0]*exp(-x/[1] + [2]*pow(x,0.25)) * x*x", 0.5*10e-5, 1);

    expect2.SetParameters(a, b, c);

    double norm = expect2.Integral(0.5*10e-5, 1);
    return norm;
}


//################# this gives the variation of a single parameter #################
double DeltaAlpha2_Eigen(double a, double b, double c, const TVectorD& eigenVec, double eigenVal){
    const double delta = std::sqrt(eigenVal);

    double a_shift = a + delta * eigenVec(0);
    double b_shift = b + delta * eigenVec(1);
    double c_shift = c + delta * eigenVec(2);

    double mean_nom = CoreMeanAlpha2(a, b, c);

    double mean_shift = CoreMeanAlpha2(a_shift, b_shift, c_shift);

    return mean_shift - mean_nom; //HERE IS ISSUE. NOT SUBTRACTING RIGHT VALUES
}


//################# this propogates the error for a single vector #################
double CoreAlpha2Error(double a, double b, double c, const TMatrixDSym& cov)
{
    TMatrixDSymEigen eig(cov);
    TVectorD eigenValues  = eig.GetEigenValues();
    TMatrixD eigenVectors = eig.GetEigenVectors();

    double err2 = 0.0;

    for (int i = 0; i < 3; ++i) {
        TVectorD v(3);
        v(0) = eigenVectors(0,i);
        v(1) = eigenVectors(1,i);
        v(2) = eigenVectors(2,i);

        double d = DeltaAlpha2_Eigen(a, b, c, v, eigenValues(i));

        err2 += d * d;
    }
    return std::sqrt(err2);
}


//################# this gives the mean #################
double CoreMeanAlpha(double a, double b, double c)
{
    TF1 f("f","[0]*exp(-x/[1]+[2]*pow(x,0.25))", 0.5*10e-5, 1);
    f.SetParameters(a,b,c);

    TF1 fx("fx","[0]*exp(-x/[1]+[2]*pow(x,0.25))*x", 0.5*10e-5, 1);
    fx.SetParameters(a,b,c);

    double norm = f.Integral(0.5*10e-5, 1);
    if (norm <= 0) return 0.0;

    return fx.Integral(0.5*10e-5, 1) / norm;
}