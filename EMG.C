#include <cmath>
#include <iostream>
inline long double Sqr(long double x) { return x*x; } //lets us use Brian's Sqr notation

double erfApprox(long double x)
{
  long double t = 1.0/(1.0 + 0.5*std::abs(x));

  long double gauss = exp(-x*x);
  long double taumod = t*std::expl(-1.26551223 + 1.00002368*t + 0.37409196*t*t + 0.09678418*std::powl(t,3) - 0.18628806*std::powl(t,4) + 0.27886807*std::powl(t,5) - 1.13520398*std::powl(t,6) + 1.48851587*std::powl(t,7) - 0.82215223*std::powl(t,8) + 0.17087277*std::powl(t,9));

  if (x > 0) return 1 - taumod*gauss;
  else return taumod*gauss - 1;
}

double erfTauMod(long double x)
{
  long double t = 1.0/(1.0 + 0.5*std::abs(x));
  return t*std::expl(-1.26551223 + 1.00002368*t + 0.37409196*t*t + 0.09678418*std::powl(t,3) - 0.18628806*std::powl(t,4) + 0.27886807*std::powl(t,5) - 1.13520398*std::powl(t,6) + 1.48851587*std::powl(t,7) - 0.82215223*std::powl(t,8) + 0.17087277*std::powl(t,9));
}

double EMG(Double_t *x, Double_t *par)
{
  const double sqrt2 = std::sqrt(2);
  
  long double integral = par[0]; 
  long double tau = par[1];
  long double sigma = std::abs(par[2]);
  long double x0 = par[3];

  long double arg = x[0] - x0;
  long double argSigRatio = arg/sigma;

  long double erfc1arg = (arg/sigma + sigma/tau)/sqrt2;
  long double erfc2arg = (arg/sigma - sigma/tau)/sqrt2;

  long double norm = std::expl(Sqr(sigma)/2.0/Sqr(tau));
  long double term1 = 0.5*std::expl(arg/tau + Sqr(sigma)/2.0/Sqr(tau))*std::erfcl(erfc1arg);
  long double term2;

  if (erfc2arg > -4) {
    term2 = 0.5*std::expl(-arg/tau + Sqr(sigma)/2.0/Sqr(tau))*(1.0 + std::erfl(erfc2arg));
  }
  else {
    long double t = 1.0/(1.0 + 0.5*std::abs(erfc2arg));
    long double tauNoGauss = t*std::expl(-1.26551223 + 1.00002368*t + 0.37409196*t*t + 0.09678418*std::powl(t,3) - 0.18628806*std::powl(t,4) + 0.27886807*std::powl(t,5) - 1.13520398*std::powl(t,6) + 1.48851587*std::powl(t,7) - 0.82215223*std::powl(t,8) + 0.17087277*std::powl(t,9));

    term2 = 0.5*std::expl(-arg/tau + Sqr(sigma)/2.0/Sqr(tau) - Sqr(erfc2arg))*tauNoGauss;
  }

//  if (std::isnan(term1) || std::isnan(term2)) {
//    std::cout << "arg, erfc1arg, erfc2 arg = " << arg << ", " << erfc1arg << ", " << erfc2arg 
//	      << "norm, term1, term2 = " << norm << ", " << term1 << ", " << term2 << std::endl;
//  }
  
  // // double erfc1 = erfc1arg > 5 ? 0 : 0.5*
  // // double erfc1 = 0.5*TMath::Erfc(erfc1arg);
  // // double erfc2 = 0.5*TMath::Erfc(erfc2arg);
  //
  // double term1 = 0;
  // if (erfc1arg < 7) {
  //   term1 = std::exp(arg/tau)*0.5*TMath::Erfc(erfc1arg);
  // }
  // else {
  //   //
  //   // Here we use an approximation to the erfc from wikipedia that allows us to combine
  //   //   the arguments of the exponential to effect the cancellation between large e^{x/tau} and small erfc
  //   //
  //   const double A = 1.98;
  //   const double B = 1.135;
  //   term1 = 0.5*(1.0 - exp(-A*erfc1arg))*std::exp(arg/tau - Sqr(erfc1arg))/B/std::sqrt(TMath::Pi())/erfc1arg;
  // }
  // // We can safely truncate the erfc here because its being subtracted from 1
  // //
  // double erfc2 = erfc2arg > 7 ? 0 : 0.5*TMath::Erfc(erfc2arg);
  // double term2 = std::exp(-arg/tau)*(1.0 - erfc2);

  double value = integral*(term1 + term2)/(2.*tau);

  return value;
}

double EMGSymm(Double_t *x, Double_t *par)
{
  std::vector<double> parext = {par[0], par[1], par[2], 0};
  return EMG(x, &parext[0]);
}

double EMGSymmC(Double_t *x, Double_t *par)
{
  return EMGSymm(x, par) + par[3];
}
