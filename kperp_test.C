/*
4/6: Reordered them so that Nn0n on their own canvas. Fixed green-line-too-big problem by removing setting parameters from ChatGPT. Added 2n0n. Plotted
XnXn on linear plot. Still missing log likelihood. Forgot to figure out combining EMG and first power law into one with a single multiplicative factor
*/
// Z upper limit of 10
//motication better by CMS by not just testing numbers but hyposeses of where the disburptions in k_[pert] and alpha come from

//create plots of the parameters as function of topologies (sigma, x_0, and alpha). Replace N with something else in notation. Done
//x_0 is vary. Done
//plot Yn1n and Zn1n with a fixed x_0 and not fixed x_0. original fit
//also rapdiity plot for pairs over 150 MeV.  Two full sets in total

//get second moment of the EMG. Fix alpha to 250 and see what happens
    //fixing alpha has not been working well

//Getting that the 2nd moment is M_2 = <x^2> = x_0^2 + sigma^2 + 2*tau^2. For limited range (here), do usual <x^2> but not over +/- inf

// fix parameter plot so that it includes 0n0n. This is from removing 0n0n from main plot

#include "AtlasStyle.C"
#include "rms_functs.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>
#include "TLatex.h"
#include "TF1.h"
#include "EMG.C" //Brian's modified Gaussian
#include <cmath>
#include <vector>

#include "TMatrixDSym.h" // these are to run cov_mtrx.C for the covariance matrix
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <stdexcept>

//############################### Complete kperp function ###############################
Double_t EMG_plus_power(Double_t *x, Double_t *par)
{
    // --- EMG part ---
    Double_t emg = EMG(x, par);  // uses par[0..3]

    // --- First Power-law part ---
    Double_t A1  = par[4];
    Double_t x01 = par[5];   // scale
    Double_t n1  = par[6];

    //if (x0 <= 0) return 0;  // safety

    Double_t power1 = A1 * (std::pow(1.0 + x[0]/x01, -n1)); //par[7] is C_0, their combined normalization

    // --- Second Power-law part ---
    Double_t A2  = par[8];
    Double_t x02 = par[9];   // scale
    Double_t n2  = par[10];

    Double_t power2 = A2 * std::pow(1.0 + x[0]/x02, -n2);

    return par[7] * (emg + power1) + power2;
}

//############################### Functions for Plotting ###############################
Double_t EMG_component(Double_t *x, Double_t *par){
    // EMG uses par[0..3]
    Double_t emg = EMG(x, par);

    // Apply global normalization
    return par[7] * emg;
}
Double_t power1_component(Double_t *x, Double_t *par){
    Double_t A1  = par[4];
    Double_t x01 = par[5];
    Double_t n1  = par[6];

    if (x01 <= 0) return 0;

    Double_t power1 = A1 * std::pow(1.0 + x[0]/x01, -n1);

    // Same global normalization as EMG
    return par[7] * power1;
}
Double_t power2_component(Double_t *x, Double_t *par){
    Double_t A2  = par[8];
    Double_t x02 = par[9];
    Double_t n2  = par[10];

    if (x02 <= 0) return 0;

    return A2 * std::pow(1.0 + x[0]/x02, -n2);
}
Double_t core_component(Double_t *x, Double_t *par){
    // EMG uses par[0..3]
    Double_t emg = EMG(x, par);

    Double_t A1  = par[4];
    Double_t x01 = par[5];
    Double_t n1  = par[6];
    Double_t power1 = A1 * std::pow(1.0 + x[0]/x01, -n1);

    // Apply global normalization
    return par[7] * (emg + power1);
}

//############################### Numerical Functions ###############################
TGraphErrors *HistToGraph(TH1D *h, double *meanN, double *par_err){ //used to collect params for plot
    TGraphErrors *g = new TGraphErrors();
    for (int i = 0; i < h->GetNbinsX(); i++) {
        g->SetPoint(i, meanN[i], h->GetBinContent(i+1));
        g->SetPointError(i, par_err[i], h->GetBinError(i+1));
    }
    return g;
}


//#################################################################### Macro ####################################################################
void kperp_test(){
    bool save_plots = false;

    TFile *f = TFile::Open("dimuon_hists.root");
    static const int neutron_bin = 13;
    static const int neutron1    = 7;
    const char* bin_names[neutron_bin] = {"0n0n", "1n0n", "2n0n", "Xn0n", "Yn0n", "Zn0n", "1n1n", "Xn1n", "Yn1n", "Zn1n", "XnXn", "YnYn", "ZnZn"}; //Organized as 0n, 1n, Nn

    TH1D *kperp_hist[neutron_bin];
    TH1D *kperp_lin[neutron_bin];
    TH1D *kperp_prapid[neutron_bin];
    TH1D *kperp_nrapid[neutron_bin];
    TH1D *n_stat[neutron_bin];

    TH1D *onon_spare={nullptr};
    TH1D *yn1n_spare1={nullptr};
    TH1D *yn1n_spare2={nullptr};
    TH1D *zn1n_spare1={nullptr};
    TH1D *zn1n_spare2={nullptr};

    TGraphErrors *g_sigma  = new TGraphErrors();
    TGraphErrors *g_x0     = new TGraphErrors();
    TGraphErrors *g_c0     = new TGraphErrors();
    TGraphErrors *g_c1     = new TGraphErrors();
    TGraphErrors *g_emgrms = new TGraphErrors();

    TGraphErrors *fit_c0_0 = new TGraphErrors();
    TGraphErrors *fit_c0_N = new TGraphErrors();
    TGraphErrors *fit_c1_0 = new TGraphErrors();
    TGraphErrors *fit_c1_N = new TGraphErrors();

    double fitp_z0_2;
    double fitp_z0_3;
    double fitp_z0_7;
    double fitp_z0_8;

    double fitn_z0_2;
    double fitn_z0_3;
    double fitn_z0_7;
    double fitn_z0_8;

    double fitp_z0_2_err;
    double fitp_z0_3_err;
    double fitp_z0_7_err;
    double fitp_z0_8_err;

    double fitn_z0_2_err;
    double fitn_z0_3_err;
    double fitn_z0_7_err;
    double fitn_z0_8_err;

    TH1D *fit_params[11]; // 6 = # of params
    for (int j = 0; j < 11; j++) {
        fit_params[j] = new TH1D(Form("fit_params_%d", j),"",neutron_bin, 0, neutron_bin);
        for (int i = 0; i < neutron_bin; i++) {
            fit_params[j]->GetXaxis()->SetBinLabel(i+1, bin_names[i]);
        }
    }
    TH1D *M2_hist = new TH1D("M2_hist", "",neutron_bin, 0, neutron_bin);
    for (int i = 0; i < neutron_bin; i++) {
        M2_hist->GetXaxis()->SetBinLabel(i+1, bin_names[i]);
    }

    for(int i = 0; i < neutron_bin; i++){
        kperp_hist[i]    = (TH1D*)f->Get(Form("h_kperp_%d", i));
        kperp_lin[i]     = (TH1D*)f->Get(Form("h_k_lin_%d", i));
        kperp_prapid[i]  = (TH1D*)f->Get(Form("h_kperp_prapid_%d", i));
        kperp_nrapid[i]  = (TH1D*)f->Get(Form("h_kperp_nrapid_%d", i));
        n_stat[i]       = (TH1D*)f->Get(Form("neutron_stat_%d", i));

        kperp_hist[i]->Sumw2();
        kperp_lin[i]->Sumw2();
        kperp_prapid[i]->Sumw2();
        kperp_nrapid[i]->Sumw2();
    }

    double zn0n_pevent_count = kperp_prapid[5]->GetEntries();
    double zn0n_nevent_count = kperp_nrapid[5]->GetEntries();
    cout << "Entry counts: " << zn0n_pevent_count << " and " << zn0n_nevent_count << endl;

    vector<double> n_means;
    vector<double> n_means_err;
    for (int i = 0; i < neutron_bin; i++) {
        double mean = n_stat[i]->GetMean();
        double rms  = n_stat[i]->GetRMS();
        double n    = n_stat[i]->GetEntries();
    
        n_means.push_back(mean);
        n_means_err.push_back(rms / sqrt(n));
    }

    onon_spare = (TH1D*)f->Get("h_kperp_0")->Clone("onon_spare");
    onon_spare->SetDirectory(nullptr); // to make sure that we're not recycling the pointer
    yn1n_spare1 = (TH1D*)f->Get("h_kperp_8")->Clone("yn1n_spare1");
    yn1n_spare1->SetDirectory(nullptr); 
    yn1n_spare2 = (TH1D*)f->Get("h_kperp_8")->Clone("yn1n_spare2");
    yn1n_spare2->SetDirectory(nullptr); 
    zn1n_spare1 = (TH1D*)f->Get("h_kperp_11")->Clone("zn1n_spare1");
    zn1n_spare1->SetDirectory(nullptr);
    zn1n_spare2 = (TH1D*)f->Get("h_kperp_11")->Clone("zn1n_spare2");
    zn1n_spare2->SetDirectory(nullptr);

//############################### 0n0n Fit and Template ###############################
    TCanvas *k_onon_emg = new TCanvas("k_onon_emg", "0n0n Linear Fit", 1200, 1000);
    TF1 *f_emg = new TF1("f_emg", EMG_plus_power, 1, 1e4, 11);

    onon_spare->Scale(1.0 / onon_spare->Integral());
    onon_spare->Scale(1.0, "width");
    
    //integral, tau, sigma, x0 for EMG
    f_emg->SetParameters(2*onon_spare->Integral(), 100.0, 6, 0, 0.0001, 100, 2, 0);
    f_emg->SetParLimits(2, 1e-3, 50);
    f_emg->SetParLimits(3, 0, 50);
    f_emg->FixParameter(7,1); //C_0. Must be 1
    f_emg->FixParameter(8,0);
    f_emg->FixParameter(9,1);
    f_emg->FixParameter(10,1);

    onon_spare->Draw("e");
    onon_spare->Fit(f_emg, "QRI");
    f_emg->SetLineColor(kBlue);
    f_emg->Draw("same");

    onon_spare->SetMinimum(5*1e-6);
    onon_spare->SetMaximum(1e-1);
    gPad->SetLogy();
    gPad->SetLogx();

    const double tau_1 = f_emg->GetParameter(1);
    const double pow_scale = f_emg->GetParameter(5);
    const double pow_exp = f_emg->GetParameter(6);

    double integral_1 = f_emg->GetParameter(0);
    double sigma1 = f_emg->GetParameter(2);
    double emg3 = f_emg->GetParameter(3); // x_0
    double emg4 = f_emg->GetParameter(4); //power law scaling factor

    onon_spare->GetXaxis()->SetTitleOffset(1.4);
    onon_spare->GetYaxis()->SetTitle("(1/N) dN/dk_{#perp}");
    onon_spare->GetXaxis()->SetTitle("MeV");

    const Double_t *p  = f_emg->GetParameters();
    const Double_t *pe = f_emg->GetParErrors();
    TF1 *f_emg_p   = new TF1("emg", EMG_component, 1, 1e4, 11);
    TF1 *f_pow1  = new TF1("pow1", power1_component, 1, 1e4, 11);

    f_emg_p->SetParameters(p);
    f_pow1->SetParameters(p);
    f_emg_p->SetLineColor(kBlue);
    f_pow1->SetLineColor(kGreen+2);

    TLatex latex;
    latex.SetNDC(); // use normalized coordinates (0–1)
    latex.SetTextSize(0.05);
    latex.SetTextAlign(22); // center horizontally
    latex.DrawLatex(0.77, 0.8, Form("%s", bin_names[0]));

    f_emg_p->Draw("SAME");
    f_pow1->Draw("SAME");


    if (save_plots == true){
        k_onon_emg->SaveAs("plots/full/kperp_0n0n_6.30.png");
    }
    //cout << tau_1 << " and " << pow_scale << " and " << pow_exp << endl;
    cout << "0n0n Fit" << endl;

//############################### Yn1n & Zn1n to vary x_0 ###############################
/*
    TCanvas *yzn1n_plot = new TCanvas("yzn1n_plot", "Yn1n x_0 Vary", 1200, 1000);
    yzn1n_plot->Divide(2,2);
    
    yn1n_spare1->Scale(1.0 / yn1n_spare1->Integral());
    yn1n_spare2->Scale(1.0 / yn1n_spare2->Integral());
    zn1n_spare1->Scale(1.0 / zn1n_spare1->Integral());
    zn1n_spare2->Scale(1.0 / zn1n_spare2->Integral());

    yn1n_spare1->Scale(1.0, "width");
    yn1n_spare2->Scale(1.0, "width");
    zn1n_spare1->Scale(1.0, "width");
    zn1n_spare2->Scale(1.0, "width");

    yzn1n_plot->cd(1);
    gPad->SetLogy();
    gPad->SetLogx();
    TF1 *yn1n_fit_1 = new TF1("yn1n_fit_1", EMG_plus_power, 1, 1e4, 11);
    yn1n_fit_1->SetParameters(0,0,sigma1,emg3,0,0,0, integral_1, 1e-2, pow_scale, pow_exp);
    yn1n_fit_1->FixParameter(0,integral_1); //emg integral
    yn1n_fit_1->FixParameter(1,tau_1); //emg tau
    yn1n_fit_1->SetParLimits(2, 1e-3, 50);
    yn1n_fit_1->SetParLimits(3, 0, 50);
    //yn1n_fit_1->FixParameter(2,sigma1); //emg sigma
    //yn1n_fit_1->FixParameter(3,emg3); //emg x_0
    yn1n_fit_1->FixParameter(4,emg4/integral_1); //pow1 factor
    yn1n_fit_1->FixParameter(5,pow_scale); //pow1 fraction
    yn1n_fit_1->FixParameter(6,pow_exp); //pow1 exponent

    yn1n_spare1->Fit(yn1n_fit_1, "QRI");
    yn1n_spare1->SetMinimum(5*1e-6);
    yn1n_spare1->SetMaximum(1e-1);
    yn1n_spare1->Draw("HIST");
    yn1n_fit_1->Draw("same");
    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.05);
    latex1.SetTextAlign(22);
    latex1.DrawLatex(0.77, 0.8, "Yn1n, free x_{0}");

    yzn1n_plot->cd(2);
    gPad->SetLogy();
    gPad->SetLogx();
    TF1 *zn1n_fit_1 = new TF1("zn1n_fit_1", EMG_plus_power, 1, 1e4, 11);
    zn1n_fit_1->SetParameters(0,0,sigma1,emg3,0,0,0, integral_1, 1e-2, pow_scale, pow_exp);
    zn1n_fit_1->FixParameter(0,integral_1); //emg integral
    zn1n_fit_1->FixParameter(1,tau_1); //emg tau
    zn1n_fit_1->SetParLimits(2, 1e-3, 50);
    zn1n_fit_1->SetParLimits(3, 0, 50);
    zn1n_fit_1->FixParameter(4,emg4/integral_1); //pow1 factor
    zn1n_fit_1->FixParameter(5,pow_scale); //pow1 fraction
    zn1n_fit_1->FixParameter(6,pow_exp); //pow1 exponent
    zn1n_spare1->Fit(zn1n_fit_1, "QRI");
    zn1n_spare1->SetMinimum(5*1e-6);
    zn1n_spare1->SetMaximum(1e-1);
    zn1n_spare1->Draw("HIST");
    zn1n_fit_1->Draw("same");
    //zn1n_spare1->GetXaxis()->SetTitle("Zn1n, free x_{0}");
    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.05);
    latex2.SetTextAlign(22);
    latex2.DrawLatex(0.77, 0.8, "Zn1n, free x_{0}");

    yzn1n_plot->cd(3);
    gPad->SetLogy();
    gPad->SetLogx();
    TF1 *yn1n_fit_2 = new TF1("yn1n_fit_2", EMG_plus_power, 1, 1e4, 11);
    yn1n_fit_2->FixParameter(2, emg3);
    yn1n_fit_2->SetParameters(0,0,sigma1,emg3,0,0,0, integral_1, 1e-2, pow_scale, pow_exp);
    yn1n_fit_2->FixParameter(0,integral_1); //emg integral
    yn1n_fit_2->FixParameter(1,tau_1); //emg tau
    yn1n_fit_2->SetParLimits(2, 1e-3, 50);
    yn1n_fit_2->FixParameter(3,emg3); //emg x_0
    yn1n_fit_2->FixParameter(4,emg4/integral_1); //pow1 factor
    yn1n_fit_2->FixParameter(5,pow_scale); //pow1 fraction
    yn1n_fit_2->FixParameter(6,pow_exp); //pow1 exponent
    yn1n_spare2->Fit(yn1n_fit_2, "QRI");
    yn1n_spare2->SetMinimum(5*1e-6);
    yn1n_spare2->SetMaximum(1e-1);
    yn1n_spare2->Draw("HIST");
    yn1n_fit_2->Draw("same");
    //yn1n_spare2->GetXaxis()->SetTitle("Yn1n, fixed x_{0}");
    TLatex latex3;
    latex3.SetNDC();
    latex3.SetTextSize(0.05);
    latex3.SetTextAlign(22);
    latex3.DrawLatex(0.77, 0.8, "Yn1n, fixed x_{0}");

    yzn1n_plot->cd(4);
    gPad->SetLogy();
    gPad->SetLogx();
    TF1 *zn1n_fit_2 = new TF1("zn1n_fit_2", EMG_plus_power, 1, 1e4, 11);
    zn1n_fit_2->FixParameter(2, emg3);
    zn1n_fit_2->SetParameters(0,0,sigma1,emg3,0,0,0, integral_1, 1e-2, pow_scale, pow_exp);
    zn1n_fit_2->FixParameter(0,integral_1); //emg integral
    zn1n_fit_2->FixParameter(1,tau_1); //emg tau
    zn1n_fit_2->SetParLimits(2, 1e-3, 50);
    zn1n_fit_2->FixParameter(3,emg3); //emg x_0
    zn1n_fit_2->FixParameter(4,emg4/integral_1); //pow1 factor
    zn1n_fit_2->FixParameter(5,pow_scale); //pow1 fraction
    zn1n_fit_2->FixParameter(6,pow_exp); //pow1 exponent
    zn1n_spare2->Fit(zn1n_fit_2, "QRI");
    zn1n_spare2->SetMinimum(5*1e-6);
    zn1n_spare2->SetMaximum(1e-1);
    zn1n_spare2->Draw("HIST");
    zn1n_fit_2->Draw("same");
    TLatex latex4;
    latex4.SetNDC();
    latex4.SetTextSize(0.05);
    latex4.SetTextAlign(22);
    latex4.DrawLatex(0.77, 0.8, "Zn1n, fixed x_{0}");

    cout << "Y and Z fits" << endl;
    */
/*
//############################### XnXn on Linear Scale ###############################
    TCanvas *xnxn_lin_solo = new TCanvas("xnxn_lin_solo", "XnXn Linear", 1200, 1000);

    kperp_lin[3]->Scale(1.0 / onon_spare->Integral());
    kperp_lin[3]->Scale(1.0, "width");

    kperp_lin[3]->Draw("HIST");

    kperp_lin[3]->GetXaxis()->SetTitleOffset(1.4);
    kperp_lin[3]->GetYaxis()->SetTitle("(1/N) dN/dk_{#perp}");
    kperp_lin[3]->GetXaxis()->SetTitle("MeV");

    if (save_plots == true){
        xnxn_lin_solo->SaveAs("plots/full/kperp_xnxn_6.30.png");
    }
*/

//############################### kperp 1 ###############################
    TCanvas *perps1 = new TCanvas("perps1", "kperp1", 1200, 1000);
    perps1->Divide(3, 2);
    bool draw_bits1 = true;
    for (int i = 0; i < neutron1; i++){
        TF1 *kperp_fit1 = new TF1(Form("kperp_fit1_%d",i), EMG_plus_power, 1, 1e4, 11);

        kperp_hist[i]->Scale(1.0 / kperp_hist[i]->Integral()); // 1/N dN/dk_perp
        kperp_hist[i]->Scale(1.0, "width");

        perps1->cd(i+1); // starting at 1, so do not include +1
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetRightMargin(0.12);

        kperp_fit1->SetParameters(0,0,sigma1,emg3,0,0,0, integral_1, 1e-2, pow_scale, pow_exp);
        kperp_fit1->FixParameter(0,1); //emg integral. Must be one for correct normalization
        //kperp_fit1->FixParameter(0,integral_1); //emg integral. Testing to not be 1
        kperp_fit1->FixParameter(1,tau_1); //emg tau
        kperp_fit1->SetParLimits(2, 1e-3, 50);
        kperp_fit1->SetParLimits(3, 0, 50);
        kperp_fit1->FixParameter(4,emg4/integral_1); //pow1 factor
        //kperp_fit1->FixParameter(4,emg4); //pow1 factor. Testing to not be normalized
        kperp_fit1->FixParameter(5,pow_scale); //pow1 fraction
        kperp_fit1->FixParameter(6,pow_exp); //pow1 exponent
        kperp_fit1->SetParLimits(8, 0, 1);
        kperp_fit1->FixParameter(9,250); //power law denominator. Brian suggested fixing to 250
        kperp_fit1->FixParameter(10,2); //this is omega. Brian suggested fix at 2

        //kperp_fit1->SetParLimits(9, 0, 1e5);
        //kperp_fit1->SetParLimits(10, 0, 1e3);


        TFitResultPtr fitResult = kperp_hist[i]->Fit(kperp_fit1, "QRI S");

        kperp_hist[i]->SetMinimum(5*1e-6);
        kperp_hist[i]->SetMaximum(1e-1);
        kperp_hist[i]->Draw("e");
        kperp_fit1->Draw("same");

        double par7 = kperp_fit1->GetParameter(7);
        double par8 = kperp_fit1->GetParameter(8);
        double par7_norm = par7 / (par7 + par8);
        double par8_norm = par8 / (par7 + par8);
        //cout << bin_names[i] << ": " << par7_norm << " + " << par8_norm << " = " << par7_norm + par8_norm << endl;
        
        const Double_t *p  = kperp_fit1->GetParameters();
        const Double_t *pe = kperp_fit1->GetParErrors();
        TF1 *f_emg   = new TF1(Form("emg_%d", i), EMG_component, 1, 1e4, 11);
        TF1 *f_pow1  = new TF1(Form("pow1_%d", i), power1_component, 1, 1e4, 11);
        TF1 *f_pow2  = new TF1(Form("pow2_%d", i), power2_component, 1, 1e4, 11);
        
        TF1 *f_core   = new TF1(Form("core_%d", i), core_component, 1, 1e4, 11);

        f_emg->SetParameters(p);
        f_pow1->SetParameters(p);
        f_pow2->SetParameters(p);
        f_core->SetParameters(p);

        double tot = kperp_fit1->Integral(0, 1e4);
        double core_tot = f_core->Integral(0, 1e4);
        double tail_tot = f_pow2->Integral(0, 1e4);
        double f_core_val = core_tot / tot;
        double f_tail_val = tail_tot / tot;
        //cout << bin_names[i] << ": " << core_tot/tot << " + " << tail_tot/tot << " = " << (core_tot + tail_tot)/tot << endl;

        // ************************* Get the errors on C_0 and C_1 *************************
        TMatrixDSym cov = fitResult->GetCovarianceMatrix();
        // Numerical gradients of core_tot and tail_tot w.r.t. each parameter
        double grad_core[11] = {0};
        double grad_tail[11] = {0};
        double dp = 1e-4;

        for (int k = 0; k < 11; ++k) {
            double p_up[11], p_dn[11];
            std::copy(p, p+11, p_up);
            std::copy(p, p+11, p_dn);
            p_up[k] *= (1.0 + dp);
            p_dn[k] *= (1.0 - dp);

            f_core->SetParameters(p_up);
            double core_up = f_core->Integral(0, 1e4);
            f_core->SetParameters(p_dn);
            double core_dn = f_core->Integral(0, 1e4);
            grad_core[k] = (core_up - core_dn) / (2.0 * dp * p[k]);

            f_pow2->SetParameters(p_up);
            double tail_up = f_pow2->Integral(0, 1e4);
            f_pow2->SetParameters(p_dn);
            double tail_dn = f_pow2->Integral(0, 1e4);
            grad_tail[k] = (tail_up - tail_dn) / (2.0 * dp * p[k]);
        }

        // Gradient of the fraction f_core = core_tot/tot
        double grad_fcore[11], grad_ftail[11];
        for (int k = 0; k < 11; ++k) {
            double grad_tot_k = grad_core[k] + grad_tail[k];  // d(tot)/dp_k
            grad_fcore[k] = (grad_core[k] - f_core_val * grad_tot_k) / tot;
            grad_ftail[k] = (grad_tail[k] - f_tail_val * grad_tot_k) / tot;
        }

        // Propagate through full covariance matrix
        double var_fcore = 0, var_ftail = 0;
        for (int a = 0; a < 11; ++a) {
            for (int b = 0; b < 11; ++b) {
                var_fcore += grad_fcore[a] * cov[a][b] * grad_fcore[b];
                var_ftail += grad_ftail[a] * cov[a][b] * grad_ftail[b];
            }
        }

        double sigma_fcore = std::sqrt(var_fcore);
        double sigma_ftail = std::sqrt(var_ftail);


        for (int j = 0; j < 11; j++) {
            if (j != 7 && j != 8){
                fit_params[j]->SetBinContent(i+1, p[j]);
                fit_params[j]->SetBinError(i+1, pe[j]);
            }
            if (j == 7){
                fit_params[j]->SetBinContent(i+1, core_tot/tot);
                fit_params[j]->SetBinError(i+1, sigma_fcore);
            }
            if (j == 8){
                fit_params[j]->SetBinContent(i+1, tail_tot/tot);
                fit_params[j]->SetBinError(i+1, sigma_ftail); 
            }
        }

        TF1 *emgOnly = new TF1(Form("emgOnly_%d", i), EMG, 0, 1e4, 4); // this is to get the second moment of the EMG
        for (int i = 0; i < 4; ++i)
            emgOnly->SetParameter(i, kperp_fit1->GetParameter(i));
        double norm = emgOnly->Integral(0, 1e4);
        TF1 *x2emg = new TF1("x2emg", [&](double *x, double *p){return x[0]*x[0]*EMG(x,p);}, 0, 1e4, 4);
        for (int i = 0; i < 4; ++i)
            x2emg->SetParameter(i, kperp_fit1->GetParameter(i));
        double M2 = x2emg->Integral(0, 1e4)/norm;
        double sqrtM2 = sqrt(M2);
        M2_hist->SetBinContent(i+1, sqrtM2); //get the RMS now just by taking the root

        // ************************* Get the errors on variance *************************
        double grad_M2[11] = {0};  // only [0..3] will be nonzero
        for (int k = 0; k < 4; ++k) {  // only EMG params
            double p_up[4], p_dn[4];
            for (int j = 0; j < 4; ++j) {
                p_up[j] = kperp_fit1->GetParameter(j);
                p_dn[j] = kperp_fit1->GetParameter(j);
            }
            // guard against p[k] == 0
            double step = (std::abs(p_up[k]) > 1e-10) ? dp * std::abs(p_up[k]) : 1e-8;
            p_up[k] += step;
            p_dn[k] -= step;

            emgOnly->SetParameters(p_up);
            x2emg->SetParameters(p_up);
            double norm_up = emgOnly->Integral(0, 1e4);
            double M2_up   = x2emg->Integral(0, 1e4) / norm_up;

            emgOnly->SetParameters(p_dn);
            x2emg->SetParameters(p_dn);
            double norm_dn = emgOnly->Integral(0, 1e4);
            double M2_dn   = x2emg->Integral(0, 1e4) / norm_dn;

            grad_M2[k] = (M2_up - M2_dn) / (2.0 * step);
        }

        // Propagate — use full 11x11 cov but only [0..3] entries contribute
        double var_M2 = 0;
        for (int a = 0; a < 11; ++a)
            for (int b = 0; b < 11; ++b)
                var_M2 += grad_M2[a] * cov[a][b] * grad_M2[b];

        double sigma_M2 = std::sqrt(var_M2);
        double sigma_sqrtM2 = sigma_M2 / (2.0 * sqrtM2);
        M2_hist->SetBinError(i+1, sigma_sqrtM2);
        // *************************  *************************

        for (int i = 0; i < neutron1; i++) {
            double y2    = fit_params[2]->GetBinContent(i+1);
            double ey2   = fit_params[2]->GetBinError(i+1);
            double y3    = fit_params[3]->GetBinContent(i+1);
            double ey3   = fit_params[3]->GetBinError(i+1);


            g_sigma->SetPoint(i, n_means[i], y2);
            g_sigma->SetPointError(i, n_means_err[i], ey2);
            g_x0->SetPoint(i, n_means[i], y3);
            g_x0->SetPointError(i, n_means_err[i], ey3);
        }
        double y7    = core_tot/tot;
        double ey7   = sigma_fcore;
        double y8    = tail_tot/tot;
        double ey8   = sigma_ftail;
        double yrms  = sqrtM2;
        double eyrms = sigma_sqrtM2;

        g_c0->SetPoint(i, n_means[i], y7);
        g_c0->SetPointError(i, n_means_err[i], ey7);
        g_c1->SetPoint(i, n_means[i], y8);
        g_c1->SetPointError(i, n_means_err[i], ey8);
        g_emgrms->SetPoint(i, n_means[i], yrms);
        g_emgrms->SetPointError(i, n_means_err[i], eyrms);
        
        if (i < 6 && i != 0){
            fit_c0_0->SetPoint(i, n_means[i], y7);
            fit_c0_0->SetPointError(i, n_means_err[i], ey7);
            fit_c1_0->SetPoint(i, n_means[i], y8);
            fit_c1_0->SetPointError(i, n_means_err[i], ey8);
        }
        else if (i >= 6){
            fit_c0_N->SetPoint(i, n_means[i], y7);
            fit_c0_N->SetPointError(i, n_means_err[i], ey7);
            fit_c1_N->SetPoint(i, n_means[i], y8);
            fit_c1_N->SetPointError(i, n_means_err[i], ey8);
        }

        f_emg->SetLineColor(kBlue);
        f_pow1->SetLineColor(kGreen+2);
        f_pow2->SetLineColor(kMagenta);

        if (draw_bits1 == true){
            f_emg->Draw("same");
            f_pow1->Draw("same");
            f_pow2->Draw("same");
        }
        
        kperp_hist[i]->GetXaxis()->SetTitleOffset(1.4);
        kperp_hist[i]->GetYaxis()->SetTitle("(1/N) dN/dk_{#perp}");
        kperp_hist[i]->GetXaxis()->SetTitle("MeV");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.77, 0.8, Form("%s", bin_names[i]));

        if (i == 0){
            latex.DrawLatex(0.75, 0.73, "X #geq 1");
            latex.DrawLatex(0.75, 0.66, "1 < Y #leq 5");
            latex.DrawLatex(0.75, 0.59, "5 < Z #leq 10");
            latex.SetTextSize(0.05);
            latex.DrawLatex(0.75, 0.88, "k_{#perp}");

            TLegend *line_names = new TLegend(0.15, 0.2, 0.9, 0.4);
            line_names->SetBorderSize(0);
            line_names->SetFillStyle(0);
            line_names->SetTextSize(0.05);
            line_names->AddEntry(f_emg, "Core", "l");
            line_names->AddEntry(f_pow1, "QED Tail", "l");
            line_names->AddEntry(f_pow2, "Dis Tail", "l");
            line_names->Draw();
        }
    }
    if (save_plots == true){
        perps1->SaveAs("plots/full/kperp11_6.30.png");
    }    
    cout << "Kperp1 Fit" << endl;

  //############################### kperp 2 ###############################
    TCanvas *perps2 = new TCanvas("perps2", "kperp2", 1200, 1000);
    perps2->Divide(3, 2);
    bool draw_bits2 = true; 

    for (int i = neutron1; i < neutron_bin; i++){
        TF1 *kperp_fit2 = new TF1(Form("kperp_fit2_%d",i), EMG_plus_power, 1, 1e4, 11);

        kperp_hist[i]->Scale(1.0 / kperp_hist[i]->Integral()); // 1/N dN/dk_perp
        kperp_hist[i]->Scale(1.0, "width");

        perps2->cd(i- (neutron1-1));
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetRightMargin(0.12);

        kperp_fit2->SetParameters(0,0,sigma1,emg3,0,0,0, integral_1, 1e-2, pow_scale, pow_exp);
        kperp_fit2->FixParameter(0,1); //emg integral
        //kperp_fit2->FixParameter(0,integral_1); //emg integral. Testing to not be 1
        kperp_fit2->FixParameter(1,tau_1); //emg tau
        kperp_fit2->SetParLimits(2, 1e-3, 50);
        kperp_fit2->SetParLimits(3, 0, 50);
        kperp_fit2->FixParameter(4,emg4/integral_1); //pow1 factor
        //kperp_fit2->FixParameter(4,emg4); //pow1 factor. Testing to not be normalized
        kperp_fit2->FixParameter(5,pow_scale); //pow1 fraction
        kperp_fit2->FixParameter(6,pow_exp); //pow1 exponent
        kperp_fit2->SetParLimits(8, 0, 1);
        kperp_fit2->FixParameter(9,250); //power law denominator. Brian suggested fixing to 250
        kperp_fit2->FixParameter(10,2); //this is omega. Brian suggested fix at 2

        //kperp_fit2->SetParLimits(9, 0, 1e5);
        //kperp_fit2->SetParLimits(10, 0, 1e3);

        TFitResultPtr fitResult2 = kperp_hist[i]->Fit(kperp_fit2, "QRI S");

        kperp_hist[i]->SetMinimum(5*1e-6);
        kperp_hist[i]->SetMaximum(1e-1);
        kperp_hist[i]->Draw("e");
        kperp_fit2->Draw("same");

        const Double_t *p  = kperp_fit2->GetParameters();
        const Double_t *pe = kperp_fit2->GetParErrors();

        TF1 *f_emg   = new TF1(Form("emg_%d", i), EMG_component, 1, 1e4, 11);
        TF1 *f_pow1  = new TF1(Form("pow1_%d", i), power1_component, 1, 1e4, 11);
        TF1 *f_pow2  = new TF1(Form("pow2_%d", i), power2_component, 1, 1e4, 11);
        
        TF1 *f_core   = new TF1(Form("core_%d", i), core_component, 1, 1e4, 11);

        f_emg->SetParameters(p);
        f_pow1->SetParameters(p);
        f_pow2->SetParameters(p);
        f_core->SetParameters(p);

        double tot = kperp_fit2->Integral(0, 1e4);
        double core_tot = f_core->Integral(0, 1e4);
        double tail_tot = f_pow2->Integral(0, 1e4);
        double f_core_val = core_tot / tot;
        double f_tail_val = tail_tot / tot;
        //cout << bin_names[i] << ": " << core_tot/tot << " + " << tail_tot/tot << " = " << (core_tot + tail_tot)/tot << endl;
    

        // ************************* Get the errors on C_0 and C_1 *************************
        TMatrixDSym cov = fitResult2->GetCovarianceMatrix();
        // Numerical gradients of core_tot and tail_tot w.r.t. each parameter
        double grad_core[11] = {0};
        double grad_tail[11] = {0};
        double dp = 1e-4;

        for (int k = 0; k < 11; ++k) {
            double p_up[11], p_dn[11];
            std::copy(p, p+11, p_up);
            std::copy(p, p+11, p_dn);
            p_up[k] *= (1.0 + dp);
            p_dn[k] *= (1.0 - dp);

            f_core->SetParameters(p_up);
            double core_up = f_core->Integral(0, 1e4);
            f_core->SetParameters(p_dn);
            double core_dn = f_core->Integral(0, 1e4);
            grad_core[k] = (core_up - core_dn) / (2.0 * dp * p[k]);

            f_pow2->SetParameters(p_up);
            double tail_up = f_pow2->Integral(0, 1e4);
            f_pow2->SetParameters(p_dn);
            double tail_dn = f_pow2->Integral(0, 1e4);
            grad_tail[k] = (tail_up - tail_dn) / (2.0 * dp * p[k]);
        }


        // Gradient of the fraction f_core = core_tot/tot
        double grad_fcore[11], grad_ftail[11];
        for (int k = 0; k < 11; ++k) {
            double grad_tot_k = grad_core[k] + grad_tail[k];  // d(tot)/dp_k
            grad_fcore[k] = (grad_core[k] - f_core_val * grad_tot_k) / tot;
            grad_ftail[k] = (grad_tail[k] - f_tail_val * grad_tot_k) / tot;
        }

        // Propagate through full covariance matrix
        double var_fcore = 0, var_ftail = 0;
        for (int a = 0; a < 11; ++a) {
            for (int b = 0; b < 11; ++b) {
                var_fcore += grad_fcore[a] * cov[a][b] * grad_fcore[b];
                var_ftail += grad_ftail[a] * cov[a][b] * grad_ftail[b];
            }
        }
        double sigma_fcore = std::sqrt(var_fcore);
        double sigma_ftail = std::sqrt(var_ftail);
        // *************************  *************************



        for (int j = 0; j < 11; j++) {
            if (j != 7 && j != 8){
                fit_params[j]->SetBinContent(i+1, p[j]);
                fit_params[j]->SetBinError(i+1, pe[j]);
            }
            if (j == 7){
                fit_params[j]->SetBinContent(i+1, core_tot/tot);
                fit_params[j]->SetBinError(i+1, sigma_fcore);
            }
            if (j == 8){
                fit_params[j]->SetBinContent(i+1, tail_tot/tot);
                fit_params[j]->SetBinError(i+1, sigma_ftail); 
            }
        }

        TF1 *emgOnly = new TF1(Form("emgOnly_%d", i), EMG, 0, 1e4, 4); // this is to get the second moment of the EMG
        for (int i = 0; i < 4; ++i)
            emgOnly->SetParameter(i, kperp_fit2->GetParameter(i));
        double norm = emgOnly->Integral(0, 1e4);
        TF1 *x2emg = new TF1("x2emg", [&](double *x, double *p){return x[0]*x[0]*EMG(x,p);}, 0, 1e4, 4);
        for (int i = 0; i < 4; ++i)
            x2emg->SetParameter(i, kperp_fit2->GetParameter(i));
        double M2 = x2emg->Integral(0, 1e4)/norm;
        double sqrtM2 = sqrt(M2);
        M2_hist->SetBinContent(i+1, sqrtM2); //get the RMS now just by taking the root

        // ************************* Get the errors on variance *************************
        double grad_M2[11] = {0};  // only [0..3] will be nonzero
        for (int k = 0; k < 4; ++k) {  // only EMG params
            double p_up[4], p_dn[4];
            for (int j = 0; j < 4; ++j) {
                p_up[j] = kperp_fit2->GetParameter(j);
                p_dn[j] = kperp_fit2->GetParameter(j);
            }
            // guard against p[k] == 0
            double step = (std::abs(p_up[k]) > 1e-10) ? dp * std::abs(p_up[k]) : 1e-8;
            p_up[k] += step;
            p_dn[k] -= step;

            emgOnly->SetParameters(p_up);
            x2emg->SetParameters(p_up);
            double norm_up = emgOnly->Integral(0, 1e4);
            double M2_up   = x2emg->Integral(0, 1e4) / norm_up;

            emgOnly->SetParameters(p_dn);
            x2emg->SetParameters(p_dn);
            double norm_dn = emgOnly->Integral(0, 1e4);
            double M2_dn   = x2emg->Integral(0, 1e4) / norm_dn;

            grad_M2[k] = (M2_up - M2_dn) / (2.0 * step);
        }

        // Propagate — use full 11x11 cov but only [0..3] entries contribute
        double var_M2 = 0;
        for (int a = 0; a < 11; ++a)
            for (int b = 0; b < 11; ++b)
                var_M2 += grad_M2[a] * cov[a][b] * grad_M2[b];

        double sigma_M2 = std::sqrt(var_M2);
        double sigma_sqrtM2 = sigma_M2 / (2.0 * sqrtM2);
        M2_hist->SetBinError(i+1, sigma_sqrtM2);
        // *************************  *************************

        for (int j = 0; j < neutron_bin; j++) {
            double y2    = fit_params[2]->GetBinContent(j+1);
            double ey2   = fit_params[2]->GetBinError(j+1);
            double y3    = fit_params[3]->GetBinContent(j+1);
            double ey3   = fit_params[3]->GetBinError(j+1);

            g_sigma->SetPoint(j, n_means[j], y2);
            g_sigma->SetPointError(j, n_means_err[j], ey2);
            g_x0->SetPoint(j, n_means[j], y3);
            g_x0->SetPointError(j, n_means_err[j], ey3);
        }
        double y7    = core_tot/tot;
        double ey7   = sigma_fcore;
        double y8    = tail_tot/tot;
        double ey8   = sigma_ftail;
        double yrms  = sqrtM2;
        double eyrms = sigma_sqrtM2;
        g_c0->SetPoint(i, n_means[i], y7);
        g_c0->SetPointError(i, n_means_err[i], ey7);
        g_c1->SetPoint(i, n_means[i], y8);
        g_c1->SetPointError(i, n_means_err[i], ey8);
        g_emgrms->SetPoint(i, n_means[i], yrms);
        g_emgrms->SetPointError(i, n_means_err[i], eyrms);

        fit_c0_N->SetPoint(i, n_means[i], y7); // in kperp2, they're all not 0n events
        fit_c0_N->SetPointError(i, n_means_err[i], ey7);
        fit_c1_N->SetPoint(i, n_means[i], y8);
        fit_c1_N->SetPointError(i, n_means_err[i], ey8);

        f_emg->SetLineColor(kBlue);
        f_pow1->SetLineColor(kGreen+2);
        f_pow2->SetLineColor(kMagenta);

        if (draw_bits2 == true){
            f_emg->Draw("same");
            f_pow1->Draw("same");
            f_pow2->Draw("same");
        }

        kperp_hist[i]->GetXaxis()->SetTitleOffset(1.4);
        kperp_hist[i]->GetYaxis()->SetTitle("(1/N) dN/dk_{#perp}");
        kperp_hist[i]->GetXaxis()->SetTitle("MeV");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.77, 0.8, Form("%s", bin_names[i]));
    }
    if (save_plots == true){
        perps2->SaveAs("plots/full/kperp12_6.30.png");
    }
    cout << "Kperp2 Fit" << endl;


  //############################### kperp by signed rapidity ###############################
    TCanvas *perp_rapid = new TCanvas("perp_rapid", "kperp by rapidity", 1500, 1000);
    perp_rapid->Divide(5, 2);
    bool draw_bits3 = true; 

    for (int i = 1; i < 6; i++){ // pos rapidity
        TF1 *kperp_fit3 = new TF1(Form("kperp_fit3_%d",i), EMG_plus_power, 1, 1e4, 11);

        kperp_prapid[i]->Scale(1.0 / kperp_prapid[i]->Integral()); // 1/N dN/dk_perp
        kperp_prapid[i]->Scale(1.0, "width");

        perp_rapid->cd(i);
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetRightMargin(0.12);

        kperp_fit3->SetParameters(0,0,sigma1,emg3,0,0,0, integral_1, 1e-2, pow_scale, pow_exp);
        kperp_fit3->FixParameter(0,1); //emg integral
        kperp_fit3->FixParameter(1,tau_1); //emg tau
        kperp_fit3->SetParLimits(2, 1e-3, 50);
        kperp_fit3->SetParLimits(3, 0, 50);
        kperp_fit3->FixParameter(4,emg4/integral_1); //pow1 factor
        kperp_fit3->FixParameter(5,pow_scale); //pow1 fraction
        kperp_fit3->FixParameter(6,pow_exp); //pow1 exponent
        kperp_fit3->SetParLimits(8, 0, 1);
        kperp_fit3->FixParameter(9,250); //power law denominator. Brian suggested fixing to 250
        kperp_fit3->FixParameter(10,2); //this is omega. Brian suggested fix at 2

        TFitResultPtr fitResult3 = kperp_prapid[i]->Fit(kperp_fit3, "QRI S");

        kperp_prapid[i]->SetMinimum(5*1e-6);
        kperp_prapid[i]->SetMaximum(1e-1);
        kperp_prapid[i]->Draw("e");
        kperp_fit3->Draw("same");

        const Double_t *p  = kperp_fit3->GetParameters();
        const Double_t *pe = kperp_fit3->GetParErrors();

        if (i == 5){
            fitp_z0_2 = p[2];
            fitp_z0_3 = p[3];
            fitp_z0_7 = p[7];
            fitp_z0_8 = p[8];

            fitp_z0_2_err = pe[2];
            fitp_z0_3_err = pe[3];
            fitp_z0_7_err = pe[7];
            fitp_z0_8_err = pe[8];
        }

        TF1 *f_emg   = new TF1(Form("emg_%d", i), EMG_component, 1, 1e4, 11);
        TF1 *f_pow1  = new TF1(Form("pow1_%d", i), power1_component, 1, 1e4, 11);
        TF1 *f_pow2  = new TF1(Form("pow2_%d", i), power2_component, 1, 1e4, 11);
        
        TF1 *f_core   = new TF1(Form("core_%d", i), core_component, 1, 1e4, 11);

        f_emg->SetParameters(p);
        f_pow1->SetParameters(p);
        f_pow2->SetParameters(p);
        f_core->SetParameters(p);

        f_emg->SetLineColor(kBlue);
        f_pow1->SetLineColor(kGreen+2);
        f_pow2->SetLineColor(kMagenta);

        if (draw_bits3 == true){
            f_emg->Draw("same");
            f_pow1->Draw("same");
            f_pow2->Draw("same");
        }

        kperp_prapid[i]->GetXaxis()->SetTitleOffset(1.4);
        kperp_prapid[i]->GetYaxis()->SetTitle("(1/N) dN/dk_{#perp}");
        kperp_prapid[i]->GetXaxis()->SetTitle("MeV");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.77, 0.8, Form("%s", bin_names[i]));
        if (i==1){
            latex.DrawLatex(0.68, 0.75, "Positive Rapidity");
        }
    }
    for (int i = 1; i < 6; i++){ // neg rapidity
        TF1 *kperp_fit3 = new TF1(Form("kperp_fit3_%d",i), EMG_plus_power, 1, 1e4, 11);

        kperp_nrapid[i]->Scale(1.0 / kperp_nrapid[i]->Integral()); // 1/N dN/dk_perp
        kperp_nrapid[i]->Scale(1.0, "width");

        perp_rapid->cd(i+5);
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetRightMargin(0.12);

        kperp_fit3->SetParameters(0,0,sigma1,emg3,0,0,0, integral_1, 1e-2, pow_scale, pow_exp);
        kperp_fit3->FixParameter(0,1); //emg integral
        kperp_fit3->FixParameter(1,tau_1); //emg tau
        kperp_fit3->SetParLimits(2, 1e-3, 50);
        kperp_fit3->SetParLimits(3, 0, 50);
        kperp_fit3->FixParameter(4,emg4/integral_1); //pow1 factor
        kperp_fit3->FixParameter(5,pow_scale); //pow1 fraction
        kperp_fit3->FixParameter(6,pow_exp); //pow1 exponent
        kperp_fit3->SetParLimits(8, 0, 1);
        kperp_fit3->FixParameter(9,250); //power law denominator. Brian suggested fixing to 250
        kperp_fit3->FixParameter(10,2); //this is omega. Brian suggested fix at 2

        TFitResultPtr fitResult3 = kperp_nrapid[i]->Fit(kperp_fit3, "QRI S");

        kperp_nrapid[i]->SetMinimum(5*1e-6);
        kperp_nrapid[i]->SetMaximum(1e-1);
        kperp_nrapid[i]->Draw("e");
        kperp_fit3->Draw("same");

        const Double_t *p  = kperp_fit3->GetParameters();
        const Double_t *pe = kperp_fit3->GetParErrors();

        if (i == 5){
            fitn_z0_2 = p[2];
            fitn_z0_3 = p[3];
            fitn_z0_7 = p[7];
            fitn_z0_8 = p[8];

            fitn_z0_2_err = pe[2];
            fitn_z0_3_err = pe[3];
            fitn_z0_7_err = pe[7];
            fitn_z0_8_err = pe[8];
        }

        TF1 *f_emg   = new TF1(Form("emg_%d", i), EMG_component, 1, 1e4, 11);
        TF1 *f_pow1  = new TF1(Form("pow1_%d", i), power1_component, 1, 1e4, 11);
        TF1 *f_pow2  = new TF1(Form("pow2_%d", i), power2_component, 1, 1e4, 11);
        
        TF1 *f_core   = new TF1(Form("core_%d", i), core_component, 1, 1e4, 11);

        f_emg->SetParameters(p);
        f_pow1->SetParameters(p);
        f_pow2->SetParameters(p);
        f_core->SetParameters(p);

        f_emg->SetLineColor(kBlue);
        f_pow1->SetLineColor(kGreen+2);
        f_pow2->SetLineColor(kMagenta);

        if (draw_bits3 == true){
            f_emg->Draw("same");
            f_pow1->Draw("same");
            f_pow2->Draw("same");
        }

        kperp_nrapid[i]->GetXaxis()->SetTitleOffset(1.4);
        kperp_nrapid[i]->GetYaxis()->SetTitle("(1/N) dN/dk_{#perp}");
        kperp_nrapid[i]->GetXaxis()->SetTitle("MeV");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.77, 0.8, Form("%s", bin_names[i]));
        if (i==1){
            latex.DrawLatex(0.68, 0.75, "Negative Rapidity");
        }
    }
    if (save_plots == true){
        perp_rapid->SaveAs("plots/full/kperp_srapid_6.30.png");
    }
    cout << "Kperp Rapid Fit" << endl;



//############################### Parameters ###############################
    TCanvas *params = new TCanvas("params", "params", 1000, 800);
    params->Divide(3, 2);
    const char* par_name[6] = {"#sigma, p2", "x_{0}, p3", "C_{0}, p7", "C_{1}, p8", "#alpha, p9", "#omega, p10"};
    vector<int> varying_params = {2, 3, 7, 8}; //omega removed after being fixed
    int spot = 0;
    for (int k = 0; k < varying_params.size(); k++) {
        int j = varying_params[k];
        params->cd(k+1);

        fit_params[j]->SetMarkerStyle(20);
        fit_params[j]->SetMarkerSize(1.2);
        fit_params[j]->GetXaxis()->SetTitle(par_name[spot]);
        spot++;
        fit_params[j]->LabelsOption("v","X");

        fit_params[j]->Draw("P");
    }
    params->cd(5);
    M2_hist->SetMarkerStyle(20);
    M2_hist->SetMarkerSize(1.2);
    M2_hist->SetLineWidth(2);
    M2_hist->GetXaxis()->SetTitle("EMG RMS");
    M2_hist->LabelsOption("v","X");
    M2_hist->Draw("P");

    if (save_plots == true){
        params->SaveAs("plots/full/params_6.30.png");
    }

//############################### ZnZn Disassociative Calculation ###############################

    // Integrated, the tail function is -C_1 * alpha/(1-omega) * (1+x/alpha)^(-omega+1) //
    //Integrating k: [150, 1000] MeV and multiply by N counts takes (1/N)dN/dk -> N_diss //
    // The integral, without C_1, is 425/4 = 106.30 //

    double c_1 = fitp_z0_8; double alpha = 250.0; double omega = 2.0; int N_total = zn0n_pevent_count + zn0n_nevent_count;
    double tail_integral = c_1 * alpha/(1.0-omega) * ( pow((1.0 + 1000.0/alpha), -omega + 1.0) - pow((1.0 + 150.0/alpha), -omega + 1.0) );
        
    cout << N_total * tail_integral << endl;

//############################### Params by neutron mean ###############################
    TCanvas *par_nmean = new TCanvas("par_test", "test parameters", 1500, 1000);
    par_nmean->Divide(3,2);

    par_nmean->cd(1);
    g_sigma->SetMarkerStyle(20); g_sigma->SetMarkerColor(kRed); 
    g_sigma->SetTitle(";Mean neutron number;#sigma");
    g_sigma->Draw("AP");
    gPad->Update();
    for (int i = 0; i < neutron_bin; i++) {
        double x, y;
        g_sigma->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.7*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
        if (i < 6){
            TMarker *m = new TMarker(x, y, 20);
            m->SetMarkerColor(kGreen+2);
            m->SetMarkerSize(g_sigma->GetMarkerSize());
            m->Draw();
        }
    }

    par_nmean->cd(2);
    g_x0->SetMarkerStyle(20); g_x0->SetMarkerColor(kRed);
    g_x0->SetTitle("; ;x_{0}");
    g_x0->Draw("AP");
    gPad->Update();
    for (int i = 0; i < neutron_bin; i++) {
        double x, y;
        g_x0->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.7*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
        if (i < 6){
            TMarker *m = new TMarker(x, y, 20);
            m->SetMarkerColor(kGreen+2);
            m->SetMarkerSize(g_x0->GetMarkerSize());
            m->Draw();
        }
    }

    par_nmean->cd(3);
    g_c0->SetMarkerStyle(20); g_c0->SetMarkerColor(kRed);
    g_c0->SetTitle("; ;C_{0}");
    g_c0->Draw("AP");
    gPad->Update();
    for (int i = 0; i < neutron_bin; i++) {
        double x, y;
        g_c0->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.7*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
        if (i < 6){
            TMarker *m = new TMarker(x, y, 20);
            m->SetMarkerColor(kGreen+2);
            m->SetMarkerSize(g_c0->GetMarkerSize());
            m->Draw();
        }
    }

    par_nmean->cd(4);
    g_c1->SetMarkerStyle(20); g_c1->SetMarkerColor(kRed);
    g_c1->SetTitle("; ;C_{1}");
    g_c1->Draw("AP");
    gPad->Update();
    for (int i = 0; i < neutron_bin; i++) {
        double x, y;
        g_c1->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.7*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
        if (i < 6){
            TMarker *m = new TMarker(x, y, 20);
            m->SetMarkerColor(kGreen+2);
            m->SetMarkerSize(g_c1->GetMarkerSize());
            m->Draw();
        }
    }

    par_nmean->cd(5);
    g_emgrms->SetMarkerStyle(20); g_emgrms->SetMarkerColor(kRed);
    g_emgrms->SetTitle("; ;EMG RMS");
    g_emgrms->Draw("AP");
    gPad->Update();
    for (int i = 0; i < neutron_bin; i++) {
        double x, y;
        g_emgrms->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.7*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
        if (i < 6){
            TMarker *m = new TMarker(x, y, 20);
            m->SetMarkerColor(kGreen+2);
            m->SetMarkerSize(g_emgrms->GetMarkerSize());
            m->Draw();
        }
    }

    par_nmean->Modified();
    par_nmean->Update();
    if (save_plots == true){
        par_nmean->SaveAs("plots/full/params_nmean_6.30.png");
    }

    TCanvas *pars_fit = new TCanvas("pars_fit", "c0 and c1 fit", 1200, 1000); //both of these fits give a "lower/upper bounds outside current parameter value"
    pars_fit->Divide(1,2);

    pars_fit->cd(1);
    TF1 *c0_0_func = new TF1("c0_0_func", "[0] * pow((1 + x/[1]), -[2]) + [3]");
    TF1 *c0_N_func = new TF1("c0_N_func", "[0] * pow((1 + x/[1]), -[2]) + [3]");
    fit_c0_0->SetMinimum(0.7);
    fit_c0_0->GetXaxis()->SetLimits(0, 17);
    fit_c0_0->GetYaxis()->SetTitle("C_{0}");
    fit_c0_0->SetMarkerColor(kGreen+2);
    c0_0_func->SetLineColor(kGreen+2);
    fit_c0_N->SetMarkerColor(kRed);
    c0_N_func->SetLineColor(kRed);

    c0_0_func->SetParLimits(1, 0, 1e5);
    c0_0_func->SetParLimits(2, 1e-5, 100);
    fit_c0_0->Fit(c0_0_func, "QI");
    fit_c0_0->Draw("AP");
    gPad->Update();
    for (int i = 1; i < 6; i++) {
        double x, y;
        fit_c0_0->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.2*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
    }
    c0_N_func->SetParLimits(1, 0, 1e5);
    c0_N_func->SetParLimits(2, 1e-5, 100);
    fit_c0_N->Fit(c0_N_func, "QI");
    fit_c0_N->Draw("P SAME");
    gPad->Update();
    for (int i = 5; i < neutron_bin; i++) {
        double x, y;
        fit_c0_N->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.2*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
    }


    pars_fit->cd(2);
    TF1 *c1_0_func = new TF1("c1_0_func", "[0] * pow((1 + x/[1]), -[2]) + [3]");
    TF1 *c1_N_func = new TF1("c1_N_func", "-[0] * pow((1 + x/[1]), -[2]) + [3]");
    fit_c1_0->SetMinimum(0.01);
    fit_c1_0->GetXaxis()->SetLimits(0.0, 17.0);
    fit_c1_0->GetYaxis()->SetTitle("C_{1}");
    fit_c1_0->SetMarkerColor(kGreen+2);
    c1_0_func->SetLineColor(kGreen+2);
    fit_c1_N->SetMarkerColor(kRed);
    c1_N_func->SetLineColor(kRed);

    c1_0_func->SetParLimits(0, 0, 10);
    c1_0_func->SetParLimits(1, 0, 50);
    fit_c1_0->Fit(c1_0_func, "QI");
    fit_c1_0->Draw("AP");
    gPad->Update();
    for (int i = 1; i < 6; i++) {
        double x, y;
        fit_c1_0->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.2*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
    }

    c1_N_func->SetParLimits(0, 0, 10);
    c1_N_func->SetParLimits(1, 0, 50);
    fit_c1_N->Fit(c1_N_func, "QI");
    fit_c1_N->Draw("P SAME");
    gPad->Update();
    for (int i = 6; i < neutron_bin; i++) {
        double x, y;
        fit_c1_N->GetPoint(i, x, y);
        double xmin = gPad->GetUxmin(); double xmax = gPad->GetUxmax(); double dx = 0.03*(xmax-xmin);
        TLatex *lab = new TLatex(x + 1.2*dx, y, bin_names[i]);
        lab->SetTextSize(0.025); lab->SetTextAlign(22); lab->Draw();
    }

    pars_fit->Modified();
    pars_fit->Update();

    if (save_plots == true){
        pars_fit->SaveAs("plots/full/pars_fit_6.30.png");
    }


}