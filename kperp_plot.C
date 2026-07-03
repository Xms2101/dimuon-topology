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

//############################### Initial 0n0n function ###############################
Double_t EMG_plus_power(Double_t *x, Double_t *par) //this is the fit we want for the 0n0n. EMG + [4]*pow(1+x/[5], -[6]) + [7]
{
    // --- EMG part ---
    Double_t emg = EMG(x, par);  // uses par[0..3]

    // --- Power-law part ---
    Double_t A  = par[4];
    Double_t x0 = par[5];   // scale
    Double_t n  = par[6];

    if (x0 <= 0) return 0;  // safety

    Double_t power = A * std::pow(1.0 + x[0]/x0, -n);

    // --- Constant offset ---
    Double_t c = par[7];

    return emg + power + c;
}

//############################### Template + Normal Fit ###############################
Double_t EMG_template(Double_t *x, Double_t *par)
{
    // ===== FIXED VALUES (these are copy and pasted from the first fit) =====
    const double tau_fixed   = 17.4291;   
    //const double x0_fixed    = 143.295;   
    const double n1_fixed    = 2.34883;    
    const double scale1_fixed = 143.295;

    // ===== EMG part =====
    Double_t emg_par[4];
    emg_par[0] = par[0];        // integral
    emg_par[1] = tau_fixed;     // FIXED tau
    emg_par[2] = par[1];        // sigma (free)
    emg_par[3] = par[2];      // x0 (free. not the power law one)

    Double_t emg = EMG(x, emg_par);

    // ===== Power law 1 (fixed shape, free amplitude) =====
    Double_t A1 = par[3];
    Double_t power1 = A1 * std::pow(1.0 + x[0]/scale1_fixed, -n1_fixed);

    // ===== Power law 2 (fully free) =====
    Double_t A2     = par[4];
    Double_t scale2 = par[5];
    Double_t n2     = par[6];

    if (scale2 <= 0) return 0;

    Double_t power2 = A2 * std::pow(1.0 + x[0]/scale2, -n2);

    // ===== pedestal =====
    Double_t c = par[7];

    Double_t val = emg + power1 + power2 + c;

    if (!std::isfinite(val)) return 0;
    //if (!std::isfinite(val) || val <= 0) return 1e-300;

    return val;
}

//############################### Functions for Plotting ###############################
Double_t EMG_component(Double_t *x, Double_t *par)
{
    const double tau_fixed = 17.4291;
    const double x0_emg    = 143.295;

    Double_t emg_par[4];
    emg_par[0] = par[0];        // integral
    emg_par[1] = tau_fixed;
    emg_par[2] = par[1];        // sigma
    emg_par[3] = x0_emg;

    return EMG(x, emg_par);
}
Double_t power1_component(Double_t *x, Double_t *par)
{
    const double scale1_fixed = 143.295;
    const double n1_fixed     = 2.34883;

    return par[2] * std::pow(1.0 + x[0]/scale1_fixed, -n1_fixed);
}
Double_t power2_component(Double_t *x, Double_t *par)
{
    if (par[4] <= 0) return 1e-300;

    return par[3] * std::pow(1.0 + x[0]/par[4], -par[5]);
}
Double_t pedestal_component(Double_t *x, Double_t *par)
{
    return par[6];
}



void kperp_plot(){    
    TFile *f = TFile::Open("dimuon_hists.root");
    static const int neutron_bin = 12;
    static const int neutron1    = 6;
    const char* bin_names[neutron_bin] = {"0n0n", "1n0n", "1n1n", "Xn0n", "Xn1n", "XnXn", "0nYn", "Yn1n", "YnYn", "Zn0n", "Zn1n", "ZnZn"}; //lined up so that they're easier to compare

    TH1D *kperp_hist[neutron_bin];
    TH1D *kperp_lin[neutron_bin];
    TH1D *onon_spare={nullptr};
    TH1D *onon_lin={nullptr};

    for(int i = 0; i < neutron_bin; i++){
        kperp_hist[i]    = (TH1D*)f->Get(Form("h_kperp_%d", i));
        kperp_lin[i]     = (TH1D*)f->Get(Form("h_k_lin_%d", i));
    }
    onon_spare = (TH1D*)f->Get("h_kperp_0")->Clone("onon_spare");
    onon_spare->SetDirectory(nullptr); // to make sure that we're not recycling the pointer


//############################### 0n0n Fit and Template ###############################
    TCanvas *k_onon_emg = new TCanvas("k_onon_emg", "0n0n Linear Fit", 1200, 1000);
    TF1 *f_emg = new TF1("f_emg", EMG_plus_power, 1, 1e4, 8);

    onon_spare->Scale(1.0 / onon_spare->Integral());
    onon_spare->Scale(1.0, "width");

    f_emg->SetParLimits(1, 1e-6, 100);
    f_emg->SetParLimits(2, 1e-6, 50);

    f_emg->SetParLimits(4, 0, 1);
    f_emg->SetParLimits(5, 1e-4, 1e8);
    f_emg->SetParLimits(6, 1, 1e4);
    f_emg->SetParLimits(7, 0, 1e-6);
    
    //integral, tau, sigma, x0 for EMG
    f_emg->SetParameters(onon_spare->Integral(), 100.0, 5, onon_spare->GetMean(), 0.01, 1, 1, 1, 1, onon_spare->GetBinContent(onon_spare->GetNbinsX()));

    onon_spare->Draw("HIST");
    onon_spare->Fit(f_emg, "RI");
    f_emg->Draw("same");

    gPad->SetLogy();
    gPad->SetLogx();

    const double tau_1 = f_emg->GetParameter(1);
    const double pow_scale = f_emg->GetParameter(5);
    const double pow_exp = f_emg->GetParameter(6);

    double emg0 = f_emg->GetParameter(0);
    double emg2 = f_emg->GetParameter(2);
    double emg3 = f_emg->GetParameter(3);
    double emg4 = f_emg->GetParameter(4);
    double emg7 = f_emg->GetParameter(7);


/*
//############################### kperp 1 ###############################
    //this is a modified Guassian + a power law + pedestal 
    //TGaxis::SetMaxDigits(3);
    TCanvas *perps1 = new TCanvas("perps1", "kperp1", 1200, 1000);
    perps1->Divide(3, 2);
    bool draw_bits1 = false;

    for (int i = 0; i < neutron1; i++){
        TF1 *kperp_fit1 = new TF1(Form("kperp_fit1_%d",i), EMG_template, 1, 1e4, 7);

        kperp_hist[i]->Scale(1.0 / kperp_hist[i]->Integral()); // 1/N dN/dk_perp
        kperp_hist[i]->Scale(1.0, "width");

        perps1->cd(i+1);
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetRightMargin(0.12);

        kperp_fit1->SetParameters(1.0, kperp_hist[i]->GetRMS(), 0.01, 0.01, 100.0, 3.0, 1e-8);

        kperp_fit1->SetParLimits(1, 1e-4, 100);
        kperp_fit1->SetParLimits(2, 0, 10);
        kperp_fit1->SetParLimits(3, 0, 10);
        kperp_fit1->SetParLimits(4, 1e-3, 1e10);
        kperp_fit1->SetParLimits(5, 0.5, 1e8);
        kperp_fit1->SetParLimits(6, 0, 1e-4);
    
        kperp_hist[i]->Draw("HIST");
        kperp_hist[i]->Fit(kperp_fit1, "QRI");
        kperp_fit1->Draw("same");


        Double_t *p = kperp_fit1->GetParameters();
        TF1 *f_emg   = new TF1(Form("emg_%d", i), EMG_component, 1, 1e4, 7);
        TF1 *f_pow1  = new TF1(Form("pow1_%d", i), power1_component, 1, 1e4, 7);
        TF1 *f_pow2  = new TF1(Form("pow2_%d", i), power2_component, 1, 1e4, 7);
        TF1 *f_ped   = new TF1(Form("ped_%d", i), pedestal_component, 1, 1e4, 7);

        f_emg->SetParameters(p);
        f_pow1->SetParameters(p);
        f_pow2->SetParameters(p);
        f_ped->SetParameters(p);

        f_emg->SetLineColor(kBlue);
        f_pow1->SetLineColor(kGreen+2);
        f_pow2->SetLineColor(kMagenta);
        f_ped->SetLineColor(kGray+2);

        if (draw_bits1 == true){
            f_emg->Draw("same");
            f_pow1->Draw("same");
            f_pow2->Draw("same");
            f_ped->Draw("same");
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
            latex.DrawLatex(0.77, 0.73, "X #geq 1");
            latex.DrawLatex(0.77, 0.66, "1 < Y #leq 5");
            latex.DrawLatex(0.77, 0.59, "Z > 5");
            latex.SetTextSize(0.05);
            latex.DrawLatex(0.77, 0.88, "k_{#perp}");
        }

    }
    //perps1->SaveAs("plots/full/kperp11_3.11.png");



  //############################### kperp 2 ###############################
    TCanvas *perps2 = new TCanvas("perps2", "kperp2", 1200, 1000);
    perps2->Divide(3, 2);
    //TGaxis::SetMaxDigits(3);
    bool draw_bits2 = false; 

    for (int i = 6; i < neutron_bin; i++){
        TF1 *kperp_fit2 = new TF1(Form("kperp_fit2_%d",i), EMG_template, 1, 1e4, 7);

        kperp_hist[i]->Scale(1.0 / kperp_hist[i]->Integral()); // 1/N dN/dk_perp
        kperp_hist[i]->Scale(1.0, "width");

        perps2->cd(i-5);
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetRightMargin(0.12);

        kperp_fit2->SetParameters(1.0, kperp_hist[i]->GetRMS(), 0.01, 0.01, 100.0, 3.0, 1e-8);

        kperp_fit2->SetParLimits(0, 0, 100);
        kperp_fit2->SetParLimits(1, 1e-4, 1000);
        kperp_fit2->SetParLimits(2, 1e-6, 10);
        kperp_fit2->SetParLimits(3, 0, 10);
        kperp_fit2->SetParLimits(4, 1e-3, 1e10);
        kperp_fit2->SetParLimits(5, 0.5, 1e8);
        kperp_fit2->SetParLimits(6, 0, 1e-4);

        kperp_hist[i]->Draw("HIST");
        kperp_hist[i]->Fit(kperp_fit2, "RI");
        kperp_fit2->Draw("same");

        Double_t *p = kperp_fit2->GetParameters();
        TF1 *f_emg   = new TF1(Form("emg_%d", i), EMG_component, 1, 1e4, 7);
        TF1 *f_pow1  = new TF1(Form("pow1_%d", i), power1_component, 1, 1e4, 7);
        TF1 *f_pow2  = new TF1(Form("pow2_%d", i), power2_component, 1, 1e4, 7);
        TF1 *f_ped   = new TF1(Form("ped_%d", i), pedestal_component, 1, 1e4, 7);

        f_emg->SetParameters(p);
        f_pow1->SetParameters(p);
        f_pow2->SetParameters(p);
        f_ped->SetParameters(p);

        f_emg->SetLineColor(kBlue);
        f_pow1->SetLineColor(kGreen+2);
        f_pow2->SetLineColor(kMagenta);
        f_ped->SetLineColor(kGray+2);

        if (draw_bits2 == true){
            //f_emg->Draw("same");
            //f_pow1->Draw("same");
            //f_pow2->Draw("same");
            //f_ped->Draw("same");
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
    //perps2->SaveAs("plots/full/kperp12_3.11.png");
*/
}
