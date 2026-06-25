//rapidity where sign is determined on which nucleus has the 0n (set one of those neg, other pos. Add).
    //Direction that 0n came from, set pos, other neg (only for the *n0n)


#include "AtlasStyle.C"
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


double TrueRMS(TH1D* h){
    double sumN  = 0.0;
    double sumN2 = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double n = h->GetBinLowEdge(b);
        double N = h->GetBinContent(b);
        sumN  += N;
        sumN2 += N*n*n;
    }
    if (sumN == 0) return 0.0;
    return std::sqrt(sumN2/sumN);
}
double TrueRMSError(TH1D* h){
    double sumN  = 0.0;
    double sumN2 = 0.0;
    double sumN4 = 0.0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double n = h->GetBinLowEdge(b);   // neutron count
        double N = h->GetBinContent(b);
        sumN  += N;
        sumN2 += N*n*n;
        sumN4 += N*n*n*n*n;
    }
    if (sumN <= 1) return 0.0;
    double m2 = sumN2/sumN;
    double m4 = sumN4/sumN;
    double rms = std::sqrt(m2);
    if (rms == 0) return 0.0;
    return std::sqrt((m4 - m2*m2)/sumN)/(2.0*rms);
}
auto ShiftHistogram = [](TH1D* h, const char* newname, double shift){ //shifts bins by 0.5
    int nbins = h->GetNbinsX();
    std::vector<double> edges(nbins + 1);
    for (int i = 1; i <= nbins + 1; i++) {
        edges[i-1] = h->GetXaxis()->GetBinLowEdge(i) + shift;
    }
    TH1D* h_shift = new TH1D(newname, h->GetTitle(), nbins, edges.data());
    for (int i = 1; i <= nbins; i++) {
        h_shift->SetBinContent(i, h->GetBinContent(i));
        h_shift->SetBinError(i, h->GetBinError(i));
    }
    return h_shift;
};


void rapidity_analysis(){
    bool save_plots = false;

    TFile *f = TFile::Open("dimuon_hists.root");
    static const int neutron_bin = 13;
    const char* bin_names[neutron_bin] = {"0n0n", "1n0n", "2n0n", "Xn0n", "Yn0n", "Zn0n", "1n1n", "Xn1n", "Yn1n", "Zn1n", "XnXn", "YnYn", "ZnZn"}; //Organized as 0n, 1n, Nn

    TH1D *rapid_sgn_hist[neutron_bin];
    TH1D *rapid_a150_hist[neutron_bin];
    TH1D *rapid_b150_hist[neutron_bin];
    TH1D* rapid_150_div[neutron_bin];
    TH1D* neutron_stat[neutron_bin];
    TH1D* n_raw_a;
    TH1D* n_raw_c;
    TH1D* n_a;
    TH1D* n_c;

    n_raw_a = (TH1D*)f->Get("h_n_raw_a");
    n_raw_c = (TH1D*)f->Get("h_n_raw_c");
    n_a     = (TH1D*)f->Get("h_NeutronsA");
    n_c     = (TH1D*)f->Get("h_NeutronsC");
    for(int i = 0; i < neutron_bin; i++){
        rapid_sgn_hist[i]   = (TH1D*)f->Get(Form("h_rapid_sgn_%d", i));
        rapid_a150_hist[i]  = (TH1D*)f->Get(Form("h_rapid_a150_%d", i));
        rapid_b150_hist[i]  = (TH1D*)f->Get(Form("h_rapid_b150_%d", i));
        neutron_stat[i]     = (TH1D*)f->Get(Form("neutron_stat_%d", i));

        rapid_a150_hist[i]->Sumw2();
        rapid_b150_hist[i]->Sumw2();

        rapid_150_div[i] = (TH1D*) rapid_a150_hist[i]->Clone(Form("rapid_150_div_%d", i));
        rapid_150_div[i]->SetDirectory(nullptr);
        rapid_150_div[i]->Divide(rapid_sgn_hist[i]); // not divided by entire sample for proper fraction
    }

    TGraphErrors* g_mean = new TGraphErrors(neutron_bin);
    for (int i = 0; i < neutron_bin; i++) {
        double mean = neutron_stat[i]->GetMean();
        double rms  = neutron_stat[i]->GetRMS();
        double n    = neutron_stat[i]->GetEntries();

        g_mean->SetPoint(i, i, mean);
        g_mean->SetPointError(i, 0, rms / sqrt(n));
    }
    TGraphErrors* g_rms = new TGraphErrors(neutron_bin);
    for (int i = 0; i < neutron_bin; i++) {
        double rms = TrueRMS(neutron_stat[i]);
        double rms_err = TrueRMSError(neutron_stat[i]);

        g_rms->SetPoint(i, i, rms);
        g_rms->SetPointError(i, 0, rms_err);
        //cout << bin_names[i] << " RMS = " << rms << endl;
    }

//############################### Signed Rapidity ###############################
    TCanvas *rapids = new TCanvas("rapids", "Signed Rapidity", 1200, 1000); // histogram only filled from 1-6
    rapids->Divide(3,2);
    for (int i = 1; i < 6; i++){
        rapids->cd(i);

        rapid_sgn_hist[i]->SetMinimum(rapid_sgn_hist[i]->GetMinimum()/1.05);
        rapid_sgn_hist[i]->SetMaximum(1.05*rapid_sgn_hist[i]->GetMaximum());
        rapid_sgn_hist[i]->Draw("HIST");

        rapid_sgn_hist[i]->GetYaxis()->SetTitle("counts");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.82, 0.8, Form("%s", bin_names[i]));
        rapid_sgn_hist[i]->GetXaxis()->SetTitle("signed rapidity");

        TLine *v_line= new TLine(0, rapid_sgn_hist[i]->GetMinimum(), 0, rapid_sgn_hist[i]->GetMaximum());
        v_line->SetLineColor(kRed);
        v_line->SetLineWidth(3);
        v_line->Draw("same");
    }
    if (save_plots == true){
        rapids->SaveAs("plots/full/rapidity_sgn_6.25.png");
    }

    TCanvas *rapids_high = new TCanvas("rapids_high", "rapidity, kperp>150 MeV", 1200, 1000);
    rapids_high->Divide(3,2);
    for (int i = 1; i < 6; i++){
        rapids_high->cd(i);

        rapid_a150_hist[i]->SetMinimum(rapid_a150_hist[i]->GetMinimum()/1.05);
        rapid_a150_hist[i]->SetMaximum(1.05*rapid_a150_hist[i]->GetMaximum());
        rapid_a150_hist[i]->Draw("HIST");

        TLine *v_line= new TLine(0, rapid_a150_hist[i]->GetMinimum(), 0, rapid_a150_hist[i]->GetMaximum());
        v_line->SetLineColor(kRed);
        v_line->SetLineWidth(3);
        v_line->Draw("same");

        rapid_a150_hist[0]->GetYaxis()->SetTitle("counts");
        rapid_a150_hist[i]->GetXaxis()->SetTitle("signed rapidity");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.87, 0.8, Form("%s", bin_names[i]));
    }
    if (save_plots == true){
        rapids_high->SaveAs("plots/full/rapidity_150a_sgn_6.25.png");
    }

    TCanvas *rapids_low = new TCanvas("rapids_low", "rapidity, kperp<150 MeV", 1200, 1000); // histogram only filled from 8-12
    rapids_low->Divide(3,2);
    for (int i = 1; i < 6; i++){
        rapids_low->cd(i);
        rapid_b150_hist[i]->SetMinimum(rapid_b150_hist[i]->GetMinimum()/1.05);
        rapid_b150_hist[i]->SetMaximum(1.05*rapid_b150_hist[i]->GetMaximum());
        rapid_b150_hist[i]->Draw("HIST");

        TLine *v_line= new TLine(0, rapid_b150_hist[i]->GetMinimum(), 0, rapid_b150_hist[i]->GetMaximum());
        v_line->SetLineColor(kRed);
        v_line->SetLineWidth(3);
        v_line->Draw("same");

        rapid_b150_hist[i]->GetYaxis()->SetTitle("counts");
        rapid_b150_hist[i]->GetXaxis()->SetTitle("signed rapidity");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.87, 0.8, Form("%s", bin_names[i]));
    }
    if (save_plots == true){
        rapids_low->SaveAs("plots/full/rapidity_150b_sgn_6.25.png");
    }

    TCanvas *divided = new TCanvas("divided", "divided", 1200, 1000);
    divided->Divide(3,2);
    for (int i = 1; i < 6; i++){
        divided->cd(i);

        rapid_150_div[i]->SetMinimum(rapid_150_div[i]->GetMinimum()/1.05);
        rapid_150_div[i]->SetMaximum(1.05*rapid_150_div[i]->GetMaximum());
        rapid_150_div[i]->Draw("HIST");

        TLine *v_line= new TLine(0, rapid_150_div[i]->GetMinimum(), 0, rapid_150_div[i]->GetMaximum());
        v_line->SetLineColor(kRed);
        v_line->SetLineWidth(3);
        v_line->Draw("same");

        rapid_150_div[i]->GetYaxis()->SetTitle("hist above 150 MeV / hist below 150 MeV");
        rapid_150_div[i]->GetXaxis()->SetTitle("signed rapidity");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.68, 0.8, Form("%s", bin_names[i]));
    }
    if (save_plots == true){
        divided->SaveAs("plots/full/rapidity_150_division_6.25.png");
    }

//############################### Neutron Statistics ###############################
    TCanvas *neutrons = new TCanvas("neutrons", "neutrons", 1200, 1000);
    neutrons->Divide(5,3);
    for (int i = 0; i < neutron_bin; i++){
        neutrons->cd(i+1);
        neutron_stat[i]->Draw("HIST");

        neutron_stat[i]->GetXaxis()->SetTitle("Total Neutrons");
        neutron_stat[i]->GetYaxis()->SetTitle("Count");
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.65, 0.8, Form("%s", bin_names[i]));
    }
    if (save_plots == true){
        neutrons->SaveAs("plots/full/neutron_count_6.25.png");
    }

    TCanvas *n_mean = new TCanvas("n_mean", "n_mean", 1200, 1000);
    TH1F *hframe_mean = new TH1F("hframe_mean", "", neutron_bin, -0.5, neutron_bin-0.5);
    for (int i = 0; i < neutron_bin; i++) {
        hframe_mean->GetXaxis()->SetBinLabel(i+1, bin_names[i]);
    }
    hframe_mean->GetYaxis()->SetTitle("Mean Neutrons");
    hframe_mean->SetMinimum(0.0);
    hframe_mean->SetMaximum(1.2*g_mean->GetY()[TMath::LocMax(neutron_bin,g_mean->GetY())]);

    hframe_mean->Draw();
    g_mean->Draw("P SAME");
    if (save_plots == true){
        n_mean->SaveAs("plots/full/neutrons_mean_6.25.png");
    }


    TCanvas *n_rms = new TCanvas("n_rms", "n_rms", 1200, 1000);
    TH1F *hframe_rms = new TH1F("hframe_rms", "", neutron_bin, -0.5, neutron_bin-0.5);
    for (int i = 0; i < neutron_bin; i++) {
        hframe_rms->GetXaxis()->SetBinLabel(i+1, bin_names[i]);
    }
    hframe_rms->GetYaxis()->SetTitle("Neutron Count RMS");
    hframe_rms->SetMaximum(1.2*g_rms->GetY()[TMath::LocMax(neutron_bin,g_rms->GetY())]);

    hframe_rms->Draw();
    g_rms->Draw("P SAME");
    if (save_plots == true){
        n_rms->SaveAs("plots/full/neutron_rms_6.25.png");
    }

//############################### Neutron Statistics ###############################
    TCanvas *multiplied = new TCanvas("multiplied", "Multiplied", 1200, 1000);
    TH1D *rapid_prod = (TH1D*)rapid_150_div[5]->Clone("rapid_prod");
    rapid_prod->Multiply(rapid_sgn_hist[5]);

    rapid_prod->GetXaxis()->SetTitle("Signed Rapidity");
    rapid_prod->GetYaxis()->SetTitle("Dissasociative^{*} count");

    rapid_prod->Draw("HIST");
    if (save_plots == true){
        multiplied->SaveAs("plots/full/diss_rapid_6.25.png");
    }


    TCanvas *n_raw = new TCanvas("n_raw", "Raw Neutron Energies", 1200, 1000);
    TH1D* n_a_shift = ShiftHistogram(n_a, "n_a_shift", -0.5);
    TH1D* n_c_shift = ShiftHistogram(n_c, "n_c_shift", -0.5);

    n_raw_a->SetLineColor(kRed);
    n_raw_c->SetLineColor(kBlue);
    gPad->SetLogy();

    n_a_shift->SetLineColor(kRed);
    n_a_shift->SetLineWidth(4);
    n_a_shift->SetLineStyle(9);
    n_c_shift->SetLineColor(kBlue);
    n_c_shift->SetLineWidth(4);
    n_c_shift->SetLineStyle(9);

    n_raw_a->GetXaxis()->SetTitle("Neutrons");
    n_raw_a->GetYaxis()->SetTitle("Count");

    n_raw->SetGridx();      // enable vertical grid lines
    //n_raw->SetGridy(0); 

    n_raw_a->Draw("HIST");
    n_raw_c->Draw("HIST SAME");
    n_a_shift->Draw("HIST SAME");
    n_c_shift->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.65, 0.72, 0.90, 0.90);
    leg->SetTextSize(0.04);
    leg->AddEntry(n_raw_a, "A Raw", "l");
    leg->AddEntry(n_raw_c, "C Raw", "l");
    leg->AddEntry(n_a, "A Rounded", "l");
    leg->AddEntry(n_c, "C Rounded", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    n_raw->Update();

    if (save_plots == true){
        n_raw->SaveAs("plots/full/neutron_energies_6.25.png");
    }

int start_b = rapid_prod->FindBin(-2.4);
int min_b = rapid_prod->FindBin(0);
int max_b = rapid_prod->FindBin(2.4);
cout << rapid_prod->Integral(min_b, max_b) << endl;

}