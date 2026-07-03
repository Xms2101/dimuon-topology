#include <math.h>
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
#include <cmath>
#include <vector>


void part_acop_plot(){
    TFile *f = TFile::Open("dimuon_hists.root");
    static const int neutron_bin = 6;
    const char* bin_names[neutron_bin] = {"0n0n", "0n1n", "0nXn", "1n1n", "1nXn", "XnXn"}; //as from the CMS paper

    TF1 * core_fit = new TF1("core_fit","[0]*exp(-x/[1]+[2]*pow(x,0.25))", 0.5*10e-5, 1); //lower bound taken from current def of histogram
    TF1 * tail_fit = new TF1("tail_fit", "[0]*pow((1 + x*[1]/[2]), -[2])", 0.5*10e-5, 1);


    TH2D *a_n_hist[neutron_bin];

    TH1D *acop_hist[neutron_bin];
    TH1D *asym_hist[neutron_bin];
    TH1D *kperp_hist[neutron_bin];
    TH1D *dphi_hist[neutron_bin];
    TH1D *ptbar_hist[neutron_bin];
    TH1D *mass_hist[neutron_bin];
    TH1D *rapidity_hist[neutron_bin];
    TH1D *Egamma_a_hist[neutron_bin];
    TH1D *Egamma_b_hist[neutron_bin];
    TH1D *pt1_hist[neutron_bin];
    TH1D *pt2_hist[neutron_bin];
    TH1D *eta1_hist[neutron_bin];
    TH1D *eta2_hist[neutron_bin];
    TH1D *phi1_hist[neutron_bin];
    TH1D *phi2_hist[neutron_bin];
    TH1D *d01_hist[neutron_bin];
    TH1D *d02_hist[neutron_bin];
    TH1D *z01_hist[neutron_bin];
    TH1D *z02_hist[neutron_bin];
    TH1D *Zvtx_hist[neutron_bin];
    TH1D *NeutronsA_hist[neutron_bin];
    TH1D *NeutronsB_hist[neutron_bin];
    TH1D *qual1_hist[neutron_bin];
    TH1D *qual2_hist[neutron_bin];

    TH1D *h_acop_n1 = (TH1D*)f->Get("h_acop_n1");
    TH1D *h_acop_n1_c = (TH1D*)f->Get("h_acop_n1");
    TH1D *h_acop_n2 = (TH1D*)f->Get("h_acop_n2");
    TH1D *h_acop_n2_c = (TH1D*)f->Get("h_acop_n2");
    TH1D *h_acop_n3 = (TH1D*)f->Get("h_acop_n3");
    TH1D *h_acop_n3_c = (TH1D*)f->Get("h_acop_n3");
    TH1D *h_acop_n4 = (TH1D*)f->Get("h_acop_n4");
    TH1D *h_acop_n4_c = (TH1D*)f->Get("h_acop_n4");
    TH1D *h_acop_n5 = (TH1D*)f->Get("h_acop_n5");
    TH1D *h_acop_n5_c = (TH1D*)f->Get("h_acop_n5");
    TH1D *h_acop_n6 = (TH1D*)f->Get("h_acop_n6");
    TH1D *h_acop_n6_c = (TH1D*)f->Get("h_acop_n6");



    for(int i = 0; i < neutron_bin; i++){
        a_n_hist[i]=      (TH2D*)f->Get(Form("h_a_n_%d", i));

        acop_hist[i]=      (TH1D*)f->Get(Form("h_acop_%d", i)); 
        asym_hist[i]=      (TH1D*)f->Get(Form("h_asym_%d", i)); 
        kperp_hist[i]=     (TH1D*)f->Get(Form("h_kperp_%d", i)); 
        dphi_hist[i]=      (TH1D*)f->Get(Form("h_dphi_%d", i)); 
        ptbar_hist[i]=     (TH1D*)f->Get(Form("h_ptbar_%d", i)); 
        mass_hist[i]=      (TH1D*)f->Get(Form("h_mass_%d", i)); 
        rapidity_hist[i]=  (TH1D*)f->Get(Form("h_rapidity_%d", i)); 
        Egamma_a_hist[i]=  (TH1D*)f->Get(Form("h_Egamma_a_%d", i)); 
        Egamma_b_hist[i]=  (TH1D*)f->Get(Form("h_Egamma_b_%d", i)); 
        pt1_hist[i]=       (TH1D*)f->Get(Form("h_pt1_%d", i)); 
        pt2_hist[i]=       (TH1D*)f->Get(Form("h_pt2_%d", i)); 
        eta1_hist[i]=      (TH1D*)f->Get(Form("h_eta1_%d", i)); 
        eta2_hist[i]=      (TH1D*)f->Get(Form("h_eta2_%d", i)); 
        phi1_hist[i]=      (TH1D*)f->Get(Form("h_phi1_%d", i)); 
        phi2_hist[i]=      (TH1D*)f->Get(Form("h_phi2_%d", i)); 
        d01_hist[i]=       (TH1D*)f->Get(Form("h_d01_%d", i)); 
        d02_hist[i]=       (TH1D*)f->Get(Form("h_d02_%d", i)); 
        z01_hist[i]=       (TH1D*)f->Get(Form("h_z01_%d", i)); 
        z02_hist[i]=       (TH1D*)f->Get(Form("h_z02_%d", i)); 
        Zvtx_hist[i]=      (TH1D*)f->Get(Form("h_Zvtx_%d", i));
        qual1_hist[i]=     (TH1D*)f->Get(Form("h_qual1_%d", i));
        qual2_hist[i]=     (TH1D*)f->Get(Form("h_qual2_%d", i));
        //NeutronsA_hist[i]= (TH1D*)f->Get(Form("h_NeutronsA_%d", i));
        //NeutronsB_hist[i]= (TH1D*)f->Get(Form("h_NeutronsB_%d", i));
    }


    //Make the differential plot *****************
    TCanvas *acop = new TCanvas("acops", "acops", 1200, 1000);
    acop->Divide(3,2);
    for (int i = 0; i < neutron_bin; i++){
        acop->cd(i+1);
        TString name = Form("h_acop_%d", i);
        TH1D* h_acop = a_n_hist[i]->ProjectionX(name);

        // Optional: convert to differential (divide by bin width)
        for (int b = 1; b <= h_acop->GetNbinsX(); b++) {
            double c = h_acop->GetBinContent(b);
            double e = h_acop->GetBinError(b);
            double w = h_acop->GetBinWidth(b);
            h_acop->SetBinContent(b, c / w);
            h_acop->SetBinError(b, e / w);
        }

        // Optional: normalize to 1 if you want (1/N)dN/dα
        //h_acop->Scale(1.0 / h_acop->Integral("width"));

        // Label axes
        h_acop->GetXaxis()->SetTitle("Acoplanarity #alpha");
        h_acop->GetYaxis()->SetTitle("(1/N) dN/d#alpha");

        // Store or draw
        h_acop->SetTitle(Form("Neutron bin %d", i));
        h_acop->Draw("ELP");
    }

/*
    TCanvas *cops = new TCanvas("cops", "acop", 1200, 1000);
    cops->Divide(3,2);

    for (int b = 1; b <= h_acop_n1->GetNbinsX(); b++) {
        double c = h_acop_n1->GetBinContent(b);
        double e = h_acop_n1->GetBinError(b);
        double w = h_acop_n1->GetBinWidth(b);
        h_acop_n1->SetBinContent(b, c / w);
        h_acop_n1_c->SetBinContent(b, c / w);
        h_acop_n1->SetBinError(b, e / w);
        h_acop_n1_c->SetBinError(b, e / w);
    }

    for (int b = 1; b <= h_acop_n2->GetNbinsX(); b++) {
        double c = h_acop_n2->GetBinContent(b);
        double e = h_acop_n2->GetBinError(b);
        double w = h_acop_n2->GetBinWidth(b);
        h_acop_n2->SetBinContent(b, c / w);
        h_acop_n2_c->SetBinContent(b, c / w);
        h_acop_n2->SetBinError(b, e / w);
        h_acop_n2_c->SetBinError(b, e / w);
    }

    for (int b = 1; b <= h_acop_n3->GetNbinsX(); b++) {
        double c = h_acop_n3->GetBinContent(b);
        double e = h_acop_n3->GetBinError(b);
        double w = h_acop_n3->GetBinWidth(b);
        h_acop_n3->SetBinContent(b, c / w);
        h_acop_n3_c->SetBinContent(b, c / w);
        h_acop_n3->SetBinError(b, e / w);
        h_acop_n3_c->SetBinError(b, e / w);
    }

    for (int b = 1; b <= h_acop_n4->GetNbinsX(); b++) {
        double c = h_acop_n4->GetBinContent(b);
        double e = h_acop_n4->GetBinError(b);
        double w = h_acop_n4->GetBinWidth(b);
        h_acop_n4->SetBinContent(b, c / w);
        h_acop_n4_c->SetBinContent(b, c / w);
        h_acop_n4->SetBinError(b, e / w);
        h_acop_n4_c->SetBinError(b, e / w);
    }

    for (int b = 1; b <= h_acop_n5->GetNbinsX(); b++) {
        double c = h_acop_n5->GetBinContent(b);
        double e = h_acop_n5->GetBinError(b);
        double w = h_acop_n5->GetBinWidth(b);
        h_acop_n5->SetBinContent(b, c / w);
        h_acop_n5_c->SetBinContent(b, c / w);
        h_acop_n5->SetBinError(b, e / w);
        h_acop_n5_c->SetBinError(b, e / w);
    }

    for (int b = 1; b <= h_acop_n6->GetNbinsX(); b++) {
        double c = h_acop_n6->GetBinContent(b);
        double e = h_acop_n6->GetBinError(b);
        double w = h_acop_n6->GetBinWidth(b);
        h_acop_n6->SetBinContent(b, c / w);
        h_acop_n6_c->SetBinContent(b, c / w);
        h_acop_n6->SetBinError(b, e / w);
        h_acop_n6_c->SetBinError(b, e / w);
    }

    h_acop_n1->Scale(1.0 / h_acop_n1->Integral("width"));
    h_acop_n1_c->Scale(1.0 / h_acop_n1_c->Integral("width"));

    h_acop_n2->Scale(1.0 / h_acop_n2->Integral("width"));
    h_acop_n2_c->Scale(1.0 / h_acop_n2_c->Integral("width"));

    h_acop_n3->Scale(1.0 / h_acop_n3->Integral("width"));
    h_acop_n3_c->Scale(1.0 / h_acop_n3_c->Integral("width"));

    h_acop_n4->Scale(1.0 / h_acop_n4->Integral("width"));
    h_acop_n4_c->Scale(1.0 / h_acop_n4_c->Integral("width"));

    h_acop_n5->Scale(1.0 / h_acop_n5->Integral("width"));
    h_acop_n5_c->Scale(1.0 / h_acop_n5_c->Integral("width"));

    h_acop_n6->Scale(1.0 / h_acop_n6->Integral("width"));
    h_acop_n6_c->Scale(1.0 / h_acop_n6_c->Integral("width"));
*/


////////////************************ FIRST 0N0N ************************////////////
/* 
    cops->cd(1);

    gStyle->SetOptTitle(0);

    h_acop_n1->SetLineWidth(2);
    h_acop_n1_c->SetLineWidth(2);
    
    h_acop_n1->SetMinimum(10e-4);
    h_acop_n1->SetMaximum(800.);

    h_acop_n1->SetMarkerStyle(20);
    h_acop_n1_c->SetMarkerStyle(20);
    h_acop_n1->SetMarkerSize(0.6);
    h_acop_n1_c->SetMarkerSize(0.6);
    h_acop_n1->SetMarkerColor(kBlack);
    h_acop_n1_c->SetMarkerColor(kBlack);
    h_acop_n1->SetLineColor(kBlack);
    h_acop_n1_c->SetLineColor(kBlack);

    tail_fit->SetParameters(1, 1, 1, 2);
    tail_fit->SetParLimits(0, 0, 1e5);
    tail_fit->SetParLimits(1, 0, 1e3);
    tail_fit->SetParLimits(2, 1e-6, 1e3);
    tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1, 2);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
    core_fit->SetParLimits(3, 0, 10);
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n1->Fit(tail_fit, "QRI");
    h_acop_n1_c->Fit(core_fit, "QRI+");

    float c_1_1 = core_fit->GetParameter(0);
    float c_2_1 = core_fit->GetParameter(1);
    float c_3_1 = core_fit->GetParameter(2);

    float t_1_1 = tail_fit->GetParameter(0);
    float t_2_1 = tail_fit->GetParameter(1);
    float t_3_1 = tail_fit->GetParameter(2);

    TF1 *core_funct1 = new TF1("core_funct", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_1, c_2_1, c_3_1), 0.5*10e-5, 1);
    TF1 *tail_funct1 = new TF1("tail_funct", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_1, t_2_1, t_3_1, t_3_1), 0.5*10e-5, 1);
    TF1 *sum_funct1 = new TF1("sum_funct", Form("%f * exp(-x / %f + %f * pow(x, 0.25)) + %f * pow((1 + x* %f/%f), -%f)", c_1_1, c_2_1, c_3_1, t_1_1, t_2_1, t_3_1, t_3_1), 0.5*10e-5, 1);

    core_funct1->SetLineColor(kRed);
    core_funct1->SetLineWidth(3);
    core_funct1->SetLineStyle(8);
    tail_funct1->SetLineColor(kBlue);
    tail_funct1->SetLineWidth(3);
    tail_funct1->SetLineStyle(2);

    sum_funct1->SetLineColor(kMagenta);


    h_acop_n1->Draw("ELP");
    h_acop_n1_c->Draw("SAME");
    core_funct1->Draw("SAME");
    tail_funct1->Draw("SAME");
    sum_funct1->Draw("SAME");

    gPad->SetLogx();
    gPad->SetLogy();

    TLatex latex1;
    latex1.SetNDC(); // use normalized coordinates (0–1)
    latex1.SetTextSize(0.05);
    latex1.SetTextAlign(22); // center horizontally
    latex1.DrawLatex(0.8, 0.8, Form("%s", bin_names[0]));
*/


}