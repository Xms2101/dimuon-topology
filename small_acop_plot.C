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


void small_acop_plot(){
    TFile *f = TFile::Open("dimuon_hists.root");
    static const int neutron_bin = 9;
    const char* bin_names[neutron_bin] = {"0n0n", "0n1n", "0nYn", "1n1n", "1nYn", "YnYn", "0nXn", "1nXn", "XnXn"}; //as from the CMS paper

    TF1 * core_fit = new TF1("core_fit","[0]*exp(-x/[1]+[2]*pow(x,0.25))", 0, 0.008); //lower bound taken from current def of histogram
    TF1 * tail_fit = new TF1("tail_fit", "[0]*pow((1 + x*[1]/[2]), -[2])", 0, 0.008);
    TF1 * sum_fit  = new TF1("sum_fit", "[0]*exp(-x/[1]+[2]*pow(x,0.25))+[3]*pow((1 + x*[4]/[5]), -[6])", 0, 0.008);



    TH2D *a_n_hist[neutron_bin];
    TH1D *kperp_hist[neutron_bin];
    TH1D *mass_hist[neutron_bin];
    TH1D *rapidity_hist[neutron_bin];

    TH1D *h_acop_n1 = (TH1D*)f->Get("h_acop_s_n1");
    TH1D *h_acop_n1_c = (TH1D*)h_acop_n1->Clone("h_acop_n1_c");
    TH1D *h_acop_n2 = (TH1D*)f->Get("h_acop_s_n2");
    TH1D *h_acop_n2_c = (TH1D*)h_acop_n2->Clone("h_acop_n2_c");
    TH1D *h_acop_n3 = (TH1D*)f->Get("h_acop_s_n3");
    TH1D *h_acop_n3_c = (TH1D*)h_acop_n3->Clone("h_acop_n3_c");
    TH1D *h_acop_n4 = (TH1D*)f->Get("h_acop_s_n4");
    TH1D *h_acop_n4_c = (TH1D*)h_acop_n4->Clone("h_acop_n4_c");
    TH1D *h_acop_n5 = (TH1D*)f->Get("h_acop_s_n5");
    TH1D *h_acop_n5_c = (TH1D*)h_acop_n5->Clone("h_acop_n5_c");
    TH1D *h_acop_n6 = (TH1D*)f->Get("h_acop_s_n6");
    TH1D *h_acop_n6_c = (TH1D*)h_acop_n6->Clone("h_acop_n6_c");
    TH1D *h_acop_n7 = (TH1D*)f->Get("h_acop_s_n7");
    TH1D *h_acop_n7_c = (TH1D*)h_acop_n7->Clone("h_acop_n7_c");
    TH1D *h_acop_n8 = (TH1D*)f->Get("h_acop_s_n8");
    TH1D *h_acop_n8_c = (TH1D*)h_acop_n8->Clone("h_acop_n8_c");
    TH1D *h_acop_n9 = (TH1D*)f->Get("h_acop_s_n9");
    TH1D *h_acop_n9_c = (TH1D*)h_acop_n9->Clone("h_acop_n9_c");

    for(int i = 0; i < neutron_bin; i++){
        a_n_hist[i] = (TH2D*)f->Get(Form("h_a_n_s_%d", i));
        kperp_hist[i] = (TH1D*)f->Get(Form("h_kperp_s_%d", i));
        mass_hist[i]  = (TH1D*)f->Get(Form("h_mass_s_%d", i));
        rapidity_hist[i] = (TH1D*)f->Get(Form("h_rapidity_s_%d", i));
    }

    //Make the differential plot *****************
    for (int i = 0; i < neutron_bin; i++){
        
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
        h_acop->GetXaxis()->SetTitle("#alpha");
        h_acop->GetYaxis()->SetTitle("(1/N) dN/d#alpha");

        // Store or draw
        a_n_hist[i]->SetTitle(Form("Neutron bin %d", i));
    }


    TCanvas *cops = new TCanvas("cops", "acop", 1200, 1000);
    cops->Divide(3,3);




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

    for (int b = 1; b <= h_acop_n7->GetNbinsX(); b++) {
        double c = h_acop_n7->GetBinContent(b);
        double e = h_acop_n7->GetBinError(b);
        double w = h_acop_n7->GetBinWidth(b);
        h_acop_n7->SetBinContent(b, c / w);
        h_acop_n7_c->SetBinContent(b, c / w);
        h_acop_n7->SetBinError(b, e / w);
        h_acop_n7_c->SetBinError(b, e / w);
    }

    for (int b = 1; b <= h_acop_n8->GetNbinsX(); b++) {
        double c = h_acop_n8->GetBinContent(b);
        double e = h_acop_n8->GetBinError(b);
        double w = h_acop_n8->GetBinWidth(b);
        h_acop_n8->SetBinContent(b, c / w);
        h_acop_n8_c->SetBinContent(b, c / w);
        h_acop_n8->SetBinError(b, e / w);
        h_acop_n8_c->SetBinError(b, e / w);
    }

    for (int b = 1; b <= h_acop_n9->GetNbinsX(); b++) {
        double c = h_acop_n9->GetBinContent(b);
        double e = h_acop_n9->GetBinError(b);
        double w = h_acop_n9->GetBinWidth(b);
        h_acop_n9->SetBinContent(b, c / w);
        h_acop_n9_c->SetBinContent(b, c / w);
        h_acop_n9->SetBinError(b, e / w);
        h_acop_n9_c->SetBinError(b, e / w);
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

    h_acop_n7->Scale(1.0 / h_acop_n7->Integral("width"));
    h_acop_n7_c->Scale(1.0 / h_acop_n7_c->Integral("width"));

    h_acop_n8->Scale(1.0 / h_acop_n8->Integral("width"));
    h_acop_n8_c->Scale(1.0 / h_acop_n8_c->Integral("width"));

    h_acop_n9->Scale(1.0 / h_acop_n9->Integral("width"));
    h_acop_n9_c->Scale(1.0 / h_acop_n9_c->Integral("width"));



////////////************************ FIRST 0N0N ************************////////////
    cops->cd(1);

    gStyle->SetOptTitle(0);

    h_acop_n1->SetLineWidth(2);
    h_acop_n1_c->SetLineWidth(2);
    
    h_acop_n1->SetMinimum(-50);
    h_acop_n1->SetMaximum(800.);

    h_acop_n1->SetMarkerStyle(20);
    h_acop_n1_c->SetMarkerStyle(20);
    h_acop_n1->SetMarkerSize(0.6);
    h_acop_n1_c->SetMarkerSize(0.6);
    h_acop_n1->SetMarkerColor(kBlack);
    h_acop_n1_c->SetMarkerColor(kBlack);
    h_acop_n1->SetLineColor(kBlack);
    h_acop_n1_c->SetLineColor(kBlack);

    tail_fit->SetParameters(1, 1, 1);
    tail_fit->SetParLimits(0, 0, 1e5);
    tail_fit->SetParLimits(1, 0, 1e3);
    tail_fit->SetParLimits(2, 1e-6, 1e3);
    tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n1->Fit(tail_fit, "QRI");
    h_acop_n1_c->Fit(core_fit, "QRI+");

    float c_1_1 = core_fit->GetParameter(0);
    float c_2_1 = core_fit->GetParameter(1);
    float c_3_1 = core_fit->GetParameter(2);

    float t_1_1 = tail_fit->GetParameter(0);
    float t_2_1 = tail_fit->GetParameter(1);
    float t_3_1 = tail_fit->GetParameter(2);

    TF1 *core_funct1 = new TF1("core_funct1", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_1, c_2_1, c_3_1), 0, 0.008);
    TF1 *tail_funct1 = new TF1("tail_funct1", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_1, t_2_1, t_3_1, t_3_1), 0, 0.008);

    core_funct1->SetLineColor(kRed);
    core_funct1->SetLineWidth(3);
    core_funct1->SetLineStyle(10);
    tail_funct1->SetLineColor(kBlue);
    tail_funct1->SetLineWidth(3);
    tail_funct1->SetLineStyle(2);

    h_acop_n1->GetXaxis()->SetTitle("#alpha");
    h_acop_n1->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n1->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n1->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n1->Draw("ELP");
    h_acop_n1_c->Draw("SAME");
    core_funct1->Draw("SAME");
    tail_funct1->Draw("SAME");
    //sum_funct1->Draw("SAME");

    TLatex latex1;
    latex1.SetNDC(); // use normalized coordinates (0–1)
    latex1.SetTextSize(0.05);
    latex1.SetTextAlign(22); // center horizontally
    latex1.DrawLatex(0.8, 0.8, Form("%s", bin_names[0]));
    latex1.DrawLatex(0.8, 0.73, "X #geq 1");
    latex1.DrawLatex(0.8, 0.66, "Y > 1");

    TLegend *leg1 = new TLegend(0.45, 0.8, 1.1, 0.65);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.05);
    leg1->AddEntry(core_funct1, "Core", "l");
    leg1->AddEntry(tail_funct1, "Tail", "l");
    leg1->Draw();



////////////************************ SECOND 0N1N ************************////////////
    cops->cd(2);

    gStyle->SetOptTitle(0);

    h_acop_n2->SetLineWidth(2);
    h_acop_n2_c->SetLineWidth(2);
    
    h_acop_n2->SetMinimum(-50);
    h_acop_n2->SetMaximum(800.);

    h_acop_n2->SetMarkerStyle(20);
    h_acop_n2_c->SetMarkerStyle(20);
    h_acop_n2->SetMarkerSize(0.6);
    h_acop_n2_c->SetMarkerSize(0.6);
    h_acop_n2->SetMarkerColor(kBlack);
    h_acop_n2_c->SetMarkerColor(kBlack);
    h_acop_n2->SetLineColor(kBlack);
    h_acop_n2_c->SetLineColor(kBlack);

    tail_fit->SetParameters(1, 1, 1);
    tail_fit->SetParLimits(0, 0, 1e5);
    tail_fit->SetParLimits(1, 0, 1e3);
    tail_fit->SetParLimits(2, 1e-6, 1e3);
    tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n2->Fit(tail_fit, "QRI");
    h_acop_n2_c->Fit(core_fit, "QRI+");

    float c_1_2 = core_fit->GetParameter(0);
    float c_2_2 = core_fit->GetParameter(1);
    float c_3_2 = core_fit->GetParameter(2);

    float t_1_2 = tail_fit->GetParameter(0);
    float t_2_2 = tail_fit->GetParameter(1);
    float t_3_2 = tail_fit->GetParameter(2);

    TF1 *core_funct2 = new TF1("core_funct2", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_2, c_2_2, c_3_2), 0, 0.008);
    TF1 *tail_funct2 = new TF1("tail_funct2", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_2, t_2_2, t_3_2, t_3_2), 0, 0.008);

    core_funct2->SetLineColor(kRed);
    core_funct2->SetLineWidth(3);
    core_funct2->SetLineStyle(10);
    tail_funct2->SetLineColor(kBlue);
    tail_funct2->SetLineWidth(3);
    tail_funct2->SetLineStyle(2);

    h_acop_n2->GetXaxis()->SetTitle("#alpha");
    h_acop_n2->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n2->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n2->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n2->Draw("ELP");
    h_acop_n2_c->Draw("SAME");
    core_funct2->Draw("SAME");
    tail_funct2->Draw("SAME");
    //sum_funct2->Draw("SAME");

    TLatex latex2;
    latex2.SetNDC(); // use normalized coordinates (0–1)
    latex2.SetTextSize(0.05);
    latex2.SetTextAlign(22); // center horizontally
    latex2.DrawLatex(0.8, 0.8, Form("%s", bin_names[1]));


////////////************************ THIRD 0NYN ************************//////////// put him in 4th spot for clarity
    cops->cd(4);

    gStyle->SetOptTitle(0);

    h_acop_n3->SetLineWidth(2);
    h_acop_n3_c->SetLineWidth(2);
    
    h_acop_n3->SetMinimum(-50);
    h_acop_n3->SetMaximum(800.);

    h_acop_n3->SetMarkerStyle(20);
    h_acop_n3_c->SetMarkerStyle(20);
    h_acop_n3->SetMarkerSize(0.6);
    h_acop_n3_c->SetMarkerSize(0.6);
    h_acop_n3->SetMarkerColor(kBlack);
    h_acop_n3_c->SetMarkerColor(kBlack);
    h_acop_n3->SetLineColor(kBlack);
    h_acop_n3_c->SetLineColor(kBlack);

    tail_fit->SetParameters(1, 1, 1);
    tail_fit->SetParLimits(0, 0, 1e5);
    tail_fit->SetParLimits(1, 0, 1e3);
    tail_fit->SetParLimits(2, 1e-6, 1e3);
    tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
     
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n3->Fit(tail_fit, "QRI");
    h_acop_n3_c->Fit(core_fit, "QRI+");

    float c_1_3 = core_fit->GetParameter(0);
    float c_2_3 = core_fit->GetParameter(1);
    float c_3_3 = core_fit->GetParameter(2);

    float t_1_3 = tail_fit->GetParameter(0);
    float t_2_3 = tail_fit->GetParameter(1);
    float t_3_3 = tail_fit->GetParameter(2);

    TF1 *core_funct3 = new TF1("core_funct3", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_3, c_2_3, c_3_3), 0, 0.008);
    TF1 *tail_funct3 = new TF1("tail_funct3", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_3, t_2_3, t_3_3, t_3_3), 0, 0.008);

    core_funct3->SetLineColor(kRed);
    core_funct3->SetLineWidth(3);
    core_funct3->SetLineStyle(10);
    tail_funct3->SetLineColor(kBlue);
    tail_funct3->SetLineWidth(3);
    tail_funct3->SetLineStyle(2);

    h_acop_n3->GetXaxis()->SetTitle("#alpha");
    h_acop_n3->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n3->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n3->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n3->Draw("ELP");
    h_acop_n3_c->Draw("SAME");
    core_funct3->Draw("SAME");
    tail_funct3->Draw("SAME");
    //sum_funct3->Draw("SAME");

    TLatex latex3;
    latex3.SetNDC(); // use normalized coordinates (0–1)
    latex3.SetTextSize(0.05);
    latex3.SetTextAlign(22); // center horizontally
    latex3.DrawLatex(0.8, 0.8, Form("%s", bin_names[2]));



////////////************************ FOURTH 1N1N ************************//////////// put him in 3rd spot for clarity
    cops->cd(3);

    gStyle->SetOptTitle(0);

    h_acop_n4->SetLineWidth(2);
    h_acop_n4_c->SetLineWidth(2);
    
    h_acop_n4->SetMinimum(-50);
    h_acop_n4->SetMaximum(800.);

    h_acop_n4->SetMarkerStyle(20);
    h_acop_n4_c->SetMarkerStyle(20);
    h_acop_n4->SetMarkerSize(0.6);
    h_acop_n4_c->SetMarkerSize(0.6);
    h_acop_n4->SetMarkerColor(kBlack);
    h_acop_n4_c->SetMarkerColor(kBlack);
    h_acop_n4->SetLineColor(kBlack);
    h_acop_n4_c->SetLineColor(kBlack);

    tail_fit->SetParameters(1, 1, 1);
    tail_fit->SetParLimits(0, 0, 1e5);
    tail_fit->SetParLimits(1, 0, 1e3);
    tail_fit->SetParLimits(2, 1e-6, 1e3);
    tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
     
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n4->Fit(tail_fit, "QRI");
    h_acop_n4_c->Fit(core_fit, "QRI+");

    float c_1_4 = core_fit->GetParameter(0);
    float c_2_4 = core_fit->GetParameter(1);
    float c_3_4 = core_fit->GetParameter(2);

    float t_1_4 = tail_fit->GetParameter(0);
    float t_2_4 = tail_fit->GetParameter(1);
    float t_3_4 = tail_fit->GetParameter(2);

    TF1 *core_funct4 = new TF1("core_funct4", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_4, c_2_4, c_3_4), 0, 0.008);
    TF1 *tail_funct4 = new TF1("tail_funct4", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_4, t_2_4, t_3_4, t_3_4), 0, 0.008);

    core_funct4->SetLineColor(kRed);
    core_funct4->SetLineWidth(3);
    core_funct4->SetLineStyle(10);
    tail_funct4->SetLineColor(kBlue);
    tail_funct4->SetLineWidth(3);
    tail_funct4->SetLineStyle(2);

    h_acop_n4->GetXaxis()->SetTitle("#alpha");
    h_acop_n4->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n4->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n4->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n4->Draw("ELP");
    h_acop_n4_c->Draw("SAME");
    core_funct4->Draw("SAME");
    tail_funct4->Draw("SAME");
    //sum_funct4->Draw("SAME");

    TLatex latex4;
    latex4.SetNDC(); // use normalized coordinates (0–1)
    latex4.SetTextSize(0.05);
    latex4.SetTextAlign(22); // center horizontally
    latex4.DrawLatex(0.8, 0.8, Form("%s", bin_names[3]));




////////////************************ FIFTH 1NYN ************************////////////
    cops->cd(5);

    gStyle->SetOptTitle(0);

    h_acop_n5->SetLineWidth(2);
    h_acop_n5_c->SetLineWidth(2);
    
    h_acop_n5->SetMinimum(-50);
    h_acop_n5->SetMaximum(800.);

    h_acop_n5->SetMarkerStyle(20);
    h_acop_n5_c->SetMarkerStyle(20);
    h_acop_n5->SetMarkerSize(0.6);
    h_acop_n5_c->SetMarkerSize(0.6);
    h_acop_n5->SetMarkerColor(kBlack);
    h_acop_n5_c->SetMarkerColor(kBlack);
    h_acop_n5->SetLineColor(kBlack);
    h_acop_n5_c->SetLineColor(kBlack);

    tail_fit->SetParameters(1, 1, 1);
    tail_fit->SetParLimits(0, 0, 1e5);
    tail_fit->SetParLimits(1, 0, 1e3);
    tail_fit->SetParLimits(2, 1e-6, 1e3);
    tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
     
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n5->Fit(tail_fit, "QRI");
    h_acop_n5_c->Fit(core_fit, "QRI+");

    float c_1_5 = core_fit->GetParameter(0);
    float c_2_5 = core_fit->GetParameter(1);
    float c_3_5 = core_fit->GetParameter(2);

    float t_1_5 = tail_fit->GetParameter(0);
    float t_2_5 = tail_fit->GetParameter(1);
    float t_3_5 = tail_fit->GetParameter(2);

    TF1 *core_funct5 = new TF1("core_funct5", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_5, c_2_5, c_3_5), 0, 0.008);
    TF1 *tail_funct5 = new TF1("tail_funct5", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_5, t_2_5, t_3_5, t_3_5), 0, 0.008);

    core_funct5->SetLineColor(kRed);
    core_funct5->SetLineWidth(3);
    core_funct5->SetLineStyle(10);
    tail_funct5->SetLineColor(kBlue);
    tail_funct5->SetLineWidth(3);
    tail_funct5->SetLineStyle(2);

    h_acop_n5->GetXaxis()->SetTitle("#alpha");
    h_acop_n5->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n5->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n5->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n5->Draw("ELP");
    h_acop_n5_c->Draw("SAME");
    core_funct5->Draw("SAME");
    tail_funct5->Draw("SAME");
    //sum_funct5->Draw("SAME");

    TLatex latex5;
    latex5.SetNDC(); // use normalized coordinates (0–1)
    latex5.SetTextSize(0.05);
    latex5.SetTextAlign(22); // center horizontally
    latex5.DrawLatex(0.8, 0.8, Form("%s", bin_names[4]));



////////////************************ SIXTH YNYN ************************////////////
    cops->cd(6);

    gStyle->SetOptTitle(0);

    h_acop_n6->SetLineWidth(2);
    h_acop_n6_c->SetLineWidth(2);
    
    h_acop_n6->SetMinimum(-50);
    h_acop_n6->SetMaximum(800.);

    h_acop_n6->SetMarkerStyle(20);
    h_acop_n6_c->SetMarkerStyle(20);
    h_acop_n6->SetMarkerSize(0.6);
    h_acop_n6_c->SetMarkerSize(0.6);
    h_acop_n6->SetMarkerColor(kBlack);
    h_acop_n6_c->SetMarkerColor(kBlack);
    h_acop_n6->SetLineColor(kBlack);
    h_acop_n6_c->SetLineColor(kBlack);

    tail_fit->SetParameters(2, 5, 20);
    tail_fit->SetParLimits(500, 1000, 10);
    //tail_fit->SetParLimits(1, 0, 1e3);
    //tail_fit->SetParLimits(2, 1e-6, 1e3);
    //tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
     
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n6->Fit(tail_fit, "QRI");
    h_acop_n6_c->Fit(core_fit, "QRI+");

    float c_1_6 = core_fit->GetParameter(0);
    float c_2_6 = core_fit->GetParameter(1);
    float c_3_6 = core_fit->GetParameter(2);

    float t_1_6 = tail_fit->GetParameter(0);
    float t_2_6 = tail_fit->GetParameter(1);
    float t_3_6 = tail_fit->GetParameter(2);

    TF1 *core_funct6 = new TF1("core_funct6", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_6, c_2_6, c_3_6), 0, 0.008);
    TF1 *tail_funct6 = new TF1("tail_funct6", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_1, t_2_6, t_3_6, t_3_6), 0, 0.008);

    core_funct6->SetLineColor(kRed);
    core_funct6->SetLineWidth(3);
    core_funct6->SetLineStyle(10);
    tail_funct6->SetLineColor(kBlue);
    tail_funct6->SetLineWidth(3);
    tail_funct6->SetLineStyle(2);

    h_acop_n6->GetXaxis()->SetTitle("#alpha");
    h_acop_n6->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n6->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n6->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n6->Draw("ELP");
    h_acop_n6_c->Draw("SAME");
    core_funct6->Draw("SAME");
    tail_funct6->Draw("SAME");
    //sum_funct6->Draw("SAME");

    TLatex latex6;
    latex6.SetNDC(); // use normalized coordinates (0–1)
    latex6.SetTextSize(0.05);
    latex6.SetTextAlign(22); // center horizontally
    latex6.DrawLatex(0.8, 0.8, Form("%s", bin_names[5]));


////////////************************ SEVENTH 0NXN ************************////////////
    cops->cd(7);

    gStyle->SetOptTitle(0);

    h_acop_n7->SetLineWidth(2);
    h_acop_n7_c->SetLineWidth(2);
    
    h_acop_n7->SetMinimum(-50);
    h_acop_n7->SetMaximum(800.);

    h_acop_n7->SetMarkerStyle(20);
    h_acop_n7_c->SetMarkerStyle(20);
    h_acop_n7->SetMarkerSize(0.6);
    h_acop_n7_c->SetMarkerSize(0.6);
    h_acop_n7->SetMarkerColor(kBlack);
    h_acop_n7_c->SetMarkerColor(kBlack);
    h_acop_n7->SetLineColor(kBlack);
    h_acop_n7_c->SetLineColor(kBlack);

    tail_fit->SetParameters(2, 5, 20);
    tail_fit->SetParLimits(500, 1000, 10);
    //tail_fit->SetParLimits(1, 0, 1e3);
    //tail_fit->SetParLimits(2, 1e-6, 1e3);
    //tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
     
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n7->Fit(tail_fit, "QRI");
    h_acop_n7_c->Fit(core_fit, "QRI+");

    float c_1_7 = core_fit->GetParameter(0);
    float c_2_7 = core_fit->GetParameter(1);
    float c_3_7 = core_fit->GetParameter(2);

    float t_1_7 = tail_fit->GetParameter(0);
    float t_2_7 = tail_fit->GetParameter(1);
    float t_3_7 = tail_fit->GetParameter(2);

    TF1 *core_funct7 = new TF1("core_funct7", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_7, c_2_7, c_3_7), 0, 0.008);
    TF1 *tail_funct7 = new TF1("tail_funct7", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_1, t_2_7, t_3_7, t_3_7), 0, 0.008);

    core_funct7->SetLineColor(kRed);
    core_funct7->SetLineWidth(3);
    core_funct7->SetLineStyle(10);
    tail_funct7->SetLineColor(kBlue);
    tail_funct7->SetLineWidth(3);
    tail_funct7->SetLineStyle(2);

    h_acop_n7->GetXaxis()->SetTitle("#alpha");
    h_acop_n7->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n7->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n7->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n7->Draw("ELP");
    h_acop_n7_c->Draw("SAME");
    core_funct7->Draw("SAME");
    tail_funct7->Draw("SAME");
    //sum_funct7->Draw("SAME");

    TLatex latex7;
    latex7.SetNDC(); // use normalized coordinates (0–1)
    latex7.SetTextSize(0.05);
    latex7.SetTextAlign(22); // center horizontally
    latex7.DrawLatex(0.8, 0.8, Form("%s", bin_names[6]));



////////////************************ EIGTH 1NXN ************************////////////
    cops->cd(8);

    gStyle->SetOptTitle(0);

    h_acop_n8->SetLineWidth(2);
    h_acop_n8_c->SetLineWidth(2);
    
    h_acop_n8->SetMinimum(-50);
    h_acop_n8->SetMaximum(800.);

    h_acop_n8->SetMarkerStyle(20);
    h_acop_n8_c->SetMarkerStyle(20);
    h_acop_n8->SetMarkerSize(0.6);
    h_acop_n8_c->SetMarkerSize(0.6);
    h_acop_n8->SetMarkerColor(kBlack);
    h_acop_n8_c->SetMarkerColor(kBlack);
    h_acop_n8->SetLineColor(kBlack);
    h_acop_n8_c->SetLineColor(kBlack);

    tail_fit->SetParameters(2, 5, 20);
    tail_fit->SetParLimits(500, 1000, 10);
    //tail_fit->SetParLimits(1, 0, 1e3);
    //tail_fit->SetParLimits(2, 1e-6, 1e3);
    //tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
     
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n8->Fit(tail_fit, "QRI");
    h_acop_n8_c->Fit(core_fit, "QRI+");

    float c_1_8 = core_fit->GetParameter(0);
    float c_2_8 = core_fit->GetParameter(1);
    float c_3_8 = core_fit->GetParameter(2);

    float t_1_8 = tail_fit->GetParameter(0);
    float t_2_8 = tail_fit->GetParameter(1);
    float t_3_8 = tail_fit->GetParameter(2);

    TF1 *core_funct8 = new TF1("core_funct8", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_8, c_2_8, c_3_8), 0, 0.008);
    TF1 *tail_funct8 = new TF1("tail_funct8", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_1, t_2_8, t_3_8, t_3_8), 0, 0.008);

    core_funct8->SetLineColor(kRed);
    core_funct8->SetLineWidth(3);
    core_funct8->SetLineStyle(10);
    tail_funct8->SetLineColor(kBlue);
    tail_funct8->SetLineWidth(3);
    tail_funct8->SetLineStyle(2);

    h_acop_n8->GetXaxis()->SetTitle("#alpha");
    h_acop_n8->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n8->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n8->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n8->Draw("ELP");
    h_acop_n8_c->Draw("SAME");
    core_funct8->Draw("SAME");
    tail_funct8->Draw("SAME");
    //sum_funct8->Draw("SAME");

    TLatex latex8;
    latex8.SetNDC(); // use normalized coordinates (0–1)
    latex8.SetTextSize(0.05);
    latex8.SetTextAlign(22); // center horizontally
    latex8.DrawLatex(0.8, 0.8, Form("%s", bin_names[7]));



////////////************************ NINTH XNXN ************************////////////
    cops->cd(9);

    gStyle->SetOptTitle(0);

    h_acop_n9->SetLineWidth(2);
    h_acop_n9_c->SetLineWidth(2);
    
    h_acop_n9->SetMinimum(-50);
    h_acop_n9->SetMaximum(800.);

    h_acop_n9->SetMarkerStyle(20);
    h_acop_n9_c->SetMarkerStyle(20);
    h_acop_n9->SetMarkerSize(0.6);
    h_acop_n9_c->SetMarkerSize(0.6);
    h_acop_n9->SetMarkerColor(kBlack);
    h_acop_n9_c->SetMarkerColor(kBlack);
    h_acop_n9->SetLineColor(kBlack);
    h_acop_n9_c->SetLineColor(kBlack);

    tail_fit->SetParameters(2, 5, 20);
    tail_fit->SetParLimits(500, 1000, 10);
    //tail_fit->SetParLimits(1, 0, 1e3);
    //tail_fit->SetParLimits(2, 1e-6, 1e3);
    //tail_fit->SetParLimits(3, 0, 10);
    tail_fit->SetRange(0.75*10e-3, 1);

    core_fit->SetParameters(1, 1, 1);
    core_fit->SetParLimits(0, 0, 1e5);
    core_fit->SetParLimits(1, 0, 1e3);
    core_fit->SetParLimits(2, 1e-6, 1e3);
     
    core_fit->SetRange(0.5*10e-5, 0.75*10e-3);

    h_acop_n9->Fit(tail_fit, "QRI");
    h_acop_n9_c->Fit(core_fit, "QRI+");

    float c_1_9 = core_fit->GetParameter(0);
    float c_2_9 = core_fit->GetParameter(1);
    float c_3_9 = core_fit->GetParameter(2);

    float t_1_9 = tail_fit->GetParameter(0);
    float t_2_9 = tail_fit->GetParameter(1);
    float t_3_9 = tail_fit->GetParameter(2);

    TF1 *core_funct9 = new TF1("core_funct9", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", c_1_9, c_2_9, c_3_9), 0, 0.008);
    TF1 *tail_funct9 = new TF1("tail_funct9", Form("%f * pow((1 + x* %f/%f), -%f)", t_1_1, t_2_9, t_3_9, t_3_9), 0, 0.008);

    core_funct9->SetLineColor(kRed);
    core_funct9->SetLineWidth(3);
    core_funct9->SetLineStyle(10);
    tail_funct9->SetLineColor(kBlue);
    tail_funct9->SetLineWidth(3);
    tail_funct9->SetLineStyle(2);

    h_acop_n9->GetXaxis()->SetTitle("#alpha");
    h_acop_n9->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
    h_acop_n9->GetXaxis()->SetTitleOffset(1.4);
    h_acop_n9->GetYaxis()->SetTitleOffset(1.4);

    h_acop_n9->Draw("ELP");
    h_acop_n9_c->Draw("SAME");
    core_funct9->Draw("SAME");
    tail_funct9->Draw("SAME");
    //sum_funct9->Draw("SAME");

    TLatex latex9;
    latex9.SetNDC(); // use normalized coordinates (0–1)
    latex9.SetTextSize(0.05);
    latex9.SetTextAlign(22); // center horizontally
    latex9.DrawLatex(0.8, 0.8, Form("%s", bin_names[8]));


    cops->SaveAs("plots/small/plot2_12.26.png");



    TCanvas *perps = new TCanvas("perps", "kperp", 1200, 1000);
    perps->Divide(3, 3);
    for (int i = 0; i < neutron_bin; i++){
        perps->cd(i+1);
        gPad->SetLogy();
        kperp_hist[i]->Draw();

        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.8, 0.8, Form("%s", bin_names[i]));
        if (i == 0){
            latex.DrawLatex(0.8, 0.73, "X #geq 1");
            latex.DrawLatex(0.8, 0.66, "Y > 1");

            latex.SetTextSize(0.06);
            latex.DrawLatex(0.8, 0.88, "k_{#perp}  < 0.008");
        }
    }
    perps->SaveAs("plots/small/kperp2_12.26.png");



// plot wanted here for upsilon concern. Masses below 60 MeV
    TCanvas *masses = new TCanvas("masses", "mass", 1200, 1000);
    masses->Divide(3,3);
    for (int i = 0; i < neutron_bin; i++){
        masses->cd(i+1);
        gPad->SetLogy();
        mass_hist[i]->Draw();

        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.8, 0.8, Form("%s", bin_names[i]));

        if (i == 0){
            latex.DrawLatex(0.8, 0.73, "X #geq 1");
            latex.DrawLatex(0.8, 0.66, "Y > 1");

            latex.SetTextSize(0.06);
            latex.DrawLatex(0.8, 0.88, "m_{#mu} < 60 MeV");
        }
    }
    masses->SaveAs("plots/small/masses2_12.26.png");





// rms(alpha) plot. rms^2 = STD^2 + mean^2 #######################
    TH1D* alpha_rms = new TH1D("alpha_rms",";Neutron bin;RMS",neutron_bin, 0.5, neutron_bin + 0.5);

    double std1 = h_acop_n1->GetRMS();
    double std2 = h_acop_n2->GetRMS();
    double std3 = h_acop_n3->GetRMS();
    double std4 = h_acop_n4->GetRMS();
    double std5 = h_acop_n5->GetRMS();
    double std6 = h_acop_n6->GetRMS();
    double std7 = h_acop_n7->GetRMS();
    double std8 = h_acop_n8->GetRMS();
    double std9 = h_acop_n9->GetRMS();

    double std1_err = h_acop_n1->GetRMSError();
    double std2_err = h_acop_n2->GetRMSError();
    double std3_err = h_acop_n3->GetRMSError();
    double std4_err = h_acop_n4->GetRMSError();
    double std5_err = h_acop_n5->GetRMSError();
    double std6_err = h_acop_n6->GetRMSError();
    double std7_err = h_acop_n7->GetRMSError();
    double std8_err = h_acop_n8->GetRMSError();
    double std9_err = h_acop_n9->GetRMSError();

    double mean1 = h_acop_n1->GetMean();
    double mean2 = h_acop_n2->GetMean();
    double mean3 = h_acop_n3->GetMean();
    double mean4 = h_acop_n4->GetMean();
    double mean5 = h_acop_n5->GetMean();
    double mean6 = h_acop_n6->GetMean();
    double mean7 = h_acop_n7->GetMean();
    double mean8 = h_acop_n8->GetMean();
    double mean9 = h_acop_n9->GetMean();

    double mean1_err = h_acop_n1->GetMeanError();
    double mean2_err = h_acop_n2->GetMeanError();
    double mean3_err = h_acop_n3->GetMeanError();
    double mean4_err = h_acop_n4->GetMeanError();
    double mean5_err = h_acop_n5->GetMeanError();
    double mean6_err = h_acop_n6->GetMeanError();
    double mean7_err = h_acop_n7->GetMeanError();
    double mean8_err = h_acop_n8->GetMeanError();
    double mean9_err = h_acop_n9->GetMeanError();

    double rms1 = sqrt(pow(std1, 2) + pow(mean1, 2));
    double rms2 = sqrt(pow(std2, 2) + pow(mean2, 2));
    double rms3 = sqrt(pow(std3, 2) + pow(mean3, 2));
    double rms4 = sqrt(pow(std4, 2) + pow(mean4, 2));
    double rms5 = sqrt(pow(std5, 2) + pow(mean5, 2));
    double rms6 = sqrt(pow(std6, 2) + pow(mean6, 2));
    double rms7 = sqrt(pow(std7, 2) + pow(mean7, 2));
    double rms8 = sqrt(pow(std8, 2) + pow(mean8, 2));
    double rms9 = sqrt(pow(std9, 2) + pow(mean9, 2));

    double rms1_err = sqrt( pow(std1 * std1_err / rms1, 2) + pow(mean1 * mean1_err / rms1, 2));
    double rms2_err = sqrt( pow(std2 * std2_err / rms2, 2) + pow(mean2 * mean2_err / rms2, 2));
    double rms3_err = sqrt( pow(std3 * std3_err / rms3, 2) + pow(mean3 * mean3_err / rms3, 2));
    double rms4_err = sqrt( pow(std4 * std4_err / rms4, 2) + pow(mean4 * mean4_err / rms4, 2));
    double rms5_err = sqrt( pow(std5 * std5_err / rms5, 2) + pow(mean5 * mean5_err / rms5, 2));
    double rms6_err = sqrt( pow(std6 * std6_err / rms6, 2) + pow(mean6 * mean6_err / rms6, 2));
    double rms7_err = sqrt( pow(std7 * std7_err / rms7, 2) + pow(mean7 * mean7_err / rms7, 2));
    double rms8_err = sqrt( pow(std8 * std8_err / rms8, 2) + pow(mean8 * mean8_err / rms8, 2));
    double rms9_err = sqrt( pow(std9 * std9_err / rms9, 2) + pow(mean9 * mean9_err / rms9, 2));

    alpha_rms->SetBinContent(1, rms1);
    alpha_rms->SetBinContent(2, rms2);
    alpha_rms->SetBinContent(3, rms3);
    alpha_rms->SetBinContent(4, rms4);
    alpha_rms->SetBinContent(5, rms5);
    alpha_rms->SetBinContent(6, rms6);
    alpha_rms->SetBinContent(7, rms7);
    alpha_rms->SetBinContent(8, rms8);
    alpha_rms->SetBinContent(9, rms9);

    alpha_rms->SetBinError  (1, rms1_err);
    alpha_rms->SetBinError  (2, rms2_err);
    alpha_rms->SetBinError  (3, rms3_err);
    alpha_rms->SetBinError  (4, rms4_err);
    alpha_rms->SetBinError  (5, rms5_err);
    alpha_rms->SetBinError  (6, rms6_err);
    alpha_rms->SetBinError  (7, rms7_err);
    alpha_rms->SetBinError  (8, rms8_err);
    alpha_rms->SetBinError  (9, rms9_err);

    alpha_rms->GetXaxis()->SetBinLabel(1, bin_names[0]);
    alpha_rms->GetXaxis()->SetBinLabel(2, bin_names[1]);
    alpha_rms->GetXaxis()->SetBinLabel(3, bin_names[2]);
    alpha_rms->GetXaxis()->SetBinLabel(4, bin_names[3]);
    alpha_rms->GetXaxis()->SetBinLabel(5, bin_names[4]);
    alpha_rms->GetXaxis()->SetBinLabel(6, bin_names[5]);
    alpha_rms->GetXaxis()->SetBinLabel(7, bin_names[6]);
    alpha_rms->GetXaxis()->SetBinLabel(8, bin_names[7]);
    alpha_rms->GetXaxis()->SetBinLabel(9, bin_names[8]);

    TCanvas *rms = new TCanvas("rms", "alpha RMS", 1200, 1000);
    alpha_rms->SetStats(0);
    alpha_rms->SetFillColor(kCyan + 1);
    alpha_rms->Draw("E1 BAR");

    rms->SaveAs("plots/small/alpha_rms2_12.26.png");


//rapidity plots
    TCanvas *rapid = new TCanvas("rapid", "rapidity", 1200, 1000);
    rapid->Divide(3,3);
    for (int i = 0; i < neutron_bin; i++){
        rapid->cd(i+1);
        rapidity_hist[i]->Draw();
    
        TLatex latex;
        latex.SetNDC(); // use normalized coordinates (0–1)
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center horizontally
        latex.DrawLatex(0.85, 0.8, Form("%s", bin_names[i]));

        if (i == 0){
            latex.DrawLatex(0.85, 0.73, "X #geq 1");
            latex.DrawLatex(0.85, 0.66, "Y > 1");

            latex.SetTextSize(0.06);
            latex.DrawLatex(0.85, 0.88, "Rapidity");
        }
    }
    rapid->SaveAs("plots/small/rapidity2_12.26.png");

}