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


void acop_try(){
    TFile *f = TFile::Open("dimuon_hists.root");
    static const int neutron_bin = 12;
    static const int neutron1    = 6;
    static const int neutron2    = 12;
    const char* bin_names[neutron_bin] = {"0n0n", "0n1n", "1n1n", "Xn0n", "Xn1n", "XnXn", "0nYn", "Yn1n", "YnYn", "Zn0n", "Zn1n", "ZnZn"}; //lined up so that they're easier to compare

    TF1 * core_fit = new TF1("core_fit","[0]*exp(-x/[1]+[2]*pow(x,0.25))", 0.5*10e-5, 1); //lower bound taken from current def of histogram
    TF1 * tail_fit = new TF1("tail_fit", "[0]*pow((1 + x*[1]/[2]), -[2])", 0.5*10e-5, 1);
    TF1 * sum_fit  = new TF1("sum_fit", "[0]*exp(-x/[1]+[2]*pow(x,0.25))+[3]*pow((1 + x*[4]/[5]), -[6])", 0.5*10e-5, 1);


    TH2D *a_n_hist[neutron_bin];
    TH1D *kperp_hist[neutron_bin];
    TH1D *mass_hist[neutron_bin];
    TH1D *rapidity_hist[neutron_bin];
    TH1D *acop_hist[neutron_bin];
    TH1D *acop_hist_c[neutron_bin];

    for(int i = 0; i < neutron_bin; i++){
        acop_hist[i]     = (TH1D*)f->Get(Form("h_acop_n%d", i+1));
        acop_hist_c[i]   = (TH1D*)acop_hist[i]->Clone(Form("h_acop_n%d_c",i+1));

        a_n_hist[i]      = (TH2D*)f->Get(Form("h_a_n_%d", i));
        kperp_hist[i]    = (TH1D*)f->Get(Form("h_kperp_%d", i));
        mass_hist[i]     = (TH1D*)f->Get(Form("h_mass_%d", i));
        rapidity_hist[i] = (TH1D*)f->Get(Form("h_rapidity_%d", i));
    }

    TCanvas *acop1 = new TCanvas("acop1", "acop", 1200, 1000);
    acop1->Divide(3,2);
    for (int i = 0; i < neutron1; i++){ // for clarity, only plot half of them on each
        for (int b = 1; b <= acop_hist[i]->GetNbinsX(); b++) {
            double c = acop_hist[i]->GetBinContent(b);
            double e = acop_hist[i]->GetBinError(b);
            double w = acop_hist[i]->GetBinWidth(b);
            acop_hist[i]->SetBinContent(b, c / w);
            acop_hist_c[i]->SetBinContent(b, c / w);
            acop_hist[i]->SetBinError(b, e / w);
            acop_hist_c[i]->SetBinError(b, e / w);
        }
        acop_hist[i]->Scale(1.0 / acop_hist[i]->Integral("width"));
        acop_hist_c[i]->Scale(1.0 / acop_hist_c[i]->Integral("width"));

        acop1->cd(i+1);

        gStyle->SetOptTitle(0);

        acop_hist[i]->SetLineWidth(2);
        acop_hist_c[i]->SetLineWidth(2);
        
        acop_hist[i]->SetMinimum(10e-4);
        acop_hist[i]->SetMaximum(800.);

        acop_hist[i]->SetMarkerStyle(20);
        acop_hist_c[i]->SetMarkerStyle(20);
        acop_hist[i]->SetMarkerSize(0.6);
        acop_hist_c[i]->SetMarkerSize(0.6);
        acop_hist[i]->SetMarkerColor(kBlack);
        acop_hist_c[i]->SetMarkerColor(kBlack);
        acop_hist[i]->SetLineColor(kBlack);
        acop_hist_c[i]->SetLineColor(kBlack);

        sum_fit->SetParameters(1, 1, 1, 1, 1, 1, 1);
        sum_fit->SetParLimits(0, 0, 1e5);
        sum_fit->SetParLimits(1, 0, 1e3);
        sum_fit->SetParLimits(2, 1e-6, 1e3);
        sum_fit->SetParLimits(3, 0, 10);
        sum_fit->SetParLimits(4, 0, 1e5);
        sum_fit->SetParLimits(5, 0, 1e3);
        sum_fit->SetParLimits(6, 1e-6, 1e3);
        sum_fit->SetParLimits(7, 0, 10);

        acop_hist[i]->Fit(sum_fit, "QRI");

        float t_1_1 = tail_fit->GetParameter(0); //cut when done

        float s_1_1 = sum_fit->GetParameter(0);
        float s_2_1 = sum_fit->GetParameter(1);
        float s_3_1 = sum_fit->GetParameter(2);
        float s_4_1 = sum_fit->GetParameter(3);
        float s_5_1 = sum_fit->GetParameter(4);
        float s_6_1 = sum_fit->GetParameter(5);
        float s_7_1 = sum_fit->GetParameter(6);

        TF1 *core_funct1 = new TF1("core_funct1", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", s_1_1, s_2_1, s_3_1), 0.5*10e-5, 1);
        TF1 *tail_funct1 = new TF1("tail_funct1", Form("%f * pow((1 + x* %f/%f), -%f)", s_4_1, s_5_1, s_6_1, s_7_1), 0.5*10e-5, 1);
        TF1 *sum_funct1 = new TF1("sum_funct1", Form("%f * exp(-x / %f + %f * pow(x, 0.25)) + %f * pow((1 + x* %f/%f), -%f)", s_1_1, s_2_1, s_3_1, s_4_1, s_5_1, s_6_1, s_7_1), 0.5*10e-5, 1);

        core_funct1->SetLineColor(kRed);
        core_funct1->SetLineWidth(3);
        core_funct1->SetLineStyle(10);
        tail_funct1->SetLineColor(kBlue);
        tail_funct1->SetLineWidth(3);
        tail_funct1->SetLineStyle(2);

        sum_funct1->SetLineColor(kMagenta);

        acop_hist[i]->GetXaxis()->SetTitle("#alpha");
        acop_hist[i]->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
        acop_hist[i]->GetXaxis()->SetTitleOffset(1.4);
        acop_hist[i]->GetYaxis()->SetTitleOffset(1.4);


        acop_hist[i]->Draw("ELP");
        core_funct1->Draw("SAME");
        tail_funct1->Draw("SAME");
        sum_funct1->Draw("SAME");

        gPad->SetLogx();
        gPad->SetLogy();

        TLatex latex1;
        latex1.SetNDC(); // use normalized coordinates (0–1)
        latex1.SetTextSize(0.05);
        latex1.SetTextAlign(22); // center horizontally
        latex1.DrawLatex(0.8, 0.8, Form("%s", bin_names[i]));
        latex1.DrawLatex(0.8, 0.75, "X #geq 1");
        latex1.DrawLatex(0.8, 0.68, "1 < Y #leq 5");
        latex1.DrawLatex(0.8, 0.61, "Z > 5");

        TLegend *leg1 = new TLegend(0.2, 0.2, 1.1, 0.4);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->SetTextSize(0.05);
        leg1->AddEntry(core_funct1, "Core", "l");
        leg1->AddEntry(tail_funct1, "Tail", "l");
        leg1->AddEntry(sum_funct1, "Sum", "l");
        leg1->Draw();

    }


TCanvas *acop2 = new TCanvas("acop2", "acop", 1200, 1000);
    acop2->Divide(3,2);
    for (int i = 6; i < neutron_bin; i++){ // second set of plots
        for (int b = 1; b <= acop_hist[i]->GetNbinsX(); b++) {
            double c = acop_hist[i]->GetBinContent(b);
            double e = acop_hist[i]->GetBinError(b);
            double w = acop_hist[i]->GetBinWidth(b);
            acop_hist[i]->SetBinContent(b, c / w);
            acop_hist_c[i]->SetBinContent(b, c / w);
            acop_hist[i]->SetBinError(b, e / w);
            acop_hist_c[i]->SetBinError(b, e / w);
        }
        acop_hist[i]->Scale(1.0 / acop_hist[i]->Integral("width"));
        acop_hist_c[i]->Scale(1.0 / acop_hist_c[i]->Integral("width"));

        acop2->cd(i-5);

        gStyle->SetOptTitle(0);

        acop_hist[i]->SetLineWidth(2);
        acop_hist_c[i]->SetLineWidth(2);
        
        acop_hist[i]->SetMinimum(10e-4);
        acop_hist[i]->SetMaximum(800.);

        acop_hist[i]->SetMarkerStyle(20);
        acop_hist_c[i]->SetMarkerStyle(20);
        acop_hist[i]->SetMarkerSize(0.6);
        acop_hist_c[i]->SetMarkerSize(0.6);
        acop_hist[i]->SetMarkerColor(kBlack);
        acop_hist_c[i]->SetMarkerColor(kBlack);
        acop_hist[i]->SetLineColor(kBlack);
        acop_hist_c[i]->SetLineColor(kBlack);

        sum_fit->SetParameters(1, 1, 1, 1, 1, 1, 1);
        sum_fit->SetParLimits(0, 0, 1e5);
        sum_fit->SetParLimits(1, 0, 1e3);
        sum_fit->SetParLimits(2, 1e-6, 1e3);
        sum_fit->SetParLimits(3, 0, 10);
        sum_fit->SetParLimits(4, 0, 1e5);
        sum_fit->SetParLimits(5, 0, 1e3);
        sum_fit->SetParLimits(6, 1e-6, 1e3);
        sum_fit->SetParLimits(7, 0, 10);

        acop_hist[i]->Fit(sum_fit, "QRI");

        float t_1_1 = tail_fit->GetParameter(0); //cut when done

        float s_1_1 = sum_fit->GetParameter(0);
        float s_2_1 = sum_fit->GetParameter(1);
        float s_3_1 = sum_fit->GetParameter(2);
        float s_4_1 = sum_fit->GetParameter(3);
        float s_5_1 = sum_fit->GetParameter(4);
        float s_6_1 = sum_fit->GetParameter(5);
        float s_7_1 = sum_fit->GetParameter(6);

        TF1 *core_funct1 = new TF1("core_funct1", Form("%f * exp(-x / %f + %f * pow(x, 0.25))", s_1_1, s_2_1, s_3_1), 0.5*10e-5, 1);
        TF1 *tail_funct1 = new TF1("tail_funct1", Form("%f * pow((1 + x* %f/%f), -%f)", s_4_1, s_5_1, s_6_1, s_7_1), 0.5*10e-5, 1);
        TF1 *sum_funct1 = new TF1("sum_funct1", Form("%f * exp(-x / %f + %f * pow(x, 0.25)) + %f * pow((1 + x* %f/%f), -%f)", s_1_1, s_2_1, s_3_1, s_4_1, s_5_1, s_6_1, s_7_1), 0.5*10e-5, 1);

        core_funct1->SetLineColor(kRed);
        core_funct1->SetLineWidth(3);
        core_funct1->SetLineStyle(10);
        tail_funct1->SetLineColor(kBlue);
        tail_funct1->SetLineWidth(3);
        tail_funct1->SetLineStyle(2);

        sum_funct1->SetLineColor(kMagenta);

        acop_hist[i]->GetXaxis()->SetTitle("#alpha");
        acop_hist[i]->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
        acop_hist[i]->GetXaxis()->SetTitleOffset(1.4);
        acop_hist[i]->GetYaxis()->SetTitleOffset(1.4);


        acop_hist[i]->Draw("ELP");
        core_funct1->Draw("SAME");
        tail_funct1->Draw("SAME");
        sum_funct1->Draw("SAME");

        gPad->SetLogx();
        gPad->SetLogy();

        TLatex latex1;
        latex1.SetNDC(); // use normalized coordinates (0–1)
        latex1.SetTextSize(0.05);
        latex1.SetTextAlign(22); // center horizontally
        latex1.DrawLatex(0.8, 0.8, Form("%s", bin_names[i]));
        latex1.DrawLatex(0.8, 0.75, "X #geq 1");
        latex1.DrawLatex(0.8, 0.68, "1 < Y #leq 5");
        latex1.DrawLatex(0.8, 0.61, "Z > 5");

        TLegend *leg1 = new TLegend(0.2, 0.2, 1.1, 0.4);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->SetTextSize(0.05);
        leg1->AddEntry(core_funct1, "Core", "l");
        leg1->AddEntry(tail_funct1, "Tail", "l");
        leg1->AddEntry(sum_funct1, "Sum", "l");
        leg1->Draw();

    }
}
