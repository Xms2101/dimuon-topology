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
#include <cmath>
#include <vector>

#include "TMatrixDSym.h" // these are to run cov_mtrx.C for the covariance matrix
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <stdexcept>

void full_acop_plot(){
    TFile *f = TFile::Open("dimuon_hists.root");
    static const int neutron_bin = 13;
    static const int neutron1    = 7;
    const char* bin_names[neutron_bin] = {"0n0n", "1n0n", "2n0n", "Xn0n", "Yn0n", "Zn0n", "1n1n", "Xn1n", "Yn1n", "Zn1n", "XnXn", "YnYn", "ZnZn"}; //Organized as 0n, 1n, Nn

    bool save_plots = false;

    TF1 * core_fit = new TF1("core_fit","[0]*exp(-x/[1]+[2]*pow(x,0.25))", 0.5*10e-5, 1); //lower bound taken from current def of histogram
    TF1 * tail_fit = new TF1("tail_fit", "[0]*pow((1 + x*[1]/[2]), -[2])", 0.5*10e-5, 1);
    TF1 * sum_fit  = new TF1("sum_fit", "[0]*exp(-x/[1]+[2]*pow(x,0.25))+[3]*pow((1 + x*([4]/[5])), -[6])", 0.5*10e-5, 1);

    TH2D *a_n_hist[neutron_bin];
    TH1D *kperp_hist[neutron_bin];
    TH1D *mass_hist[neutron_bin];
    TH1D *rapidity_hist[neutron_bin];
    TH1D *acop_hist[neutron_bin];
    TH1D *acop_hist_c[neutron_bin];

    TH1D *core_expectation = new TH1D("core_expectation", "", neutron_bin, 0, neutron_bin);

    TH1D *qed_tail={nullptr};
    TH1D *dis_tail={nullptr};

    vector<double> core_1(neutron_bin);
    vector<double> core_2(neutron_bin);
    vector<double> core_3(neutron_bin);

    vector<double> core_1_err(neutron_bin);
    vector<double> core_2_err(neutron_bin);
    vector<double> core_3_err(neutron_bin);


    for(int i = 0; i < neutron_bin; i++){
        acop_hist[i]     = (TH1D*)f->Get(Form("h_acop_n%d", i+1));
        acop_hist_c[i]   = (TH1D*)acop_hist[i]->Clone(Form("h_acop_n%d_c",i+1));

        a_n_hist[i]      = (TH2D*)f->Get(Form("h_a_n_%d", i));
        kperp_hist[i]    = (TH1D*)f->Get(Form("h_kperp_%d", i));
        mass_hist[i]     = (TH1D*)f->Get(Form("h_mass_%d", i));
        rapidity_hist[i] = (TH1D*)f->Get(Form("h_rapidity_%d", i));

        core_expectation->GetXaxis()->SetBinLabel(i+1, bin_names[i]);
    }
    qed_tail = (TH1D*)f->Get("h_acop_n1")->Clone("qed_tail"); //0n0n
    qed_tail->SetDirectory(nullptr); // to make sure that we're not recycling the pointer
    dis_tail = (TH1D*)f->Get("h_acop_n5")->Clone("dis_tail"); //Zn0n
    dis_tail->SetDirectory(nullptr);


//############################### create template fit ###############################
    TF1 *sum_fit_qed = (TF1*)sum_fit->Clone("sum_fit_qed");
    TF1 *sum_fit_dis = (TF1*)sum_fit->Clone("sum_fit_dis");

    qed_tail->SetMinimum(10e-4);
    qed_tail->SetMaximum(800.);
    dis_tail->SetMinimum(10e-4);
    dis_tail->SetMaximum(800.);

    for (int b = 1; b <= qed_tail->GetNbinsX(); b++){
            double c = qed_tail->GetBinContent(b);
            double e = qed_tail->GetBinError(b);
            double w = qed_tail->GetBinWidth(b);
            qed_tail->SetBinContent(b, c / w);
            qed_tail->SetBinError(b, e / w);
        }
        qed_tail->Scale(1.0 / qed_tail->Integral("width"));
    for (int b = 1; b <= dis_tail->GetNbinsX(); b++){
            double c = dis_tail->GetBinContent(b);
            double e = dis_tail->GetBinError(b);
            double w = dis_tail->GetBinWidth(b);
            dis_tail->SetBinContent(b, c / w);
            dis_tail->SetBinError(b, e / w);
        }
        dis_tail->Scale(1.0 / dis_tail->Integral("width"));

    sum_fit_qed->SetParameters(qed_tail->GetBinContent(1), 2e-3, 0.3, 0.02, 30.0, 4.5, 1);
    sum_fit_qed->SetParLimits(0, 1e2, 1e5);
    sum_fit_qed->SetParLimits(1, 0, 1e3);
    sum_fit_qed->SetParLimits(2, 1e-6, 1e3);
    sum_fit_qed->SetParLimits(3, 0, 1e3);
    sum_fit_qed->SetParLimits(4, 0, 1e5);
    sum_fit_qed->SetParLimits(5, 1e-6, 1e3);
    sum_fit_qed->SetParLimits(6, 1e-6, 1e3);
    qed_tail->Fit(sum_fit_qed, "QLR");

    double qed_tail_3 = sum_fit_qed->GetParameter(3);
    double qed_tail_4 = sum_fit_qed->GetParameter(4);
    double qed_tail_5 = sum_fit_qed->GetParameter(5);
    double qed_tail_6 = sum_fit_qed->GetParameter(6);

    sum_fit_dis->SetParameters(dis_tail->GetBinContent(1), 1, 1, 1, 1, 1, 1);
    sum_fit_dis->SetParLimits(0, 1e2, 1e5);
    sum_fit_dis->SetParLimits(1, 0, 1e3);
    sum_fit_dis->SetParLimits(2, 1e-6, 1e3);
    sum_fit_dis->SetParLimits(3, 0, 1e3);
    sum_fit_dis->SetParLimits(4, 0, 1e5);
    sum_fit_dis->SetParLimits(5, 1e-6, 1e3);
    sum_fit_dis->SetParLimits(6, 1e-6, 1e3);
    dis_tail->Fit(sum_fit_dis, "QLR");

    double dis_tail_3 = sum_fit_dis->GetParameter(3);
    double dis_tail_4 = sum_fit_dis->GetParameter(4);
    double dis_tail_5 = sum_fit_dis->GetParameter(5);
    double dis_tail_6 = sum_fit_dis->GetParameter(6);


//############################### acop ###############################
    TCanvas *acop1 = new TCanvas("acop1", "acop 1", 800, 800);
    acop1->Divide(3,3);
    for (int i = 0; i < neutron1; i++){ // for clarity, only plot half of them on each
        // this is our template fit. We assume that the tail is composed of mostly-QED effects and mostly-dissociative effects
        TF1 *template_fit = new TF1(Form("template_fit1_%d", i), Form("[0]*exp(-x/[1]+[2]*pow(x,0.25)) + [3]*pow((1+x*%f/%f), -%f) + [4]*pow((1+x*%f/%f), -%f)", qed_tail_4, qed_tail_5, qed_tail_6, dis_tail_4, dis_tail_5, dis_tail_6), 0.5*10e-5, 1);

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
        acop_hist[i]->SetMaximum(5000.);

        acop_hist[i]->SetMarkerStyle(20); acop_hist[i]->SetMarkerSize(0.6); acop_hist[i]->SetMarkerColor(kBlack); acop_hist[i]->SetLineColor(kBlack);
        acop_hist_c[i]->SetMarkerStyle(20); acop_hist_c[i]->SetMarkerSize(0.6); acop_hist_c[i]->SetMarkerColor(kBlack); acop_hist_c[i]->SetLineColor(kBlack);
        
        template_fit->SetParameters(acop_hist[i]->GetBinContent(1), 1, 1, qed_tail_3, dis_tail_3);
        //template_fit->SetParLimits(0, 1e2, 1e5);
        //template_fit->SetParLimits(1, 1e-5, 1);
        //template_fit->SetParLimits(2, 1e-6, 1e3);
        //template_fit->SetParLimits(3, 1e-4, 1e2);
        //template_fit->SetParLimits(4, 1e-4, 1e2);

        TFitResultPtr acop_fit_result = acop_hist[i]->Fit(template_fit, "QLRS");

        double s_1 = template_fit->GetParameter(0);
        double s_2 = template_fit->GetParameter(1);
        double s_3 = template_fit->GetParameter(2);
        double s_4 = template_fit->GetParameter(3);
        double s_5 = template_fit->GetParameter(4);

        core_1[i] = s_1;
        core_2[i] = s_2;
        core_3[i] = s_3;

        core_1_err[i] = template_fit->GetParError(0);
        core_2_err[i] = template_fit->GetParError(1);
        core_3_err[i] = template_fit->GetParError(2);

        TF1 *alpha_x_int = new TF1(Form("alpha_int_%d", i), "x * [0]*exp(-x/[1]+[2]*pow(x,0.25))", 0.5*10e-5, 1); // get alpha expectation values for core
        TF1 *alpha_int_norm = new TF1(Form("alpha_int_%d", i), "[0]*exp(-x/[1]+[2]*pow(x,0.25))", 0.5*10e-5, 1);
        alpha_x_int->FixParameter(0, s_1); alpha_x_int->FixParameter(1, s_2); alpha_x_int->FixParameter(2, s_3);
        alpha_int_norm->FixParameter(0, s_1); alpha_int_norm->FixParameter(1, s_2); alpha_int_norm->FixParameter(2, s_3);
        double alpha_expect = alpha_x_int->Integral(0.5*10e-5, 1) / alpha_int_norm->Integral(0.5*10e-5, 1);

        // ************************* Covariant Error Calculation *************************
        TMatrixDSym cov = acop_fit_result->GetCovarianceMatrix();
        double grad[3];
        grad[0] = 0.0;      // normalization cancels exactly

        for (int k = 1; k < 3; k++) {
            double p  = template_fit->GetParameter(k);
            double dp = template_fit->GetParError(k);
            template_fit->SetParameter(k, p + dp);
                for (int j = 0; j < 3; j++) {
                    alpha_x_int->SetParameter(j, template_fit->GetParameter(j));
                    alpha_int_norm->SetParameter(j, template_fit->GetParameter(j));
                }
            double plus = alpha_x_int->Integral(0.5e-5,1) / alpha_int_norm->Integral(0.5e-5,1);

            template_fit->SetParameter(k, p - dp);
                for (int j = 0; j < 3; j++) {
                    alpha_x_int->SetParameter(j, template_fit->GetParameter(j));
                    alpha_int_norm->SetParameter(j, template_fit->GetParameter(j));
                }
            double minus = alpha_x_int->Integral(0.5e-5,1) / alpha_int_norm->Integral(0.5e-5,1);
            template_fit->SetParameter(k, p);
            grad[k] = (plus - minus)/(2.0*dp);
        }
        double err2 = 0.0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                err2 += grad[a] * cov(a,b) * grad[b];
        double alpha_expect_err = std::sqrt(err2);
        // *************************  *************************

        core_expectation->SetBinContent(i+1, alpha_expect);
        core_expectation->SetBinError(i+1, alpha_expect_err);


        TF1 *sum_funct1 = new TF1(Form("sum_funct1_%d", i), Form("%f * exp(-x / %f + %f * pow(x, 0.25)) +%f*pow((1 + x* %f/%f), -%f) + %f*pow((1 + x* %f/%f), -%f)", s_1, s_2, s_3, s_4, qed_tail_4, qed_tail_5, qed_tail_6, s_5, dis_tail_4, dis_tail_5, dis_tail_6), 0.5*10e-5, 1);
        TF1 *core_funct = new TF1(Form("tcore_funct1_%d", i), Form("%f*exp(-x/%f+%f*pow(x,0.25))", s_1, s_2, s_3), 0.5*10e-5, 1);
        TF1 *tail_funct1 = new TF1(Form("tail_funct1_1_%d", i), Form("%f*pow((1+x*%f/%f), -%f)", s_4, qed_tail_4, qed_tail_5, qed_tail_6), 0.5*10e-5, 1);
        TF1 *tail_funct2 = new TF1(Form("tail_funct1_2_%d", i), Form("%f*pow((1+x*%f/%f), -%f)", s_5, dis_tail_4, dis_tail_5, dis_tail_6), 0.5*10e-5, 1);


        core_funct->SetLineColor(kRed); core_funct->SetLineWidth(3); core_funct->SetLineStyle(10);
        tail_funct1->SetLineColor(kBlue); tail_funct1->SetLineWidth(3); tail_funct1->SetLineStyle(2);
        tail_funct2->SetLineColor(kGreen); tail_funct2->SetLineWidth(3); tail_funct2->SetLineStyle(2);

        sum_funct1->SetLineColor(kMagenta);

        acop_hist[i]->GetXaxis()->SetTitle("#alpha");
        acop_hist[i]->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
        acop_hist[i]->GetXaxis()->SetTitleOffset(1.4);
        acop_hist[i]->GetYaxis()->SetTitleOffset(1.4);


        acop_hist[i]->Draw("ELP");
        core_funct->Draw("SAME");
        tail_funct1->Draw("SAME");
        tail_funct2->Draw("SAME");
        sum_funct1->Draw("SAME");

        gPad->SetLogx();
        gPad->SetLogy();

        if (i == 0){
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
            leg1->AddEntry(core_funct, "Core", "l");
            leg1->AddEntry(tail_funct1, "QED Tail", "l");
            leg1->AddEntry(tail_funct2, "Dis Tail", "l");
            leg1->AddEntry(sum_funct1, "Sum", "l");
            leg1->Draw();
        }
        TLatex latex1;
        latex1.SetNDC(); // use normalized coordinates (0–1)
        latex1.SetTextSize(0.05);
        latex1.SetTextAlign(22); // center horizontally
        latex1.DrawLatex(0.8, 0.8, Form("%s", bin_names[i]));


    }
    if (save_plots == true){
        acop1->SaveAs("plots/full/acop1_6.30.png");
    }
    TCanvas *acop2 = new TCanvas("acop2", "acop 2", 1200, 1000);
    acop2->Divide(3,2);
    for (int i = neutron1; i < neutron_bin; i++){ // second set of plots
        // this is our template fit. We assume that the tail is composed of mostly-QED effects and mostly-dissociative effects
        TF1 *template_fit = new TF1(Form("template_fit2_%d", i), Form("[0]*exp(-x/[1]+[2]*pow(x,0.25)) + [3]*pow((1+x*%f/%f), -%f) + [4]*pow((1+x*%f/%f), -%f)", qed_tail_4, qed_tail_5, qed_tail_6, dis_tail_4, dis_tail_5, dis_tail_6), 0.5*10e-5, 1);

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

        acop2->cd(i-(neutron1-1));

        gStyle->SetOptTitle(0);

        acop_hist[i]->SetLineWidth(2);
        acop_hist_c[i]->SetLineWidth(2);
        acop_hist[i]->SetMinimum(10e-4);
        acop_hist[i]->SetMaximum(5000.);

        acop_hist[i]->SetMarkerStyle(20); acop_hist[i]->SetMarkerSize(0.6); acop_hist[i]->SetMarkerColor(kBlack); acop_hist[i]->SetLineColor(kBlack);
        acop_hist_c[i]->SetMarkerStyle(20); acop_hist_c[i]->SetMarkerSize(0.6); acop_hist_c[i]->SetMarkerColor(kBlack); acop_hist_c[i]->SetLineColor(kBlack);

        template_fit->SetParameters(acop_hist[i]->GetBinContent(1), 1, 1, qed_tail_3, dis_tail_3);
        template_fit->SetParLimits(0, 1e2, 1e5);
        template_fit->SetParLimits(1, 1e-5, 1);
        template_fit->SetParLimits(2, 1e-6, 1e3);
        template_fit->SetParLimits(3, 1e-4, 5e1);
        template_fit->SetParLimits(4, 1e-4, 5e1);
        // p3 and p4 gave issues because of how small p0 is in ZnZn. Notice that acop1 and acop2 plots have diff limits bc of this

        TFitResultPtr acop_fit_result = acop_hist[i]->Fit(template_fit, "QLRS");

        double s_1 = template_fit->GetParameter(0);
        double s_2 = template_fit->GetParameter(1);
        double s_3 = template_fit->GetParameter(2);
        double s_4 = template_fit->GetParameter(3);
        double s_5 = template_fit->GetParameter(4);

        core_1[i] = s_1;
        core_2[i] = s_2;
        core_3[i] = s_3;

        core_1_err[i] = template_fit->GetParError(0);
        core_2_err[i] = template_fit->GetParError(1);
        core_3_err[i] = template_fit->GetParError(2);


        TF1 *alpha_x_int = new TF1("alpha_int", "x * [0]*exp(-x/[1]+[2]*pow(x,0.25))", 0.5*10e-5, 1); // get alpha expectation values for core
        TF1 *alpha_int_norm = new TF1("alpha_int", "[0]*exp(-x/[1]+[2]*pow(x,0.25))", 0.5*10e-5, 1);
        alpha_x_int->FixParameter(0, s_1); alpha_x_int->FixParameter(1, s_2); alpha_x_int->FixParameter(2, s_3);
        alpha_int_norm->FixParameter(0, s_1); alpha_int_norm->FixParameter(1, s_2); alpha_int_norm->FixParameter(2, s_3);
        double alpha_expect = alpha_x_int->Integral(0.5*10e-5, 1) / alpha_int_norm->Integral(0.5*10e-5, 1);

        // ************************* Covariant Error Calculation *************************
        TMatrixDSym cov = acop_fit_result->GetCovarianceMatrix();
        double grad[3];
        grad[0] = 0.0;      //normalization cancels exactly

        for (int k = 1; k < 3; k++) {
            double p  = template_fit->GetParameter(k);
            double dp = template_fit->GetParError(k);
            template_fit->SetParameter(k, p + dp);
                for (int j = 0; j < 3; j++) {
                    alpha_x_int->SetParameter(j, template_fit->GetParameter(j));
                    alpha_int_norm->SetParameter(j, template_fit->GetParameter(j));
                }
            double plus = alpha_x_int->Integral(0.5e-5,1) / alpha_int_norm->Integral(0.5e-5,1);

            template_fit->SetParameter(k, p - dp);
                for (int j = 0; j < 3; j++) {
                    alpha_x_int->SetParameter(j, template_fit->GetParameter(j));
                    alpha_int_norm->SetParameter(j, template_fit->GetParameter(j));
                }
            double minus = alpha_x_int->Integral(0.5e-5,1) / alpha_int_norm->Integral(0.5e-5,1);
            template_fit->SetParameter(k, p);
            grad[k] = (plus - minus)/(2.0*dp);
        }
        double err2 = 0.0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                err2 += grad[a] * cov(a,b) * grad[b];
        double alpha_expect_err = std::sqrt(err2);
        // *************************  *************************
        core_expectation->SetBinContent(i+1, alpha_expect);
        core_expectation->SetBinError(i+1, alpha_expect_err);


        TF1 *sum_funct1 = new TF1(Form("sum_funct2_%d", i), Form("%f * exp(-x / %f + %f * pow(x, 0.25)) +%f*pow((1 + x* %f/%f), -%f) + %f*pow((1 + x* %f/%f), -%f)", s_1, s_2, s_3, s_4, qed_tail_4, qed_tail_5, qed_tail_6, s_5, dis_tail_4, dis_tail_5, dis_tail_6), 0.5*10e-5, 1);
        TF1 *core_funct = new TF1(Form("tcore_funct2_%d", i), Form("%f*exp(-x/%f+%f*pow(x,0.25))", s_1, s_2, s_3), 0.5*10e-5, 1);
        TF1 *tail_funct1 = new TF1(Form("tail_funct2_1_%d", i), Form("%f*pow((1+x*%f/%f), -%f)", s_4, qed_tail_4, qed_tail_5, qed_tail_6), 0.5*10e-5, 1);
        TF1 *tail_funct2 = new TF1(Form("tail_funct2_2_%d", i), Form("%f*pow((1+x*%f/%f), -%f)", s_5, dis_tail_4, dis_tail_5, dis_tail_6), 0.5*10e-5, 1);


        core_funct->SetLineColor(kRed); core_funct->SetLineWidth(3); core_funct->SetLineStyle(10);
        tail_funct1->SetLineColor(kBlue); tail_funct1->SetLineWidth(3); tail_funct1->SetLineStyle(2);
        tail_funct2->SetLineColor(kGreen); tail_funct2->SetLineWidth(3); tail_funct2->SetLineStyle(2);

        sum_funct1->SetLineColor(kMagenta);

        acop_hist[i]->GetXaxis()->SetTitle("#alpha");
        acop_hist[i]->GetYaxis()->SetTitle("(1/N) dN/d#alpha");
        acop_hist[i]->GetXaxis()->SetTitleOffset(1.4);
        acop_hist[i]->GetYaxis()->SetTitleOffset(1.4);


        acop_hist[i]->Draw("ELP");
        core_funct->Draw("SAME");
        tail_funct1->Draw("SAME");
        tail_funct2->Draw("SAME");
        sum_funct1->Draw("SAME");

        gPad->SetLogx();
        gPad->SetLogy();

        TLatex latex1;
        latex1.SetNDC(); // use normalized coordinates (0–1)
        latex1.SetTextSize(0.05);
        latex1.SetTextAlign(22); // center horizontally
        latex1.DrawLatex(0.8, 0.8, Form("%s", bin_names[i]));
    }
    if (save_plots == true){
        acop2->SaveAs("plots/full/acop2_6.30.png");   
    }

//############################### RMS ###############################
/*    TH1D* alpha_rms = new TH1D("alpha_rms",";Neutron bin;RMS",neutron_bin, 0.5, neutron_bin + 0.5);

    vector<double> std(neutron_bin);
    vector<double> std_err(neutron_bin);
    vector<double> mean(neutron_bin);
    vector<double> mean_err(neutron_bin);
    vector<double> rms(neutron_bin);
    vector<double> rms_err_1p(neutron_bin);
    vector<double> rms_err_1m(neutron_bin);
    vector<double> rms_err_2p(neutron_bin);
    vector<double> rms_err_2m(neutron_bin);
    vector<double> rms_err_3p(neutron_bin);
    vector<double> rms_err_3m(neutron_bin);
    vector<double> rms_err(neutron_bin);

    vector<double> x_expect(neutron_bin);

    int c; //for indexing
    for(c = 0; c < neutron_bin; c++){ //x^2 expecation value of the TF1 divide by just integral of TF1 (do for both x^2 and x)
        std[c]      = acop_hist[c]->GetRMS(); //actually gives standard deviation
        std_err[c]  = acop_hist[c]->GetRMSError();
        mean[c]     = acop_hist[c]->GetMean();
        mean_err[c] = acop_hist[c]->GetMeanError();

        //diagonalize covar matrix, get eigenvalues/vectors

        //<x^2> given as integral of core(x)*x^2 / integral of core(x) w.r.t. x. RMS = sqrt(<x^2>)
        TF1 *x_base = new TF1(Form("x_base_%d", c), "[0] * exp(-x / [1] + [2] * pow(x, 0.25))", 0.5*10e-5, 1);
        TF1 *x_sqr = new TF1(Form("x_sqr_%d", c), "x*x * [0] * exp(-x / [1] + [2] * pow(x, 0.25))", 0.5*10e-5, 1);

        x_base->SetParameters(core_1[c], core_2[c], core_3[c]);
        x_sqr->SetParameters(core_1[c], core_2[c], core_3[c]);
        x_expect[c] = x_sqr->Integral(0.5*10e-5, 1) / x_base->Integral(0.5*10e-5, 1);
        rms[c] = sqrt(x_expect[c]);

        // begin varying parameters
        double del[3] = {0.0, 0.0, 0.0}; //for doing our delta terms        
        for (int j = 0; j < 3; j++){
            double p[3]  = {core_1[c], core_2[c], core_3[c]};
            double dp[3] = {core_1_err[c], core_2_err[c], core_3_err[c]};

            //pos var
            p[j] += dp[j];
            x_base->SetParameters(p);
            x_sqr->SetParameters(p);
            double x_sqr_p = x_sqr->Integral(0.5*10e-5, 1) / x_base->Integral(0.5*10e-5, 1) - x_expect[c];

            //neg var
            p[j] -= 2.0 * dp[j]; //first to subtract off previous action, second to give result
            x_base->SetParameters(p);
            x_sqr->SetParameters(p);
            double x_sqr_m = x_sqr->Integral(0.5*10e-5, 1) / x_base->Integral(0.5*10e-5, 1) - x_expect[c];

            del[j] = 0.5 * (abs(x_sqr_p) + abs(x_sqr_m));
        }

        double x_sqr_err = sqrt(pow(del[0], 2) + pow(del[1], 2) + pow(del[2], 2));
        rms_err[c] = 0.5 * x_sqr_err / rms[c]; //this is from standard error prop

        alpha_rms->SetBinContent(c+1, rms[c]);
        alpha_rms->SetBinError  (c+1, rms_err[c]);
        alpha_rms->GetXaxis()->SetBinLabel(c+1, bin_names[c]);
    }
    TCanvas *rms_plots = new TCanvas("rms", "alpha RMS", 1200, 1000);
    alpha_rms->SetStats(0);
    alpha_rms->SetFillColor(kCyan + 1);
    alpha_rms->Draw("E1 BAR");

    if (save_plots == true){
        rms_plots->SaveAs("plots/full/alpha_rms1_6.30.png");
    }
*/
/*    TMatrixDSym cov;
    TVectorD eigenValues;
    TMatrixD eigenVectors;
    covar_evectors(core_1, core_2, core_3, cov, eigenValues, eigenVectors);
    
    TH1D* alpha_rms = new TH1D("alpha_rms",";Neutron bin;RMS",neutron_bin, 0.5, neutron_bin + 0.5);

    for(int i = 0; i < neutron_bin; i++){
        double s1 = core_1[i];
        double s2 = core_2[i];
        double s3 = core_3[i];

        double rms = sqrt(rms2(s1, s2, s3)); //Soumya gave that rms^2 = <x^2> / <x>

        double rms_err = rms_err_f(s1, s2, s3, cov);

        alpha_rms->SetBinContent(i + 1, rms);
        alpha_rms->SetBinError(i + 1, 0.5 * rms_err / rms);
        alpha_rms->GetXaxis()->SetBinLabel(i + 1, bin_names[i]);
    }

    TCanvas *rms_plots = new TCanvas("rms", "alpha RMS", 1200, 1000);
    gStyle->SetOptStat(0);

    alpha_rms->SetFillColor(kCyan + 1);

    alpha_rms->Draw("E1 BAR");
    if (save_plots == true){
        rms_plots->SaveAs("plots/full/alpha_rms1_6.30.png");
    }
*/

//############################### Core expectation value (plus approx CMS result) ###############################
    TCanvas *a_exp = new TCanvas("a_exp", "alpha expectation values", 1200, 1000);
    core_expectation->Draw("P");

    TMarker *onon_point = new TMarker(0.5, 1.16e-3, 8);   onon_point->SetMarkerColor(kRed); onon_point->Draw("SAME");
    TMarker *onln_point = new TMarker(1.5, 1.28e-3, 8);   onln_point->SetMarkerColor(kRed); onln_point->Draw("SAME");
    TMarker *onxn_point = new TMarker(3.5, 1.3e-3, 8);    onxn_point->SetMarkerColor(kRed); onxn_point->Draw("SAME");
    TMarker *lnln_point = new TMarker(6.5, 1.38e-3, 8);   lnln_point->SetMarkerColor(kRed); lnln_point->Draw("SAME");
    TMarker *lnxn_point = new TMarker(7.5, 1.425e-3, 8);  lnxn_point->SetMarkerColor(kRed); lnxn_point->Draw("SAME");
    TMarker *xnxn_point = new TMarker(10.5, 1.485e-3, 8); xnxn_point->SetMarkerColor(kRed); xnxn_point->Draw("SAME");

    if (save_plots == true){
        a_exp->SaveAs("plots/full/alpha_expectation_6.30.png");
    }
}