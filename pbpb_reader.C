#define pbpb_reader_cxx
#include "pbpb_reader.h"
#include "log_bin.cxx"
#include "TH1D.h"
#include "Riostream.h"
#include <math.h>

pbpb_reader::pbpb_reader(){
    InitTree();
    TFile *output_file = new TFile("dimuon_hists.root", "recreate");
    InitHists();

    Long64_t nentries = fChain->GetEntries(); //runs over all events
    cout << nentries << " original events" << endl;
    int kept = 0;
    int count = 0;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        if (ch1 == ch2 || mass < 8 || mass > 60 || abs(asym) > 0.06 || abs(rapidity) > 2.4 || abs(eta1) > 2.4 || abs(eta2) > 2.4 || num_trk_good != 0 || abs(d01) > 0.1 || abs(d02) > 0.1 || abs(z01) > 0.2 || abs(z02) > 0.2) continue;
        //NOTE (1/25): all num_trk_good == 0 from upstream. None that do not equal strictly 0

        h_n_raw_a->Fill(NeutronsA);
        h_n_raw_c->Fill(NeutronsC);

        int neutronsA_count = round(NeutronsA);
        int neutronsC_count = round(NeutronsC);

        if (neutronsA_count > 10 || neutronsC_count > 10) continue;
        kept++;

        h_NeutronsA->Fill(neutronsA_count);
        h_NeutronsC->Fill(neutronsC_count);

        //wrap angles back from 2pi
        if (dphi> (1.5)*(M_PI)) {
            dphi-=2*M_PI;
        }
        else if(dphi<-(0.5)*(M_PI)) {
            dphi+=2*M_PI;
        }

//sort results by number of outgoing neutrons. They are all "if" as they are NOT mutually exclusive conditions
        //NOTE: These were updated to no longer seperate Nn0n cases, so acop bins aren't correct (6/4/2026)
    if ((neutronsA_count == 0 && neutronsC_count == 0)){ //(0,0)
        h_a_n[0]     ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[0]   ->Fill(kperp);
        h_k_lin[0]   ->Fill(kperp);
        h_mass[0]    ->Fill(mass);
        h_rapidity[0]->Fill(rapidity);
        h_acop[0]    ->Fill(acop);
        neutron_stat[0] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[0]->Fill(acop);
            h_a_n_s[0] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[0] ->Fill(kperp);
            h_mass_s[0]    ->Fill(mass);
            h_rapidity_s[0]->Fill(rapidity);
        }
    }
    if ((neutronsA_count == 1 && neutronsC_count == 0) || (neutronsA_count == 0 && neutronsC_count == 1)){ //(1,0)
        h_a_n[1]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[1]   ->Fill(kperp);
        h_k_lin[1]   ->Fill(kperp);
        h_mass[1]    ->Fill(mass);
        h_rapidity[1]->Fill(rapidity);
        h_acop[1]    ->Fill(acop);
        neutron_stat[1] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[1]->Fill(acop);
            h_a_n_s[1] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[1] ->Fill(kperp);
            h_mass_s[1]    ->Fill(mass);
            h_rapidity_s[1]->Fill(rapidity);
        }
        if(neutronsA_count == 0){
            h_rapid_sgn[1]->Fill(rapidity);
            h_kperp_prapid[1]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[1]->Fill(rapidity);
            }
            if(kperp < 150){
                h_rapid_b150[1]->Fill(rapidity);
            }
        }
        if(neutronsC_count == 0){
            h_rapid_sgn[1]->Fill(-rapidity);
            h_kperp_nrapid[1]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[1]->Fill(-rapidity);
            }
            if(kperp < 150){
                h_rapid_b150[1]->Fill(-rapidity);
            }
        }
    }
    if ((neutronsA_count == 2 && neutronsC_count == 0) || (neutronsA_count == 0 && neutronsC_count == 2)){ //(2,0)
        h_a_n[2]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[2]   ->Fill(kperp);
        h_k_lin[2]   ->Fill(kperp);
        h_mass[2]    ->Fill(mass);
        h_rapidity[2]->Fill(rapidity);
        h_acop[2]    ->Fill(acop);
        neutron_stat[2] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[2]->Fill(acop);
            h_a_n_s[2] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[2] ->Fill(kperp);
            h_mass_s[2]    ->Fill(mass);
            h_rapidity_s[2]->Fill(rapidity);
        } 
        if(neutronsA_count == 0){
            h_rapid_sgn[2]->Fill(rapidity);
            h_kperp_prapid[2]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[2]->Fill(rapidity);
            }
            if(kperp < 150){
                h_rapid_b150[2]->Fill(rapidity);
            }
        }
        if(neutronsC_count == 0){
            h_rapid_sgn[2]->Fill(-rapidity);
            h_kperp_nrapid[2]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[2]->Fill(-rapidity);
            }
            if(kperp < 150){
                h_rapid_b150[2]->Fill(-rapidity);
            }
        }
    }
    if ((neutronsA_count == 0 && neutronsC_count > 0) || (neutronsA_count > 0 && neutronsC_count == 0)){ //(X,0)
        h_a_n[3]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[3]   ->Fill(kperp);
        h_k_lin[3]   ->Fill(kperp);
        h_mass[3]    ->Fill(mass);
        h_rapidity[3]->Fill(rapidity);
        h_acop[3]    ->Fill(acop);
        neutron_stat[3] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[3]->Fill(acop);
            h_a_n_s[3] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[3] ->Fill(kperp);
            h_mass_s[3]    ->Fill(mass);
            h_rapidity_s[3]->Fill(rapidity);
        }  
        if(neutronsA_count == 0){
            h_rapid_sgn[3]->Fill(rapidity);
            h_kperp_prapid[3]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[3]->Fill(rapidity);
            } 
            if(kperp < 150){
                h_rapid_b150[3]->Fill(rapidity);
            }
        }
        if(neutronsC_count == 0){
            h_rapid_sgn[3]->Fill(-rapidity);
            h_kperp_nrapid[3]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[3]->Fill(-rapidity);
            }
            if(kperp < 150){
                h_rapid_b150[3]->Fill(-rapidity);
            }
        }
    }
    if ((neutronsA_count > 1 && neutronsC_count == 0) || (neutronsA_count == 0 && neutronsC_count > 1)  && (neutronsA_count <= 5) && (neutronsC_count <= 5)){ //(Y,0)
        h_a_n[4]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[4]   ->Fill(kperp);
        h_k_lin[4]   ->Fill(kperp);
        h_mass[4]    ->Fill(mass);
        h_rapidity[4]->Fill(rapidity);
        h_acop[4]    ->Fill(acop);
        neutron_stat[4] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[4]->Fill(acop);
            h_a_n_s[4] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[4] ->Fill(kperp);
            h_mass_s[4]    ->Fill(mass);
            h_rapidity_s[4]->Fill(rapidity);
        }
        if(neutronsA_count == 0){
            h_rapid_sgn[4]->Fill(rapidity);
            h_kperp_prapid[4]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[4]->Fill(rapidity);
            } 
            if(kperp < 150){
                h_rapid_b150[4]->Fill(rapidity);
            }
        }
        if(neutronsC_count == 0){
            h_rapid_sgn[4]->Fill(-rapidity);
            h_kperp_nrapid[4]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[4]->Fill(-rapidity);
            }
            if(kperp < 150){
                h_rapid_b150[4]->Fill(-rapidity);
            }
        }
    }
    if ((neutronsA_count > 5 && neutronsC_count == 0) || (neutronsA_count == 0 && neutronsC_count > 5)){ //(Z,0)
        h_a_n[5]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[5]   ->Fill(kperp);
        h_k_lin[5]   ->Fill(kperp);
        h_mass[5]    ->Fill(mass);
        h_rapidity[5]->Fill(rapidity);
        h_acop[5]    ->Fill(acop);
        neutron_stat[5] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[5]->Fill(acop);
            h_a_n_s[5] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[5] ->Fill(kperp);
            h_mass_s[5]    ->Fill(mass);
            h_rapidity_s[5]->Fill(rapidity);
        } 
        if(neutronsA_count == 0){
            h_rapid_sgn[5]->Fill(rapidity);
            h_kperp_prapid[5]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[5]->Fill(rapidity);
            } 
            if(kperp < 150){
                h_rapid_b150[5]->Fill(rapidity);
            }
        }
        if(neutronsC_count == 0){
            h_rapid_sgn[5]->Fill(-rapidity);
            h_kperp_nrapid[5]->Fill(kperp);
            if(kperp > 150){
                h_rapid_a150[5]->Fill(-rapidity);
            }
            if(kperp < 150){
                h_rapid_b150[5]->Fill(-rapidity);
            }
        }
    }
    if ((neutronsA_count == 1 && neutronsC_count == 1)){ //(1,1)
        h_a_n[6]     ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[6]   ->Fill(kperp);
        h_k_lin[6]   ->Fill(kperp);
        h_mass[6]    ->Fill(mass);
        h_rapidity[6]->Fill(rapidity);
        h_acop[6]    ->Fill(acop);
        neutron_stat[6] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[6]->Fill(acop);
            h_a_n_s[6] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[6] ->Fill(kperp);
            h_mass_s[6]    ->Fill(mass);
            h_rapidity_s[6]->Fill(rapidity);
        }
    }    
    if ((neutronsA_count > 0 && neutronsC_count == 1) || (neutronsA_count == 1 && neutronsC_count > 0)){ //(X,1)
        h_a_n[7]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[7]   ->Fill(kperp);
        h_k_lin[7]   ->Fill(kperp);
        h_mass[7]    ->Fill(mass);
        h_rapidity[7]->Fill(rapidity);
        h_acop[7]    ->Fill(acop);
        neutron_stat[7] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[7]->Fill(acop);
            h_a_n_s[7] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[7] ->Fill(kperp);
            h_mass_s[7]    ->Fill(mass);
            h_rapidity_s[7]->Fill(rapidity);
        }
    }
    if ((neutronsA_count > 1 && neutronsC_count == 1) || (neutronsA_count == 1 && neutronsC_count > 1)  && (neutronsA_count <= 5) && (neutronsC_count <= 5)){ //(Y,1)
        h_a_n[8]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[8]   ->Fill(kperp);
        h_k_lin[8]   ->Fill(kperp);
        h_mass[8]    ->Fill(mass);
        h_rapidity[8]->Fill(rapidity);
        h_acop[8]    ->Fill(acop);
        neutron_stat[8] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[8]->Fill(acop);
            h_a_n_s[8] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[8] ->Fill(kperp);
            h_mass_s[8]    ->Fill(mass);
            h_rapidity_s[8]->Fill(rapidity);
        }
    }
    if ((neutronsA_count > 5 && neutronsC_count == 1) || (neutronsA_count == 1 && neutronsC_count > 5)){ //(Z,1)
        h_a_n[9]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[9]   ->Fill(kperp);
        h_k_lin[9]   ->Fill(kperp);
        h_mass[9]    ->Fill(mass);
        h_rapidity[9]->Fill(rapidity);
        h_acop[9]    ->Fill(acop);
        neutron_stat[9] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[9]->Fill(acop);
            h_a_n_s[9] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[9] ->Fill(kperp);
            h_mass_s[9]    ->Fill(mass);
            h_rapidity_s[9]->Fill(rapidity);
        }
    }
    if (neutronsA_count > 0 && neutronsC_count > 0){ //(X,X)
        h_a_n[10]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[10]   ->Fill(kperp);
        h_k_lin[10]   ->Fill(kperp);
        h_mass[10]    ->Fill(mass);
        h_rapidity[10]->Fill(rapidity);
        h_acop[10]    ->Fill(acop);
        neutron_stat[10] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[10]->Fill(acop);
            h_a_n_s[10] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[10] ->Fill(kperp);
            h_mass_s[10]    ->Fill(mass);
            h_rapidity_s[10]->Fill(rapidity);
        }
    }
    if (neutronsA_count > 1 && neutronsC_count > 1  && (neutronsA_count <= 5) && (neutronsC_count <= 5)){ //(Y,Y)
        h_a_n[11]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[11]   ->Fill(kperp);
        h_k_lin[11]   ->Fill(kperp);
        h_mass[11]    ->Fill(mass);
        h_rapidity[11]->Fill(rapidity);
        h_acop[11]    ->Fill(acop);
        neutron_stat[11] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[11]->Fill(acop);
            h_a_n_s[11] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[11] ->Fill(kperp);
            h_mass_s[11]    ->Fill(mass);
            h_rapidity_s[11]->Fill(rapidity);
        }
    }
    if (neutronsA_count > 5 && neutronsC_count > 5){ //(Z,Z)
        h_a_n[12]       ->Fill(acop, neutronsA_count + neutronsC_count);
        h_kperp[12]   ->Fill(kperp);
        h_k_lin[12]   ->Fill(kperp);
        h_mass[12]    ->Fill(mass);
        h_rapidity[12]->Fill(rapidity);
        h_acop[12]    ->Fill(acop);
        neutron_stat[12] ->Fill(neutronsA_count + neutronsC_count);
        if (abs(acop) < 0.008){
            h_acop_s[12]->Fill(acop);
            h_a_n_s[12] ->Fill(acop, neutronsA_count + neutronsC_count);
            h_kperp_s[12] ->Fill(kperp);
            h_mass_s[12]    ->Fill(mass);
            h_rapidity_s[12]->Fill(rapidity);
        }    
    }

    }

    output_file->cd();

    h_NeutronsA  ->Write();
    h_NeutronsC  ->Write();
    h_n_raw_a    ->Write();
    h_n_raw_c    ->Write();

    for(int i = 0; i<neutron_bin; i++){
        h_a_n[i]         ->Write();
        h_a_n_s[i]       ->Write();
        h_kperp[i]       ->Write();
        h_k_lin[i]       ->Write();
        h_kperp_s[i]     ->Write();
        h_mass_s[i]      ->Write();
        h_rapidity_s[i]  ->Write();

        h_acop[i]        ->Write();
        h_asym[i]        ->Write();
        h_dphi[i]        ->Write();
        h_ptbar[i]       ->Write();
        h_mass[i]        ->Write();
        h_rapidity[i]    ->Write();
        h_Egamma_a[i]    ->Write();
        h_Egamma_b[i]    ->Write();
        h_pt1[i]         ->Write();
        h_pt2[i]         ->Write();
        h_eta1[i]        ->Write();
        h_eta2[i]        ->Write();
        h_phi1[i]        ->Write();
        h_phi2[i]        ->Write();
        h_d01[i]         ->Write();
        h_d02[i]         ->Write();
        h_z01[i]         ->Write();
        h_z02[i]         ->Write();
        h_Zvtx[i]        ->Write();
        h_qual1[i]       ->Write();
        h_qual2[i]       ->Write();
        h_rapid_sgn[i]   ->Write();
        h_rapid_a150[i]  ->Write();
        h_rapid_b150[i]  ->Write();
        neutron_stat[i]  ->Write();
        h_kperp_prapid[i]->Write();
        h_kperp_nrapid[i]->Write();
    }
    output_file->Close();
    cout << kept << " by end" << endl;
}


void pbpb_reader::InitHists() { //not running this will cause a segmentation error
    vector<double> acop_bins1;
    vector<double> acop_s_bins1;
    vector<double> acop_bins2;
    vector<double> acop_s_bins2;
    vector<double> acop_bins3;
    vector<double> acop_s_bins3;
    vector<double> acop_bins4;
    vector<double> acop_s_bins4;
    vector<double> acop_bins5;
    vector<double> acop_s_bins5;
    vector<double> acop_bins6;
    vector<double> acop_s_bins6;
    vector<double> acop_bins7;
    vector<double> acop_s_bins7;
    vector<double> acop_bins8;
    vector<double> acop_s_bins8;
    vector<double> acop_bins9;
    vector<double> acop_s_bins9;
    vector<double> acop_bins10;
    vector<double> acop_s_bins10;
    vector<double> acop_bins11;
    vector<double> acop_s_bins11;
    vector<double> acop_bins12;
    vector<double> acop_s_bins12;
    vector<double> acop_bins13;

    vector<double> k_log_bins;
    vector<double> kperp_bins;
    fillLogBinningArray(kperp_bins, 80, 100.0, 10000.0);

    //function takes (std::vector<double>& bins, int nBins, double low, double high)
    fillLogBinningArray(acop_bins1,  100, 0.5*10e-5, 1); //0n0n
    fillLogBinningArray(acop_bins2,  100, 0.5*10e-5, 1); //1n0n
    fillLogBinningArray(acop_bins3,  85, 0.5*10e-5, 1);  //1n1n
    fillLogBinningArray(acop_bins4,  100, 0.5*10e-5, 1); //Xn0n
    fillLogBinningArray(acop_bins5,  80, 0.5*10e-5, 1);  //Xn1n
    fillLogBinningArray(acop_bins6,  100, 0.5*10e-5, 1); //XnXn
    fillLogBinningArray(acop_bins7,  100, 0.5*10e-5, 1); //Yn0n
    fillLogBinningArray(acop_bins8,  70, 0.5*10e-5, 1);  //Yn1n
    fillLogBinningArray(acop_bins9,  100, 0.5*10e-5, 1); //YnYn
    fillLogBinningArray(acop_bins10, 100, 0.5*10e-5, 1); //Zn0n
    fillLogBinningArray(acop_bins11, 100, 0.5*10e-5, 1); //Zn1n
    fillLogBinningArray(acop_bins12, 100, 0.5*10e-5, 1); //ZnZn
    fillLogBinningArray(acop_bins13, 100, 0.5*10e-5, 1); //


    //############################### kperp bin creation ###############################
    double Emin = 1;      // Min at 1 MeV since log cannot start at 0
    double ElogMax = 100.0;    // log until 100 MeV
    double Emax = 10000.0;      // 10 GeV max

    int nLogBins = 30;
    int nLinBins = 90;
    std::vector<double> k_bins;
    // --- Logarithmic region ---
    double logMin = std::log10(Emin);
    double logMax = std::log10(ElogMax);
    double logStep = (logMax - logMin) / nLogBins;
    for (int i = 0; i <= nLogBins; i++) {
        k_bins.push_back(std::pow(10, logMin + i * logStep));
    }
    // --- Linear region ---
    double linStep = (Emax - ElogMax) / nLinBins;
    for (int i = 1; i <= nLinBins; i++) {
        k_bins.push_back(ElogMax + i * linStep);
    }

    // --- Again, but all log ---
    int nLogBins2 = 50;
    double logMin2 = log10(Emin);
    double logMax2 = log10(10000.0);
    double logStep2 = (logMax2 - logMin2) / nLogBins2;
    for (int i = 0; i <= nLogBins2; i++){
        k_log_bins.push_back(pow(10, logMin2 + i * logStep2));
    }

    h_NeutronsA= new TH1D("h_NeutronsA", "count;NeutronsA", 11, 0, 11);  
    h_NeutronsC= new TH1D("h_NeutronsC", "count;NeutronsB", 11, 0, 11);
    h_n_raw_a=   new TH1D("h_n_raw_a", "count;NeutronsA", 105, 0, 10.5);
    h_n_raw_c=   new TH1D("h_n_raw_c", "count;NeutronsC", 105, 0, 10.5);
    //neutron_stat = new TH1D("neutron_stat", "Neutrons by top", 100, 0, 10);

    h_acop[0]=    new TH1D("h_acop_n1", "acop n0n0", acop_bins1.size()-1, acop_bins1.data());
    h_acop_s[0]=  new TH1D("h_acop_s_n1", "acop n0n0", 80, 0, 0.008);
    h_acop[1]=    new TH1D("h_acop_n2", "acop n0n1", acop_bins2.size()-1, acop_bins2.data());
    h_acop_s[1]=  new TH1D("h_acop_s_n2", "acop n0n1", 80, 0, 0.008);
    h_acop[2]=    new TH1D("h_acop_n3", "acop n0nY", acop_bins3.size()-1, acop_bins3.data());
    h_acop_s[2]=  new TH1D("h_acop_s_n3", "acop n0nY", 80, 0, 0.008);
    h_acop[3]=    new TH1D("h_acop_n4", "acop n1n1", acop_bins4.size()-1, acop_bins4.data());
    h_acop_s[3]=  new TH1D("h_acop_s_n4", "acop n1n1", 80, 0, 0.008);
    h_acop[4]=    new TH1D("h_acop_n5", "acop n1nY", acop_bins5.size()-1, acop_bins5.data());
    h_acop_s[4]=  new TH1D("h_acop_s_n5", "acop n1nY", 80, 0, 0.008);
    h_acop[5]=    new TH1D("h_acop_n6", "acop nYnY", acop_bins6.size()-1, acop_bins6.data());
    h_acop_s[5]=  new TH1D("h_acop_s_n6", "acop nYnY", 80, 0, 0.008);
    h_acop[6]=    new TH1D("h_acop_n7", "acop n0nX", acop_bins7.size()-1, acop_bins7.data());
    h_acop_s[6]=  new TH1D("h_acop_s_n7", "acop n0nX", 80, 0, 0.008);
    h_acop[7]=    new TH1D("h_acop_n8", "acop n1nX", acop_bins8.size()-1, acop_bins8.data());
    h_acop_s[7]=  new TH1D("h_acop_s_n8", "acop n1nX", 80, 0, 0.008);
    h_acop[8]=    new TH1D("h_acop_n9", "acop nXnX", acop_bins9.size()-1, acop_bins9.data());
    h_acop_s[8]=  new TH1D("h_acop_s_n9", "acop nXnX", 80, 0, 0.008);
    h_acop[9]=    new TH1D("h_acop_n10", "acop nXnX", acop_bins10.size()-1, acop_bins10.data());
    h_acop_s[9]=  new TH1D("h_acop_s_n10", "acop nXnX", 80, 0, 0.008);
    h_acop[10]=   new TH1D("h_acop_n11", "acop nXnX", acop_bins11.size()-1, acop_bins11.data());
    h_acop_s[10]= new TH1D("h_acop_s_n11", "acop nXnX", 80, 0, 0.008);
    h_acop[11]=   new TH1D("h_acop_n12", "acop nXnX", acop_bins12.size()-1, acop_bins12.data());
    h_acop_s[11]= new TH1D("h_acop_s_n12", "acop nXnX", 80, 0, 0.008);
    h_acop[12]=   new TH1D("h_acop_n13", "acop nXnX", acop_bins13.size()-1, acop_bins13.data());
    h_acop_s[12]= new TH1D("h_acop_s_n13", "acop nXnX", 80, 0, 0.008);

    for(int i = 0; i < neutron_bin; i++){
        h_a_n[i]        = new TH2D(Form("h_a_n_%d", i), "acop;#alpha;neutrons;count", 100, 0, 8*10e-3, 10, 0, 10);
        h_a_n_s[i]      = new TH2D(Form("h_a_n_s_%d", i), "acop;#alpha;neutrons;count", 100, 0, 8*10e-3, 10, 0, 10);
        //h_kperp[i]      = new TH1D(Form("h_kperp_%d", i), "count;k_{#perp}", k_bins.size()-1, &k_bins[0]); //small log region
        //h_kperp[i]      = new TH1D(Form("h_kperp_%d", i), "count;k_{#perp}", 30, 0, 1e4); //all linear region
        h_kperp[i]      = new TH1D(Form("h_kperp_%d", i), "count;k_{#perp}", k_log_bins.size()-1, &k_log_bins[0]); //all log
        h_k_lin[i]      = new TH1D(Form("h_k_lin_%d", i), "k_{#perp};MeV", 5000, 0, 1e4);
        h_kperp_s[i]    = new TH1D(Form("h_kperp_s_%d", i), "count;k_{#perp}", 50, 0, 1e4);
        h_mass_s[i]     = new TH1D(Form("h_mass_s_%d", i), "count;mass", 13, 8, 60);  
        h_rapidity_s[i] = new TH1D(Form("h_rapidity_s_%d", i), "count;y", 24, -2.4, 2.4);
        
        h_rapid_sgn[i]    = new TH1D(Form("h_rapid_sgn_%d", i), "count;rapidity", 12, -2.4, 2.4);
        h_rapid_a150[i]   = new TH1D(Form("h_rapid_a150_%d", i), "count;rapidity", 12, -2.4, 2.4);
        h_rapid_b150[i]   = new TH1D(Form("h_rapid_b150_%d", i), "count;rapidity", 12, -2.4, 2.4);
        neutron_stat[i]   = new TH1D(Form("neutron_stat_%d", i), "Neutrons by top", 20, 0, 20);
        h_kperp_prapid[i] = new TH1D(Form("h_kperp_prapid_%d", i), "count;k_{#perp}", k_log_bins.size()-1, &k_log_bins[0]);
        h_kperp_nrapid[i] = new TH1D(Form("h_kperp_nrapid_%d", i), "count;k_{#perp}", k_log_bins.size()-1, &k_log_bins[0]);
        
        h_asym[i]     = new TH1D(Form("h_asym_%d", i), "count;asym", 100, 0, 1);    
        h_dphi[i]     = new TH1D(Form("h_dphi_%d", i), "count;dphi", 100, -0.5*M_PI, 1.5*M_PI);  
        h_ptbar[i]    = new TH1D(Form("h_ptbar_%d", i), "count;ptbar", 100, 0, 30);  
        h_mass[i]     = new TH1D(Form("h_mass_%d", i), "count;mass", 13, 8, 60);  
        h_rapidity[i] = new TH1D(Form("h_rapidity_%d", i), "count;y", 24, -2.4, 2.4);  
        h_Egamma_a[i] = new TH1D(Form("h_Egamma_a_%d", i), "count;Egamma_a", 100, 0, 150);  
        h_Egamma_b[i] = new TH1D(Form("h_Egamma_b_%d", i), "count;Egamma_b", 100, 0, 150);  
        h_pt1[i]      = new TH1D(Form("h_pt1_%d", i), "count;pt1", 100, 0, 100);  
        h_pt2[i]      = new TH1D(Form("h_pt2_%d", i), "count;pt2", 100, 0, 100);  
        h_eta1[i]     = new TH1D(Form("h_eta1_%d", i), "count;eta1", 100, -M_PI, M_PI);  
        h_eta2[i]     = new TH1D(Form("h_eta2_%d", i), "count;eta2", 100, -2.4, 2.4);  
        h_phi1[i]     = new TH1D(Form("h_phi1_%d", i), "count;phi1", 100, -M_PI, M_PI);  
        h_phi2[i]     = new TH1D(Form("h_phi2_%d", i), "count;phi2", 100, -M_PI, M_PI);  
        h_d01[i]      = new TH1D(Form("h_d01_%d", i), "count;d01", 100, -1, 1);  
        h_d02[i]      = new TH1D(Form("h_d02_%d", i), "count;d02", 100, -1, 1);
        h_z01[i]      = new TH1D(Form("h_z01_%d", i), "count;z01", 100, -200, 200);  
        h_z02[i]      = new TH1D(Form("h_z02_%d", i), "count;z02", 100, -200, 200);  
        h_Zvtx[i]     = new TH1D(Form("h_Zvtx_%d", i), "count;Zvtx", 100, -250, 250);
        h_qual1[i]    = new TH1D(Form("h_qual1_%d", i), "count;qual1", 100, 100, 600);
        h_qual2[i]    = new TH1D(Form("h_qual2_%d", i), "count;qual2", 100, 100, 600);
    }
}