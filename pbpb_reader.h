#ifndef pbpb_reader_h
#define pbpb_reader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include "vector"
#include <TH2D.h>
#include <TH3D.h>

class pbpb_reader {
    private:
    TChain            *fChain;

    UInt_t RunNumber = 0;
    static const int neutron_bin = 13; //(0,0), (1,1)
                                       //(X,1), (X,X)
                                       //(Y,1), (Y,Y)
                                       //(Z,1), (Z,Z)
                                       //(1,0), (2,0), (X,0), (Y,0), (Z,0)
    //X >= 1, 1<Y<=5, 5<5<=10. Seperate (n,0) cases

    //Decloration of leaf types
    float       acop = 0;
    float       asym = 0;
    float       kperp = 0;
    float       dphi = 0;
    float       ptbar = 0;
    float       mass = 0;   
    float       rapidity = 0;
    float       Egamma_a = 0;
    float       Egamma_b = 0;
    float       pt1 = 0;
    float       pt2 = 0;
    float       eta1 = 0;
    float       eta2 = 0;
    float       phi1 = 0;
    float       phi2 = 0;
    int         ch1 = 0;
    int         ch2 = 0;
    float       d01 = 0;
    float       d02 = 0;
    float       z01 = 0;
    float       z02 = 0;
    int         qual1 = 0;
    int         qual2 = 0;
    float       NeutronsA = 0;
    float       NeutronsC = 0;
    float       Zvtx = 0;
    int num_trk_good = 0;

    void InitTree ();//Initialize the branches in the Tree
    void InitHists();//Initialize the histograms

    public:
    //put any histograms here that will be used later
    TH1D* h_acop[neutron_bin]={nullptr};
    TH1D* h_acop_s[neutron_bin]={nullptr};

    TH1D* h_asym[neutron_bin]={nullptr};
    TH1D* h_dphi[neutron_bin]={nullptr};
    TH1D* h_ptbar[neutron_bin]={nullptr};
    TH1D* h_mass[neutron_bin]={nullptr};
    TH1D* h_rapidity[neutron_bin]={nullptr};
    TH1D* h_Egamma_a[neutron_bin]={nullptr};
    TH1D* h_Egamma_b[neutron_bin]={nullptr};
    TH1D* h_pt1[neutron_bin]  ={nullptr};
    TH1D* h_pt2[neutron_bin]  ={nullptr};
    TH1D* h_eta1[neutron_bin] ={nullptr};
    TH1D* h_eta2[neutron_bin] ={nullptr};
    TH1D* h_phi1[neutron_bin] ={nullptr};
    TH1D* h_phi2[neutron_bin] ={nullptr};
    TH1D* h_ch1[neutron_bin]  ={nullptr};
    TH1D* h_ch2[neutron_bin]  ={nullptr};
    TH1D* h_d01[neutron_bin]  ={nullptr};
    TH1D* h_d02[neutron_bin]  ={nullptr};
    TH1D* h_z01[neutron_bin]  ={nullptr};
    TH1D* h_z02[neutron_bin]  ={nullptr};
    TH1D* h_NeutronsA         ={nullptr};
    TH1D* h_NeutronsC         ={nullptr};
    TH1D* h_Zvtx[neutron_bin] ={nullptr};
    TH1D* h_qual1[neutron_bin]={nullptr};
    TH1D* h_qual2[neutron_bin]={nullptr};

    TH2D* h_a_n[neutron_bin]={nullptr};
    TH2D* h_a_n_s[neutron_bin]={nullptr};
    TH1D* h_kperp[neutron_bin]={nullptr};
    TH1D* h_kperp_s[neutron_bin]={nullptr};
    TH1D* h_mass_s[neutron_bin]={nullptr};
    TH1D* h_rapidity_s[neutron_bin]={nullptr};

    TH1D* h_k_lin[neutron_bin]={nullptr};

    TH1D* h_rapid_sgn[neutron_bin]={nullptr};
    TH1D* h_rapid_a150[neutron_bin]={nullptr}; //above 150
    TH1D* h_rapid_b150[neutron_bin]={nullptr};

    TH1D* h_kperp_prapid[neutron_bin]={nullptr}; //kperp for pos rapidity
    TH1D* h_kperp_nrapid[neutron_bin]={nullptr}; //kperp for neg rapidity

    TH1D* neutron_stat[neutron_bin]={nullptr}; //Soumya implied a TProfile, but I think a histogram might actually be easier?
    TH1D* h_n_raw_a={nullptr};
    TH1D* h_n_raw_c={nullptr};

    pbpb_reader();
};

void pbpb_reader::InitTree(){

    fChain = new TChain("HeavyIonD3PD");
    fChain->Add("data/*"); //all data
    //fChain->Add("data/Data_Run2018.root"); //just one data file for now

    //36491 bad pairs in 2018 of 179757. 20.3%
    //41548 bad pairs in 2023 of 197457. 21%
    //76420 bad pairs in 2024 of 271373. 28.2%


    fChain->SetMakeClass(1);

    //define all the branches and where they map onto
    fChain->SetBranchAddress("acop", &acop);
    fChain->SetBranchAddress("asym", &asym);
    fChain->SetBranchAddress("kperp", &kperp);
    fChain->SetBranchAddress("dphi", &dphi);
    fChain->SetBranchAddress("ptbar", &ptbar);
    fChain->SetBranchAddress("mass", &mass);
    fChain->SetBranchAddress("rapidity", &rapidity);
    fChain->SetBranchAddress("Egamma_a", &Egamma_a);
    fChain->SetBranchAddress("Egamma_b", &Egamma_b);
    fChain->SetBranchAddress("pt1", &pt1);
    fChain->SetBranchAddress("pt2", &pt2);
    fChain->SetBranchAddress("eta1", &eta1);
    fChain->SetBranchAddress("eta2", &eta2);
    fChain->SetBranchAddress("phi1", &phi1);
    fChain->SetBranchAddress("phi2", &phi2);
    fChain->SetBranchAddress("ch1", &ch1);
    fChain->SetBranchAddress("ch2", &ch2);
    fChain->SetBranchAddress("d01", &d01);
    fChain->SetBranchAddress("d02", &d02);
    fChain->SetBranchAddress("z01", &z01);
    fChain->SetBranchAddress("z02", &z02);
    //fChain->SetBranchAddress("RunNumber", &RunNumber);
    fChain->SetBranchAddress("NeutronsA", &NeutronsA);
    fChain->SetBranchAddress("NeutronsC", &NeutronsC);
    fChain->SetBranchAddress("Zvtx", &Zvtx);
    fChain->SetBranchAddress("qual1", &qual1);
    fChain->SetBranchAddress("qual2", &qual2);

    //Turn off all branches (for any residue) and then turn on the ones we want
    fChain->SetBranchStatus("*",0);

    fChain->SetBranchStatus("acop", 1);
    fChain->SetBranchStatus("asym", 1);
    fChain->SetBranchStatus("kperp", 1);
    fChain->SetBranchStatus("dphi", 1);
    fChain->SetBranchStatus("ptbar", 1);
    fChain->SetBranchStatus("mass", 1);
    fChain->SetBranchStatus("rapidity", 1);
    fChain->SetBranchStatus("Egamma_a", 1);
    fChain->SetBranchStatus("Egamma_b", 1);
    fChain->SetBranchStatus("pt1", 1);
    fChain->SetBranchStatus("pt2", 1);
    fChain->SetBranchStatus("eta1", 1);
    fChain->SetBranchStatus("eta2", 1);
    fChain->SetBranchStatus("phi1", 1);
    fChain->SetBranchStatus("phi2", 1);
    fChain->SetBranchStatus("ch1", 1);
    fChain->SetBranchStatus("ch2", 1);
    fChain->SetBranchStatus("d01", 1);
    fChain->SetBranchStatus("d02", 1);
    fChain->SetBranchStatus("z01", 1);
    fChain->SetBranchStatus("z02", 1);
    fChain->SetBranchStatus("RunNumber", 1);
    fChain->SetBranchStatus("NeutronsA", 1);
    fChain->SetBranchStatus("NeutronsC", 1);
    fChain->SetBranchStatus("Zvtx", 1);
    fChain->SetBranchStatus("qual1", 1);
    fChain->SetBranchStatus("qual2", 1);

}

#endif