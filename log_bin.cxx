void fillLogBinningArray(std::vector<double>& bins, int nBins, double low, double high) {
    bins.clear();

    double logLow = std::log10(low);
    double logHigh = std::log10(high);
    double logStep = (logHigh - logLow) / nBins;

    for (int i = 0; i <= nBins; ++i) {
        bins.push_back(std::pow(10, logLow + i * logStep));
    }
}

void test(){
    // // Create logarithmic pair-pT bin edges
    std::vector<double> pT_bins_60;

    fillLogBinningArray(pT_bins_60,  20, 8.0, 60.0);  // 20 log bins from 8 to 60 GeV
    
    std::cout << "Print out bins for pT 8-60 GeV" << std::endl;
    for (auto bin : pT_bins_60){
        std::cout << bin << ", ";
    }
    std::cout << std::endl;
}
