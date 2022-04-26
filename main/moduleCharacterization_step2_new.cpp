#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
//#include "interface/Na22SpectrumAnalyzerSingleBar.h"
#include "interface/Na22SpectrumAnalyzerSingleBar_TOFHIR2.h"
//#include "interface/Na22SpectrumAnalyzerModule_TOFHIR2.h"
#include "interface/Co60SpectrumAnalyzer_2Peaks.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TSpectrum.h"



Double_t langaufun(Double_t *x, Double_t *par)
{
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0];
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}






//---- find energy bins
void GetEnergyBins(TH1F *h, std::vector<float> *r, std::map<int, float> & b){

  for(unsigned int i = 1; i < r->size(); i++){
    TH1F *binHisto = new TH1F ( "binHisto", "binHisto", h -> FindBin(r->at(i)) - h->FindBin(r-> at(i-1)), r-> at(i-1), r->at(i));
    int j = 1;
    for (int bin = h->FindBin(r->at(i-1)) ; bin < h -> FindBin(r->at(i))+1 ; bin++){
      binHisto -> SetBinContent( j, h->GetBinContent(bin));
      j++;
    }
    b[i] = binHisto -> GetMean();
    binHisto -> Delete();
  }
}



// ============  ********* MAIN  ************ =============== //
int main(int argc, char** argv)
{
  setTDRStyle();
  
  typedef std::numeric_limits<double> dbl;
        std::cout.precision(dbl::max_digits10);
  if( argc < 2 )
  {
    std::cout << ">>> moduleCharacterization_step2::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  

  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  
  //--- get parameters
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  //system(Form("rm -r %s", plotDir.c_str())); // questo non va bene se stiamo lavorando in parallelo
  system(Form("mkdir -p %s",plotDir.c_str()));
  //system(Form("mkdir -p %s/qfine/",plotDir.c_str()));
  system(Form("mkdir -p %s/tot/",plotDir.c_str()));
  system(Form("mkdir -p %s/totRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/t1fine/",plotDir.c_str()));
  //system(Form("mkdir -p %s/CTR/",plotDir.c_str()));
  system(Form("mkdir -p %s/test_plots/",plotDir.c_str()));
  system(Form("mkdir -p %s/amplitudeWalkCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_amplitudeWalkCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_amplitudeWalkPosCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_amplitudeWalkPhasePosCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatioCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/totRatioCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/phaseCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/positionCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_energyRatioCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_totRatioCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_energyRatioPhaseCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_totRatioPhaseCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_energyRatioPhasePosCorr/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_totRatioPhasePosCorr/",plotDir.c_str()));
  

  
  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  LRLabels.push_back("L-R");
  
  
  std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
  std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  
  std::map<float,int> map_energyMins;
  std::map<float,int> map_energyMaxs;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii)
    {
      map_energyMins[Vov[ii]] = energyMins[ii];
      map_energyMaxs[Vov[ii]] = energyMaxs[ii];
    }
  
  int useTrackInfo = opts.GetOpt<int>("Input.useTrackInfo");
  
  //--- open files
  std::string step1FileName= opts.GetOpt<std::string>("Input.step1FileName");
  TFile* inFile = TFile::Open(step1FileName.c_str(),"READ");

/*
// -- read minimum energy for each bar from file
  std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
  std::map < std::pair<int, float>, float> minE; 
  std::cout << minEnergiesFileName <<std::endl;
  if (minEnergiesFileName != "") {
    std::ifstream minEnergiesFile;
    minEnergiesFile.open(minEnergiesFileName);
    std::string line;
    int bar;
    float ov;
    float value;
    while ( minEnergiesFile.good() ){
      getline(minEnergiesFile, line);
      std::istringstream ss(line);
      ss >> bar >> ov >> value; 
      minE[std::make_pair(bar,ov)] = value; 
      std::cout<< bar <<  "   " << ov << "  " << minE[std::make_pair(bar,ov)] <<std::endl;
    }
  }
*/

  std::map<std::string,TTree*> trees;
  
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;  
  
  TList* list = inFile -> GetListOfKeys();
  TIter next(list);
  TObject* object = 0;
  while( (object = next()) )
  {
    std::string name(object->GetName());
    std::vector<std::string> tokens = GetTokens(name,'_');
    std::size_t found;
     
    found = name.find("data_");
    //tree
    if( found!=std::string::npos )
    {
      std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
    }
    found = name.find("h1_energy_b");
    if( found!=std::string::npos )
    {
     //Vov e th
      std::string stepLabel = tokens[3]+"_"+tokens[4];
      VovLabels[tokens[3]] += 1;
      thLabels[tokens[4]] += 1;
      stepLabels.push_back(stepLabel);
      std::string string_Vov = tokens[3];
      string_Vov.erase(0,3);
      map_Vovs[stepLabel] = atof(string_Vov.c_str());
      std::string string_th = tokens[4];
      string_th.erase(0,2);
      map_ths[stepLabel] = atof(string_th.c_str());
    }
  }
  std::sort(stepLabels.begin(),stepLabels.end());
  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep2");
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  outFile->cd();

  std::map<int,TTree*> outTrees;
  std::map<double,TTree*> outTrees2;
  float energy511LR;
  float energy1275LR;
  float energy1786LR;
  float energy511R;
  float energy1275R;
  float energy1786R;
  float energy511L;
  float energy1275L;
  float energy1786L;
  float theIndex;
  


  std::map<double,TH1F*> h1_energyRatio;
  std::map<double,TH1F*> h1_totRatio;
  std::map<double,TH1F*> h1_t1fineMean;
  std::map<double,TH1F*> h1_qT1Mean;
  std::map<int,TH1F*> h1_totRatioRef;
  
  std::map<double,TH1F*> h1_deltaT_raw;
  std::map<double,TH1F*> h1_deltaT;
  std::map<double,TH1F*> h1_deltaT_P;
  std::map<double,TH1F*> h1_deltaT_tConverted;
  //std::map<double,TH1F*> h1_tMCP;

  std::map<double,TProfile*> p1_deltaT_vs_energyRatio;
  std::map<double,TProfile*> p1_deltaT_vs_totRatio;
  std::map<double,TProfile*> p1_tL_tMCP_vs_energyL;
  std::map<double,TProfile*> p1_tR_tMCP_vs_energyR;
  std::map<double,TH2F*> h2_tL_tMCP_vs_energyL;
  std::map<double,TH2F*> h2_tR_tMCP_vs_energyR;
  
  std::map<double,TProfile*> p1_deltaT_energyRatioCorr_vs_t1fineMean;
  std::map<double,TProfile*> p1_deltaT_totRatioCorr_vs_t1fineMean;
  std::map<double,TProfile*> p1_deltaT_energyRatioCorr_vs_posX;
  std::map<double,TProfile*> p1_deltaT_totRatioCorr_vs_posX;
  std::map<double,TProfile*> p1_deltaT_ampWalkEnergyCorr_vs_posX;
  std::map<double,TProfile*> p1_deltaT_ampWalkTotCorr_vs_posX;
  std::map<double,TProfile*> p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean;
  std::map<double,TProfile*> p1_deltaT_ampWalkTotCorr_vs_t1fineMean;
  
  std::map<double,TH1F*> h1_tR_tMCP;
  std::map<double,TH1F*> h1_tL_tMCP;
  std::map<double,TProfile*> p1_deltaT_totRatioCorr_vs_totRatio;
  std::map<double,TH2F*> h2_deltaT_vs_totRatio;
  std::map<double,TH2F*> h2_deltaT_totRatioCorr_vs_totRatio;

  std::map<double,TH1F*> h1_deltaT_ampWalkEnergyCorr;
  //std::map<double,TH1F*> h1_deltaT_ampWalkTotCorr;
  std::map<double,TH1F*> h1_deltaT_ampWalkEnergyPosCorr;
  //std::map<double,TH1F*> h1_deltaT_ampWalkTotPhaseCorr;
  std::map<double,TH1F*> h1_deltaT_ampWalkEnergyPhasePosCorr;
  //std::map<double,TH1F*> h1_deltaT_ampWalkTotPhasePosCorr;

  std::map<double,TH1F*> h1_deltaT_energyRatioCorr;
  std::map<double,TH1F*> h1_deltaT_totRatioCorr;
  std::map<double,TH1F*> h1_deltaT_energyRatioPhaseCorr;
  std::map<double,TH1F*> h1_deltaT_totRatioPhaseCorr;
  std::map<double,TH1F*> h1_deltaT_energyRatioPhasePosCorr;
  std::map<double,TH1F*> h1_deltaT_totRatioPhasePosCorr;

  
  std::map<std::string, std::map<int, std::vector<float>*> > ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int, std::vector<float>*> > mpv_energy;
  std::map<std::string, std::map<int, std::map<std::string,std::pair<float,float> > > > peaks;	//peaks[LRlabel][index][energyPeak]
  std::map<std::string, std::map<int, std::map<int,float> > > energyBin; // energyBin[LRlabel][index]

  std::map<int,TF1*>  f_langaus; // f_langaus[index]
  std::map<int,TF1*>  f_gaus; // f_gaus[index]
  std::map<int,TF1*>  f_landau; // f_gaus[index]
  std::map<float,int>  Vov_LandauMin; //Vov_LandauMin[Vov]
  Vov_LandauMin[1.0] = 0;
  Vov_LandauMin[1.5] = 50;
  Vov_LandauMin[2.0] = 100;
  Vov_LandauMin[2.5] = 180;
  Vov_LandauMin[3.0] = 250;
  Vov_LandauMin[3.5] = 300;
  Vov_LandauMin[6.0] = 300;
  
  
  //--- get plot settings
  TCanvas* c;
  TCanvas* c2;
  float* vals = new float[6];
  TLatex* latex;
  TLatex* latexRef;
  TH1F* histo, histo1;
  TH2F* histo2;
  TProfile* prof, prof1;

  
  //------------------
  //--- draw 1st plots
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  std::string Na22 = "Na22";
  std::string Na22SingleBar = "Na22SingleBar";
  std::string Co60 = "Co60";
  std::string Co60SumPeak = "Co60SumPeak";	
  std::string Laser = "Laser";	
  std::string TB = "TB";
  std::string keepAll = "keepAll";
  std::vector<int> barList = opts.GetOpt<std::vector<int> >("Plots.barList");// list of bars to be analyzed read from cfg
  

  for(auto stepLabel : stepLabels) {
      
    float Vov = map_Vovs[stepLabel];
    float vth1 = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.2f",Vov));
    std::string thLabel(Form("th%02.0f",vth1));
      
    //--------------------------------------------------------
    // --- loop over bars
    for(int iBar = 0; iBar < 16; ++iBar) {

      bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
      if (!barFound) continue;

      int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );

	  
      
      outTrees[index] = new TTree(Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth1),Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth1));
      outTrees[index] -> Branch("energyPeak511LR",&energy511LR);
      outTrees[index] -> Branch("energyPeak1275LR",&energy1275LR);
      outTrees[index] -> Branch("energyPeak1786LR",&energy1786LR);
      outTrees[index] -> Branch("energyPeak511L",&energy511L);
      outTrees[index] -> Branch("energyPeak1275L",&energy1275L);
      outTrees[index] -> Branch("energyPeak1786L",&energy1786L);
      outTrees[index] -> Branch("energyPeak511R",&energy511R);
      outTrees[index] -> Branch("energyPeak1275R",&energy1275R);
      outTrees[index] -> Branch("energyPeak1786R",&energy1786R);
      outTrees[index] -> Branch("indexID",&theIndex);
      
	  	      

      // -- loop over L, R, LR
      for(auto LRLabel : LRLabels ) {
	
	//label histo
        std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));
		//std::string labelRef(Form("%s_%s",LRLabel.c_str(),stepLabel.c_str()));

        latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth1)));
		//latexRef = new TLatex(0.40,0.85,Form("#splitline{reference bar%s}{V_{OV} = %.2f V, th. = %d DAC}",LRLabel.c_str(),Vov,int(vth1)));
        if (LRLabel == "L-R") { 
          latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	}

        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);

        // -- qfine and tot only per L, R
        if (LRLabel == "R" || LRLabel == "L") {

          /*histo = (TH1F*)( inFile->Get(Form("h1_qfine_%s",label.c_str())) );
          if( !histo ) continue;
          if( histo->GetEntries() < 100 ) continue;
          
          c = new TCanvas(Form("c_qfine_%s",label.c_str()),Form("c_qfine_%s",label.c_str()));
          gPad -> SetLogy();
                  
          histo -> SetTitle(";Q_{fine} [ADC];entries");
          histo -> SetLineColor(kRed);
          histo -> Draw();
          histo -> GetXaxis() -> SetRangeUser(0.,600.);
          // TLine* line_qfineAcc1 = new TLine(cut_qfineAcc[chID][Vov],histo->GetMinimum(),cut_qfineAcc[chID][Vov],histo->GetMaximum());
          // line_qfineAcc1 -> SetLineColor(kBlack);
          // line_qfineAcc1 -> Draw("same");
          latex -> Draw("same");      
          histo -> Write();
          c -> Print(Form("%s/qfine/c_qfine__%s.png",plotDir.c_str(),label.c_str()));
          c -> Print(Form("%s/qfine/c_qfine__%s.pdf",plotDir.c_str(),label.c_str()));
          delete c;*/
        
	  // -- draw ToT
          c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
          gPad -> SetLogy();
        
          histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
		  if( !histo ) continue;
          histo -> SetTitle(";ToT [ns];entries");
          histo -> SetLineColor(kRed);
          histo -> Draw();
          // TLine* line_totAcc1 = new TLine(cut_totAcc[chID][Vov],histo->GetMinimum(),cut_totAcc[chID][Vov],histo->GetMaximum());
          // line_totAcc1 -> SetLineColor(kBlack);
          // line_totAcc1 -> Draw("same");
          latex -> Draw("same");      
          histo -> Write();
          c -> Print(Form("%s/tot/c_tot__%s.png",plotDir.c_str(),label.c_str()));
          c -> Print(Form("%s/tot/c_tot__%s.pdf",plotDir.c_str(),label.c_str()));
          delete c;
	}
        
        // -- draw energy
        c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
        gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );      
        if( !histo ) continue;
        histo -> SetTitle(";energy [a.u.];entries");
        histo -> SetLineColor(kRed);
        histo -> SetLineWidth(2);
        histo -> Draw();
	

        // --- look for peaks and define energy ranges
        if( source.compare(Na22) && source.compare(Na22SingleBar) && source.compare(Co60) && source.compare(Co60SumPeak) && source.compare(Laser) && source.compare(TB) && source.compare(keepAll) )
	  {
	    std::cout << " Source not found !!! " << std::endl;
	    return(0);
	  }
	
        ranges[LRLabel][index] = new std::vector<float>;
		mpv_energy[LRLabel][index] = new std::vector<float>;


        // --- Na22 or Co60 spectrum
        if(!source.compare(Na22)  || !source.compare(Na22SingleBar) ||  !source.compare(Co60) )
        {
          std::string firstPeak = "";
          if (!source.compare(Na22)) {
            peaks[LRLabel][index] = Na22SpectrumAnalyzer(histo,ranges[LRLabel][index]);
            firstPeak = "0.511 MeV";
          }
          if (!source.compare(Na22SingleBar)) {
	    //peaks[LRLabel][index] = Na22SpectrumAnalyzerSingleBar(histo,ranges[LRLabel][index]);//
	    // !!!!! fare in modo di avere un singolo codice Na22SpectrumAnalyzer, eventualmente con delle flag...
            peaks[LRLabel][index] = Na22SpectrumAnalyzerSingleBar_TOFHIR2(histo,ranges[LRLabel][index]);//run330, barInArray
	    //peaks[LRLabel][index] = Na22SpectrumAnalyzerModule_TOFHIR2(histo,ranges[LRLabel][index]);//run555, bar read out by two modules, "wirelessbar"
            firstPeak = "0.511 MeV";	
          }
          if (!source.compare(Co60)) {
            peaks[LRLabel][index] = Co60SpectrumAnalyzer_2Peaks(histo,ranges[LRLabel][index]);
            firstPeak = "1.173 MeV";
          }
          
          if (peaks[LRLabel][index][firstPeak].first > 0){
            //histo -> GetXaxis() -> SetRangeUser(0.,5.*peaks[LRLabel][index][firstPeak].first);
	    histo -> GetXaxis() -> SetRangeUser(0.,1000.);
          }
	  
          if (peaks[LRLabel][index][firstPeak].first== -9999){
            histo -> GetXaxis() -> SetRangeUser(0., map_energyMaxs[Vov]); 
            peaks[LRLabel].erase(index);
            ranges[LRLabel].erase(index);
            if (!source.compare(Na22) || !source.compare(Na22SingleBar) ) peaks[LRLabel][index]["0.511 MeV"].first = -10;
            if (!source.compare(Na22) || !source.compare(Na22SingleBar) ) peaks[LRLabel][index]["1.275 MeV"].first = -10;
            if (!source.compare(Na22SingleBar))                           peaks[LRLabel][index]["1.786 MeV"].first = -10;
            if (!source.compare(Co60))                                    peaks[LRLabel][index]["1.173 MeV"].first = -10;
            if (!source.compare(Co60))                                    peaks[LRLabel][index]["1.332 MeV"].first = -10;
            if (!source.compare(Co60))                                    peaks[LRLabel][index]["2.505 MeV"].first = -10;
          }
          
          if(peaks[LRLabel][index][firstPeak].first != -10){
            GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
          }          
        } // end if Na22/Co60
          
        
        // -- if Co60 sum peak or laser, we don't use the spectrum analyzers - just gaussian fit to the energy peak.
        if( !source.compare(Co60SumPeak) || !source.compare(Laser) )
        { 
          histo -> GetXaxis() -> SetRangeUser(map_energyMins[Vov],map_energyMaxs[Vov]);
          float max = FindXMaximum(histo,map_energyMins[Vov],map_energyMaxs[Vov]);
          
          TF1* fitFunc = new TF1 ( Form("func_%d",index), "gaus(0)", max-2*histo->GetRMS(), max+2*histo->GetRMS() );
          fitFunc -> SetLineColor(kBlack);
          fitFunc -> SetLineWidth(2);
          histo -> Fit( fitFunc, "NQR");
	  fitFunc -> SetRange(fitFunc->GetParameter(1) - fitFunc->GetParameter(2)*5, fitFunc->GetParameter(1) + fitFunc->GetParameter(2)*5 );
	  histo -> Fit( fitFunc, "QRS+");
          
          // -- questa condizione er solo per il laser... puo' andare bene tenerla anche per 60Co?
          if ((fitFunc->GetParameter(1) - fitFunc->GetParameter(2)*5)>0) 
            ranges[LRLabel][index] -> push_back( fitFunc->GetParameter(1) - fitFunc->GetParameter(2)*5);
          else 
            ranges[LRLabel][index] -> push_back(0);
          ranges[LRLabel][index] -> push_back( fitFunc->GetParameter(1) + fitFunc->GetParameter(2)*5);

          // use energy511XX branch to save the peak value 
          if (LRLabel == "L") energy511L = fitFunc->GetParameter(1);
	  if (LRLabel == "R") energy511R = fitFunc->GetParameter(1);
	  if (LRLabel == "L-R") energy511LR = fitFunc->GetParameter(1);					

          for(auto range: (*ranges[LRLabel][index])){
            float yval = std::max(10., histo->GetBinContent(histo->FindBin(range)));
            TLine* line = new TLine(range,3.,range, yval);
            line -> SetLineWidth(1);
            line -> SetLineStyle(7);
            line -> Draw("same");
	  }
	  
	  GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	}// end Co60SumPeak or laser	
	
        // -- keep all events within energyMin and energyMax
        if( !source.compare(keepAll) )
        { 
	  ranges[LRLabel][index] -> push_back( map_energyMins[Vov] );
          ranges[LRLabel][index] -> push_back( map_energyMins[Vov] );
	  
          for(auto range: (*ranges[LRLabel][index]))
	    {
	      float yval = std::max(10., histo->GetBinContent(histo->FindBin(range)));
	      TLine* line = new TLine(range,3.,range, yval);
	      line -> SetLineWidth(1);
	      line -> SetLineStyle(7);
	      line -> Draw("same");
	    }
	  
	  GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	}// end keepAll
	

	// -- if MIP peak, we don't use the spectrum analyzers - just langaus fit to the energy peak.
        if(!source.compare(TB)){ 
	  /*
	  //TF1* f_langaus = new TF1("f_langaus", langaufun, 300.,1000.,4);
	    f_langaus[index] = new TF1(Form("fit_energy_bar%02d_Vov%.2f_vth1_%02.0f",iBar,Vov,vth1),langaufun,Vov_LandauMin[Vov],1000.,4);
	    f_langaus[index] -> SetParameters(f_landau->GetParameter(2),f_landau->GetParameter(1),histo->Integral(histo->FindBin(300.),histo->FindBin(1000.))*histo->GetBinWidth(1),10.);
	    f_langaus[index] -> SetLineColor(kBlack);
	    f_langaus[index] -> SetLineWidth(2);
	    histo -> Fit(f_langaus[index],"QNRS+");
	    f_langaus[index] -> Draw("same");
	    
	    ranges[LRLabel][index] -> push_back( 0.80*f_langaus[index]->GetParameter(1));
	    ranges[LRLabel][index] -> push_back( histo -> GetBinCenter(500) );
	    
	    // use energy511XX branch to save the peak value 
	    if (LRLabel == "L") energy511L = f_langaus[index]->GetParameter(1);					
	    if (LRLabel == "R") energy511R = f_langaus[index]->GetParameter(1);					
	    if (LRLabel == "L-R") energy511LR = f_langaus[index]->GetParameter(1);	*/		
	 
	  if( opts.GetOpt<int>("Channels.array") == 1){
	    histo->GetXaxis()->SetRangeUser(25., 650.);
	    if (Vov>3.0) histo->GetXaxis()->SetRangeUser(120., 650.);
	    //histo->GetXaxis()->SetRangeUser(50., 650.);
	    //histo->GetXaxis()->SetRangeUser(30 + Vov*10,650);
	    //histo->GetXaxis()->SetRangeUser(-130. + Vov*115., 650);
	    //if (iBar==0 || iBar==14) histo->GetXaxis()->SetRangeUser(-160. + Vov*115., 650); 
	    //if (iBar==0) histo->GetXaxis()->SetRangeUser(-160. + Vov*115., 650); 
	  }
	  if( opts.GetOpt<int>("Channels.array") == 0){
	    histo->GetXaxis()->SetRangeUser(20., 650.);
	    //histo->GetXaxis()->SetRangeUser(TMath::Max(-130. + Vov*115,40.), 650);
	  }
	  float max = histo->GetBinCenter(histo->GetMaximumBin());
	  histo->GetXaxis()->SetRangeUser(0,1000);
	  
	  f_gaus[index] = new TF1(Form("fit_energy_bar%02d%s_Vov%.2f_vth1_%02.0f",iBar,LRLabel.c_str(),Vov,vth1), "gaus", max-50, max+50);
	  f_gaus[index]->SetParameters(histo->GetMaximumBin(), max, 70);
	  histo->Fit(f_gaus[index], "QRS");
	  f_gaus[index]->SetRange(f_gaus[index]->GetParameter(1)-f_gaus[index]->GetParameter(2), f_gaus[index]->GetParameter(1)+f_gaus[index]->GetParameter(2));
	  histo->Fit(f_gaus[index], "QRS");
	  f_gaus[index] -> SetLineColor(kBlue);
	  f_gaus[index] -> SetLineWidth(2);
	  f_gaus[index] -> SetLineStyle(2);
	  //f_gaus[index] -> Draw("same");
	  
	  //ranges[LRLabel][index] -> push_back( 0.80*f_gaus[index]->GetParameter(1));
	  //ranges[LRLabel][index] -> push_back( histo -> GetBinCenter(700) ); // to avoid sturation
	  
	  f_landau[index] = new TF1(Form("f_landau_bar%02d%s_Vov%.2f_vth1_%02.0f", iBar,LRLabel.c_str(),Vov,vth1),"[0]*TMath::Landau(x,[1],[2])", 0,1000.);
	  f_landau[index] -> SetRange(max * 0.75, max * 1.5);    //0.75
	  f_landau[index] -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max, Vov*5.);
	  histo -> Fit(f_landau[index],"QRS");
	  mpv_energy[LRLabel][index] -> push_back( f_landau[index]->GetParameter(1) ); 

	  f_landau[index] -> SetLineColor(kBlack);
	  f_landau[index] -> SetLineWidth(2);
	  f_landau[index] -> Draw("same");
	  
	  if (step1FileName.find("2E14")!= std::string::npos && iBar == 14) ranges[LRLabel][index] -> push_back(5); //
	  if (iBar == 7 && LRLabel == "L" ) ranges[LRLabel][index] -> push_back(20);
	  else if (iBar == 5 && LRLabel == "L" ) ranges[LRLabel][index] -> push_back(1);
	  else if (iBar == 5 && LRLabel == "R" ) ranges[LRLabel][index] -> push_back(40);
	  else if (iBar == 13 && LRLabel == "R" ) ranges[LRLabel][index] -> push_back(15);
	  else{
	    if (f_landau[index]->GetParameter(1)> 10) 
	      ranges[LRLabel][index] -> push_back( 0.75*f_landau[index]->GetParameter(1));   //0.75
	    else
	      ranges[LRLabel][index] -> push_back( 20 ); 
	  }
	  ranges[LRLabel][index] -> push_back( 940. );  //700.

          // use energy511XX branch to save the peak value 
          if (LRLabel == "L") energy511L = f_landau[index]->GetParameter(1);					
          if (LRLabel == "R") energy511R = f_landau[index]->GetParameter(1);					
          if (LRLabel == "L-R") energy511LR = f_landau[index]->GetParameter(1);

	  for(auto range: (*ranges[LRLabel][index])){
	    TLine* line = new TLine(range,0.,range, histo->GetMaximum());
            line -> SetLineWidth(2);
            line -> SetLineStyle(7);
            line -> Draw("same");
	  }
	  
	  GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	}// end MIP (TB)	

	
	
	// -- draw and print energy plots 
	latex -> Draw("same"); 
	outFile -> cd();
	histo->Write();     
	c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label.c_str()));
	c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
	delete c;
	
      }// end loop over L, R, L-R labels
      

      // --- now fill the output tree.  
      // *********** MM: Ma serve ? non basterebbe leggere i valori di posizione dei picchi e numero di eventi dagli istogrammi che vengono comunque salvati nel file di output? 
      
      theIndex = index;
      
      if(!source.compare(Na22)  || !source.compare(Na22SingleBar) ){
	energy511R  = peaks["R"][index]["0.511 MeV"].first;		
	energy511L  = peaks["L"][index]["0.511 MeV"].first;		
	energy511LR = peaks["L-R"][index]["0.511 MeV"].first;
	
	energy1275R = peaks["R"][index]["1.275 MeV"].first;
	energy1275L = peaks["L"][index]["1.275 MeV"].first;
	energy1275LR = peaks["L-R"][index]["1.275 MeV"].first;
	
	if( !source.compare(Na22SingleBar) ) {
	  energy1786R = peaks["R"][index]["1.786 MeV"].first;
	  energy1786L = peaks["L"][index]["1.786 MeV"].first;
	  energy1786LR = peaks["L-R"][index]["1.786 MeV"].first;
	}       
      }
      
      if(!source.compare(Co60)){
	
	energy511R  = peaks["R"][index]["1.173 MeV"].first;
	energy511L  = peaks["L"][index]["1.173 MeV"].first;
	energy511LR = peaks["L-R"][index]["1.173 MeV"].first;
	
	energy1275R = peaks["R"][index]["1.332 MeV"].first;
	energy1275L = peaks["L"][index]["1.332 MeV"].first;
	energy1275LR = peaks["L-R"][index]["1.332 MeV"].first;
	
	energy1786R = peaks["R"][index]["2.505 MeV"].first;
	energy1786L = peaks["L"][index]["2.505 MeV"].first;
	energy1786LR = peaks["L-R"][index]["2.505 MeV"].first;
	
      }

      if (!source.compare(TB)){
	energy1275LR = -10; 
	energy1786LR = -10;  
	energy1275R  = -10; 
	energy1786R  = -10;
	energy1275R  = -10;
	energy1786R  = -10;
      }
      
      outTrees[index] -> Fill();
      
    }// -- end loop over bars
    
  } // -- end loop over stepLabels
  
  // ---  end 1st plots
  
  
  //------------------------
  //--- 2nd loop over events
  std::map<int,std::map<int,bool> > accept;
  
  for(auto mapIt : trees)
    {
      
      ModuleEventClass* anEvent = new ModuleEventClass();
      
      mapIt.second -> SetBranchAddress("event",&anEvent);
      
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) std::cout << ">>> 2nd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	  mapIt.second -> GetEntry(entry);
	  
	  bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
	  if (!barFound) continue;

	  
	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	
	  	  
	  accept[index1][entry] = false;
	  
	  if(!ranges["L-R"][index1] ) continue;
	  
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  
	  if( energyBinAverage < 1 ) continue;
	  
	  accept[index1][entry] = true;

	  double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

	  if(anEvent->amp_MCP < 500) continue;
	  
	  if( h1_energyRatio[index2] == NULL )
	    {

	      std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
	      std::string labelL(Form("bar%02dL_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
		  std::string labelR(Form("bar%02dR_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
		  std::string labelLR(Form("bar%02dL-R_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
	      
          h1_energyRatio[index2] = new TH1F(Form("h1_energyRatio_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);
	      h1_totRatio[index2] = new TH1F(Form("h1_totRatio_%s",labelLR_energyBin.c_str()),"",2000,0.,5.);
		  //if(anEvent->barID == 8) h1_totRatioRef[index] = new TH1F(Form("h1_totRatio_%s",labelLR_Ref.c_str()),"",2000,0.,5.);
	      h1_t1fineMean[index2] = new TH1F(Form("h1_t1fineMean_%s",labelLR_energyBin.c_str()),"",1000,0.,1000.);
		  h1_qT1Mean[index1] = new TH1F(Form("h1_qT1Mean_%s",labelLR.c_str()),"",1000,0.4,1.6);
	      h1_deltaT_raw[index2] = new TH1F(Form("h1_deltaT_raw_%s",labelLR_energyBin.c_str()),"",2000,-24000.,24000.);
		  h1_tL_tMCP[index1] = new TH1F(Form("h1_tL-tMCP_%s",labelL.c_str()),"",800,-4.,4.);
		  h1_tR_tMCP[index1] = new TH1F(Form("h1_tR-tMCP_%s",labelR.c_str()),"",800,-4.,4.);
	    }
	  
	  if (fabs(anEvent->timeR-anEvent->timeL)<10000) {
	    
	    if ((anEvent->energyR / anEvent->energyL >0) & (anEvent->energyR / anEvent->energyL <5)){
	      h1_energyRatio[index2] -> Fill( anEvent->energyR / anEvent->energyL );						     
	      h1_totRatio[index2] -> Fill( anEvent->totR / anEvent->totL );
		  //if(anEvent->barID == 8) h1_totRatioRef[index] -> Fill(anEvent->totR_ref / anEvent->totL_ref);
	      h1_deltaT_raw[index2] -> Fill( anEvent->timeR-anEvent->timeL );
		  h1_tL_tMCP[index1] -> Fill(anEvent->phi_L - anEvent->phi_MCP_C);
		  h1_tR_tMCP[index1] -> Fill(anEvent->phi_R - anEvent->phi_MCP_C);

	    }
	    
	    h1_t1fineMean[index2] -> Fill( 0.5 * (anEvent->t1fineR + anEvent->t1fineL) );
		h1_qT1Mean[index1] -> Fill( 0.5 * (anEvent->qT1R + anEvent->qT1L) );
	  }
	      

	} //-- end loop over entries
      
      std::cout << std::endl;

    }
  
  
  
  
  //------------------
  //--- draw 2nd plots
  std::map<double,float> CTRMeans;
  std::map<double,float> CTRSigmas;
  
  std::map<double,TF1*> fitFunc_energyRatio;
  std::map<double,TF1*> fitFunc_totRatio;

  std::map<double,TF1*> fitFunc_tL_tMCP;
  std::map<double,TF1*> fitFunc_tR_tMCP;

  for(auto mapIt : h1_deltaT_raw)
    {
      double index = mapIt.first;
      
      FindSmallestInterval(vals,h1_deltaT_raw[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans[index] = mean;
      CTRSigmas[index] = effSigma;
    }
  
  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];

      for(int iBar = 0; iBar < 16; ++iBar) 
	{
	  bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
          if (!barFound) continue;      
	  
	  std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
	  
	  
	  int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
	  
	  if( !ranges["L-R"][index1] ) continue;
	  
	  int nEnergyBins = ranges["L-R"][index1]->size()-1;
	  
	  for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
	    {
	      //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
	      double index2( 10000000*iEnergyBin+index1 );

		  std::string labelL(Form("bar%02dL_%s",iBar,stepLabel.c_str()));
	  std::string labelR(Form("bar%02dR_%s",iBar,stepLabel.c_str()));

	     //draw tL - tMCP

		if(!h1_tL_tMCP[index1]) continue;

		c = new TCanvas(Form("c_tL-tMCP_%s",labelL.c_str()),Form("c_tL-tMCP_%s",labelL.c_str()));
	      
	      histo = h1_tL_tMCP[index1];
	      
	      histo -> SetTitle(Form("; #phi_{L} - #phi_{MCP} [clock units];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kBlue);
	      histo -> SetMarkerColor(kBlue);
	      histo -> Draw("");

		  histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-1.5*histo->GetRMS(),histo->GetMean()+1.5*histo->GetRMS());    

		  fitFunc_tL_tMCP[index1] = new TF1(Form("fitFunc_tL_tavg_%s",labelL.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
	      histo -> Fit(fitFunc_tL_tMCP[index1],"QNRS");
	      histo -> Fit(fitFunc_tL_tMCP[index1],"QNSR","",fitFunc_tL_tMCP[index1]->GetParameter(1)-2.5*fitFunc_tL_tMCP[index1]->GetParameter(2),fitFunc_tL_tMCP[index1]->GetParameter(1)+2.5*fitFunc_tL_tMCP[index1]->GetParameter(2));
	      histo -> Fit(fitFunc_tL_tMCP[index1],"QSR+","",fitFunc_tL_tMCP[index1]->GetParameter(1)-2.*fitFunc_tL_tMCP[index1]->GetParameter(2),fitFunc_tL_tMCP[index1]->GetParameter(1)+2.*fitFunc_tL_tMCP[index1]->GetParameter(2));
	      
	      fitFunc_tL_tMCP[index1] -> SetLineColor(kRed);
	      fitFunc_tL_tMCP[index1] -> SetLineWidth(2);
	      fitFunc_tL_tMCP[index1] -> Draw("same");

		  outFile -> cd();
	      histo -> Write();

		  c -> Print(Form("%s/test_plots/c_tL-tMCP_%s.pdf",plotDir.c_str(),labelL.c_str()));
	      c -> Print(Form("%s/test_plots/c_tL-tMCP_%s.png",plotDir.c_str(),labelL.c_str()));
	      delete c;

         //draw tR - tMCP

		if(!h1_tR_tMCP[index1]) continue;

		c = new TCanvas(Form("c_tR-tMCP_%s",labelR.c_str()),Form("c_tR-tMCP_%s",labelR.c_str()));
	      
	      histo = h1_tR_tMCP[index1];
	      
	      histo -> SetTitle(Form("; #phi_{R} - #phi_{MCP} [clock units];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kBlue);
	      histo -> SetMarkerColor(kBlue);
	      histo -> Draw("");

		  histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-1.5*histo->GetRMS(),histo->GetMean()+1.5*histo->GetRMS());    

		  fitFunc_tR_tMCP[index1] = new TF1(Form("fitFunc_tR_tavg_%s",labelR.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
	      histo -> Fit(fitFunc_tR_tMCP[index1],"QNRS");
	      histo -> Fit(fitFunc_tR_tMCP[index1],"QNSR","",fitFunc_tR_tMCP[index1]->GetParameter(1)-2.5*fitFunc_tR_tMCP[index1]->GetParameter(2),fitFunc_tR_tMCP[index1]->GetParameter(1)+2.5*fitFunc_tR_tMCP[index1]->GetParameter(2));
	      histo -> Fit(fitFunc_tR_tMCP[index1],"QSR+","",fitFunc_tR_tMCP[index1]->GetParameter(1)-2.*fitFunc_tR_tMCP[index1]->GetParameter(2),fitFunc_tR_tMCP[index1]->GetParameter(1)+2.*fitFunc_tR_tMCP[index1]->GetParameter(2));
	      
	      fitFunc_tR_tMCP[index1] -> SetLineColor(kRed);
	      fitFunc_tR_tMCP[index1] -> SetLineWidth(2);
	      fitFunc_tR_tMCP[index1] -> Draw("same");

		  outFile -> cd();
	      histo -> Write();

		  c -> Print(Form("%s/test_plots/c_tR-tMCP_%s.pdf",plotDir.c_str(),labelR.c_str()));
	      c -> Print(Form("%s/test_plots/c_tR-tMCP_%s.png",plotDir.c_str(),labelR.c_str()));
	      delete c;



	      if (!h1_energyRatio[index2]) continue;
	      
	      std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));
	      
	      // -- draw energy ratio 
	      c = new TCanvas(Form("c_energyRatio_%s",labelLR_energyBin.c_str()),Form("c_energyRatio_%s",labelLR_energyBin.c_str()));
	      histo = h1_energyRatio[index2];
	      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	      histo -> SetTitle(Form(";energy_{right} / energy_{left};entries"));
	      histo -> SetLineColor(kRed);
	      histo -> SetLineWidth(2);
	      histo -> Draw();
	      histo -> Write();
	      
	      fitFunc_energyRatio[index2] = new TF1(Form("fitFunc_energyRatio_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
	      histo -> Fit(fitFunc_energyRatio[index2],"QNRS");
	      histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));
	      histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));
	      
	      fitFunc_energyRatio[index2] -> SetLineColor(kBlack);
	      fitFunc_energyRatio[index2] -> SetLineWidth(2);
	      fitFunc_energyRatio[index2] -> Draw("same");
	      
	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");
	      
	      c -> Print(Form("%s/energyRatio/c_energyRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/energyRatio/c_energyRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete c;


	      // -- draw tot ratio 
	      c = new TCanvas(Form("c_totRatio_%s",labelLR_energyBin.c_str()),Form("c_totRatio_%s",labelLR_energyBin.c_str()));
	      histo = h1_totRatio[index2];
	      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	      histo -> SetTitle(Form(";tot_{right} / tot_{left};entries"));
	      histo -> SetLineColor(kRed);
	      histo -> SetLineWidth(2);
	      histo -> Draw();
	      histo -> Write();
	      
	      fitFunc_totRatio[index2] = new TF1(Form("fitFunc_totRatio_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
	      histo -> Fit(fitFunc_totRatio[index2],"QNRS");
	      histo -> Fit(fitFunc_totRatio[index2],"QSR+","",fitFunc_totRatio[index2]->GetParameter(1)-2.*fitFunc_totRatio[index2]->GetParameter(2),fitFunc_totRatio[index2]->GetParameter(1)+2.*fitFunc_totRatio[index2]->GetParameter(2));
	      histo -> Fit(fitFunc_totRatio[index2],"QSR+","",fitFunc_totRatio[index2]->GetParameter(1)-2.*fitFunc_totRatio[index2]->GetParameter(2),fitFunc_totRatio[index2]->GetParameter(1)+2.*fitFunc_totRatio[index2]->GetParameter(2));
	      
	      fitFunc_totRatio[index2] -> SetLineColor(kBlack);
	      fitFunc_totRatio[index2] -> SetLineWidth(2);
	      fitFunc_totRatio[index2] -> Draw("same");
	      
	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");
	      
	      c -> Print(Form("%s/totRatio/c_totRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/totRatio/c_totRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete c;

	      // -- draw t1fine (average between left and right)
	      c = new TCanvas(Form("c_t1fineMean_%s",labelLR_energyBin.c_str()),Form("c_t1fineMean_%s",labelLR_energyBin.c_str()));
	      histo = h1_t1fineMean[index2];
	      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	      histo -> SetTitle(Form(";(t1fine_{right}+t1fine_{left})/2;entries"));
	      histo -> SetLineColor(kRed);
	      histo -> SetLineWidth(2);
	      histo -> Draw();
	      histo -> Write();
	      	      
	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");
	      
	      c -> Print(Form("%s/t1fine/c_t1fineMean__%s.png",plotDir.c_str(),labelLR.c_str()));
	      c -> Print(Form("%s/t1fine/c_t1fineMean__%s.pdf",plotDir.c_str(),labelLR.c_str()));
	      delete c;

		  // -- draw qT1 (average between left and right)
	      c = new TCanvas(Form("c_qT1Mean_%s",labelLR_energyBin.c_str()),Form("c_qT1Mean_%s",labelLR_energyBin.c_str()));
	      histo = h1_qT1Mean[index1];
	      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	      histo -> SetTitle(Form(";(qT1_{right}+qT1_{left})/2;entries"));
	      histo -> SetLineColor(kRed);
	      histo -> SetLineWidth(2);
	      histo -> Draw();
	      histo -> Write();
	      	      
	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");
	      
	      c -> Print(Form("%s/t1fine/c_qT1Mean__%s.png",plotDir.c_str(),labelLR.c_str()));
	      c -> Print(Form("%s/t1fine/c_qT1Mean__%s.pdf",plotDir.c_str(),labelLR.c_str()));
	      delete c;
	      
	      
	    } // --- end loop over energy bins
	}// --- end loop ober bars
    }// --- end loop over stepLabels


  //------------------------
  //--- 3rd loop over events
  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);
      
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) std::cout << ">>> 3rd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	  mapIt.second -> GetEntry(entry);
	  
	  bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;
	  
	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	  //int index( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + 99 ); 
	  
	  if( !accept[index1][entry] ) continue;

	  if(anEvent->amp_MCP < 500) continue;
	  
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  
	  double  index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	  
	  float energyRatioMean = fitFunc_energyRatio[index2]->GetParameter(1);
	  float energyRatioSigma = fitFunc_energyRatio[index2]->GetParameter(2);

	  
          float totRatioMean = fitFunc_totRatio[index2]->GetParameter(1);
	  float totRatioSigma = fitFunc_totRatio[index2]->GetParameter(2);
	  

	  if( fabs(anEvent->totR/anEvent->totL-totRatioMean) > 3.*totRatioSigma  ||  (anEvent->totR/anEvent->totL)>5 || (anEvent->totR/anEvent->totL)<0 )  
	    {
	      accept[index1][entry] = false;
	      continue;
	    }
	  
	  	  
	  float energyMean = 0.5*(anEvent->energyR + anEvent->energyL );			
	  if( !source.compare(TB) && energyMean < ranges["L-R"][index1]->at(0) )
	    {
	      accept[index1][entry] = false;
	      continue;
	    }	  
	  
	   long long deltaT = anEvent->timeR - anEvent->timeL;


	  if( h1_deltaT[index2] == NULL )
	    {
	      std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
		  //std::string labelLR_Ref(Form("refbarL-R_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1));
		  std::string labelL(Form("bar%02dL_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
		  std::string labelR(Form("bar%02dR_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
	      
	      h1_deltaT[index2] = new TH1F(Form("h1_deltaT_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
		  //h1_tavg[index1] = new TH1F(Form("h1_tavg_%s",labelLR_energyBin.c_str()),"",50e6,100e12,610e12);
	      p1_deltaT_vs_energyRatio[index2] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()),"",50,energyRatioMean-3.*energyRatioSigma,energyRatioMean+3.*energyRatioSigma);
	      p1_deltaT_vs_totRatio[index2] = new TProfile(Form("p1_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,totRatioMean-3.*totRatioSigma, totRatioMean+3.*totRatioSigma);
	      h2_deltaT_vs_totRatio[index2] = new TH2F(Form("h2_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,totRatioMean-3.*totRatioSigma, totRatioMean+3.*totRatioSigma, 2000, -12000., 12000.);

		  
		  p1_tL_tMCP_vs_energyL[index1] = new TProfile(Form("p1_tL_tMCP_vs_energyL_%s",labelL.c_str()),"",40,0,1000);
		  p1_tR_tMCP_vs_energyR[index1] = new TProfile(Form("p1_tR_tMCP_vs_energyR_%s",labelR.c_str()),"",40,0,1000);
		  //p1_tavgL_vs_tot[index1] = new TProfile(Form("p1_tL-tavg_vs_tot_%s",labelL.c_str()),"",180,4,40);
		  //p1_tavgR_vs_tot[index1] = new TProfile(Form("p1_tR-tavg_vs_tot_%s",labelR.c_str()),"",180,4,40);
		  //if(anEvent->barID == 8) p1_tavg_vs_totRatio[index] = new TProfile(Form("p1_tavg_vs_totRatio_%s",labelLR_Ref.c_str()),"",50,totRatioMean_ref-3.*totRatioSigma_ref, totRatioMean_ref+3.*totRatioSigma_ref);

		  
		  h2_tL_tMCP_vs_energyL[index1] = new TH2F(Form("h2_tL_tMCP_vs_energyL_%s",labelL.c_str()),"",35, 0.7, ranges["L"][index1]->at(1)/mpv_energy["L"][index1]->at(0), 100, fitFunc_tL_tMCP[index1]->GetParameter(1) - 5*fitFunc_tL_tMCP[index1]->GetParameter(2), fitFunc_tL_tMCP[index1]->GetParameter(1) + 5*fitFunc_tL_tMCP[index1]->GetParameter(2));
		  h2_tR_tMCP_vs_energyR[index1] = new TH2F(Form("h2_tR_tMCP_vs_energyR_%s",labelR.c_str()),"",35, 0.7, ranges["R"][index1]->at(1)/mpv_energy["R"][index1]->at(0), 100, fitFunc_tR_tMCP[index1]->GetParameter(1) - 5*fitFunc_tR_tMCP[index1]->GetParameter(2), fitFunc_tR_tMCP[index1]->GetParameter(1) + 5*fitFunc_tR_tMCP[index1]->GetParameter(2));
		  //h2_tavgL_vs_tot[index1] = new TH2F(Form("h2_tavgL_vs_tot_%s",labelL.c_str()),"",144, 4, 40, 400, 500, 5000);
		  //h2_tavgR_vs_tot[index1] = new TH2F(Form("h2_tavgR_vs_tot_%s",labelR.c_str()),"",144, 4, 40, 400, 500, 5000);

	      // p1_deltaT_vs_energyRatio[index2] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()),"",10,energyRatioMean-3.*energyRatioSigma,energyRatioMean+3.*energyRatioSigma);
	      // p1_deltaT_vs_totRatio[index2] = new TProfile(Form("p1_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),"",10,totRatioMean-3.*totRatioSigma, totRatioMean+3.*totRatioSigma);
	      // h2_deltaT_vs_totRatio[index2] = new TH2F(Form("h2_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),"",10,totRatioMean-3.*totRatioSigma, totRatioMean+3.*totRatioSigma, 1000, -10000., 10000.);
	    }
	  
	  if(fabs(deltaT)>10000) continue;
	  h1_deltaT[index2] -> Fill( deltaT ); 

          float timeLow = CTRMeans[index2] - 3.* CTRSigmas[index2];
	  float timeHig = CTRMeans[index2] + 3.* CTRSigmas[index2];
	  
	  if( ( deltaT > timeLow ) && ( deltaT < timeHig ) ){
	    p1_deltaT_vs_energyRatio[index2] -> Fill( anEvent->energyR/anEvent->energyL,deltaT );
            p1_deltaT_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL,deltaT );
            h2_deltaT_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL,deltaT );

			
			p1_tL_tMCP_vs_energyL[index1] -> Fill (anEvent->energyL, anEvent->phi_L - anEvent->phi_MCP_C);
			p1_tR_tMCP_vs_energyR[index1] -> Fill (anEvent->energyR, anEvent->phi_R - anEvent->phi_MCP_C);
			//p1_tavgL_vs_tot[index1] -> Fill (anEvent->totL, anEvent->timeL - tavg_totCorr);
			//p1_tavgR_vs_tot[index1] -> Fill (anEvent->totR, anEvent->timeR - tavg_totCorr);


			h2_tL_tMCP_vs_energyL[index1] -> Fill (anEvent->energyL/mpv_energy["L"][index1]->at(0), anEvent->phi_L - anEvent->phi_MCP_C);
			h2_tR_tMCP_vs_energyR[index1] -> Fill (anEvent->energyR/mpv_energy["R"][index1]->at(0), anEvent->phi_R - anEvent->phi_MCP_C);
			//h2_tavgL_vs_tot[index1] -> Fill (anEvent->totL, anEvent->timeL - tavg_totCorr);
			//h2_tavgR_vs_tot[index1] -> Fill (anEvent->totR, anEvent->timeR - tavg_totCorr);

          }
	}

      std::cout << std::endl;
    
    }
  

  //------------------
  //--- draw 3rd plots
  std::map<double,TF1*> fitFunc_energyRatioCorr;
  std::map<double,TF1*> fitFunc_totRatioCorr;
  std::map<double,TF1*> fitFunc_ampWalkCorrEnergyL;
  std::map<double,TF1*> fitFunc_ampWalkCorrEnergyR;
  std::map<double,TF1*> fitFunc_tL_tavg;
  std::map<double,TF1*> fitFunc_tR_tavg;
  
  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];

	  int index( (10000*int(Vov*100.)) + (100*vth1) + 99 ); 
	  std::string labelLR_Ref(Form("refbarL-R_%s",stepLabel.c_str()));
	
      for(int iBar = 0; iBar < 16; ++iBar){
	
        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;                                                                                       
	if (!barFound) continue;  

	std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
	std::string labelL(Form("bar%02dL_%s",iBar,stepLabel.c_str()));
	std::string labelR(Form("bar%02dR_%s",iBar,stepLabel.c_str()));
	
	int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
	

	if( !ranges["L-R"][index1] ) continue;
	
	int nEnergyBins = ranges["L-R"][index1]->size()-1;
	
	for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
	  {
	    
	//-------------------------------------------------------- amplitude walk correction plots ---------------------------------------------------------------
	//draw tL - tMCP vs energyL
	if(!p1_tL_tMCP_vs_energyL[index1]) continue;
	c = new TCanvas(Form("c_tL-tMCP_vs_energyL_%s",labelL.c_str()),Form("c_tL-tMCP_vs_energyL_%s",labelL.c_str()));
	    
	prof = p1_tL_tMCP_vs_energyL[index1];
	prof -> SetTitle(Form(";energy [a.u.]; #phi_{L} - #phi_{MCP} [clock units]"));
	prof -> GetYaxis() -> SetRangeUser(prof->GetMean(2)-0.7*prof->GetRMS(2), prof->GetMean(2)+0.7*prof->GetRMS(2));
	prof -> Draw("");
	    
	latex = new TLatex(0.40,0.85,Form("#splitline{bar %02dL}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	latex -> SetNDC();
	latex -> SetTextFont(42);
	latex -> SetTextSize(0.04);
	latex -> SetTextColor(kRed);
	latex -> Draw("same");

	fitFunc_ampWalkCorrEnergyL[index1] = new TF1(Form("fitFunc_ampWalkCorrEnergyL_%s",labelL.c_str()),"pol3",0,1000);
	//fitFunc_ampWalkCorrEnergyL[index1] -> SetRange(ranges["L"][index1]->at(0), ranges["L"][index1]->at(1));
	prof -> Fit(fitFunc_ampWalkCorrEnergyL[index1],"QRS+");
	fitFunc_ampWalkCorrEnergyL[index1] -> SetLineColor(kRed);
	fitFunc_ampWalkCorrEnergyL[index1] -> SetLineWidth(2);
	fitFunc_ampWalkCorrEnergyL[index1] -> Draw("same");
		
	outFile->cd();
	prof->Write();

	c -> Print(Form("%s/amplitudeWalkCorr/c_tL-tMCP_vs_energyL_%s.png",plotDir.c_str(),labelL.c_str()));
	c -> Print(Form("%s/amplitudeWalkCorr/c_tL-tMCP_vs_energyL_%s.pdf",plotDir.c_str(),labelL.c_str()));
	delete c;

	if(!h2_tL_tMCP_vs_energyL[index1]) continue;
	c = new TCanvas(Form("c_h2_tL-tMCP_vs_energyL_%s",labelL.c_str()),Form("c_h2_tL-tMCP_vs_energyL_%s",labelL.c_str()));
	    
	histo2 = h2_tL_tMCP_vs_energyL[index1];
	histo2 -> SetTitle(Form(";energy [a.u.]; #phi_{L} - #phi_{MCP} [clock units]"));
	histo2 -> GetYaxis() -> SetRangeUser(histo2->GetMean(2)-3.5*histo2->GetRMS(2), histo2->GetMean(2)+3.5*histo2->GetRMS(2));
	histo2 -> Draw("colz");
	    
	latex = new TLatex(0.40,0.85,Form("#splitline{bar %02dL}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	latex -> SetNDC();
	latex -> SetTextFont(42);
	latex -> SetTextSize(0.04);
	latex -> SetTextColor(kRed);
	latex -> Draw("same");

	fitFunc_ampWalkCorrEnergyL[index1] = new TF1(Form("fitFunc_ampWalkCorrEnergyL_%s",labelL.c_str()),"pol3",0.8,ranges["L"][index1]->at(1)/mpv_energy["L"][index1]->at(0));
	//fitFunc_ampWalkCorrEnergyL[index1] -> SetRange( ranges["L"][index1]->at(0), ranges["L"][index1]->at(1));
	histo2 -> Fit(fitFunc_ampWalkCorrEnergyL[index1],"QRS+");
	fitFunc_ampWalkCorrEnergyL[index1] -> SetLineColor(kRed);
	fitFunc_ampWalkCorrEnergyL[index1] -> SetLineWidth(2);
	fitFunc_ampWalkCorrEnergyL[index1] -> Draw("same");
		
	outFile->cd();
	histo2->Write();

	c -> Print(Form("%s/amplitudeWalkCorr/c_h2_tL-tMCP_vs_energyL_%s.png",plotDir.c_str(),labelL.c_str()));
	c -> Print(Form("%s/amplitudeWalkCorr/c_h2_tL-tMCP_vs_energyL_%s.pdf",plotDir.c_str(),labelL.c_str()));
	delete c;
	
	//draw tR - tMCP vs energyR
	if(!p1_tR_tMCP_vs_energyR[index1]) continue;
	c = new TCanvas(Form("c_tR-tMCP_vs_energyR_%s",labelR.c_str()),Form("c_tR-tMCP_vs_energyR_%s",labelL.c_str()));
	    
	prof = p1_tR_tMCP_vs_energyR[index1];
	prof -> SetTitle(Form(";energy [a.u.]; #phi_{R} - #phi_{MCP} [clock units]"));
	prof -> GetYaxis() -> SetRangeUser(prof->GetMean(2)-0.7*prof->GetRMS(2), prof->GetMean(2)+0.7*prof->GetRMS(2));
	prof -> Draw("");
	    
	latex = new TLatex(0.40,0.85,Form("#splitline{bar %02dL}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	latex -> SetNDC();
	latex -> SetTextFont(42);
	latex -> SetTextSize(0.04);
	latex -> SetTextColor(kRed);
	latex -> Draw("same");

	fitFunc_ampWalkCorrEnergyR[index1] = new TF1(Form("fitFunc_ampWalkCorrEnergyR_%s",labelR.c_str()),"pol3",0,1000);
	fitFunc_ampWalkCorrEnergyR[index1] -> SetRange( ranges["R"][index1]->at(0), ranges["R"][index1]->at(1));
	prof -> Fit(fitFunc_ampWalkCorrEnergyR[index1],"QRS+");
	fitFunc_ampWalkCorrEnergyR[index1] -> SetLineColor(kRed);
	fitFunc_ampWalkCorrEnergyR[index1] -> SetLineWidth(2);
	fitFunc_ampWalkCorrEnergyR[index1] -> Draw("same");
		
	c -> Print(Form("%s/amplitudeWalkCorr/c_tR-tMCP_vs_energyR_%s.png",plotDir.c_str(),labelL.c_str()));
	c -> Print(Form("%s/amplitudeWalkCorr/c_tR-tMCP_vs_energyR_%s.pdf",plotDir.c_str(),labelL.c_str()));
	delete c;

	outFile->cd();
	prof->Write();


	if(!h2_tR_tMCP_vs_energyR[index1]) continue;
	c = new TCanvas(Form("c_h2_tR-tMCP_vs_energyR_%s",labelL.c_str()),Form("c_h2_tR-tMCP_vs_energyR_%s",labelL.c_str()));
	    
	histo2 = h2_tR_tMCP_vs_energyR[index1];
	histo2 -> SetTitle(Form(";energy [a.u.]; #phi_{R} - #phi_{MCP} [clock units]"));
	histo2 -> GetYaxis() -> SetRangeUser(histo2->GetMean(2)-3.5*histo2->GetRMS(2), histo2->GetMean(2)+3.5*histo2->GetRMS(2));
	histo2 -> Draw("colz");
	    
	latex = new TLatex(0.40,0.85,Form("#splitline{bar %02dL}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	latex -> SetNDC();
	latex -> SetTextFont(42);
	latex -> SetTextSize(0.04);
	latex -> SetTextColor(kRed);
	latex -> Draw("same");

	fitFunc_ampWalkCorrEnergyR[index1] = new TF1(Form("fitFunc_ampWalkCorrEnergyR_%s",labelL.c_str()),"pol3",0.8,ranges["R"][index1]->at(1)/mpv_energy["R"][index1]->at(0));
	//fitFunc_ampWalkCorrEnergyR[index1] -> SetRange( ranges["R"][index1]->at(0), ranges["R"][index1]->at(1));
	histo2 -> Fit(fitFunc_ampWalkCorrEnergyR[index1],"QRS+");
	fitFunc_ampWalkCorrEnergyR[index1] -> SetLineColor(kRed);
	fitFunc_ampWalkCorrEnergyR[index1] -> SetLineWidth(2);
	fitFunc_ampWalkCorrEnergyR[index1] -> Draw("same");
		
	outFile->cd();
	histo2->Write();

	c -> Print(Form("%s/amplitudeWalkCorr/c_h2_tR-tMCP_vs_energyR_%s.png",plotDir.c_str(),labelL.c_str()));
	c -> Print(Form("%s/amplitudeWalkCorr/c_h2_tR-tMCP_vs_energyR_%s.pdf",plotDir.c_str(),labelL.c_str()));
	delete c;

//--------------------------- end amplitude correction plots ----------------------------------------------------------

	    double  index2( 10000000*iEnergyBin+index1 );
	    if(!p1_deltaT_vs_energyRatio[index2]) continue;
	    
	    std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));
	    
	    // -- draw deltaT vs energy ratio
	    c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()),Form("c_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()));
	    
	    prof = p1_deltaT_vs_energyRatio[index2];
	    prof -> SetTitle(Form(";energy_{right} / energy_{left};#Deltat [ps]"));
	    prof -> GetYaxis() -> SetRangeUser(CTRMeans[index2]-3.*CTRSigmas[index2],CTRMeans[index2]+3.*CTRSigmas[index2]);
	    prof -> Draw("");
	    
	    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	    latex -> SetNDC();
	    latex -> SetTextFont(42);
	    latex -> SetTextSize(0.04);
	    latex -> SetTextColor(kRed);
	    latex -> Draw("same");
	    
	    float fitXMin = fitFunc_energyRatio[index2]->GetParameter(1) - 2.*fitFunc_energyRatio[index2]->GetParameter(2);
	    float fitXMax = fitFunc_energyRatio[index2]->GetParameter(1) + 2.*fitFunc_energyRatio[index2]->GetParameter(2);
	    
	    fitFunc_energyRatioCorr[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"pol4",fitXMin,fitXMax);
	    prof -> Fit(fitFunc_energyRatioCorr[index2],"QRS+");
	    fitFunc_energyRatioCorr[index2] -> SetLineColor(kRed);
	    fitFunc_energyRatioCorr[index2] -> SetLineWidth(2);
	    fitFunc_energyRatioCorr[index2] -> Draw("same");
	    
	    c -> Print(Form("%s/energyRatioCorr/c_deltaT_vs_energyRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	    c -> Print(Form("%s/energyRatioCorr/c_deltaT_vs_energyRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	    delete c;


	    // -- draw deltaT vs tot ratio
            if(!p1_deltaT_vs_totRatio[index2]) continue;

	    c = new TCanvas(Form("c_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),Form("c_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()));
	    
	    prof = p1_deltaT_vs_totRatio[index2];
	    prof -> SetTitle(Form(";ToT_{right} / ToT_{left};#Deltat [ps]"));
	    prof -> GetYaxis() -> SetRangeUser(CTRMeans[index2]-3.*CTRSigmas[index2],CTRMeans[index2]+3.*CTRSigmas[index2]);
	    prof -> Draw("");
	    
	    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	    latex -> SetNDC();
	    latex -> SetTextFont(42);
	    latex -> SetTextSize(0.04);
	    latex -> SetTextColor(kRed);
	    latex -> Draw("same");
	    
        fitXMin = fitFunc_totRatio[index2]->GetParameter(1) - 3.*fitFunc_totRatio[index2]->GetParameter(2);
	    fitXMax = fitFunc_totRatio[index2]->GetParameter(1) + 3.*fitFunc_totRatio[index2]->GetParameter(2);
	    
	    fitFunc_totRatioCorr[index2] = new TF1(Form("fitFunc_totRatioCorr_%s",labelLR_energyBin.c_str()),"pol3",fitXMin,fitXMax);
	    prof -> Fit(fitFunc_totRatioCorr[index2],"QRS+");
	    fitFunc_totRatioCorr[index2] -> SetLineColor(kRed);
	    fitFunc_totRatioCorr[index2] -> SetLineWidth(2);
	    fitFunc_totRatioCorr[index2] -> Draw("same");
	    
	    c -> Print(Form("%s/totRatioCorr/c_deltaT_vs_totRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	    c -> Print(Form("%s/totRatioCorr/c_deltaT_vs_totRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	    delete c;
	  }
      }
    }
  
  
  
  //------------------------
  //--- 4th loop over events
  
  gStyle->SetOptFit(1111);
  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);
      
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) std::cout << ">>> 4th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	  mapIt.second -> GetEntry(entry);
	  
	  bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	  
	  if( !accept[index1][entry] ) continue;
	  
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  double  index2( 10000000*energyBinAverage+index1 );   

	  if(anEvent->amp_MCP < 500) continue;  
	  
	  long long deltaT = anEvent->timeR - anEvent->timeL;
	  
	  float t1fineMean = 0.5 * ( anEvent->t1fineR + anEvent->t1fineL );
	  

	  if( !fitFunc_energyRatioCorr[index2] ) continue;
	  
	  float energyRatioCorr = fitFunc_energyRatioCorr[index2]->Eval(anEvent->energyR/anEvent->energyL) - 
                                  fitFunc_energyRatioCorr[index2]->Eval(fitFunc_energyRatio[index2]->GetParameter(1));

          
	  if( !fitFunc_totRatioCorr[index2] ) continue;
	  
	  float totRatioCorr = fitFunc_totRatioCorr[index2]->Eval(anEvent->totR/anEvent->totL) - 
	                       fitFunc_totRatioCorr[index2]->Eval(fitFunc_totRatio[index2]->GetParameter(1));

	 if(!fitFunc_ampWalkCorrEnergyL[index1]) continue;
	  float ampWalkCorrEnergyL = anEvent->phi_L - anEvent->phi_MCP_C - fitFunc_ampWalkCorrEnergyL[index1]->Eval(anEvent->energyL/mpv_energy["L"][index1]->at(0));
	
	  if(!fitFunc_ampWalkCorrEnergyR[index1]) continue;
	  float ampWalkCorrEnergyR = anEvent->phi_R - anEvent->phi_MCP_C - fitFunc_ampWalkCorrEnergyR[index1]->Eval(anEvent->energyR/mpv_energy["R"][index1]->at(0));

	  float tL_tR = (ampWalkCorrEnergyL - ampWalkCorrEnergyR)*6250 - floor(anEvent->timeR/6250)*6250 + floor(anEvent->timeL/6250)*6250 ;
//---------------------------------------------------------------------------------------------
	  
	  
	  if( h1_deltaT_energyRatioCorr[index2] == NULL )
	    {
	      std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.02f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
		  std::string labelLR_(Form("bar%02dL-R_Vov%.02f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
	      h1_deltaT_energyRatioCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
	      h1_deltaT_totRatioCorr[index2]    = new TH1F(Form("h1_deltaT_totRatioCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
              p1_deltaT_totRatioCorr_vs_totRatio[index2] = new TProfile(Form("p1_deltaT_totRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,fitFunc_totRatio[index2]->GetParameter(1)-3.*fitFunc_totRatio[index2]->GetParameter(2), fitFunc_totRatio[index2]->GetParameter(1)+3.*fitFunc_totRatio[index2]->GetParameter(2));
              h2_deltaT_totRatioCorr_vs_totRatio[index2] = new TH2F(Form("h2_deltaT_totRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,fitFunc_totRatio[index2]->GetParameter(1)-3.*fitFunc_totRatio[index2]->GetParameter(2), fitFunc_totRatio[index2]->GetParameter(1)+3.*fitFunc_totRatio[index2]->GetParameter(2), 2000, -12000., 12000. );
	      p1_deltaT_energyRatioCorr_vs_t1fineMean[index2] = new TProfile(Form("p1_deltaT_energyRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()),"",50,0,1000);
	      p1_deltaT_totRatioCorr_vs_t1fineMean[index2] = new TProfile(Form("p1_deltaT_totRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()),"",50,0,1000);

		  h1_deltaT_ampWalkEnergyCorr[index1] = new TH1F(Form("h1_deltaT_ampWalkEnergyCorr_%s",labelLR_.c_str()),"",2400,-2.,2.);
		  h1_deltaT_P[index1] = new TH1F(Form("h1_deltaT_P_%s",labelLR_.c_str()),"",2400,-2.,2.);
		  h1_deltaT_tConverted[index1]    = new TH1F(Form("h1_deltaT_tConverted_%s",labelLR_.c_str()),"",2000,-12000.,12000.);
		  p1_deltaT_ampWalkEnergyCorr_vs_posX[index1] = new TProfile(Form("p1_deltaT_ampWalkEnergyCorr_vs_posX_%s",labelLR_.c_str()),"",30,-10,20);
		  //h1_deltaT_ampWalkTotCorr[index1] = new TH1F(Form("h1_deltaT_ampWalkTotCorr_%s",labelLR_.c_str()),"",2000,-12000.,12000.);
	
	      //p1_deltaT_ampWalkTotCorr_vs_t1fineMean[index1] = new TProfile(Form("p1_deltaT_ampWalkTotCorr_vs_t1fineMean_%s",labelLR_.c_str()),"",50,0,1000);

	    }
	  
	  if (fabs(deltaT - energyRatioCorr)>10000 ) continue;
	  if (fabs(deltaT - totRatioCorr)>10000 ) continue;
	  h1_deltaT_energyRatioCorr[index2] -> Fill( deltaT  - energyRatioCorr );
	  h1_deltaT_totRatioCorr[index2] -> Fill( deltaT  - totRatioCorr );
	  p1_deltaT_totRatioCorr_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL, deltaT  - totRatioCorr );
	  if (fabs(deltaT - energyRatioCorr)< 2*h1_deltaT[index2]->GetRMS() ){
	    p1_deltaT_energyRatioCorr_vs_t1fineMean[index2] -> Fill( t1fineMean, deltaT - energyRatioCorr );
	    p1_deltaT_totRatioCorr_vs_t1fineMean[index2] -> Fill( t1fineMean, deltaT - totRatioCorr );

		h1_deltaT_ampWalkEnergyCorr[index1] -> Fill(ampWalkCorrEnergyL - ampWalkCorrEnergyR);
		//h1_deltaT_ampWalkTotCorr[index1] -> Fill(ampWalkCorrTotL - ampWalkCorrTotR);
		//h1_deltaT_ampWalkEnergyCorr[index1] -> Fill(tLampWalkCorr_energy - tRampWalkCorr_energy);
		h1_deltaT_P[index1] -> Fill (anEvent->phi_R - anEvent->phi_L);
		h1_deltaT_tConverted[index1]->Fill(tL_tR);
		//h1_deltaT_ampWalkTotCorr[index1] -> Fill(tLampWalkCorr_tot - tRampWalkCorr_tot);
		
		if (useTrackInfo && anEvent->nhits>0 && anEvent->x>-100){
	     p1_deltaT_ampWalkEnergyCorr_vs_posX[index1] ->Fill( anEvent->x, tL_tR);
	 } 
		
	  }
 	}		 
      std::cout << std::endl;
    }
  
  double  theIndex2;
  float enBin;
  float theDeltaTMean;
  float theEntriesCTR;
  float timeRes;
  float timeResToT;
  float timeResPhaseCorr;
  float errTimeRes;
  float errTimeResToT;
  float errTimeResPhaseCorr;
  float timeResEffSigma;
  
  //------------------
  //--- draw 4th plots
  std::map<double,TF1*> fitFunc3_posCorr;

  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];
      
      for(int iBar = 0; iBar < 16; ++iBar)
	{
	  bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	  if (!barFound) continue;

	  int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
	  std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));

	  if( !ranges["L-R"][index1] ) continue;
	    	  
	  int nEnergyBins = ranges["L-R"][index1]->size()-1;

	  for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
	    {
	      double  index2( 10000000*iEnergyBin + index1 );
	      
	      outTrees2[index2] = new TTree(Form("dataRes_bar%02dL-R_Vov%.2f_th%02.0f_enBin%02d",iBar,Vov,vth1,iEnergyBin),Form("dataRes_bar%02dL-R_Vov%.2f_th%02.0f_enBin%02d",iBar,Vov,vth1,iEnergyBin));
	      outTrees2[index2] -> Branch("energyBin", &enBin);
	      outTrees2[index2] -> Branch("timeResolution",&timeRes);
	      outTrees2[index2] -> Branch("timeResolutionToT",&timeResToT);
              outTrees2[index2] -> Branch("effSigma",&timeResEffSigma);
	      outTrees2[index2] -> Branch("errTimeResolution",&errTimeRes);
              outTrees2[index2] -> Branch("errTimeResolutionToT",&errTimeResToT);
	      outTrees2[index2] -> Branch("indexID2",&theIndex2);
	      outTrees2[index2] -> Branch("DeltaTMean",&theDeltaTMean);
	      outTrees2[index2] -> Branch("EntriesCTR",&theEntriesCTR);

		  	//------------ draw ampwalk corrected deltaT ---------------


	  if(!h1_deltaT_ampWalkEnergyCorr[index1]) continue;

	  c = new TCanvas(Form("c_deltaT_ampWalkEnergyCorr_%s",labelLR.c_str()),Form("c_deltaT_ampWalkEnergyCorr_%s",labelLR.c_str()));
	      
	      histo = h1_deltaT_ampWalkEnergyCorr[index1];
	      
	      histo -> SetTitle(Form(";amplitude walk corrected #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kBlue);
	      histo -> SetMarkerColor(kBlue);
	      histo -> Draw("");
	      
	      
	      FindSmallestInterval(vals,histo,0.68);
	      float min = vals[4];
	      float max = vals[5];
	      float delta = max-min;
	      float sigma = 0.5*delta;
	      float effSigma = sigma;

	      TF1* fitFunc = new TF1(Form("fitFunc_ampWalkEnergyCorr_%s",labelLR.c_str()),"gaus",-1, 1);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      //fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL");
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kBlue+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same"); 

		  latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.4f }{#sigma_{corr.}^{gaus} = %.4f }",effSigma,fabs(fitFunc->GetParameter(2))));
		  latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kBlue);
	      latex -> Draw("same");  

		  histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));      
	      
	      outFile -> cd();
	      histo -> Write();

		  // -- raw delta T (phase)
		  if(!h1_deltaT_P[index1]) continue;
          histo = h1_deltaT_P[index1];
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kRed);
	      histo -> SetMarkerColor(kRed);

	      c->cd();
	      histo -> Draw("same");
	      
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      
	      fitFunc = new TF1(Form("fitFunc_%s",labelLR.c_str()),"gaus",-1,1);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      //fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL");
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kRed+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");
	      
	      latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{raw}^{eff} = %.4f }{#sigma_{raw}^{gaus} = %.4f }",effSigma,fabs(fitFunc->GetParameter(2))));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
              c->cd();
	      latex -> Draw("same");

		  c -> Print(Form("%s/CTR_amplitudeWalkCorr/c_deltaT_ampWalkEnergyCorr__%s.pdf",plotDir.c_str(),labelLR.c_str()));
	      c -> Print(Form("%s/CTR_amplitudeWalkCorr/c_deltaT_ampWalkEnergyCorr__%s.png",plotDir.c_str(),labelLR.c_str()));

		  outFile -> cd();
	      histo -> Write();
		
	      delete c;

		  // draw same plots converted to ps
	      c2 = new TCanvas(Form("c_deltaT_tConverted_%s",labelLR.c_str()),Form("c_deltaT_tConverted_%s",labelLR.c_str()));
	      
	      histo = h1_deltaT_tConverted[index1];
              histo -> SetTitle(Form(";walk-corrected #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kBlue);
	      histo -> SetMarkerColor(kBlue);
              histo -> Draw("");
	      	      
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      	      
	      float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      fitFunc = new TF1(Form("fitFunc_tConverted_%s",labelLR.c_str()),"gaus",-10000, 10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL");
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.0*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.0*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL");
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+");
	      
	      fitFunc -> SetLineColor(kBlue+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");         
	      
	      outFile -> cd();
	      histo -> Write();

              outTrees2[index2]->Fill();

	      histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));


              latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
    
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlue);
              latex -> Draw("same");
              

	  
	      // -- raw delta T
              histo = h1_deltaT[index2];
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kRed);
	      histo -> SetMarkerColor(kRed);

	      c2->cd();
	      histo -> Draw("same");
	       
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      
	      fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      fitFunc = new TF1(Form("fitFunc_%s",labelLR.c_str()),"gaus",-10000,10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kRed+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");
	      
	      latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{raw}^{eff} = %.0f ps}{#sigma_{raw}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);

              c2->cd();
	      latex -> Draw("same");

	      outFile -> cd();
	      histo -> Write();

	      c2 -> Print(Form("%s/CTR_amplitudeWalkCorr/c_deltaT_tConverted__%s.pdf",plotDir.c_str(),labelLR.c_str()));
	      c2 -> Print(Form("%s/CTR_amplitudeWalkCorr/c_deltaT_tConverted__%s.png",plotDir.c_str(),labelLR.c_str()));
	      delete c2;


	// -- draw deltaT vs position
	      if (useTrackInfo) {
		c = new TCanvas(Form("c_deltaT_ampWalkEnergyPosCorr_vs_posX_%s",labelLR.c_str()),Form("c_deltaT_ampWalkEnergyPosCorr_vs_posX_%s",labelLR.c_str()));
		
		prof = p1_deltaT_ampWalkEnergyCorr_vs_posX[index1];
		prof -> SetTitle(Form("; x [mm] ;#Delta t [ps]"));
		prof -> GetYaxis() -> SetRangeUser(prof->GetMean(2)-2*prof->GetRMS(2), prof->GetMean(2)+2*prof->GetRMS(2));
		prof -> Draw("");
		
		latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		latex -> Draw("same");
		
		fitFunc3_posCorr[index1] = new TF1(Form("fitFunc3_posCorr_%s",labelLR.c_str()),"pol1",-50,50);
		fitFunc3_posCorr[index1] -> SetRange( prof -> GetMean()-1.7*prof->GetRMS(), prof -> GetMean()+1.7*prof->GetRMS());
		prof -> Fit(fitFunc3_posCorr[index1],"QRS+");
		fitFunc3_posCorr[index1] -> SetLineColor(kRed);
		fitFunc3_posCorr[index1] -> SetLineWidth(2);
		fitFunc3_posCorr[index1] -> Draw("same");
		
		c -> Print(Form("%s/positionCorr/c_deltaT_ampWalkEnergyCorr_vs_posX__%s.png",plotDir.c_str(),labelLR.c_str()));
		c -> Print(Form("%s/positionCorr/c_deltaT_ampWalkEnergyCorr_vs_posX_%s.pdf",plotDir.c_str(),labelLR.c_str()));
		delete c;
		  }

		
//------------------------------------------------------------------------------------------------------------------------------	      
	      
	      if(!h1_deltaT_energyRatioCorr[index2]) continue;
	      
	      std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));
	      
	      
	      // -- energy corr deltaT
	      c = new TCanvas(Form("c_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()));
	      
	      histo = h1_deltaT_energyRatioCorr[index2];
	      
	      histo -> SetTitle(Form(";energy-corrected #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kBlue);
	      histo -> SetMarkerColor(kBlue);
	      histo -> Draw("");
	      
	      
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      
	      
	      fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      fitFunc = new TF1(Form("fitFunc_energyCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      // non mi prende il range di fit se non faccio SetRange a mano...
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kBlue+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");         
	      
	      outFile -> cd();
	      histo -> Write();
	      
	      enBin = energyBin["L-R"][index1][iEnergyBin];
	      theIndex2 = index2;
	      timeRes = fabs(fitFunc->GetParameter(2));
	      timeResEffSigma = effSigma;
	      
	      theDeltaTMean = fitFunc->GetParameter(1);
	      theEntriesCTR = fitFunc->GetParameter(0);
	      errTimeRes = fabs(fitFunc->GetParError(2));
	      
	      int contr = 0;
	      if(!source.compare("Na22") && fitFunc->GetParameter(0)<20){
		timeRes = 0;
		effSigma = 0;
		contr = 1;
	      }
	      
              if (!source.compare("Laser")) histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-10.*fitFunc->GetParameter(2),
										fitFunc->GetParameter(1)+10.*fitFunc->GetParameter(2));
	      
	      histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));
	      
	      
	      latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
	      if(contr==1) latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,timeRes));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kBlue);
	      latex -> Draw("same");
	      
	      // -- totRatio corr deltaT
	      c2 = new TCanvas(Form("c_deltaT_totRatioCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioCorr_%s",labelLR_energyBin.c_str()));
	      
	      histo = h1_deltaT_totRatioCorr[index2];
              histo -> SetTitle(Form(";ToT-corrected #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kBlue);
	      histo -> SetMarkerColor(kBlue);
              histo -> Draw("");
	      	      
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      	      
	      fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      fitFunc = new TF1(Form("fitFunc_totCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL");
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.0*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.0*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL");
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+");
	      
	      fitFunc -> SetLineColor(kBlue+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");         
	      
	      outFile -> cd();
	      histo -> Write();

              timeResToT = fabs(fitFunc->GetParameter(2));
              errTimeResToT = fabs(fitFunc->GetParError(2));

              outTrees2[index2]->Fill();

	      if (!source.compare("Laser")) histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-10.*fitFunc->GetParameter(2),
										fitFunc->GetParameter(1)+10.*fitFunc->GetParameter(2));
	      
	      histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));


              latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
              if(contr==1) latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,timeRes));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlue);
              latex -> Draw("same");
              

	  
	      // -- raw delta T
              histo = h1_deltaT[index2];
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kRed);
	      histo -> SetMarkerColor(kRed);

	      c->cd();
	      histo -> Draw("same");

	      c2->cd();
	      histo -> Draw("same");
	      
	      
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      
	      fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      fitFunc = new TF1(Form("fitFunc_%s",labelLR_energyBin.c_str()),"gaus",-10000,10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kRed+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");
	      
	      latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{raw}^{eff} = %.0f ps}{#sigma_{raw}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
              c->cd();
	      latex -> Draw("same");
              c2->cd();
	      latex -> Draw("same");

	      outFile -> cd();
	      histo -> Write();

	      c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT_energyRatioCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT_energyRatioCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete c;

	      c2 -> Print(Form("%s/CTR_totRatioCorr/c_deltaT_energyRatioCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c2 -> Print(Form("%s/CTR_totRatioCorr/c_deltaT_energyRatioCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete c2;


	      // -- draw deltaT vs t1fine
	      if(!p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]) continue;
	      c = new TCanvas(Form("c_deltaT_energyRatioCorr__vs_t1fineMean_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()));
	      
	      prof = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2];
	      prof -> SetTitle(Form(";t1fineMean [ps];#Deltat [ps]"));
	      prof -> GetYaxis()->SetRangeUser(prof -> GetMean(2) -500., prof -> GetMean(2)+500);
	      prof -> Draw("pl");
	      
	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");
	      
	      c -> Print(Form("%s/phaseCorr/c_deltaT_energyRatioCorr_vs_t1fineMean__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/phaseCorr/c_deltaT_energyRatioCorr_vs_t1fineMean__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete c;

	      // -- draw deltaT vs t1fine
	      if(!p1_deltaT_totRatioCorr_vs_t1fineMean[index2]) continue;
	      c = new TCanvas(Form("c_deltaT_totRatioCorr__vs_t1fineMean_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()));
	      
	      prof = p1_deltaT_totRatioCorr_vs_t1fineMean[index2];
	      prof -> GetYaxis()->SetRangeUser(prof -> GetMean(2) -500., prof -> GetMean(2)+500);
	      prof -> SetTitle(Form(";t1fineMean [ps];#Deltat [ps]"));
	      prof -> Draw("pl");
	      
	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");
	      
	      c -> Print(Form("%s/phaseCorr/c_deltaT_totRatioCorr_vs_t1fineMean__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/phaseCorr/c_deltaT_totRatioCorr_vs_t1fineMean__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete c;
	      
	    }
	}
    }
  
  //------------------------
  //--- 5th loop over events
  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);
      
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
       {
	 
	 if( entry%100000 == 0 ) std::cout << ">>> 5th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	 mapIt.second -> GetEntry(entry);
	 
	 bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;                                                                             
	 if (!barFound) continue;

	 int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	 if( !accept[index1][entry] ) continue;
	 
	 int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	 double  index2( 10000000*energyBinAverage+index1 );     

	 if(anEvent->amp_MCP < 500) continue;
	 
	 long long deltaT = anEvent->timeR - anEvent->timeL;
	 	 
	 if( !fitFunc_energyRatioCorr[index2] )	continue;
	 if( !p1_deltaT_energyRatioCorr_vs_t1fineMean[index2] )	continue;
	 if( !p1_deltaT_ampWalkEnergyCorr_vs_posX[index1] )	continue;

	 float energyRatioCorr = fitFunc_energyRatioCorr[index2]->Eval(anEvent->energyR/anEvent->energyL) -
	                         fitFunc_energyRatioCorr[index2]->Eval(fitFunc_energyRatio[index2]->GetParameter(1));
	 
	 float t1fineMean = 0.5* ( anEvent->t1fineR + anEvent->t1fineL );
	 int t1fineBin = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	 int t1fineBin2 = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	 float t1fineCorr = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );

	 //----------------------------------------------------------------------------------------------------------------------
	  
	  if(!fitFunc_ampWalkCorrEnergyL[index1]) continue;
	  float ampWalkCorrEnergyL = anEvent->phi_L - anEvent->phi_MCP_C - fitFunc_ampWalkCorrEnergyL[index1]->Eval(anEvent->energyL/mpv_energy["L"][index1]->at(0));
	
	  if(!fitFunc_ampWalkCorrEnergyR[index1]) continue;
	  float ampWalkCorrEnergyR = anEvent->phi_R - anEvent->phi_MCP_C - fitFunc_ampWalkCorrEnergyR[index1]->Eval(anEvent->energyR/mpv_energy["R"][index1]->at(0));

	  float tL_tR = (ampWalkCorrEnergyL - ampWalkCorrEnergyR)*6250 - floor(anEvent->timeR/6250)*6250 + floor(anEvent->timeL/6250)*6250 ;

	  float posCorr1 = fitFunc3_posCorr[index1]->Eval(anEvent->x) - fitFunc3_posCorr[index1]->Eval( p1_deltaT_ampWalkEnergyCorr_vs_posX[index1]->GetMean());

	  float qT1Mean = 0.5* ( anEvent->qT1R + anEvent->qT1L );	  
	 //-------------------------------------------------------------------------------------------------------------------------	 
	 if( h1_deltaT_energyRatioPhaseCorr[index2] == NULL )
	   {
	     std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
		 std::string labelLR_(Form("bar%02dL-R_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
	     h1_deltaT_energyRatioPhaseCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioPhaseCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
	     p1_deltaT_energyRatioCorr_vs_posX[index2] = new TProfile(Form("p1_deltaT_energyRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()),"",15,-10,20);

		 h1_deltaT_ampWalkEnergyPosCorr[index1] = new TH1F(Form("h1_deltaT_ampWalkEnergyPosCorr_%s",labelLR_.c_str()),"",2000,-12000.,12000.);
		 p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1] = new TProfile(Form("p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean_%s",labelLR_.c_str()),"",50,0.4,1.6);
           }
	 
	 if (fabs(deltaT - energyRatioCorr)>10000 ) continue;
	 h1_deltaT_energyRatioPhaseCorr[index2] -> Fill( deltaT  - energyRatioCorr - t1fineCorr );
	 h1_deltaT_ampWalkEnergyPosCorr[index1] -> Fill( tL_tR - posCorr1);
	 p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1] -> Fill(qT1Mean, tL_tR - posCorr1);
	 if (useTrackInfo && anEvent->nhits>0 && anEvent->x>-100){
		 p1_deltaT_energyRatioCorr_vs_posX[index2] ->Fill( anEvent->x, deltaT - energyRatioCorr - t1fineCorr);
	     //p1_deltaT_ampWalkEnergyCorr_vs_posX[index1] ->Fill( anEvent->x, tL_tR - t1fineCorr1);
	 } 

         // totRatio corr
	 if( !fitFunc_totRatioCorr[index2] )	continue;
	 if( !p1_deltaT_totRatioCorr_vs_t1fineMean[index2] )	continue;

	  if( !fitFunc_totRatioCorr[index2] ) continue;
	  
	  float totRatioCorr = fitFunc_totRatioCorr[index2]->Eval(anEvent->totR/anEvent->totL) - 
	                       fitFunc_totRatioCorr[index2]->Eval(fitFunc_totRatio[index2]->GetParameter(1));

	 t1fineBin = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	 t1fineBin2 = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	 t1fineCorr = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );

	 if( h1_deltaT_totRatioPhaseCorr[index2] == NULL )
	   {
	     std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
		 std::string labelLR_(Form("bar%02dL-R_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));

	     h1_deltaT_totRatioPhaseCorr[index2] = new TH1F(Form("h1_deltaT_totRatioPhaseCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
	     p1_deltaT_totRatioCorr_vs_posX[index2] = new TProfile(Form("p1_deltaT_totRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()),"",50,-50,50);
		 
		 //h1_deltaT_ampWalkTotPhaseCorr[index1] = new TH1F(Form("h1_deltaT_ampWalkTotPhaseCorr_%s",labelLR_.c_str()),"",2000,-12000.,12000.);
		 //p1_deltaT_ampWalkTotCorr_vs_posX[index1] = new TProfile(Form("p1_deltaT_ampWalkTotCorr_vs_posX_%s",labelLR_.c_str()),"",50,-50,50);
           }
	 
	 if (fabs(deltaT - totRatioCorr)>10000 ) continue;
	 h1_deltaT_totRatioPhaseCorr[index2] -> Fill( deltaT  - totRatioCorr - t1fineCorr );
	 //h1_deltaT_ampWalkTotPhaseCorr[index1] -> Fill( tLampWalkCorr_tot - tRampWalkCorr_tot - t1fineCorr1);
	 if (useTrackInfo && anEvent->nhits>0 && anEvent->x>-100){
		 p1_deltaT_totRatioCorr_vs_posX[index2] ->Fill( anEvent->x, deltaT  - totRatioCorr - t1fineCorr);
		 //p1_deltaT_ampWalkTotCorr_vs_posX[index1] ->Fill( anEvent->x, tLampWalkCorr_tot - tRampWalkCorr_tot- t1fineCorr1);

	 } 
       }
      std::cout << std::endl;
    }
  
  //------------------
  //--- draw 5th plots
  std::map<double,TF1*> fitFunc1_posCorr;
  std::map<double,TF1*> fitFunc2_posCorr;
  std::map<double,TF1*> fitFunc4_posCorr;

  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];
      
      for(int iBar = 0; iBar < 16; ++iBar)
	{
	  bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	  if (!barFound) continue; 

	  std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));

	  int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );

	  if( !ranges["L-R"][index1] ) continue;

	//-------------- draw amplitude walk and position corrected delta T -----------------------------
	  if(!h1_deltaT_ampWalkEnergyPosCorr[index1]) continue;
	  
	  c = new TCanvas(Form("c_deltaT_ampWalkEnergyPosCorr_%s",labelLR.c_str()),Form("c_deltaT_ampWalkEnergyPosCorr_%s",labelLR.c_str()));
	    
	      histo = h1_deltaT_ampWalkEnergyPosCorr[index1];
	      
	      histo -> SetTitle(Form("; #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kGreen+1);
	      histo -> SetMarkerColor(kGreen+1);
	      histo -> Draw("");
	      
	      
	      FindSmallestInterval(vals,histo,0.68);
	      float min = vals[4];
	      float max = vals[5];
	      float delta = max-min;
	      float sigma = 0.5*delta;
	      float effSigma = sigma;

		  float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
		  float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      TF1* fitFunc = new TF1(Form("fitFunc_energyWalkPosCorr_%s",labelLR.c_str()),"gaus",-1000, 1000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL","",fitXMin,fitXMax);
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kGreen+2);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");         
	      
	      histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));
	      
	      latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corrPh.}^{eff} = %.0f ps}{#sigma_{corrPh.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
	      //if(contr==1) latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corrPh.}^{eff} = %.0f ps}{#sigma_{corrPh.}^{gaus} = %.0f ps}",effSigma,timeResPhaseCorr));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kGreen+2);
	      latex -> Draw("same");
	      
	      outFile -> cd();
	      histo -> Write();

		  c -> Print(Form("%s/CTR_amplitudeWalkPosCorr/c_deltaT_ampWalkEnergyPosCorr__%s.pdf",plotDir.c_str(),labelLR.c_str()));
	      c -> Print(Form("%s/CTR_amplitudeWalkPosCorr/c_deltaT_ampWalkEnergyPosCorr__%s.png",plotDir.c_str(),labelLR.c_str()));
	      delete c;

		    //draw ampwalk corrected deltaT vs qT1Mean
		  if(!p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1]) continue;
	      c = new TCanvas(Form("c_deltaT_ampWalkEnergyCorr_vs_t1fineMean_%s",labelLR.c_str()),Form("c_deltaT_ampWalkEnergyCorr_vs_t1fineMean_%s",labelLR.c_str()));
	      
	      prof = p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1];
	      prof -> SetTitle(Form(";t1fineMean [clock units];#Deltat [ps]"));
	      prof -> GetYaxis()->SetRangeUser(prof -> GetMean(2) -500., prof -> GetMean(2)+500);
	      prof -> Draw("pl");
	      
	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");
	      
	      c -> Print(Form("%s/phaseCorr/c_deltaT_ampWalkEnergyCorr_vs_t1fineMean__%s.png",plotDir.c_str(),labelLR.c_str()));
	      c -> Print(Form("%s/phaseCorr/c_deltaT_ampWalkEnergyCorr_vs_t1fineMean__%s.pdf",plotDir.c_str(),labelLR.c_str()));
	      delete c;

	//---------------------------------------------------------------------------------------------------------------------------------------

	  
	  int nEnergyBins = ranges["L-R"][index1]->size()-1;
	  
	  for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
	    {
	      double  index2( 10000000*iEnergyBin + index1 );
	      	      
	      if(!h1_deltaT_energyRatioPhaseCorr[index2]) continue;
	      
	      std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));
	      
	      std::cout << labelLR_energyBin.c_str()<<std::endl;
	      

	      c = new TCanvas(Form("c_deltaT_energyRatioPhaseCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioPhaseCorr_%s",labelLR_energyBin.c_str()));
	      
	      // -- energy and phase corr deltaT
	      histo = h1_deltaT_energyRatioPhaseCorr[index2];
	      
	      histo -> SetTitle(Form(";phase-corrected #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kGreen+1);
	      histo -> SetMarkerColor(kGreen+1);
	      
	      histo -> Draw("");
	      
	      
	      FindSmallestInterval(vals,histo,0.68);
	      float min = vals[4];
	      float max = vals[5];
	      float delta = max-min;
	      float sigma = 0.5*delta;
	      float effSigma = sigma;
	      
	      
	      float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      TF1* fitFunc = new TF1(Form("fitFunc_phaseCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kGreen+2);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");         
	      
	      
	      outFile -> cd();
	      histo -> Write();
	      
	      timeResPhaseCorr = fabs(fitFunc->GetParameter(2));
	      errTimeResPhaseCorr = fabs(fitFunc->GetParError(2));
	      
	      int contr = 0;
	      if(!source.compare("Na22") && fitFunc->GetParameter(0)<20){
		timeResPhaseCorr = 0;
		effSigma = 0;
		contr = 1;
	      }
	      
	      	      
	      if (!source.compare("Laser")) histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-10.*fitFunc->GetParameter(2),
										fitFunc->GetParameter(1)+10.*fitFunc->GetParameter(2));
	      
	      histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));
	      
	      
	      latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corrPh.}^{eff} = %.0f ps}{#sigma_{corrPh.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
	      if(contr==1) latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corrPh.}^{eff} = %.0f ps}{#sigma_{corrPh.}^{gaus} = %.0f ps}",effSigma,timeResPhaseCorr));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kGreen+2);
	      latex -> Draw("same");
	      
	      // -- energy corr deltaT
	      histo = h1_deltaT_energyRatioCorr[index2];
	      
	      histo -> SetTitle(Form(";energy-corrected #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kBlue);
	      histo -> SetMarkerColor(kBlue);
	      
	      histo -> Draw("same");
	      
	      
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      
	      
	      fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      fitFunc = new TF1(Form("fitFunc_energyCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kBlue+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");         
	      
	      
	      contr = 0;
	      if(!source.compare("Na22") && fitFunc->GetParameter(0)<20){
		timeRes = 0;
		effSigma = 0;
		contr = 1;
	      }
	      
	      
	      
	      latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{corrEn.}^{eff} = %.0f ps}{#sigma_{corrEn.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
	      if(contr==1) latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{corrEn.}^{eff} = %.0f ps}{#sigma_{corrEn.}^{gaus} = %.0f ps}",effSigma,timeRes));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kBlue);
	      latex -> Draw("same");
	      	      
	      
	      c -> Print(Form("%s/CTR_energyRatioPhaseCorr/c_deltaT_energyRatioPhaseCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/CTR_energyRatioPhaseCorr/c_deltaT_energyRatioPhaseCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete c;



              // -- tot and phase corr deltaT
	      if (!h1_deltaT_totRatioPhaseCorr[index2]) continue;

	      c = new TCanvas(Form("c_deltaT_totRatioPhaseCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioPhaseCorr_%s",labelLR_energyBin.c_str()));
	      
              histo = h1_deltaT_totRatioPhaseCorr[index2];
	      
	      histo -> SetTitle(Form(";phase-corrected #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kGreen+1);
	      histo -> SetMarkerColor(kGreen+1);
	      
	      histo -> Draw("");
	      
	      
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      	      
	      fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      fitFunc = new TF1(Form("fitFunc_phaseCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kGreen+2);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");         
	      
	      
	      outFile -> cd();
	      histo -> Write();
	      
              if (!source.compare("Laser")) histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-10.*fitFunc->GetParameter(2),
										fitFunc->GetParameter(1)+10.*fitFunc->GetParameter(2));
	      
	      histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	      histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));
	      
	      
	      latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corrPh.}^{eff} = %.0f ps}{#sigma_{corrPh.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
              latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kGreen+2);
	      latex -> Draw("same");
	      
              // -- tot corr deltaT
	      histo = h1_deltaT_totRatioCorr[index2];
	      
	      histo -> SetTitle(Form(";tot-corrected #Deltat [ps];entries"));
	      histo -> SetLineWidth(2);
	      histo -> SetLineColor(kBlue);
	      histo -> SetMarkerColor(kBlue);
	      
	      histo -> Draw("same");
	      
	      
	      FindSmallestInterval(vals,histo,0.68);
	      min = vals[4];
	      max = vals[5];
	      delta = max-min;
	      sigma = 0.5*delta;
	      effSigma = sigma;
	      
	      
	      fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
	      fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
	      
	      fitFunc = new TF1(Form("fitFunc_totCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
	      fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
	      fitFunc -> SetRange(fitXMin, fitXMax);
	      histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
	      fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
	      
	      fitFunc -> SetLineColor(kBlue+1);
	      fitFunc -> SetLineWidth(3);
	      fitFunc -> Draw("same");         
	      	      
	      latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{corrToT.}^{eff} = %.0f ps}{#sigma_{corrToT.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
              latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kBlue);
	      latex -> Draw("same");
	      	      
	      
	      c -> Print(Form("%s/CTR_totRatioPhaseCorr/c_deltaT_totRatioPhaseCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/CTR_totRatioPhaseCorr/c_deltaT_totRatioPhaseCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete c;


	      // -- draw deltaT vs position
	      if (useTrackInfo) {
		c = new TCanvas(Form("c_deltaT_energyRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()));
		
		prof = p1_deltaT_energyRatioCorr_vs_posX[index2];
		prof -> SetTitle(Form("; x [mm] ;#Deltat [ps]"));
		prof -> GetYaxis() -> SetRangeUser(prof->GetMean(2)-3*prof->GetRMS(2), prof->GetMean(2)+3*prof->GetRMS(2));
		prof -> Draw("");
		
		latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		latex -> Draw("same");
		
		fitFunc1_posCorr[index2] = new TF1(Form("fitFunc1_posCorr_%s",labelLR_energyBin.c_str()),"pol1",-50,50);
		fitFunc1_posCorr[index2] -> SetRange( prof -> GetMean()-3*prof->GetRMS(), prof -> GetMean()+3*prof->GetRMS());
		prof -> Fit(fitFunc1_posCorr[index2],"QRS+");
		fitFunc1_posCorr[index2] -> SetLineColor(kRed);
		fitFunc1_posCorr[index2] -> SetLineWidth(2);
		fitFunc1_posCorr[index2] -> Draw("same");
		
		c -> Print(Form("%s/positionCorr/c_deltaT_energyRatioCorr_vs_posX__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
		c -> Print(Form("%s/positionCorr/c_deltaT_energyRatioCorr_vs_posX_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
		delete c;



		
		c = new TCanvas(Form("c_deltaT_totRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()));
		
		prof = p1_deltaT_totRatioCorr_vs_posX[index2];
		prof -> SetTitle(Form("; x [mm] ;#Deltat [ps]"));
		prof -> GetYaxis() -> SetRangeUser(prof->GetMean(2)-3*prof->GetRMS(2), prof->GetMean(2)+3*prof->GetRMS(2));
		prof -> Draw("");
		
		latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		latex -> Draw("same");
		
		fitFunc2_posCorr[index2] = new TF1(Form("fitFunc2_posCorr_%s",labelLR_energyBin.c_str()),"pol1",-50,50);
		fitFunc2_posCorr[index2] -> SetRange( prof -> GetMean()-2*prof->GetRMS(), prof -> GetMean()+2*prof->GetRMS());
		prof -> Fit(fitFunc2_posCorr[index2],"QRS+");
		fitFunc2_posCorr[index2] -> SetLineColor(kRed);
		fitFunc2_posCorr[index2] -> SetLineWidth(2);
		fitFunc2_posCorr[index2] -> Draw("same");
		
		c -> Print(Form("%s/positionCorr/c_deltaT_totRatioCorr_vs_posX__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
		c -> Print(Form("%s/positionCorr/c_deltaT_totRatioCorr_vs_posX_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
		delete c;
	      }
	    }
	}
    }// end draw 5th plots


  //------------------------
  //--- 6th loop over events
  if (useTrackInfo){
    for(auto mapIt : trees)
      {
	ModuleEventClass* anEvent = new ModuleEventClass();
	mapIt.second -> SetBranchAddress("event",&anEvent);
	
	int nEntries = mapIt.second->GetEntries();
	for(int entry = 0; entry < nEntries; ++entry)
	  {
	    
	    if( entry%100000 == 0 ) std::cout << ">>> 6th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	    mapIt.second -> GetEntry(entry);
	    
	    bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;                                                                             
	    if (!barFound) continue;
	    
	    int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	    if( !accept[index1][entry] ) continue;
	    
	    int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	    double  index2( 10000000*energyBinAverage+index1 );     
	    
		if(anEvent->amp_MCP < 500) continue;

	    long long deltaT = anEvent->timeR - anEvent->timeL;
	    
	    if( !fitFunc_energyRatioCorr[index2] )	continue;
	    if( !p1_deltaT_energyRatioCorr_vs_t1fineMean[index2] )	continue;
	    if( !fitFunc1_posCorr[index2] )	continue;

	 	if( !p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1] )	continue;
		if( !fitFunc3_posCorr[index1] )	continue;
	    
	    float energyRatioCorr = fitFunc_energyRatioCorr[index2]->Eval(anEvent->energyR/anEvent->energyL) -
	                            fitFunc_energyRatioCorr[index2]->Eval(fitFunc_energyRatio[index2]->GetParameter(1));

	  
	  if(!fitFunc_ampWalkCorrEnergyL[index1]) continue;
	  float ampWalkCorrEnergyL = anEvent->phi_L - anEvent->phi_MCP_C - fitFunc_ampWalkCorrEnergyL[index1]->Eval(anEvent->energyL/mpv_energy["L"][index1]->at(0));
	
	  if(!fitFunc_ampWalkCorrEnergyR[index1]) continue;
	  float ampWalkCorrEnergyR = anEvent->phi_R - anEvent->phi_MCP_C - fitFunc_ampWalkCorrEnergyR[index1]->Eval(anEvent->energyR/mpv_energy["R"][index1]->at(0));

	  float tL_tR = (ampWalkCorrEnergyL - ampWalkCorrEnergyR)*6250 - floor(anEvent->timeR/6250)*6250 + floor(anEvent->timeL/6250)*6250 ;

		
		float qT1Mean = 0.5* ( anEvent->qT1R + anEvent->qT1L );
	   int qT1Bin = p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1]->FindBin(qT1Mean);
	   int qT1Bin2 = p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1]->FindBin(h1_qT1Mean[index1]->GetMean());
	   float qT1Corr = p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1]->GetBinContent(qT1Bin) - p1_deltaT_ampWalkEnergyCorr_vs_t1fineMean[index1]->GetBinContent( qT1Bin2 );
	   
	    float posCorr1 = fitFunc3_posCorr[index1]->Eval(anEvent->x) - fitFunc3_posCorr[index1]->Eval( p1_deltaT_ampWalkEnergyCorr_vs_posX[index1]->GetMean());

	    float t1fineMean = 0.5* ( anEvent->t1fineR + anEvent->t1fineL );
	    int t1fineBin = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	    int t1fineBin2 = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	    float t1fineCorr = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );
	    
	    float posCorr = fitFunc1_posCorr[index2]->Eval(anEvent->x) - fitFunc1_posCorr[index2]->Eval( p1_deltaT_energyRatioCorr_vs_posX[index2]->GetMean());
	    
	    if( h1_deltaT_energyRatioPhasePosCorr[index2] == NULL )
	      {
		std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
		std::string labelLR_(Form("bar%02dL-R_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));

		h1_deltaT_energyRatioPhasePosCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
		h1_deltaT_ampWalkEnergyPhasePosCorr[index1] = new TH1F(Form("h1_deltaT_ampWalkEnergyPhasePosCorr_%s",labelLR_.c_str()),"",2000,-12000.,12000.);
	      }
	    
	    if (fabs(deltaT - energyRatioCorr)>10000 ) continue;
	    h1_deltaT_energyRatioPhasePosCorr[index2] -> Fill( deltaT  - energyRatioCorr - t1fineCorr - posCorr);
		h1_deltaT_ampWalkEnergyPhasePosCorr[index1] -> Fill( tL_tR - posCorr1 - qT1Corr);
	    
	    // totRatio corr
	    if( !fitFunc_totRatioCorr[index2] )	continue;
	    if( !p1_deltaT_totRatioCorr_vs_t1fineMean[index2] )	continue;
	    if( !fitFunc2_posCorr[index2] )     continue;

	    
	    float totRatioCorr = fitFunc_totRatioCorr[index2]->Eval(anEvent->totR/anEvent->totL) -
	      fitFunc_totRatioCorr[index2]->Eval(fitFunc_totRatio[index2]->GetParameter(1));
	    
	    t1fineBin = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	    t1fineBin2 = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	    t1fineCorr = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );
	    
	    posCorr = fitFunc2_posCorr[index2]->Eval(anEvent->x) - fitFunc2_posCorr[index2]->Eval( p1_deltaT_totRatioCorr_vs_posX[index2]->GetMean());

	    if( h1_deltaT_totRatioPhasePosCorr[index2] == NULL )
	      {
		std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
		std::string labelLR_(Form("bar%02dL-R_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
		h1_deltaT_totRatioPhasePosCorr[index2] = new TH1F(Form("h1_deltaT_totRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
		//h1_deltaT_ampWalkTotPhasePosCorr[index1] = new TH1F(Form("h1_deltaT_ampWalkTotPhasePosCorr_%s",labelLR_.c_str()),"",2000,-12000.,12000.);
	      }
	    
	    if (fabs(deltaT - totRatioCorr)>10000 ) continue;
	    h1_deltaT_totRatioPhasePosCorr[index2] -> Fill( deltaT  - totRatioCorr - t1fineCorr - posCorr);
		//h1_deltaT_ampWalkTotPhasePosCorr[index1] -> Fill(tLampWalkCorr_tot - tRampWalkCorr_tot - t1fineCorr1 - posCorr1);
	  }
	std::cout << std::endl;
      }
  
    
    //------------------
    // draw 6th plots

    for(auto stepLabel : stepLabels)
      {
	float Vov = map_Vovs[stepLabel];
	float vth1 = map_ths[stepLabel];
	
	for(int iBar = 0; iBar < 16; ++iBar)
	  {
	    bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	    if (!barFound) continue;
	    
	    std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
	    
	    int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );

		if( !ranges["L-R"][index1] ) continue;

		//----------------------------------------------------------------------------------------------------------------------------
		if(!h1_deltaT_ampWalkEnergyPhasePosCorr[index1]) continue;
		
		c = new TCanvas(Form("c_deltaT_ampWalkEnergyPhasePosCorr_%s",labelLR.c_str()),Form("c_deltaT_ampWalkEnergyPhasePosCorr_%s",labelLR.c_str()));

		// -- energy, position and phase corr deltaT
		histo = h1_deltaT_ampWalkEnergyPhasePosCorr[index1];

		histo -> SetTitle(Form(";corrected #Deltat [ps];entries"));
		histo -> SetLineWidth(2);
		histo -> SetLineColor(kGreen+1);
		histo -> SetMarkerColor(kGreen+1);
		histo -> Draw("");

		FindSmallestInterval(vals,histo,0.68);
		float min = vals[4];
		float max = vals[5];
		float delta = max-min;
		float sigma = 0.5*delta;
		float effSigma = sigma;

		float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
		float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;

		TF1* fitFunc = new TF1(Form("fitFunc_ampWalkEnergyPhasePosCorr_%s",labelLR.c_str()),"gaus",-1000, 1000);
		fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
		fitFunc -> SetRange(fitXMin, fitXMax);
		histo -> Fit(fitFunc,"QNRSL","",fitXMin,fitXMax);
		fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
		histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
		fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
		histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));

		fitFunc -> SetLineColor(kMagenta);
		fitFunc -> SetLineWidth(3);
		fitFunc -> Draw("same");

		outFile -> cd();
		histo -> Write();

		histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
		histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));

		latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr}^{eff} = %.0f ps}{#sigma_{corr}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kMagenta);
		latex -> Draw("same");

		c -> Print(Form("%s/CTR_amplitudeWalkPhasePosCorr/c_deltaT_ampWalkEnergyPhasePosCorr__%s.pdf",plotDir.c_str(),labelLR.c_str()));
		c -> Print(Form("%s/CTR_amplitudeWalkPhasePosCorr/c_deltaT_ampWalkEnergyPhasePosCorr__%s.png",plotDir.c_str(),labelLR.c_str()));
		delete c;

//---------------------------------------------------------------------------------------------------------

	    
	    int nEnergyBins = ranges["L-R"][index1]->size()-1;
	    
	    for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
	      {
		double  index2( 10000000*iEnergyBin + index1 );
		
		if(!h1_deltaT_energyRatioPhasePosCorr[index2]) continue;
		
		std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));
	      
		std::cout << labelLR_energyBin.c_str()<<std::endl;


		c = new TCanvas(Form("c_deltaT_energyRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioPhasePosCorr_%s",labelLR_energyBin.c_str()));

		// -- energy phase and position corr deltaT
		histo = h1_deltaT_energyRatioPhasePosCorr[index2];

		histo -> SetTitle(Form(";corrected #Deltat [ps];entries"));
		histo -> SetLineWidth(2);
		histo -> SetLineColor(kGreen+1);
		histo -> SetMarkerColor(kGreen+1);
		histo -> Draw("");

		FindSmallestInterval(vals,histo,0.68);
		min = vals[4];
		max = vals[5];
		delta = max-min;
		sigma = 0.5*delta;
		effSigma = sigma;

		fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
		fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
		
		fitFunc = new TF1(Form("fitFunc_posCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
		fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
		fitFunc -> SetRange(fitXMin, fitXMax);
		histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
		fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
		histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
		fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
		histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));

		fitFunc -> SetLineColor(kMagenta);
		fitFunc -> SetLineWidth(3);
		fitFunc -> Draw("same");

		outFile -> cd();
		histo -> Write();

		histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
		histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));

		latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr}^{eff} = %.0f ps}{#sigma_{corr}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kMagenta);
		latex -> Draw("same");

		c -> Print(Form("%s/CTR_energyRatioPhasePosCorr/c_deltaT_energyRatioPhasePosCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
		c -> Print(Form("%s/CTR_energyRatioPhasePosCorr/c_deltaT_energyRatioPhasePosCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
		delete c;


		//tot, phase and position corr deltaT
		if(!h1_deltaT_totRatioPhasePosCorr[index2]) continue;
		c = new TCanvas(Form("c_deltaT_totRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioPhasePosCorr_%s",labelLR_energyBin.c_str()));
		histo = h1_deltaT_totRatioPhasePosCorr[index2];

		histo -> SetTitle(Form(";corrected #Deltat [ps];entries"));
		histo -> SetLineWidth(2);
		histo -> SetLineColor(kGreen+1);
		histo -> SetMarkerColor(kGreen+1);
		histo -> Draw("");

		FindSmallestInterval(vals,histo,0.68);
		min = vals[4];
		max = vals[5];
		delta = max-min;
		sigma = 0.5*delta;
		effSigma = sigma;

		fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
		fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
		
		//TF1* fitFunc = new TF1(Form("fitFunc_posCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
		fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
		fitFunc -> SetRange(fitXMin, fitXMax);
		histo -> Fit(fitFunc,"QNRSL","", fitXMin, fitXMax);
		fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
		histo -> Fit(fitFunc,"QNRSL","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
		fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
		histo -> Fit(fitFunc,"QRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));

		fitFunc -> SetLineColor(kMagenta);
		fitFunc -> SetLineWidth(3);
		fitFunc -> Draw("same");

		outFile -> cd();
		histo -> Write();

		histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
		histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));

		latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr}^{eff} = %.0f ps}{#sigma_{corr}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kMagenta);
		latex -> Draw("same");

		c -> Print(Form("%s/CTR_totRatioPhasePosCorr/c_deltaT_totRatioPhasePosCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
		c -> Print(Form("%s/CTR_totRatioPhasePosCorr/c_deltaT_totRatioPhasePosCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
		delete c;

	    }

	}

    }

  }


  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}



