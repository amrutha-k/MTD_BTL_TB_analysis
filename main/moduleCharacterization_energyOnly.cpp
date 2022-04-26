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
  system(Form("mkdir -p %s",plotDir.c_str()));
  //system(Form("mkdir -p %s/qfine/",plotDir.c_str()));
  system(Form("mkdir -p %s/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/energy2D/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyError/",plotDir.c_str()));

  
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

  int chL_ext = opts.GetOpt<int>("Coincidence.chL");
	int chR_ext = opts.GetOpt<int>("Coincidence.chR");
  
  //--- open files
  std::string step1FileName= opts.GetOpt<std::string>("Input.step1FileName");
  TFile* inFile = TFile::Open(step1FileName.c_str(),"READ");
  
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
  
  std::map<double,TH1F*> h1_deltaT_raw;
  std::map<double,TH1F*> h1_deltaT;
  
  std::map<double,TH2F*> h2_energyL_vs_energyR;
  std::map<double,TH2F*> h2_energy_error;

  std::map<std::string, std::map<int, std::vector<float>*> > mpv_energy;
  std::map<std::string, std::map<int, std::vector<float>*> > ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int, std::map<std::string,std::pair<float,float> > > > peaks;//peaks[LRlabel][index][energyPeak]
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
  TH1F* histo;
  TH2F* histo2;
  TProfile* prof;

  
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

        latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth1)));
        if (LRLabel == "L-R") { 
          latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	}
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);

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
	 
	   
	  if( opts.GetOpt<int>("Channels.array") == 1){
	    histo->GetXaxis()->SetRangeUser(25., 650.);
	    if (Vov>3.0) histo->GetXaxis()->SetRangeUser(120., 650.);
	    if (iBar == 7 && LRLabel == "R") histo->GetXaxis()->SetRangeUser(20., 650.);
	    if (iBar == 5 && LRLabel == "R") histo->GetXaxis()->SetRangeUser(10., 650.);
	    if (iBar == 13 && LRLabel == "L") histo->GetXaxis()->SetRangeUser(20., 650.);
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
	  f_landau[index] -> SetRange(max * 0.7, max * 1.5);    //0.75
	  f_landau[index] -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max, Vov*5.);
	  histo -> Fit(f_landau[index],"QRS");

    mpv_energy[LRLabel][index] -> push_back(f_landau[index]->GetParameter(1));

	  f_landau[index] -> SetLineColor(kBlack);
	  f_landau[index] -> SetLineWidth(2);
	  f_landau[index] -> Draw("same");
	    
	  if (step1FileName.find("2E14")!= std::string::npos && iBar == 14) ranges[LRLabel][index] -> push_back(5); //
	  else{
	    if (f_landau[index]->GetParameter(1)> 10) 
	      ranges[LRLabel][index] -> push_back( 0.75*f_landau[index]->GetParameter(1));   
	    else
	      ranges[LRLabel][index] -> push_back( 20 ); 
	  }
	  ranges[LRLabel][index] -> push_back( 940. );  

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
    
      theIndex = index;
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
	  std::string labelLR(Form("bar%02dL-R_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));

	  if(h2_energyL_vs_energyR[index1] == NULL) h2_energyL_vs_energyR[index1] = new TH2F(Form("h2_energyL_vs_energyR_%s",labelLR.c_str()),"",100,0,1000,100,0,1000);
	  if(anEvent->energyL > 0.6*mpv_energy["L"][index1]->at(0) && anEvent->energyR > 0.6*mpv_energy["R"][index1]->at(0)  && anEvent->energyL < 940 && anEvent->energyR < 940) h2_energyL_vs_energyR[index1] -> Fill (anEvent->energyL, anEvent->energyR);

	  accept[index1][entry] = false;
	    
	  if(!ranges["L-R"][index1] ) continue;
	    
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	    
	  if( energyBinAverage < 1 ) continue;
	    
	  accept[index1][entry] = true;

	  double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	    
	  if( h1_energyRatio[index2] == NULL )
	    {

	      std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
	            
              h1_energyRatio[index2] = new TH1F(Form("h1_energyRatio_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);
	
	      h1_deltaT_raw[index2] = new TH1F(Form("h1_deltaT_raw_%s",labelLR_energyBin.c_str()),"",2000,-24000.,24000.);
	        
	    }
	    
	  if (fabs(anEvent->timeR-anEvent->timeL)<10000) {
	        
	    if ((anEvent->energyR / anEvent->energyL >0) & (anEvent->energyR / anEvent->energyL <5)){
	      h1_energyRatio[index2] -> Fill( anEvent->energyR / anEvent->energyL );     
	      h1_deltaT_raw[index2] -> Fill( anEvent->timeR-anEvent->timeL );
	    }
	        
	  }
	        

	} //-- end loop over entries
      
      std::cout << std::endl;

    }
  
  
  
  
  //------------------
  //--- draw 2nd plots
  std::map<double,float> CTRMeans;
  std::map<double,float> CTRSigmas;
  
  std::map<double,TF1*> fitFunc_energyRatio;
  std::map<double,TF1*> fitFunc_energyL_vs_energyR;

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

	  if (!h2_energyL_vs_energyR[index1]) continue;

	  // -- draw energy L vs energy R
	  c = new TCanvas(Form("c_energyL_vs_energyR_%s",labelLR.c_str()),Form("c_energyL_vs_energyR_%s",labelLR.c_str()));
	  histo2 = h2_energyL_vs_energyR[index1];
	  //histo2 -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	  //histo2 -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	  histo2 -> GetYaxis() -> SetRangeUser(histo2->GetMean(2)-3.5*histo2->GetRMS(2), histo2->GetMean(2)+3.5*histo2->GetRMS(2));
	  histo2 -> SetTitle(Form(";energy_{left} (a.u.) ; energy_{right} (a.u.)"));
	  histo2 -> Draw("colz");
	        
	    fitFunc_energyL_vs_energyR[index1] = new TF1(Form("fitFunc_energyL_vs_energyR_%s",labelLR.c_str()),"pol1",-1000,1000);
      fitFunc_energyL_vs_energyR[index1] -> SetRange(ranges["L"][index1]->at(0), ranges["R"][index1]->at(0),ranges["L"][index1]->at(0)+420, ranges["R"][index1]->at(0)+420);
      histo2 -> Fit(fitFunc_energyL_vs_energyR[index1],"QNRS");
      
      fitFunc_energyL_vs_energyR[index1] -> SetLineColor(kBlack);
      fitFunc_energyL_vs_energyR[index1] -> SetLineWidth(2);
      fitFunc_energyL_vs_energyR[index1] -> Draw("same");
      
      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");

      histo2 -> Write();
      
      c -> Print(Form("%s/energy2D/c_energyL_vs_energyR__%s.png",plotDir.c_str(),labelLR.c_str()));
      c -> Print(Form("%s/energy2D/c_energyL_vs_energyR__%s.pdf",plotDir.c_str(),labelLR.c_str()));
      delete c;
  
  if( !ranges["L-R"][index1] ) continue;
  
  int nEnergyBins = ranges["L-R"][index1]->size()-1;
  
  for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
    {
      //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
      double index2( 10000000*iEnergyBin+index1 );
      
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

      
    } // --- end loop over energy bins
}// --- end loop ober bars
    }// --- end loop over stepLabels


//------------------------
  //--- 3rd loop over events (error in energy)
  
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
	  std::string labelLR(Form("bar%02dL-R_Vov%.2f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));

    if(anEvent->energyL < 0.6*mpv_energy["L"][index1]->at(0) || anEvent->energyR < 0.6*mpv_energy["R"][index1]->at(0) || anEvent->energyL > 940 || anEvent->energyR > 940) continue;
    //if(!ranges["L-R"][index1] ) continue;
    if(chL_ext == 1 && chR_ext == 30 && anEvent->vth1 == 5){
      if(h2_energy_error[index1] == NULL) h2_energy_error[index1] = new TH2F(Form("h2_energy_error_%s",labelLR.c_str()),"",1350,0,1350,400,-200,200);
      float theta = atan(fitFunc_energyL_vs_energyR[index1]->GetParameter(1));
      float u = cos(theta)*anEvent->energyL + sin(theta)*(anEvent->energyR - fitFunc_energyL_vs_energyR[index1]->GetParameter(0));
      float v = -sin(theta)*anEvent->energyL + cos(theta)*(anEvent->energyR - fitFunc_energyL_vs_energyR[index1]->GetParameter(0));
           //std::cout<<"xA: "<< anEvent->energyL<<"\t"<<"yA"<<anEvent->energyR<<"\t"<<"theta: "<<theta<<"\t"<<"uA: "<<u<<"\t"<<"vA: "<<v<<std::endl;
      h2_energy_error[index1] -> Fill(u,v);
    }
       

	} //-- end loop over entries
      
      std::cout << std::endl;

    }
  
  
  
  
  //------------------
  //--- draw 3rd (energy error) plots
  
  //std::map<double,TF1*> fitFunc_energyL_vs_energyR;

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

	  if (!h2_energy_error[index1]) continue;

	  // -- draw energy error each bar
	  c = new TCanvas(Form("c_energy_error_%s",labelLR.c_str()),Form("c_energy_error_%s",labelLR.c_str()));
	  histo2 = h2_energy_error[index1];
	  //histo2 -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	  //histo2 -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	  //histo2 -> GetYaxis() -> SetRangeUser(histo2->GetMean(2)-3.5*histo2->GetRMS(2), histo2->GetMean(2)+3.5*histo2->GetRMS(2));
    //histo2 -> GetXaxis() -> SetRangeUser(histo2->GetMean(1)-3.5*histo2->GetRMS(1), histo2->GetMean(1)+3.5*histo2->GetRMS(1));
	  histo2 -> SetTitle(Form(";energy_{left} (a.u.) ; energy_error (a.u.)"));
	  histo2 -> Draw("colz");
      
      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");

      histo2 -> Write();
      
      c -> Print(Form("%s/energyError/c_energy_error__%s.png",plotDir.c_str(),labelLR.c_str()));
      c -> Print(Form("%s/energyError/c_energy_error__%s.pdf",plotDir.c_str(),labelLR.c_str()));
      delete c;
  
 
}// --- end loop over bars
    }// --- end loop over stepLabels



  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}



