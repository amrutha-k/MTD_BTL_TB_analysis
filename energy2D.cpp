//#include "interface/AnalysisUtils.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>

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
#include "TPad.h"
#include "TStyle.h"
#include "TClass.h"

class Events {

public:
  int barID;
  float Vov;
  int vth1;
  float energyL;
  float energyR;
  float totL;
  float totR;
  long long timeL;
  long long timeR;
  unsigned short t1fineL;
  unsigned short t1fineR;
  int nhits;
  float x;
  float y;
  float qT1L;
  float qT1R;
  float amp_MCP;
  float t_MCP;
  float t_CLK_P;
  float t_CLK_C;
  float phi_L;
  float phi_R;
  float phi_MCP;
  float phi_MCP_C;

};


using namespace std;


int main (int argc, char** argv){

  TFile *outFile = new TFile("./data/energy2D.root","RECREATE");
  TCanvas *c = new TCanvas();
  
  //vector<string> th = {"03","05","07","10","15","20","30"}; 
  vector<string> bars = {"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14"};
  
  for(auto ibar : bars){

    TCanvas *c1 = new TCanvas();
    TGraph *gr_energy_intersection = new TGraph();
    //gStyle->SetOptFit(1111);
    
    
    for (int refbar = 4; refbar<13; refbar++){

     /* TFile *f = TFile::Open(Form("/afs/cern.ch/user/a/amkrishn/Amrutha/Lab5015Analysis/data/refBar%d/moduleCharacterization_step1_run%s.root",refbar,argv[1]),"READ");
      if(!f || !f->IsOpen()) continue;

      TH2F* h2 = new TH2F(Form("h2_energyL_vs_energyR_bar%sL-R_refbar%d_Vov3.50_th05",ibar.c_str(),refbar),"",100,0,1000,100,0,1000);

      TTree *t = (TTree*)f->Get(Form("data_bar%sL-R_Vov3.50_th05",ibar.c_str()));
      long nentries = t->GetEntries();
      std::cout<<"nentries: "<<nentries<<std::endl;
      Events* anEvent = new Events();
      t -> SetBranchAddress("event",&anEvent);
      
      for(int entry = 0; entry < nentries; entry++){
        if( entry%100000 == 0 ) std::cout << ">>> 1st loop: << reading entry " << entry << " / " << nentries << " (" << 100.*entry/nentries << "%)" << "\r" << std::flush;
        t->GetEntry(entry);
        if(anEvent->energyL > 0 && anEvent->energyR > 0 && anEvent->energyL < 940 && anEvent->energyR < 940){
        h2 -> Fill (anEvent->energyL, anEvent->energyR); 
        }
        
      } 
      */

      TFile *f2 = TFile::Open(Form("/eos/user/a/amkrishn/BTL/MTDTB_OCT2021_CERN_afs_backup/data/refBar%d/moduleCharacterization_step2_run%s.root",refbar,argv[1]),"READ");
      if(!f2 || !f2->IsOpen()) continue;

      TH2F* h2 = (TH2F*)f2->Get(Form("h2_energyL_vs_energyR_bar%sL-R_Vov3.50_th05",ibar.c_str()));
      if(h2 == NULL) continue;

      TH1F *hLR = (TH1F*)f2->Get(Form("h1_energy_bar%sL-R_Vov3.50_th05",ibar.c_str()));
      if(hLR == NULL) continue;
      TF1 *fit = (TF1*)hLR->GetFunction(Form("f_landau_bar%sL-R_Vov3.50_vth1_05",ibar.c_str()));

      TH1F *hL = (TH1F*)f2->Get(Form("h1_energy_bar%sL_Vov3.50_th05",ibar.c_str()));
      if(hL == NULL) continue;
      TF1 *fitL = (TF1*)hL->GetFunction(Form("f_landau_bar%sL_Vov3.50_vth1_05",ibar.c_str()));
      
      float minX = 0.75*fitL->GetParameter(1);
      float maxX = 0.75*fitL->GetParameter(1) + 420;

      TH1F *hR = (TH1F*)f2->Get(Form("h1_energy_bar%sR_Vov3.50_th05",ibar.c_str()));
      if(hR == NULL) continue;
      TF1 *fitR = (TF1*)hR->GetFunction(Form("f_landau_bar%sR_Vov3.50_vth1_05",ibar.c_str()));
      
      float minY = 0.75*fitR->GetParameter(1);
      //if(ibar == "13") minY = 480;
      float maxY = 0.75*fitR->GetParameter(1) + 420;
      //if(ibar == "13") maxY = 700;
      
      TF1 *fitFunc = new TF1((Form("fitFunc_energyL_vs_energyR_refbar%d_bar%sL-R_Vov3.50_th05",refbar,ibar.c_str())),"pol1",-1000,1000);
      fitFunc->SetRange(minX, minY, maxX, maxY);

      if(h2 == NULL) continue;
      h2 -> GetXaxis()->SetRangeUser(0.5*fit->GetParameter(1), 1000);
      h2 -> GetYaxis()->SetRangeUser(0.5*fit->GetParameter(1), 1000);
      c->cd();
      h2 -> GetYaxis() -> SetRangeUser(h2->GetMean(2)-3.5*h2->GetRMS(2), h2->GetMean(2)+3.5*h2->GetRMS(2));
	    h2 -> SetTitle(Form(";energy_{left} (a.u.) ; energy_{right} (a.u.)"));
      h2->Fit(fitFunc,"QR");
      h2->Draw("colz");
      c->SetLogz();

      TLatex *latex = new TLatex(0.40,0.84,Form("#splitline{bar %sLR}{Vov = 3.50  th. = 5 DAC}",ibar.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");

      outFile->cd();
      h2->Write();

      c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_energyL_vs_energyR_bar%s_refbar%d_Vov3.5_th05.png",argv[1],ibar.c_str(),refbar));
      c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_energyL_vs_energyR_bar%s_refbar%d_Vov3.5_th05.pdf",argv[1],ibar.c_str(),refbar));

/*
      //---------- find the error in energy 

      if (refbar == 8){
        TGraph *gr_energy_err = new TGraph();
        vector<int> binShift = {0,-10,10,20};
        for(int i : binShift) {
          int bin = h2->FindBin(fitL->GetParameter(1)+i);
          //cout<<"bin = "<<bin<<endl;
          TH1D *proj = h2->ProjectionY(Form("h2_energyL_vs_energyR_bar%sL-R_Vov3.50_th05_py%d",ibar.c_str(),i),bin,bin);
          TF1* fgaus = new TF1(Form("fgaus_energyL_vs_energyR_bar%sL-R_Vov3.50_th05_py%d",ibar.c_str(),i),"gaus",proj->GetMean()-2*proj->GetRMS(),proj->GetMean()+2*proj->GetRMS());
          proj->Fit(fgaus,"QR");
          float err = 100*(fgaus->GetParameter(2) / fgaus->GetParameter(1));

          gr_energy_err->SetPoint(gr_energy_err->GetN(),gr_energy_err->GetN()+1,err);

        }
        c->Clear();
        c->cd();
        gr_energy_err->SetTitle(";slice no.;energy error (%)}");
        gr_energy_err->SetMarkerStyle(20);
        gr_energy_err->Draw("AP");
        c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_energy_error_bar%s_refbar%d_Vov3.5_th05.png",argv[1],ibar.c_str(),refbar));
        c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_energy_error_bar%s_refbar%d_Vov3.5_th05.pdf",argv[1],ibar.c_str(),refbar));
        
       /* TFile *f = TFile::Open(Form("/afs/cern.ch/user/a/amkrishn/Amrutha/Lab5015Analysis/data/refBar%d/moduleCharacterization_step1_run%s.root",refbar,argv[1]),"READ");
        if(!f || !f->IsOpen()) continue;
        t = (TTree*)f->Get(Form("data_bar%sL-R_Vov3.50_th05",ibar.c_str()));
        long nentries = t->GetEntries();
        std::cout<<"nentries: "<<nentries<<std::endl;
        //Events* anEvent = new Events();
        t -> SetBranchAddress("event",&anEvent);

        TH2F* h2_err = new TH2F(Form("h2_energy_err_bar%sL-R_refbar%d_Vov3.50_th05",ibar.c_str(),refbar),"",100,0,1000,100,0,1000);
        for(int ientry = 0; ientry < nentries; ientry++){
        if( ientry%100000 == 0 ) std::cout << ">>> 2nd loop: << reading entry " << ientry << " / " << nentries << " (" << 100.*ientry/nentries << "%)" << "\r" << std::flush;
        t->GetEntry(ientry);
        if(anEvent->energyL > 0 && anEvent->energyR > 0 && anEvent->energyL < 940 && anEvent->energyR < 940){
          float theta = atan(fitFunc->GetParameter(1));
          float u = cos(theta)*anEvent->energyL + sin(theta)*anEvent->energyR ;
          float v = -sin(theta)*anEvent->energyL + cos(theta)*anEvent->energyR ;

          h2_err -> Fill (u, v); 
          }
        }

        if(h2_err == NULL) continue;
        c->Clear();
        c->cd();
        h2_err -> GetYaxis() -> SetRangeUser(h2->GetMean(2)-3.5*h2->GetRMS(2), h2->GetMean(2)+3.5*h2->GetRMS(2));
	      h2_err -> SetTitle(Form(";energy_{left} (a.u.) ; energy_{right} (a.u.)"));
        h2_err->Fit(fitFunc,"QR");
        h2_err->Draw("colz");
        c->SetLogz();

        latex = new TLatex(0.40,0.84,Form("#splitline{bar %sLR}{Vov = 3.50  th. = 5 DAC}",ibar.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");

        outFile->cd();
        h2_err->Write();

        c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_energy_err_bar%s_refbar%d_Vov3.5_th05.png",argv[1],ibar.c_str(),refbar));
        c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_energy_err_bar%s_refbar%d_Vov3.5_th05.pdf",argv[1],ibar.c_str(),refbar));

      }
      
*/
      float p0_ref, p1_ref, x, y;
      c1->cd();
      fitFunc->SetRange(-500, 600);
      if(refbar == 4) {
        fitFunc->Draw();
        p0_ref = fitFunc->GetParameter(0);
        p1_ref = fitFunc->GetParameter(1);
        //cout<<"bar"<<ibar<<"\t"<<"refbar"<<refbar"\t"<<p0_ref<<","<<p1_ref
      }
      else {
        fitFunc->Draw("same");
        x = (fitFunc->GetParameter(0) - p0_ref) / (p1_ref - fitFunc->GetParameter(1)) ;
        y = p1_ref*x + p0_ref;

        if(x > 1000 || x < -1000) continue;
        if(y > 1000 || y < -1000) continue;
        gr_energy_intersection->SetPoint(gr_energy_intersection->GetN(), x, y);
      }

      f2->Close();

  }
  //gStyle->SetPalette(kBird);
  gPad->SetTitle(";energy_{left} (a.u.) ; energy_{right} (a.u.)");
  c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_fit_energyL_vs_energyR_bar%s_Vov3.5_th05.png",argv[1],ibar.c_str()));
  c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_fit_energyL_vs_energyR_bar%s_Vov3.5_th05.pdf",argv[1],ibar.c_str()));
  c1->Write();

  c1->Clear();
  c1->cd();
  gr_energy_intersection->SetTitle(";offset energy_L (a.u); offset energy_R (a.u)");
  gr_energy_intersection->SetMarkerStyle(kOpenSquare);
  gr_energy_intersection->Draw("AP");

  outFile->cd();
  gr_energy_intersection->Write();

  c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_energy_intersection_bar%s_Vov3.5_th05.png",argv[1],ibar.c_str()));
  c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/%s/SummaryPlots/c_energy_intersection_bar%s_Vov3.5_th05.pdf",argv[1],ibar.c_str()));
  
}

}
