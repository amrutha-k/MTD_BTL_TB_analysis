#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;


struct pulse_t {
  Int_t event;
  Short_t tc[2];
  Float_t channel[18][1024];
  Float_t times[2][1024];
} pulse;

double getmin(TGraph *gr) {
  double mean = 0;
  for(int i=0; i!=250; ++i)
    mean += gr->GetPointY(i);
  mean /= 250;

  double min = 11111111111111111;
  for(int i=380; i!=750; ++i) {
    if(gr->GetPointY(i)<min)
      min = gr->GetPointY(i);
  }
  return  -( min - mean );
}


int main() {
  // TFile *file = new TFile(Form("merged/Run_%d.root",run));
  TFile *fout = new TFile("/afs/cern.ch/user/a/amkrishn/public/BTL/uvlaser_40022.root","RECREATE");
  TFile *file = new TFile("/eos/cms/store/group/dpg_mtd/comm_mtd/TB/benchtest_2021_07/merged/Run_40022.root");
  TTree *tree = (TTree*) file->Get("pulse");
  //tree->Print();
  tree->SetBranchAddress("channel",&pulse.channel);
  tree->SetBranchAddress("times",&pulse.times);

  Long_t nevts = tree->GetEntries();
  int nbins = 200;
  TH1D *amp_sipm0 = new TH1D("amp_sipm0","SIPM0; amp (mV)",nbins,0,500);
  TH1D *amp_sipm1 = new TH1D("amp_sipm1","SIPM1;amp (mV)",nbins,0,500);
  //TH2D *hWF1 = new TH2D("hWF1","hWF1;ns;mV",100,0,200,100,-500,+500 );
  //TH2D *hWF2 = new TH2D("hWF2","hWF2;ns;mV",100,0,200,100,-500,+500 );
  Double_t source[nbins];
  
  for(int i=0; i!=nevts; ++i) {
    tree->GetEntry(i);
    TGraph *gr1 = new TGraph( 1024, pulse.times[0], pulse.channel[2] );
    TGraph *gr2 = new TGraph( 1024, pulse.times[0], pulse.channel[4] );
    if((i%500)==0)
      cout << "Events processed: " << i << endl;
    /*for(int k=0; k!=1024; ++k) {
      double x = gr1->GetPointX(k);
      double y = gr1->GetPointY(k);
      hWF1->Fill( x, y );
      x = gr2->GetPointX(k);
      y = gr2->GetPointY(k);
      hWF2->Fill( x, y );
      
    }*/
    
    /// DO STUFF HERE
    amp_sipm0->Fill( getmin( gr1 ) );
    amp_sipm1->Fill( getmin( gr2 ) );

    /*Double_t fpeaks(Double_t *x, Double_t *par) {
      Double_t result = par[0] + par[1]*x[0];
      for(Int_t p=0;p<npeaks;p++) {
  Double_t norm = par[3*p+2]; //height
  Double_t mean = par[3*p+3];
  Double_t sigma = par[3*p+4];
  norm /= sigma * (TMath::Sqrt(TMath::TwoPi()));
  result += norm*TMath::Gaus(x[0],mean,sigma);
      }
      return result;
      }*/

    delete gr1;
    delete gr2;

   }

    /*TCanvas *main = new TCanvas();
    main->Divide(2,1);
    main->cd(1);
    amp_sipm0->Draw();
    main->cd(2);
    amp_sipm1->Draw();
    main->SaveAs("/eos/user/a/amkrishn/BTL/40008_amp.png");*/

    tree->GetEntry(100);
    TGraph *pulse_sipm0 = new TGraph( 1024, pulse.times[0], pulse.channel[2] );
    TGraph *pulse_sipm1 = new TGraph( 1024, pulse.times[0], pulse.channel[4] );
    pulse_sipm0->GetXaxis()->SetTitle("time (ns)");
    pulse_sipm1->GetXaxis()->SetTitle("time (ns)");
    pulse_sipm0->GetYaxis()->SetTitle("amp (mV)");
    pulse_sipm1->GetYaxis()->SetTitle("amp (mV)");
    pulse_sipm0->SetTitle("SIPM 0");
    pulse_sipm1->SetTitle("SIPM 1");


    TCanvas *c = new TCanvas();
    //c->Divide(2,1);
    //c->cd(1);
    pulse_sipm0->Draw();
    TCanvas *c1 = new TCanvas();
    pulse_sipm1->Draw();

    fout->cd();
    amp_sipm0->Write();
    amp_sipm1->Write();
    c->Write();
    c1->Write();

    fout->Close();

    //c->SaveAs("/eos/user/a/amkrishn/BTL/40008_pulse.png");

    delete pulse_sipm0, pulse_sipm1, amp_sipm0, amp_sipm1,c,c1;

    return 0;
  }