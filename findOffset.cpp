#include <iostream>
#include <fstream>
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

using namespace std;


int main (int argc, char** argv){

  TFile *outFile = new TFile("./data/energyOffset.root","RECREATE");
  
  vector<string> bars = {"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14"};
  vector<string> vov = {"1.75", "2.00", "2.50", "3.50"};

  for(auto ibar : bars){
    
    TGraphErrors* gr_energy_offset_L = new TGraphErrors();
    TGraphErrors* gr_energy_offset_R = new TGraphErrors();
    TGraphErrors* gr_energy_offset1_L = new TGraphErrors();
    TGraphErrors* gr_energy_offset1_R = new TGraphErrors();
    
    TFile *f = TFile::Open("/afs/cern.ch/user/a/amkrishn/Amrutha/Lab5015Analysis/data/conf_23/moduleCharacterization_step2_run4894.root","READ");
    if(!f || !f->IsOpen()) continue;

    for (auto ov : vov){

      //redoing Andrea's plots
      float gain = (3.69 + 9.76*stof(ov))*pow(10,4);
      float pde =  0.384*(1 - exp(-0.583*stof(ov)));
      float ecf = 1 + 1.78*pow(10,-4)*stof(ov)*3.69;

      float gain_ref = (3.69 + 9.76*3.5)*pow(10,4);
      float pde_ref =  0.384*(1 - exp(-0.583*3.5));
      float ecf_ref = 1 + 1.78*pow(10,-4)*3.5*3.69;

      float y = (4.2*gain*pde*ecf) / (gain_ref*pde_ref*ecf_ref) ;

      // ------ left SiPM
      TH1F* h_ = (TH1F*)f->Get(Form("h1_energy_bar%sL_Vov%s_th05",ibar.c_str(),ov.c_str()));
	    if(h_ == NULL) continue;
	    TF1* fit_ = (TF1*)h_->GetFunction(Form("f_landau_bar%sL_Vov%s_vth1_05",ibar.c_str(),ov.c_str()));
	    

      if ((ibar == "10" || ibar == "12") && (ov == "2.00" || ov == "2.50" || ov == "1.75")){
        h_->GetXaxis()->SetRangeUser(100., 650.);
        float max = h_->GetBinCenter(h_->GetMaximumBin());
				TF1* f_landau1 = new TF1("f_landau1","landau", 10, 900);
        f_landau1->SetRange(max * 0.7, max * 1.5);
				h_->Fit(f_landau1,"QR");
				fit_ = f_landau1;
			}
			if(fit_->GetParError(1) > 40) continue;

	    gr_energy_offset_L->SetPoint(gr_energy_offset_L->GetN(), stof(ov), fit_->GetParameter(1));
	    gr_energy_offset_L->SetPointError(gr_energy_offset_L->GetN()-1, 0, fit_->GetParError(1));

      gr_energy_offset1_L->SetPoint(gr_energy_offset1_L->GetN(), fit_->GetParameter(1), y);
	    gr_energy_offset1_L->SetPointError(gr_energy_offset1_L->GetN()-1, fit_->GetParError(1), 0);

      // ------ right SiPM
	    h_ = (TH1F*)f->Get(Form("h1_energy_bar%sR_Vov%s_th05",ibar.c_str(), ov.c_str()));
	    if(h_ == NULL) continue;
	    fit_ = (TF1*)h_->GetFunction(Form("f_landau_bar%sR_Vov%s_vth1_05",ibar.c_str(), ov.c_str()));
	    if(fit_->GetParError(1) > 40) continue;

	    gr_energy_offset_R->SetPoint(gr_energy_offset_R->GetN(), stof(ov), fit_->GetParameter(1));
	    gr_energy_offset_R->SetPointError(gr_energy_offset_R->GetN()-1, 0, fit_->GetParError(1));

      gr_energy_offset1_R->SetPoint(gr_energy_offset1_R->GetN(), fit_->GetParameter(1), y);
	    gr_energy_offset1_R->SetPointError(gr_energy_offset1_R->GetN()-1, fit_->GetParError(1), 0);

	    delete h_, fit_;  


  }


    // ------- draw all

    //gStyle->SetOptFit(1111);
    
    TCanvas *c = new TCanvas();
    
    gr_energy_offset_L->SetMarkerStyle(20);
    gr_energy_offset_L->SetMarkerColor(kBlue);
    gr_energy_offset_L->SetTitle("; Vov ; energy MPV (a.u.)");
    //gr_energy_offset_L->GetXaxis()->SetLimits(0,14);
    TF1* f_pol1 = new TF1("f_pol1","pol1", 1.75, 3.5);
    //f_pol1->SetRange(1.75, 3.5);
    gr_energy_offset_L->Fit(f_pol1,"QR");
    gr_energy_offset_L->Draw("AP");

    float offset = f_pol1->GetParameter("p0");

    TLatex *latexL = new TLatex(0.40,0.84,Form("#splitline{bar %sL  th. = 5 DAC}{offset = %.1f}",ibar.c_str(),offset));
    latexL -> SetNDC();
    latexL -> SetTextFont(42);
    latexL -> SetTextSize(0.04);
    latexL -> SetTextColor(kRed);
    latexL->Draw("same");

    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_energy_offset_L_bar%s_th05.png",ibar.c_str()));
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_energy_offset_L_bar%s_th05.pdf",ibar.c_str()));

    outFile->cd();
    gr_energy_offset_L->Write(Form("energy_offset_L_bar%s_th05_conf23",ibar.c_str()));

    c->Clear();
    gr_energy_offset_R->SetMarkerStyle(20);
    gr_energy_offset_R->SetMarkerColor(kBlue);
    gr_energy_offset_R->SetTitle("; Vov ; energy MPV (a.u.)");
    gr_energy_offset_R->Fit(f_pol1,"QR");
    gr_energy_offset_R->Draw("AP");

    offset = f_pol1->GetParameter("p0");

    TLatex *latexR = new TLatex(0.40,0.84,Form("#splitline{bar %sR  th. = 5 DAC}{offset = %.1f}",ibar.c_str(),offset));
    latexR -> SetNDC();
    latexR -> SetTextFont(42);
    latexR -> SetTextSize(0.04);
    latexR -> SetTextColor(kRed);
    latexR->Draw("same");

    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_energy_offset_R_bar%s_th05.png",ibar.c_str()));
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_energy_offset_R_bar%s_th05.pdf",ibar.c_str()));

    outFile->cd();
    gr_energy_offset_R->Write(Form("energy_offset_R_bar%s_th05_conf23",ibar.c_str()));

    //------------- plot Andrea's method

    c->Clear();

    TF1* fitL_pol1 = new TF1("fitL_pol1","[0]+[1]*x", 20, 600);
    fitL_pol1->SetParameters(gr_energy_offset1_L->GetPointY(0),(gr_energy_offset1_L->GetPointY(3) - gr_energy_offset1_L->GetPointY(0)) / (gr_energy_offset1_L->GetPointX(3) - gr_energy_offset1_L->GetPointX(0)));
    //cout<<gr_energy_offset1_L->GetPointX(0)<<endl;
    gr_energy_offset1_L->SetMarkerStyle(20);
    //gr_energy_offset1_L->SetMarkerColor(kBlue);
    gr_energy_offset1_L->SetTitle("; energy MPV (a.u.); E*G*PDE*ECF / [G*PDE*ECF]_{3.5V}");
    fitL_pol1->SetRange(gr_energy_offset1_L->GetPointX(0), gr_energy_offset1_L->GetPointX(3));
    gr_energy_offset1_L->Fit(fitL_pol1,"QR");
    gr_energy_offset1_L->Draw("AP");

    offset = - fitL_pol1->GetParameter("p0")/fitL_pol1->GetParameter("p1");

    latexL = new TLatex(0.40,0.84,Form("#splitline{bar %sL  th. = 5 DAC}{offset = %.1f}",ibar.c_str(),offset));
    latexL -> SetNDC();
    latexL -> SetTextFont(42);
    latexL -> SetTextSize(0.04);
    latexL -> SetTextColor(kRed);
    latexL->Draw("same");

    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_energy_offset2_L_bar%s_th05.png",ibar.c_str()));
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_energy_offset2_L_bar%s_th05.pdf",ibar.c_str()));

    outFile->cd();
    gr_energy_offset1_L->Write(Form("energy_offset1_L_bar%s_th05_conf23",ibar.c_str()));

    c->Clear();
    TF1* fitR_pol1 = new TF1("fitR_pol1","pol1", -20, 600);
    gr_energy_offset1_R->SetMarkerStyle(20);
    //gr_energy_offset1_R->SetMarkerColor(kBlue);
    //cout<<gr_energy_offset1_R->GetPointX(0)<<"\t"<<gr_energy_offset1_R->GetPointX(3)<<endl;
    gr_energy_offset1_R->SetTitle("; energy MPV (a.u.); E*G*PDE*ECF / [G*PDE*ECF]_{3.5V}");
    fitR_pol1->SetParameters(1, (gr_energy_offset1_R->GetPointY(3) - gr_energy_offset1_R->GetPointY(0)) / (gr_energy_offset1_R->GetPointX(3) - gr_energy_offset1_R->GetPointX(0)));
    fitR_pol1->SetRange(gr_energy_offset1_R->GetPointX(0), gr_energy_offset1_R->GetPointX(3));
    gr_energy_offset1_R->Fit(fitR_pol1,"QR");
    gr_energy_offset1_R->Draw("AP");

    offset = - fitR_pol1->GetParameter("p0")/fitR_pol1->GetParameter("p1");

    latexR = new TLatex(0.40,0.85,Form("#splitline{bar %sR  th. = 5 DAC}{offset = %.1f}",ibar.c_str(), offset));
    latexR -> SetNDC();
    latexR -> SetTextFont(42);
    latexR -> SetTextSize(0.04);
    latexR -> SetTextColor(kRed);
    latexR->Draw("same");


    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_energy_offset2_R_bar%s_th05.png",ibar.c_str()));
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_energy_offset2_R_bar%s_th05.pdf",ibar.c_str()));

    outFile->cd();
    gr_energy_offset1_R->Write(Form("energy_offset1_R_bar%s_th05_conf23",ibar.c_str()));
    
    delete c, gr_energy_offset_L, gr_energy_offset_R, f_pol1, gr_energy_offset1_L, gr_energy_offset1_R, fitL_pol1,fitR_pol1;

  }

  TGraph* gr_pde_gain = new TGraph();
    //gain and pde calculation
    for(float iov = 0; iov < 5; iov = iov+0.25){
      float gain = (3.69 + 9.76*iov)*pow(10,4);
      float pde =  0.384*(1 - exp(-0.583*iov));
      float ecf = 1 + 1.78*pow(10,-4)*iov*3.69;

      gr_pde_gain->SetPoint(gr_pde_gain->GetN(), iov, gain*pde);
    }


    TF1* fit_pol3 = new TF1("fit_pol3","pol3", 0, 4);
    TCanvas *c1 = new TCanvas();
    gr_pde_gain->SetMarkerStyle(21);
    gr_pde_gain->SetTitle("; Vov (V); Gain x PDE");
    gr_pde_gain->Fit(fit_pol3,"QR");
    gr_pde_gain->Draw("AP");

    c1->SaveAs("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_pdeGain_vs_ov.png");
    c1->SaveAs("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/conf_23/SummaryPlots/c_pdeGain_vs_ov.pdf");

delete c1, gr_pde_gain, fit_pol3;
  
}
