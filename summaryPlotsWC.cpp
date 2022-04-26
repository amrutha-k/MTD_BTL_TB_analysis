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

using namespace std;


int main (int argc, char** argv){

	TFile *outFile = new TFile("./data/SummaryPlotsWC.root","RECREATE");
	
	vector<string> th = {"03","05","07","10","15","20","30"}; 
	vector<string> bars = {"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15"};

	
	for(auto ibar : bars){
		
		TGraphErrors* gr_tRes_vs_x = new TGraphErrors();
		//TTree *t = (TTree*)f->Get("data");
		for (int slice = -3; slice<14; slice=slice+2){

			TFile *f = TFile::Open(Form("/afs/cern.ch/user/a/amkrishn/Amrutha/Lab5015Analysis/data/%dto%d/moduleCharacterization_step2_run%s.root",slice,(slice+2),argv[1]),"READ");
			if(!f || !f->IsOpen()) continue;

			string bestThr = "05";
			int bestRes = 500;

			for(auto ith : th){

				//const char* name = Form("h1_deltaT_energyRatioPhaseCorr_bar%sL-R_Vov3.50_th%s_energyBin01",ibar.c_str(),ith.c_str());
				//if(!(f->GetListOfKeys()->Contains(name))) continue;  

				//if(!(gDirectory->FindObject(name))) continue;
				//if((ibar == "15") && (ith == "03")) continue;   //doesn't exist
				TH1F* h = (TH1F*)f->Get(Form("h1_deltaT_totRatioPhaseCorr_bar%sL-R_Vov3.50_th%s_energyBin01",ibar.c_str(),ith.c_str()));
				if(h == NULL) continue;
				TF1 *fit = (TF1*)h->GetFunction(Form("fitFunc_phaseCorr_bar%sL-R_Vov3.50_th%s_energyBin01",ibar.c_str(),ith.c_str()));

				if (0.5*fit->GetParameter(2) < 35 ) continue;
				if (fit->GetParError(2) > 10) continue;

				if (fit->GetParameter(2) < bestRes){
					bestThr = ith;
					bestRes = fit->GetParameter(2);
				}  
			delete h, fit;
			}

		if (bestRes == 500) continue;
		//if (refbar == 0 && ibar == "03")   bestThr = "07";
	
		TH1F* h_bestThr = (TH1F*)f->Get(Form("h1_deltaT_totRatioPhaseCorr_bar%sL-R_Vov3.50_th%s_energyBin01",ibar.c_str(),bestThr.c_str()));
		TF1 *fit_bestThr = (TF1*)h_bestThr->GetFunction(Form("fitFunc_phaseCorr_bar%sL-R_Vov3.50_th%s_energyBin01",ibar.c_str(),bestThr.c_str()));

		gr_tRes_vs_x->SetPoint(gr_tRes_vs_x->GetN(), slice+1, fit_bestThr->GetParameter(2)/2);
		gr_tRes_vs_x->SetPointError(gr_tRes_vs_x->GetN()-1, 1, fit_bestThr->GetParError(2)/2);

		delete h_bestThr, fit_bestThr;

		}
	
			TCanvas *c = new TCanvas();
			gr_tRes_vs_x->SetMarkerStyle(20);
			gr_tRes_vs_x->SetTitle(";x (mm) ; #sigma_{bar}(ns)");
			//gr_tRes_vs_refBarNo->SetTitleStyle(1001);
			gr_tRes_vs_x->GetXaxis()->SetLimits(-4,14);
			gr_tRes_vs_x->Draw("AP");
			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/test/SummaryPlotsWC/c_tRes_vs_x_bar%s_Vov3.5.png",argv[1],ibar.c_str(),argv[1]));
			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/test/SummaryPlotsWC/c_tRes_vs_x_bar%s_Vov3.5.pdf",argv[1],ibar.c_str(),argv[1]));

			outFile->cd();
			gr_tRes_vs_x->Write(Form("tRes_vs_x_bar%s_Vov3.50_run%s",ibar.c_str(),argv[1]));

 			delete c, gr_tRes_vs_x; 

	}	

}

