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

	TFile *outFile = new TFile("./data/SummaryPlots.root","RECREATE");
	
	vector<string> th = {"03","05","07","10","15","20","30"}; 
	vector<string> bars = {"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14"};

	

	for(auto ibar : bars){
		
		//TGraphErrors* gr_tRes_vs_refBarNo = new TGraphErrors();
		TGraphErrors* gr_mpvL_vs_refBarNo_1 = new TGraphErrors();
		TGraphErrors* gr_mpvR_vs_refBarNo_1 = new TGraphErrors();
		TGraphErrors* gr_mpvLR_vs_refBarNo_1 = new TGraphErrors();

		TGraphErrors* gr_mpvL_vs_refBarNo_2 = new TGraphErrors();
		TGraphErrors* gr_mpvR_vs_refBarNo_2 = new TGraphErrors();
		TGraphErrors* gr_mpvLR_vs_refBarNo_2 = new TGraphErrors();

		TGraphErrors* gr_mpvL_vs_refBarNo_3 = new TGraphErrors();
		TGraphErrors* gr_mpvR_vs_refBarNo_3 = new TGraphErrors();
		TGraphErrors* gr_mpvLR_vs_refBarNo_3 = new TGraphErrors();

		TFile *f_offset = TFile::Open("/afs/cern.ch/user/a/amkrishn/Amrutha/Lab5015Analysis/data/energyOffset.root","READ");

		//----- offset using Andrea's method 
		TGraphErrors *gr_offsetL = (TGraphErrors*)f_offset->Get(Form("energy_offset1_L_bar%s_th05_conf23",ibar.c_str()));
		TF1* fitFuncL = (TF1*)gr_offsetL->GetFunction("fitL_pol1");
		float offsetL_1 = - fitFuncL->GetParameter(0) / fitFuncL->GetParameter(1);
		TGraphErrors *gr_offsetR = (TGraphErrors*)f_offset->Get(Form("energy_offset1_R_bar%s_th05_conf23",ibar.c_str()));
		TF1* fitFuncR = (TF1*)gr_offsetR->GetFunction("fitR_pol1");
		float offsetR_1 = - fitFuncR->GetParameter(0) / fitFuncR->GetParameter(1);

		

		//----- offset using mpv vs OV
		gr_offsetL = (TGraphErrors*)f_offset->Get(Form("energy_offset_L_bar%s_th05_conf23",ibar.c_str()));
		fitFuncL = (TF1*)gr_offsetL->GetFunction("f_pol1");
		float offsetL_2 = fitFuncL->GetParameter(0);
		gr_offsetR = (TGraphErrors*)f_offset->Get(Form("energy_offset_R_bar%s_th05_conf23",ibar.c_str()));
		fitFuncR = (TF1*)gr_offsetR->GetFunction("f_pol1");
		float offsetR_2 = fitFuncR->GetParameter(0);


		//----- offset using Sasha's method
		map<string,vector<float> > map_bar_offset3 = {
			{"00",{0,0}},
			{"01",{-175,-68}},
			{"02",{-35,-98}},
			{"03",{-407,-228}},
			{"04",{-368,-445}},
			{"05",{-121.4,-231.5}},
			{"06",{-334,-231.5}},
			{"07",{-155,-342}},
			{"08",{-39,-214}},
			{"09",{-204,-104}},
			{"10",{-68,-115.5}},
			{"11",{-172.6,-348.6}},
			{"12",{-90.2,-48}},
			{"13",{0,0}},
			{"14",{0,0}},

		};
	/*	map_bar_offset3["00"] = ;
		map_bar_offset3["01"] = ;
		map_bar_offset3["02"] = -35;
		map_bar_offset3["03"] = ;
		map_bar_offset3["04"] = ;
		map_bar_offset3["05"] = ;
		map_bar_offset3["06"] = ;
		map_bar_offset3["07"] = ;
		map_bar_offset3["08"] = ;
		map_bar_offset3["09"] = ;
		map_bar_offset3["10"] = ;
		map_bar_offset3["11"] = ;
		map_bar_offset3["12"] = ;
		map_bar_offset3["13"] = ;
		map_bar_offset3["14"] = ;
*/
		
		for (int refbar = 0; refbar<15; refbar++){

			TFile *f = TFile::Open(Form("/eos/user/a/amkrishn/BTL/MTDTB_OCT2021_CERN_afs_backup/data/refBar%d/moduleCharacterization_step2_run%s.root",refbar,argv[1]),"READ");
			if(!f || !f->IsOpen()) continue;

/*			string bestThr = "05";
			int bestRes = 500;

			// ---------- find the best threshold

			for(auto ith : th){

				TH1F* h = (TH1F*)f->Get(Form("h1_deltaT_totRatioPhaseCorr_bar%sL-R_Vov3.50_th%s_energyBin01",ibar.c_str(),ith.c_str()));
				if(h == NULL) continue;
				TF1 *fit = (TF1*)h->GetFunction(Form("fitFunc_phaseCorr_bar%sL-R_Vov3.50_th%s_energyBin01",ibar.c_str(),ith.c_str()));

				if (0.5*fit->GetParameter(2) < 37 ) continue;
				if (fit->GetParError(2) > 14) continue;

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

		gr_tRes_vs_refBarNo->SetPoint(gr_tRes_vs_refBarNo->GetN(), refbar, fit_bestThr->GetParameter(2)/2);
		gr_tRes_vs_refBarNo->SetPointError(gr_tRes_vs_refBarNo->GetN()-1, 0, fit_bestThr->GetParError(2)/2);

		// ------------ mpv vs impact point

		
		delete h_bestThr, fit_bestThr;

*/

		if(refbar >= 4 && refbar <= 12){
			TH1F* h_ = (TH1F*)f->Get(Form("h1_energy_bar%sL-R_Vov3.50_th05",ibar.c_str()));
			if(h_ == NULL) continue;
			TF1 *fit_ = (TF1*)h_->GetFunction(Form("f_landau_bar%sL-R_Vov3.50_vth1_05",ibar.c_str()));
		//	if(fit_->GetParError(1) > 40) continue;
		//	if(fit_->GetParameter(1) < 10) continue;
		//	gr_mpvLR_vs_refBarNo_1->SetPoint(gr_mpvLR_vs_refBarNo_1->GetN(), refbar, fit_->GetParameter(1));
		//	gr_mpvLR_vs_refBarNo_1->SetPointError(gr_mpvLR_vs_refBarNo_1->GetN()-1, 0, fit_->GetParError(1));

			h_ = (TH1F*)f->Get(Form("h1_energy_bar%sL_Vov3.50_th05",ibar.c_str()));
			if(h_ == NULL) continue;
			fit_ = (TF1*)h_->GetFunction(Form("f_landau_bar%sL_Vov3.50_vth1_05",ibar.c_str()));
			if(fit_->GetParError(1) > 40) continue;
			if(fit_->GetParameter(1) < 10) continue;
			//gr_mpvL_vs_refBarNo_1->SetPoint(gr_mpvL_vs_refBarNo_1->GetN(), refbar, fit_->GetParameter(1));
			gr_mpvL_vs_refBarNo_1->SetPoint(gr_mpvL_vs_refBarNo_1->GetN(), refbar, fit_->GetParameter(1)-offsetL_1);
			gr_mpvL_vs_refBarNo_1->SetPointError(gr_mpvL_vs_refBarNo_1->GetN()-1, 0, fit_->GetParError(1));
			gr_mpvL_vs_refBarNo_2->SetPoint(gr_mpvL_vs_refBarNo_2->GetN(), refbar, fit_->GetParameter(1)-offsetL_2);
			gr_mpvL_vs_refBarNo_2->SetPointError(gr_mpvL_vs_refBarNo_2->GetN()-1, 0, fit_->GetParError(1));
			gr_mpvL_vs_refBarNo_3->SetPoint(gr_mpvL_vs_refBarNo_3->GetN(), refbar, fit_->GetParameter(1)- map_bar_offset3[Form("%s",ibar.c_str())].at(0));
			gr_mpvL_vs_refBarNo_3->SetPointError(gr_mpvL_vs_refBarNo_3->GetN()-1, 0, fit_->GetParError(1));


			h_ = (TH1F*)f->Get(Form("h1_energy_bar%sR_Vov3.50_th05",ibar.c_str()));
			if(h_ == NULL) continue;
			fit_ = (TF1*)h_->GetFunction(Form("f_landau_bar%sR_Vov3.50_vth1_05",ibar.c_str()));
			if (ibar == "07"){
				TF1* f_landau1 = new TF1("f_landau1","landau", 10, 500);
				h_->Fit(f_landau1,"QR");
				fit_ = f_landau1;
			}
			if (ibar == "11"){
			  TF1* f_landau2 = new TF1("f_landau2","landau", 20, 500);
			  h_->Fit(f_landau2,"QR");
			  fit_ = f_landau2;                                                                                                                                                
            }
			if (ibar == "05"){
			  TF1* f_landau3 = new TF1("f_landau3","landau", 10, 50);
			  h_->Fit(f_landau3,"QR");
			  fit_ = f_landau3;   
			  cout<<f_landau3->GetParameter(1)<<endl;                                                                                                                                                  
            }

			if(fit_->GetParError(1) > 40) continue;
			if(fit_->GetParameter(1) < 10) continue;
			gr_mpvR_vs_refBarNo_1->SetPoint(gr_mpvR_vs_refBarNo_1->GetN(), refbar, fit_->GetParameter(1)-offsetR_1);
			gr_mpvR_vs_refBarNo_1->SetPointError(gr_mpvR_vs_refBarNo_1->GetN()-1, 0, fit_->GetParError(1));
			gr_mpvR_vs_refBarNo_2->SetPoint(gr_mpvR_vs_refBarNo_2->GetN(), refbar, fit_->GetParameter(1)-offsetR_2);
			gr_mpvR_vs_refBarNo_2->SetPointError(gr_mpvR_vs_refBarNo_2->GetN()-1, 0, fit_->GetParError(1));
			gr_mpvR_vs_refBarNo_3->SetPoint(gr_mpvR_vs_refBarNo_3->GetN(), refbar, fit_->GetParameter(1)- map_bar_offset3[Form("%s",ibar.c_str())].at(1));
			gr_mpvR_vs_refBarNo_3->SetPointError(gr_mpvR_vs_refBarNo_3->GetN()-1, 0, fit_->GetParError(1));

			delete h_, fit_;
		}


		}

			float mpv_refL1 = gr_mpvL_vs_refBarNo_1->GetPointY(4);
			float mpv_refR1 = gr_mpvR_vs_refBarNo_1->GetPointY(4);
			float mpv_refL2 = gr_mpvL_vs_refBarNo_2->GetPointY(4);
			float mpv_refR2 = gr_mpvR_vs_refBarNo_2->GetPointY(4);
			float mpv_refL3 = gr_mpvL_vs_refBarNo_3->GetPointY(4);
			float mpv_refR3 = gr_mpvR_vs_refBarNo_3->GetPointY(4);

			for(int i =0; i<gr_mpvL_vs_refBarNo_1->GetN(); i++){

				gr_mpvL_vs_refBarNo_1->SetPoint(i, gr_mpvL_vs_refBarNo_1->GetPointX(i), gr_mpvL_vs_refBarNo_1->GetPointY(i)/mpv_refL1);
				gr_mpvL_vs_refBarNo_1->SetPointError(i, 0, gr_mpvL_vs_refBarNo_1->GetErrorY(i)/mpv_refL1);
				gr_mpvR_vs_refBarNo_1->SetPoint(i, gr_mpvR_vs_refBarNo_1->GetPointX(i), gr_mpvR_vs_refBarNo_1->GetPointY(i)/mpv_refR1);
				gr_mpvR_vs_refBarNo_1->SetPointError(i, 0, gr_mpvR_vs_refBarNo_1->GetErrorY(i)/mpv_refR1);

				gr_mpvL_vs_refBarNo_2->SetPoint(i, gr_mpvL_vs_refBarNo_2->GetPointX(i), gr_mpvL_vs_refBarNo_2->GetPointY(i)/mpv_refL2);
				gr_mpvL_vs_refBarNo_2->SetPointError(i, 0, gr_mpvL_vs_refBarNo_2->GetErrorY(i)/mpv_refL2);
				gr_mpvR_vs_refBarNo_2->SetPoint(i, gr_mpvR_vs_refBarNo_2->GetPointX(i), gr_mpvR_vs_refBarNo_2->GetPointY(i)/mpv_refR2);
				gr_mpvR_vs_refBarNo_2->SetPointError(i, 0, gr_mpvR_vs_refBarNo_2->GetErrorY(i)/mpv_refR2);

				gr_mpvL_vs_refBarNo_3->SetPoint(i, gr_mpvL_vs_refBarNo_3->GetPointX(i), gr_mpvL_vs_refBarNo_3->GetPointY(i)/mpv_refL3);
				gr_mpvL_vs_refBarNo_3->SetPointError(i, 0, gr_mpvL_vs_refBarNo_3->GetErrorY(i)/mpv_refL3);
				gr_mpvR_vs_refBarNo_3->SetPoint(i, gr_mpvR_vs_refBarNo_3->GetPointX(i), gr_mpvR_vs_refBarNo_3->GetPointY(i)/mpv_refR3);
				gr_mpvR_vs_refBarNo_3->SetPointError(i, 0, gr_mpvR_vs_refBarNo_3->GetErrorY(i)/mpv_refR3);

				gr_mpvLR_vs_refBarNo_1->SetPoint(i, gr_mpvL_vs_refBarNo_1->GetPointX(i), 0.5* (gr_mpvL_vs_refBarNo_1->GetPointY(i)+ gr_mpvR_vs_refBarNo_1->GetPointY(i)) );
				gr_mpvLR_vs_refBarNo_2->SetPoint(i, gr_mpvL_vs_refBarNo_2->GetPointX(i), 0.5* (gr_mpvL_vs_refBarNo_2->GetPointY(i)+ gr_mpvR_vs_refBarNo_2->GetPointY(i)) );
				gr_mpvLR_vs_refBarNo_3->SetPoint(i, gr_mpvL_vs_refBarNo_3->GetPointX(i), 0.5* (gr_mpvL_vs_refBarNo_3->GetPointY(i)+ gr_mpvR_vs_refBarNo_3->GetPointY(i)) );
			}
			// ------- draw all

		 	//gStyle->SetOptFit(1111);
	
			TCanvas *c = new TCanvas();
		/*	gr_tRes_vs_refBarNo->SetMarkerStyle(20);
			gr_tRes_vs_refBarNo->SetTitle(";upstream bar # ; #sigma_{bar}(ns)");
			//gr_tRes_vs_refBarNo->SetTitleStyle(1001);
			gr_tRes_vs_refBarNo->GetXaxis()->SetLimits(0,14);
			gr_tRes_vs_refBarNo->Draw("ALP");
			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_tRes_vs_refBarNo_bar%s_Vov3.5.png",argv[1],ibar.c_str()));
			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_tRes_vs_refBarNo_bar%s_Vov3.5.pdf",argv[1],ibar.c_str()));

			outFile->cd();
			gr_tRes_vs_refBarNo->Write(Form("tRes_vs_refBar#_bar%s_Vov3.50_run%s",ibar.c_str(),argv[1]));
*/

			//----- draw Andrea's method
			
			//c->Clear();
			gr_mpvL_vs_refBarNo_1->SetMarkerStyle(20);
			gr_mpvL_vs_refBarNo_1->SetMarkerColor(kBlue);
			gr_mpvL_vs_refBarNo_1->SetLineColor(kBlue);
			gr_mpvL_vs_refBarNo_1->SetTitle(";upstream bar # ; energy MPV / MPV_{central bar}");
			gr_mpvL_vs_refBarNo_1->GetXaxis()->SetLimits(0,14);
			gr_mpvL_vs_refBarNo_1->SetMinimum(0.85);
			gr_mpvL_vs_refBarNo_1->SetMaximum(1.15);
			//TF1* f_pol2 = new TF1("f_pol2","pol2", 4, 12);
			//gr_mpvL_vs_refBarNo_1->Fit(f_pol2,"QR");
			gr_mpvL_vs_refBarNo_1->Draw("ALP");

		/*	TLatex *latexL = new TLatex(0.40,0.85,Form("#splitline{bar %s}{Vov = 3.5V  th. = 5 DAC}",ibar.c_str(),offset));
    		latexL -> SetNDC();
    		latexL -> SetTextFont(42);
    		latexL -> SetTextSize(0.04);
    		latexL -> SetTextColor(kRed);
    		latexL->Draw("same");

			//c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpvL_vs_refBarNo_bar%s_Vov3.5_th05.png",argv[1],ibar.c_str()));
			//c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpvL_vs_refBarNo_bar%s_Vov3.5_th05.pdf",argv[1],ibar.c_str()));
		*/
			outFile->cd();
			gr_mpvL_vs_refBarNo_1->Write(Form("mpvL_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));

			gr_mpvR_vs_refBarNo_1->SetMarkerStyle(20);
			gr_mpvR_vs_refBarNo_1->SetMarkerColor(kRed);
			gr_mpvR_vs_refBarNo_1->SetLineColor(kRed);
			//gr_mpvR_vs_refBarNo_1->SetTitle(";upstream bar # ; energy MPV / MPV_{central bar}");
			//gr_mpvR_vs_refBarNo_1->GetXaxis()->SetLimits(0,14);
			//gr_mpvR_vs_refBarNo_1->Fit(f_pol2,"QR");
			gr_mpvR_vs_refBarNo_1->Draw("PL same");

			//c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpvR_vs_refBarNo_bar%s_Vov3.5_th05.png",argv[1],ibar.c_str()));
			//c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpvR_vs_refBarNo_bar%s_Vov3.5_th05.pdf",argv[1],ibar.c_str()));

			outFile->cd();
			gr_mpvR_vs_refBarNo_1->Write(Form("mpvR_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));

 			gr_mpvLR_vs_refBarNo_1->SetMarkerStyle(20);
			gr_mpvLR_vs_refBarNo_1->SetMarkerColor(kBlack);
			//gr_mpvLR_vs_refBarNo_1->SetTitle(";upstream bar # ; energy MPV (a.u.)");
			//gr_tRes_vs_refBarNo->SetTitleStyle(1001);
			//gr_mpvLR_vs_refBarNo_1->GetXaxis()->SetLimits(0,14);
			gr_mpvLR_vs_refBarNo_1->Draw("LP");

			TLegend* legend = new TLegend(0.55,0.78,0.73,0.88);
          	legend->AddEntry(gr_mpvL_vs_refBarNo_1, Form("bar %sL",ibar.c_str()), "P");
			legend->AddEntry(gr_mpvR_vs_refBarNo_1, Form("bar %sR",ibar.c_str()), "P");
			legend->AddEntry(gr_mpvLR_vs_refBarNo_1, "(L+R)/2" , "P");
          	legend -> SetBorderSize(0);
          	legend->Draw();

			outFile->cd();
			gr_mpvLR_vs_refBarNo_1->Write(Form("mpvLR_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));

			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpv_vs_refBarNo_m1_bar%s_Vov3.5_th05.png",argv[1],ibar.c_str()));
			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpv_vs_refBarNo_m1_bar%s_Vov3.5_th05.pdf",argv[1],ibar.c_str()));

	 
			//----- draw mpv vs OV method

			c->Clear();
			gr_mpvL_vs_refBarNo_2->SetMarkerStyle(20);
			gr_mpvL_vs_refBarNo_2->SetMarkerColor(kBlue);
			gr_mpvL_vs_refBarNo_2->SetLineColor(kBlue);
			gr_mpvL_vs_refBarNo_2->SetTitle(";upstream bar # ; energy MPV / MPV_{central bar}");
			gr_mpvL_vs_refBarNo_2->GetXaxis()->SetLimits(0,14);
			gr_mpvL_vs_refBarNo_2->SetMinimum(0.85);
			gr_mpvL_vs_refBarNo_2->SetMaximum(1.15);
			//TF1* f_pol2 = new TF1("f_pol2","pol2", 4, 12);
			//gr_mpvL_vs_refBarNo_2->Fit(f_pol2,"QR");
			gr_mpvL_vs_refBarNo_2->Draw("ALP");

			outFile->cd();
			gr_mpvL_vs_refBarNo_2->Write(Form("mpvL_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));

			gr_mpvR_vs_refBarNo_2->SetMarkerStyle(20);
			gr_mpvR_vs_refBarNo_2->SetMarkerColor(kRed);
			gr_mpvR_vs_refBarNo_2->SetLineColor(kRed);
			gr_mpvR_vs_refBarNo_2->GetXaxis()->SetLimits(0,14);
			//gr_mpvR_vs_refBarNo_2->Fit(f_pol2,"QR");
			gr_mpvR_vs_refBarNo_2->Draw("PL same");

			outFile->cd();
			gr_mpvR_vs_refBarNo_2->Write(Form("mpvR_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));

 			gr_mpvLR_vs_refBarNo_2->SetMarkerStyle(20);
			gr_mpvLR_vs_refBarNo_2->SetMarkerColor(kBlack);
			//gr_mpvLR_vs_refBarNo_2->GetXaxis()->SetLimits(0,14);
			gr_mpvLR_vs_refBarNo_2->Draw("PL");

			legend = new TLegend(0.55,0.78,0.73,0.88);
          	legend->AddEntry(gr_mpvL_vs_refBarNo_2, Form("bar %sL",ibar.c_str()), "P");
			legend->AddEntry(gr_mpvR_vs_refBarNo_2, Form("bar %sR",ibar.c_str()), "P");
			legend->AddEntry(gr_mpvLR_vs_refBarNo_2, "(L+R)/2" , "P");
          	legend -> SetBorderSize(0);
          	legend->Draw();

			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpv_vs_refBarNo_m2_bar%s_Vov3.5_th05.png",argv[1],ibar.c_str()));
			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpv_vs_refBarNo_m2_bar%s_Vov3.5_th05.pdf",argv[1],ibar.c_str()));

			outFile->cd();
			gr_mpvLR_vs_refBarNo_2->Write(Form("mpvLR_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));


			//----- draw Sasha's method 

			c->Clear();
			gr_mpvL_vs_refBarNo_3->SetMarkerStyle(20);
			gr_mpvL_vs_refBarNo_3->SetMarkerColor(kBlue);
			gr_mpvL_vs_refBarNo_3->SetLineColor(kBlue);
			gr_mpvL_vs_refBarNo_3->SetTitle(";upstream bar # ; energy MPV / MPV_{central bar}");
			gr_mpvL_vs_refBarNo_3->GetXaxis()->SetLimits(0,14);
			gr_mpvL_vs_refBarNo_3->SetMinimum(0.85);
			gr_mpvL_vs_refBarNo_3->SetMaximum(1.15);
			//TF1* f_pol2 = new TF1("f_pol2","pol2", 4, 12);
			//gr_mpvL_vs_refBarNo_3->Fit(f_pol2,"QR");
			gr_mpvL_vs_refBarNo_3->Draw("ALP");

			outFile->cd();
			gr_mpvL_vs_refBarNo_3->Write(Form("mpvL_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));

			gr_mpvR_vs_refBarNo_3->SetMarkerStyle(20);
			gr_mpvR_vs_refBarNo_3->SetMarkerColor(kRed);
			gr_mpvR_vs_refBarNo_3->SetLineColor(kRed);
			//gr_mpvR_vs_refBarNo_3->GetXaxis()->SetLimits(0,14);
			//gr_mpvR_vs_refBarNo_3->Fit(f_pol2,"QR");
			gr_mpvR_vs_refBarNo_3->Draw("PL same");

			outFile->cd();
			gr_mpvR_vs_refBarNo_3->Write(Form("mpvR_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));

 			gr_mpvLR_vs_refBarNo_3->SetMarkerStyle(20);
			gr_mpvLR_vs_refBarNo_3->SetMarkerColor(kBlack);
			gr_mpvLR_vs_refBarNo_3->Draw("PL");

			legend = new TLegend(0.55,0.78,0.73,0.88);
          	legend->AddEntry(gr_mpvL_vs_refBarNo_3, Form("bar %sL",ibar.c_str()), "P");
			legend->AddEntry(gr_mpvR_vs_refBarNo_3, Form("bar %sR",ibar.c_str()), "P");
			legend->AddEntry(gr_mpvLR_vs_refBarNo_3, "(L+R)/2" , "P");
          	legend -> SetBorderSize(0);
          	legend->Draw();

			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpv_vs_refBarNo_m3_bar%s_Vov3.5_th05.png",argv[1],ibar.c_str()));
			c->SaveAs(Form("/eos/user/a/amkrishn/www/MTDBTL_TB_OCT2021_CERN/run%s/SummaryPlots/c_mpv_vs_refBarNo_m3_bar%s_Vov3.5_th05.pdf",argv[1],ibar.c_str()));

			outFile->cd();
			gr_mpvLR_vs_refBarNo_3->Write(Form("mpvLR_vs_refBar#_bar%s_Vov3.50_th05_run%s",ibar.c_str(),argv[1]));
			

			 
			 delete c, gr_mpvLR_vs_refBarNo_1, gr_mpvR_vs_refBarNo_1, gr_mpvL_vs_refBarNo_1, gr_mpvL_vs_refBarNo_2, gr_mpvR_vs_refBarNo_2, gr_mpvLR_vs_refBarNo_2;
			 delete gr_mpvL_vs_refBarNo_3, gr_mpvLR_vs_refBarNo_3, gr_mpvR_vs_refBarNo_3; 
			//delete gr_tRes_vs_refBarNo;

	}	

	
}
