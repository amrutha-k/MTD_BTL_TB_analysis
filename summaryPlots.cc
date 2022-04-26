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
#include "TPad.h"
#include "TMultiGraph.h"
#include "TLegend.h"

using namespace std;


int main(){


float Vbias[] = {38.5, 38.7, 39.1, 39.5, 39.9, 40.3, 40.7};   //0deg
float OV0[] = {0.86, 1.11, 1.47, 1.81, 2.13, 2.39, 2.68};
float OV1[] = {0.81, 1.07, 1.42, 1.75, 2.06, 2.32, 2.55};
	
TGraphErrors *gr_tRes_vs_vbias =  new TGraphErrors(); 
TGraphErrors *gr_MPV_vs_OV_sipm0 =  new TGraphErrors(); 
TGraphErrors *gr_MPV_vs_OV_sipm1 =  new TGraphErrors(); 
TGraphErrors *gr_tRes_vs_amp = new TGraphErrors();
TMultiGraph *mg  = new TMultiGraph();
TCanvas *c1 = new TCanvas();
TLegend* legenda = new TLegend(0.72,0.65,0.86,0.87);

TFile *fout = new TFile("summaryPlots_5e13.root","RECREATE"); 

gStyle->SetPalette(kVisibleSpectrum);

int count = 0;

for (int i = 46119; i <= 46131; i++ ){
	TFile *f = TFile::Open(Form("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/output_timeRes/Run_%d_5e13.root",i),"READ");
	if(!f || !f->IsOpen()) continue;

	//time resolution vs bias V  
	TGraphErrors *gr_tRes = (TGraphErrors*)f->Get("tRes_vs_threshold");

	double x,y;
	double min_y = 11111111111;
	int min_point = 6;
	for (int j = 0; j<15; j++)
	{
		gr_tRes->GetPoint(j,x,y);
		if (y < min_y){
			min_y = y;
        	min_point = j;
		}
	}
	
	gr_tRes_vs_vbias->SetPoint(gr_tRes_vs_vbias->GetN(),Vbias[count],1000*min_y);
	gr_tRes_vs_vbias->SetPointError(gr_tRes_vs_vbias->GetN()-1,0,1000*gr_tRes->GetErrorY(min_point));   

	//MPV vs OV
	TH1F *h_MPV_simp0 = (TH1F*)f->Get("MPV_amp0");
	TH1F *h_MPV_simp1 = (TH1F*)f->Get("MPV_amp1");
	TF1 *fit_MPV_sipm0 = (TF1*)h_MPV_simp0->GetFunction("l2");
	TF1 *fit_MPV_sipm1 = (TF1*)h_MPV_simp1->GetFunction("l3");
	
	gr_MPV_vs_OV_sipm0->SetPoint(gr_MPV_vs_OV_sipm0->GetN(),OV0[count],fit_MPV_sipm0->GetParameter(1));
	gr_MPV_vs_OV_sipm0->SetPointError(gr_MPV_vs_OV_sipm0->GetN()-1,0,fit_MPV_sipm0->GetParError(1));
	gr_MPV_vs_OV_sipm1->SetPoint(gr_MPV_vs_OV_sipm1->GetN(),OV1[count],fit_MPV_sipm1->GetParameter(1));
	gr_MPV_vs_OV_sipm1->SetPointError(gr_MPV_vs_OV_sipm1->GetN()-1,0,fit_MPV_sipm1->GetParError(1));

	//time resolution vs threshold all in one graph
	TGraphErrors *gr_tRes_vs_thr = (TGraphErrors*)f->Get("tRes_vs_threshold");
	
	legenda->AddEntry(gr_tRes_vs_thr, Form("Vbias = %.1f V",Vbias[count]), "PE");
	mg->Add(gr_tRes_vs_thr,"PL");

	if (i == 46123) gr_tRes_vs_amp = (TGraphErrors*)f->Get("gr_tRes_vs_amp_bestthr");

	count++;

	delete h_MPV_simp0, h_MPV_simp1, fit_MPV_sipm0, fit_MPV_sipm1, gr_tRes_vs_thr;
}
mg->SetTitle("0 deg; threshold (mV); time resolution (ns)");
mg->GetYaxis()->SetRangeUser(0.15,0.4);
mg->Draw("A pmc plc");
legenda -> SetBorderSize(0);
legenda->Draw("same");

c1->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_thr_0deg_5e13.png");
c1->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_thr_0deg_5e13.pdf");

TCanvas *c = new TCanvas();
c->SetGrid();
gr_tRes_vs_vbias->SetMarkerStyle(20);
gr_tRes_vs_vbias->SetLineColor(kRed);
gr_tRes_vs_vbias->SetMarkerColor(kRed);
gr_tRes_vs_vbias->SetTitle("0 deg; bias (V); time resolution (ps)");
gr_tRes_vs_vbias->Draw("ALP");

c->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_Vbias_0deg_5e13.png");
c->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_Vbias_0deg_5e13.pdf");


// ------------------------------------- 52 deg ----------------------------------------------

float OV0_52deg[] = {0.78, 1.01, 1.44, 1.71, 2.00, 2.30, 2.52};
float OV1_52deg[] = {0.82, 1.01, 1.45, 1.75, 2.05, 2.39, 2.51};
	
TGraphErrors *gr_tRes_vs_vbias_52deg =  new TGraphErrors(); 
TGraphErrors *gr_MPV_vs_OV_sipm0_52deg =  new TGraphErrors(); 
TGraphErrors *gr_MPV_vs_OV_sipm1_52deg =  new TGraphErrors(); 
TMultiGraph *mg_tRes_vs_thr_52deg  = new TMultiGraph();
TLegend* legendb = new TLegend(0.14,0.65,0.3,0.87);
TGraphErrors *gr_tRes_vs_amp_52deg =  new TGraphErrors(); 

TCanvas *c2 = new TCanvas();

count = 0;

for (int i = 46040; i <= 46049; i++ ){
	TFile *f = TFile::Open(Form("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/output_timeRes/Run_%d_5e13.root",i),"READ");
	if(!f || !f->IsOpen()) continue;

	//time resolution vs bias V  
	TGraphErrors *gr_tRes = (TGraphErrors*)f->Get("tRes_vs_threshold");

	double x,y;
	double min_y = 11111111111;
	int min_point = 6;
	for (int j = 0; j<15; j++)
	{
		gr_tRes->GetPoint(j,x,y);
		if (y < min_y){
			min_y = y;
        	min_point = j;
		}
	}
	
	gr_tRes_vs_vbias_52deg->SetPoint(gr_tRes_vs_vbias_52deg->GetN(),Vbias[count],min_y*1000);
	gr_tRes_vs_vbias_52deg->SetPointError(gr_tRes_vs_vbias_52deg->GetN()-1,0,1000*gr_tRes->GetErrorY(min_point));   

	//MPV vs OV
	TH1F *h_MPV_simp0 = (TH1F*)f->Get("MPV_amp0");
	TH1F *h_MPV_simp1 = (TH1F*)f->Get("MPV_amp1");
	TF1 *fit_MPV_sipm0 = (TF1*)h_MPV_simp0->GetFunction("l2");
	TF1 *fit_MPV_sipm1 = (TF1*)h_MPV_simp1->GetFunction("l3");
	
	gr_MPV_vs_OV_sipm0_52deg->SetPoint(gr_MPV_vs_OV_sipm0_52deg->GetN(),OV0_52deg[count],fit_MPV_sipm0->GetParameter(1));
	gr_MPV_vs_OV_sipm0_52deg->SetPointError(gr_MPV_vs_OV_sipm0_52deg->GetN()-1,0,fit_MPV_sipm0->GetParError(1));
	gr_MPV_vs_OV_sipm1_52deg->SetPoint(gr_MPV_vs_OV_sipm1_52deg->GetN(),OV1_52deg[count],fit_MPV_sipm1->GetParameter(1));
	gr_MPV_vs_OV_sipm1_52deg->SetPointError(gr_MPV_vs_OV_sipm1_52deg->GetN()-1,0,fit_MPV_sipm1->GetParError(1));

	//time resolution vs threshold all in one graph
	TGraphErrors *gr_tRes_vs_thr = (TGraphErrors*)f->Get("tRes_vs_threshold");
	legendb->AddEntry(gr_tRes_vs_thr, Form("Vbias = %.1f V",Vbias[count]), "PE");
	mg_tRes_vs_thr_52deg->Add(gr_tRes_vs_thr,"PL");

	if (i == 46043) gr_tRes_vs_amp_52deg = (TGraphErrors*)f->Get("gr_tRes_vs_amp_bestthr");

	count++;

	delete h_MPV_simp0, h_MPV_simp1, fit_MPV_sipm0, fit_MPV_sipm1, gr_tRes_vs_thr;
}
mg_tRes_vs_thr_52deg->SetTitle("52 deg; threshold (mV); time resolution (ns)");
mg_tRes_vs_thr_52deg->GetYaxis()->SetRangeUser(0.07,0.14);
mg_tRes_vs_thr_52deg->Draw("A pmc plc");

legendb->SetBorderSize(0);
legendb->Draw("same");

c2->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_thr_52deg_5e13.png");
c2->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_thr_52deg_5e13.pdf");

c2->Clear();
gr_tRes_vs_amp_52deg->SetMarkerStyle(20);
gr_tRes_vs_amp_52deg->SetLineColor(kBlue);
gr_tRes_vs_amp_52deg->SetMarkerColor(kBlue);
gr_tRes_vs_amp_52deg->SetTitle("; amp/MPV ; time resolution (ns)");
gr_tRes_vs_amp_52deg->GetYaxis()->SetRangeUser(0.04, 0.23);
gr_tRes_vs_amp_52deg->Draw("ALP");

gr_tRes_vs_amp->SetMarkerStyle(20);
gr_tRes_vs_amp->SetLineColor(kRed);
gr_tRes_vs_amp->SetMarkerColor(kRed);
gr_tRes_vs_amp->Draw("LP same");

TLegend* legend0 = new TLegend(0.72,0.75,0.86,0.85);
legend0->AddEntry(gr_tRes_vs_amp,"0 deg","PE");
legend0->AddEntry(gr_tRes_vs_amp_52deg, "52 deg","PE");
legend0->SetBorderSize(0);
legend0->Draw("same");

c2->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_amp_5e13.png");
c2->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_amp_5e13.pdf");



TGraphErrors *gr_tRes_vs_amp_52deg_shifted = new TGraphErrors();

for (int ii = 0; ii < gr_tRes_vs_amp_52deg->GetN(); ii++){
	gr_tRes_vs_amp_52deg_shifted->SetPoint(ii, 2.5*gr_tRes_vs_amp_52deg->GetPointX(ii),gr_tRes_vs_amp_52deg->GetPointY(ii));
	gr_tRes_vs_amp_52deg_shifted->SetPointError(ii, gr_tRes_vs_amp_52deg->GetErrorX(ii),gr_tRes_vs_amp_52deg->GetErrorY(ii));
}


c2->Clear();

TF1 *fitFunc = new TF1("fitFunc","[0]/pow(x,[1])",0.8,6.1);
fitFunc -> SetLineColor(kBlack);
fitFunc->SetRange(2.2,6.1);

gr_tRes_vs_amp_52deg_shifted->SetMarkerStyle(20);
gr_tRes_vs_amp_52deg_shifted->SetLineColor(kBlue);
gr_tRes_vs_amp_52deg_shifted->SetMarkerColor(kBlue);
gr_tRes_vs_amp_52deg_shifted->SetTitle("; amp/MPV ; time resolution (ns)");
gr_tRes_vs_amp_52deg_shifted->GetYaxis()->SetRangeUser(0.04, 0.23);
gr_tRes_vs_amp_52deg_shifted->GetXaxis()->SetLimits(0,7);
gr_tRes_vs_amp_52deg_shifted->Fit(fitFunc,"QNRS");
fitFunc->SetParameter(0,fitFunc->GetParameter(0));
fitFunc->SetParameter(1,fitFunc->GetParameter(1));
gr_tRes_vs_amp_52deg_shifted->Fit(fitFunc,"R");
gr_tRes_vs_amp_52deg_shifted->Draw("ALP");

cout<<"chi2 blue:"<<fitFunc->GetChisquare()<<endl;

gr_tRes_vs_amp->SetMarkerStyle(20);
gr_tRes_vs_amp->SetLineColor(kRed);
gr_tRes_vs_amp->SetMarkerColor(kRed);
fitFunc->SetRange(0.8,2.4);
gr_tRes_vs_amp->Fit(fitFunc,"QNRS");
fitFunc->SetParameter(0,fitFunc->GetParameter(0));
fitFunc->SetParameter(1,fitFunc->GetParameter(1));
gr_tRes_vs_amp->Fit(fitFunc,"R");
gr_tRes_vs_amp->Draw("LP same");

cout<<"chi2 red:"<<fitFunc->GetChisquare()<<endl;

legend0->Clear();
legend0->AddEntry(gr_tRes_vs_amp,"0 deg","PE");
legend0->AddEntry(gr_tRes_vs_amp_52deg, "52 deg","PE");
legend0->SetBorderSize(0);
legend0->Draw("same");

c2->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_ampShifted_5e13.png");
c2->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_ampShifted_5e13.pdf");

TGraphErrors *gr_tRes_vs_amp_combined = new TGraphErrors();
for(int k = 0; k<6; k++){
	gr_tRes_vs_amp_combined->SetPoint(gr_tRes_vs_amp_combined->GetN(), gr_tRes_vs_amp->GetPointX(k),gr_tRes_vs_amp->GetPointY(k));
	gr_tRes_vs_amp_combined->SetPointError(gr_tRes_vs_amp_combined->GetN()-1, gr_tRes_vs_amp->GetErrorX(k),gr_tRes_vs_amp->GetErrorY(k));
}
for(int n = 0; n<6; n++){
	gr_tRes_vs_amp_combined->SetPoint(gr_tRes_vs_amp_combined->GetN(), 2.2*gr_tRes_vs_amp_52deg->GetPointX(n),gr_tRes_vs_amp_52deg->GetPointY(n));
	gr_tRes_vs_amp_combined->SetPointError(gr_tRes_vs_amp_combined->GetN()-1, gr_tRes_vs_amp_52deg->GetErrorX(n),gr_tRes_vs_amp_52deg->GetErrorY(n));
}

gStyle->SetOptFit(1111);

c2->Clear();
gr_tRes_vs_amp_combined->SetMarkerStyle(20);
gr_tRes_vs_amp_combined->SetLineColor(kBlack);
gr_tRes_vs_amp_combined->SetMarkerColor(kBlack);
gr_tRes_vs_amp_combined->SetTitle("; amp/MPV ; time resolution (ns)");
fitFunc -> SetLineColor(kRed);
fitFunc->SetRange(0.8,5.3);
gr_tRes_vs_amp_combined->Fit(fitFunc,"QNRS");
fitFunc->SetParameter(0,fitFunc->GetParameter(0));
fitFunc->SetParameter(1,fitFunc->GetParameter(1));
gr_tRes_vs_amp_combined->Fit(fitFunc,"QR");
gr_tRes_vs_amp_combined->Draw("AP");

c2->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_ampCombined_5e13.png");
c2->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_ampCombined_5e13.pdf");


TCanvas *c3 = new TCanvas();
c3->SetGrid();
gr_tRes_vs_vbias_52deg->SetMarkerStyle(20);
gr_tRes_vs_vbias_52deg->SetLineColor(kRed);
gr_tRes_vs_vbias_52deg->SetMarkerColor(kRed);
gr_tRes_vs_vbias_52deg->SetTitle("52 deg; bias (V); time resolution (ps)");
gr_tRes_vs_vbias_52deg->Draw("ALP");
c3->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_Vbias_52deg_5e13.png");
c3->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_Vbias_52deg_5e13.pdf");

//draw time resolution vs bias voltage for 0 deg and 52 deg on the same canvas
TCanvas *c4 = new TCanvas();
c4->SetGrid();
gr_tRes_vs_vbias_52deg->SetMarkerStyle(20);
gr_tRes_vs_vbias_52deg->SetLineColor(kRed);
gr_tRes_vs_vbias_52deg->SetMarkerColor(kRed);
gr_tRes_vs_vbias_52deg->GetYaxis()->SetRangeUser(0.05,.35);
gr_tRes_vs_vbias_52deg->SetTitle("; bias (V); time resolution (ps)");
gr_tRes_vs_vbias_52deg->Draw("ALP");

gr_tRes_vs_vbias->SetMarkerStyle(20);
gr_tRes_vs_vbias->SetLineColor(kBlue);
gr_tRes_vs_vbias->SetMarkerColor(kBlue);
gr_tRes_vs_vbias->Draw("LP same");

TLegend* legend = new TLegend(0.18,0.75,0.34,0.85);
legend->AddEntry(gr_tRes_vs_vbias,"0 deg","PE");
legend->AddEntry(gr_tRes_vs_vbias_52deg, "52 deg","PE");
legend->SetBorderSize(0);
legend->Draw("same");

c4->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_Vbias_5e13.png");
c4->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/tRes_vs_Vbias_5e13.pdf");

//--------------------------------------------------------------------------------------------------------------------------

//for OV correction 
double x1, y1, x2, y2;
gr_MPV_vs_OV_sipm0->GetPoint(0,x1,y1);
gr_MPV_vs_OV_sipm0->GetPoint(1,x2,y2);
TF1 *f_OV0_corr = new TF1("f_OV0_corr","((x - [0])*([3] - [2])/([1] - [0])) + [2]",0,3); 
f_OV0_corr->SetParameter(0,x1);
f_OV0_corr->SetParameter(1,x2);
f_OV0_corr->SetParameter(2,y1);
f_OV0_corr->SetParameter(3,y2);
int iter = 0;
for(auto ov : OV0){
	float ov0_corr = ov - f_OV0_corr->GetX(gr_MPV_vs_OV_sipm0->GetPointY(iter));
	//cout<<"OV0 corrected 0 deg: "<<ov0_corr<< endl;
	iter ++;
} 

gr_MPV_vs_OV_sipm1->GetPoint(0,x1,y1);
gr_MPV_vs_OV_sipm1->GetPoint(1,x2,y2);
TF1 *f_OV1_corr = new TF1("f_OV1_corr","((x - [0])*([3] - [2])/([1] - [0])) + [2]",0,3); 
f_OV1_corr->SetParameter(0,x1);
f_OV1_corr->SetParameter(1,x2);
f_OV1_corr->SetParameter(2,y1);
f_OV1_corr->SetParameter(3,y2);
iter = 0;
for(auto ov : OV1){
	float ov1_corr = ov - f_OV1_corr->GetX(gr_MPV_vs_OV_sipm1->GetPointY(iter));
	//cout<<"OV1 corrected 0 deg: "<<ov1_corr<< endl;
	iter ++;
} 

gr_MPV_vs_OV_sipm0_52deg->GetPoint(0,x1,y1);
gr_MPV_vs_OV_sipm0_52deg->GetPoint(1,x2,y2);
TF1 *f_OV0_corr_52deg = new TF1("f_OV0_corr_52deg","((x - [0])*([3] - [2])/([1] - [0])) + [2]",0,3); 
f_OV0_corr_52deg->SetParameter(0,x1);
f_OV0_corr_52deg->SetParameter(1,x2);
f_OV0_corr_52deg->SetParameter(2,y1);
f_OV0_corr_52deg->SetParameter(3,y2);
iter = 0;
for(auto ov : OV0_52deg){
	float ov0_corr_52deg = ov - f_OV0_corr_52deg->GetX(gr_MPV_vs_OV_sipm0_52deg->GetPointY(iter));
	//cout<<"OV0 corrected 52 deg: "<<ov0_corr_52deg<< endl;
	iter ++;
} 

gr_MPV_vs_OV_sipm1_52deg->GetPoint(0,x1,y1);
gr_MPV_vs_OV_sipm1_52deg->GetPoint(1,x2,y2);
TF1 *f_OV1_corr_52deg = new TF1("f_OV1_corr_52deg","((x - [0])*([3] - [2])/([1] - [0])) + [2]",0,3); 
f_OV1_corr_52deg->SetParameter(0,x1);
f_OV1_corr_52deg->SetParameter(1,x2);
f_OV1_corr_52deg->SetParameter(2,y1);
f_OV1_corr_52deg->SetParameter(3,y2);
iter = 0;
for(auto ov : OV1_52deg){
	float ov1_corr_52deg = ov - f_OV1_corr_52deg->GetX(gr_MPV_vs_OV_sipm1_52deg->GetPointY(iter));
	//cout<<"OV1 corrected 52 deg: "<<ov1_corr_52deg<< endl;
	iter ++;
} 

//draw MPV vs OV
c->cd();
c->Clear();
c->SetGrid();
gr_MPV_vs_OV_sipm0->GetYaxis()->SetRangeUser(10,215);
gr_MPV_vs_OV_sipm0->SetMarkerStyle(20);
gr_MPV_vs_OV_sipm0->SetLineColor(kRed);
gr_MPV_vs_OV_sipm0->SetMarkerColor(kRed);
gr_MPV_vs_OV_sipm0->SetTitle("; OV (V); amp MPV (mV)");
gr_MPV_vs_OV_sipm0->Draw("ALP");
f_OV0_corr->SetLineStyle(kDashed);
f_OV0_corr->SetLineColor(kBlack);
f_OV0_corr->Draw("same");

gr_MPV_vs_OV_sipm0_52deg->SetMarkerStyle(20);
gr_MPV_vs_OV_sipm0_52deg->SetLineColor(kBlue);
gr_MPV_vs_OV_sipm0_52deg->SetMarkerColor(kBlue);
gr_MPV_vs_OV_sipm0_52deg->Draw("LP same");
f_OV0_corr_52deg->SetLineStyle(kDashed);
f_OV0_corr_52deg->SetLineColor(kBlack);
f_OV0_corr_52deg->Draw("same");

TLegend* legendc = new TLegend(0.14,0.75,0.3,0.85);
legendc->AddEntry(gr_MPV_vs_OV_sipm0,"0 deg","PE");
legendc->AddEntry(gr_MPV_vs_OV_sipm0_52deg, "52 deg","PE");
legendc->SetBorderSize(0);
legendc->Draw("same");

c->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/MPV_vs_OV_sipm0_5e13.png");
c->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/MPV_vs_OV_sipm0_5e13.pdf");

c->cd();
c->Clear();
c->SetGrid();
gr_MPV_vs_OV_sipm1->GetYaxis()->SetRangeUser(10,215);
gr_MPV_vs_OV_sipm1->SetMarkerStyle(20);
gr_MPV_vs_OV_sipm1->SetLineColor(kRed);
gr_MPV_vs_OV_sipm1->SetMarkerColor(kRed);
gr_MPV_vs_OV_sipm1->SetTitle("; OV (V); amp MPV (mV)");
gr_MPV_vs_OV_sipm1->Draw("ALP");
f_OV1_corr->SetLineStyle(kDashed);
f_OV1_corr->SetLineColor(kBlack);
f_OV1_corr->Draw("same");

gr_MPV_vs_OV_sipm1_52deg->SetMarkerStyle(20);
gr_MPV_vs_OV_sipm1_52deg->SetLineColor(kBlue);
gr_MPV_vs_OV_sipm1_52deg->SetMarkerColor(kBlue);
gr_MPV_vs_OV_sipm1_52deg->Draw("LP same");
f_OV1_corr_52deg->SetLineStyle(kDashed);
f_OV1_corr_52deg->SetLineColor(kBlack);
f_OV1_corr_52deg->Draw("same");

TLegend* legendd = new TLegend(0.14,0.75,0.3,0.85);
legendd->AddEntry(gr_MPV_vs_OV_sipm0,"0 deg","PE");
legendd->AddEntry(gr_MPV_vs_OV_sipm0_52deg, "52 deg","PE");
legendd->SetBorderSize(0);
legendd->Draw("same");

c->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/MPV_vs_OV_sipm1_5e13.png");
c->SaveAs("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/summaryPlots/MPV_vs_OV_sipm1_5e13.pdf");

fout->cd();
gr_MPV_vs_OV_sipm0->Write();
gr_MPV_vs_OV_sipm1->Write();
gr_tRes_vs_vbias->Write();
mg->Write();
gr_MPV_vs_OV_sipm0_52deg->Write();
gr_MPV_vs_OV_sipm1_52deg->Write();
gr_tRes_vs_vbias_52deg->Write();
mg_tRes_vs_thr_52deg->Write();

delete gr_MPV_vs_OV_sipm0, gr_MPV_vs_OV_sipm1, gr_tRes_vs_vbias, gr_MPV_vs_OV_sipm0_52deg, gr_MPV_vs_OV_sipm1_52deg, gr_tRes_vs_vbias_52deg, mg, mg_tRes_vs_thr_52deg;
delete gr_tRes_vs_amp_52deg, gr_tRes_vs_amp, legend, legend0, legendb, legenda, legendc, legendd;

return 0;
}
