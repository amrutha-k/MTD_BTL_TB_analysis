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

using namespace std;

struct anEvent{
  Double_t        channel[2][1024];
  Double_t        times[1024];
  Int_t           evt;
  Double_t        xtrk;
  Double_t        ytrk;
  Double_t        ampMCP;
  Double_t        timeMCP;
  Double_t        rmsMCP;
  Double_t        chi2MCP;
  Double_t        ampLo[2];
  Double_t        timeLo[2];
  Double_t        rmsLo[2];
  Double_t        chi2Lo[2];
} event;

vector<vector<double>> totVals;

void Initialize(){

totVals.clear();

}

bool cuts(double x, double y, double tmcp){

  bool x_cut = true, y_cut = true, mcp_cut = true, amp0_cut = true, amp1_cut = true;

  if(x < 25. ||  x > 29.)  x_cut = false;   //track cuts
  if(y < 7. || y > 10.)  y_cut = false;
  if(tmcp < 100 ||  tmcp > 118)  mcp_cut = false;  
  bool pass = (x_cut && y_cut && mcp_cut);
  return pass;
}


double getmax(TGraph *gr) {
  /*double mean = 0;
  for(int i=0; i!=250; ++i)
    mean += gr->GetPointY(i);
  mean /= 250;*/

  double max = -11111111111;
  int p = 0;
  for(int i=250; i!=350; ++i) {
    if(gr->GetPointY(i)>max){
      max = gr->GetPointY(i);
      p = i;
    }
  }
  //return  -( min - mean );
  return p;
}

vector<vector<double>> tot(int th){

//to find the time over threshold for various threshold values via linear interpolation of amplitude values
// just below and above the threshold, after applying delay and invert technique

      Double_t vbelow0 = -1111, nbelow0 = 200, vabove0 = 1111, nabove0 = 200; 
      Double_t vbelow1 = -1111, nbelow1 = 200, vabove1 = 1111, nabove1 = 200; 
      Double_t m0, cy0, m1, cy1, y0, y1, xb0, xa0, yb0, yb1, ya0, xb1, xa1, ya1, t1, t0; 
      
      TGraph *grWF0 = new TGraph(1024,event.times,event.channel[0]);
      TGraph *grWF1 = new TGraph(1024,event.times,event.channel[1]);

      //TGraph *gr_t_plus0 = new TGraph();
      //TGraph *gr_t_minus0 = new TGraph();
      //TGraph *gr_t_plus1 = new TGraph();
      //TGraph *gr_t_minus1 = new TGraph();
      TGraph *gr_DI0 = new TGraph();
      TGraph *gr_DI1 = new TGraph();
      //TCanvas *c1 = new TCanvas();
      //TCanvas *c2 = new TCanvas();

      for (int n = 0; n<1024; n++){
       /*gr_t_plus0->SetPoint(n,event.times[n] + 0.2,event.channel[0][n]);
       gr_t_minus0->SetPoint(n,event.times[n] - 0.2,event.channel[0][n]);
       gr_t_plus1->SetPoint(n,event.times[n] + 0.2,event.channel[1][n]);
       gr_t_minus1->SetPoint(n,event.times[n] - 0.2,event.channel[1][n]);
      */

      //}

      //for (int j=0; j<1024; j++){
        double yp_0 = grWF0->Eval(event.times[n]+0.2);
        double ym_0 = grWF0->Eval(event.times[n]-0.2);
        double yp_1 = grWF1->Eval(event.times[n]+0.2);
        double ym_1 = grWF1->Eval(event.times[n]-0.2);

        gr_DI0->SetPoint(n,event.times[n],(yp_0 - ym_0)/0.4);
        gr_DI1->SetPoint(n,event.times[n],(yp_1 - ym_1)/0.4);

      }

   /*   grWF0->SetMarkerStyle(kCircle);
      grWF0->SetMarkerColor(kCyan);
      grWF0->SetMarkerSize(0.2);
      grWF1->SetMarkerStyle(kCircle);
      grWF1->SetMarkerColor(kBlue);
      grWF1->SetMarkerSize(0.2);
      gr_DI0->SetMarkerStyle(kCircle);
      gr_DI0->SetMarkerColor(kRed);
      gr_DI0->SetMarkerSize(0.2);
      gr_DI1->SetMarkerStyle(kCircle);
      gr_DI1->SetMarkerColor(kGreen);
      gr_DI1->SetMarkerSize(0.2);  */

     /*c1->Divide(2,2);
      c1->cd(1);
      grWF0->Draw("ALP");
      c1->cd(2);
      grWF1->Draw("ALP"); 

      //c2->Divide(2,1);
      c1->cd(3);
      gr_DI0->Draw("ALP");
      //gr_DI0->GetXaxis()->SetRangeUser(20,70);
      //gr_DI0->Draw("AP");
      //c1->Update();
      c1->cd(4);
      gr_DI1->Draw("ALP");
      //gr_DI1->GetXaxis()->SetRangeUser(20,70);
      //gr_DI1->Draw("AP");
      //c1->Update();
      
      ifstream file(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46123/c_WFDI_test_th%d.png",th));
      if(!file){
       // c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46076/c_WF_test_th%d.png",th));
       c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46123/c_WFDI_test_th%d.png",th));
      } 
      //c1->Clear();
      //c2->Clear();
      //fout->cd();
      //grWF0->Write("WF0");
      //grWF1->Write("WF1");
      //gr_DI0->Write("DI0");
      //gr_DI1->Write("DI1");
      
    */
      double point0 = getmax(gr_DI0);
      double point1 = getmax(gr_DI1);

      for (int i = 240; i<point0; i++){
        y0 = gr_DI0->GetPointY(i) - th;
        if (y0 <= 0 && y0 > vbelow0) {
            vbelow0 = y0;
            nbelow0 = i; 
          }

        if (y0 > 0 && y0 < vabove0){
            vabove0 = y0;
            nabove0 = i;
          }
        }

      for (int i = 240; i<point1; i++){
        y1 = gr_DI1->GetPointY(i) - th;
        if (y1 <= 0 && y1 > vbelow1) {
                  vbelow1 = y1;
                  nbelow1 = i;
                }
        if (y1 > 0 && y1 < vabove1){
                  vabove1 = y1;
                  nabove1 = i;
                }
     }


        xb0 = gr_DI0->GetPointX(nbelow0);
        xa0 = gr_DI0->GetPointX(nabove0);
        yb0 = gr_DI0->GetPointY(nbelow0);
        ya0 = gr_DI0->GetPointY(nabove0);

        m0 = (ya0 - yb0) / (xa0 - xb0);
        cy0 = ya0 - m0*xa0; 

        t0 = (th - cy0)/m0;  //time over threshold (tot)

        xb1 = gr_DI1->GetPointX(nbelow1);
        xa1 = gr_DI1->GetPointX(nabove1);
        yb1 = gr_DI1->GetPointY(nbelow1);
        ya1 = gr_DI1->GetPointY(nabove1);

        m1 = (ya1 - yb1) / (xa1 - xb1);
        cy1 = ya1 - m1*xa1; 

        t1 = (th - cy1)/m1;

        totVals.push_back({t0,t1});

        delete grWF0, grWF1, gr_DI0, gr_DI1;

        return totVals;

}

//----------------sliced resolution-----------------------------

    vector<int> xbins_a{1,4,6,9,14,25,49};
    vector<int> xbins_p{8,17,20,23,26,53};

    vector<double> xvalues(int a, double ll, double ul)
     {
      
     if (a == 0){
      vector<double> vals_a;
      double x;
      for(int i = 0; i<xbins_a.size()-1; i++){
        //x = (25+xbins_a[i]*0.043) + (25+xbins_a[i+1]*0.043);
        x = (ll+xbins_a[i]*((ul - ll)/50)) + (ll+xbins_a[i+1]*((ul - ll)/50));
        vals_a.push_back(x/2);
      }
      return vals_a;
    }

    if (a == 1){
      vector<double> vals_p;
      double x;
      for(int i = 0; i<xbins_p.size()-1; i++){
        x = (22+xbins_p[i]*0.1246) + (22+xbins_p[i+1]*0.1246);
        vals_p.push_back(x/2);
      }
      return vals_p;
    }

  }

  vector<double> xerrors(double ll, double ul)
     {
      
      vector<double> xerr;
      double e;
      for(int i = 0; i<xbins_a.size()-1; i++){
        //x = (25+xbins_a[i]*0.043) + (25+xbins_a[i+1]*0.043);
        e = (ll+xbins_a[i+1]*((ul - ll)/50)) - (ll+xbins_a[i]*((ul - ll)/50));
        xerr.push_back(e/2);
      }
      return xerr;
    }



//---------- main code ------------------

int main(int argc, char* argv[]){

  TFile *fout = new TFile(Form("./test/Run_%s_5e13.root",argv[1]),"RECREATE");  
  TFile *infile = TFile::Open(Form("/eos/user/a/amkrishn/BTL/MTDTB_FNAL_April-May2021/Run_%s_bar5E13.root",argv[1]), "READ");
  if(!infile || !infile->IsOpen()) cout<<"file does not exist!";
  TTree *mytree = (TTree*)infile->Get("bar");


  mytree->SetBranchAddress("channel",&event.channel);
  mytree->SetBranchAddress("times",&event.times);
  mytree->SetBranchAddress("evt",&event.evt);
  mytree->SetBranchAddress("xtrk",&event.xtrk);
  mytree->SetBranchAddress("ytrk",&event.ytrk);
  mytree->SetBranchAddress("ampMCP",&event.ampMCP);
  mytree->SetBranchAddress("timeMCP",&event.timeMCP);
  mytree->SetBranchAddress("rmsMCP",&event.rmsMCP);
  mytree->SetBranchAddress("chi2MCP",&event.chi2MCP);
  mytree->SetBranchAddress("ampLo",&event.ampLo);
  mytree->SetBranchAddress("timeLo",&event.timeLo);
  mytree->SetBranchAddress("rmsLo",&event.rmsLo);
  mytree->SetBranchAddress("chi2Lo",&event.chi2Lo);

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  //--------- create histograms ---------------

  TCanvas *c = new TCanvas();
  TH2D *amp_sipm0_vs_sipm1 = new TH2D("amp_sipm0_vs_sipm1", "; amp_SiPM0 (mV); amp_SiPM1 (mV)", 1000, -100, 2000, 1000, -100, 2000);
  TH2D *amp_sipm0_vs_sipm1Corr = new TH2D("amp_sipm0_vs_sipm1Corr", "; amp_SiPM0 (mV); amp_SiPM1 (mV)", 1000, -100, 2000, 1000, -100, 2000);
  TH2D *MCP_time_vs_amp = new TH2D("MCP_time_vs_amp", "; time_MCP; amp_MCP", 1000, 0, 200, 1000, -100, 1000);
  TH2D *amp_vs_x = new TH2D("amp_vs_x", "; x (mm); amp_SiPM0 (mV)", 80, 20, 35, 50, 10, 150);
  TH2D *amp_vs_y = new TH2D("amp_vs_y", "; y (mm); amp_SiPM0 (mV)", 100, 0, 20, 50, 10, 150);
  TProfile2D* hprof2D = new TProfile2D("hprof2D", "Profile of Amp vs x and y tracks", 120, 20, 30, 120, 0, 40, 0, 300);
  TH2D *hBeam = new TH2D("hBeam", "Beam; x (mm); y (mm); events", 120, 10, 40, 120, 0, 40);
  TH1D *MPV_amp0 = new TH1D("MPV_amp0", "SiPM0; amplitude (mV); counts", 400, 50, 700);
  TH1D *MPV_amp1 = new TH1D("MPV_amp1", "SiPM1; amplitude (mV); counts", 400, 50, 700);
  TH1D *MCP_amp = new TH1D("MCP_amp", "MCP ; amplitude (mV); counts", 200, 30, 600);
  TProfile2D* MCP_efficiency = new TProfile2D("MCP_efficiency", "MCP hit position; x (mm); y (mm)", 60, 0, 40, 60, 0, 30, 0, 1.1);
  TProfile2D* MCP_hitsProfile = new TProfile2D("MCP_hitsProfile", "MCP hits profile; x (mm); y (mm)", 60, 20, 30, 60, 12, 18, 5, 250);

  TGraphErrors *tRes_vs_threshold = new TGraphErrors();


  Double_t tdiff[45000], delta_t, y, tdiffPosCorr, y_amp, t0_AmpCorrected[45000], t1_AmpCorrected[45000], tdiffAmpPosCorr, y1; 
  vector<vector<double>> t;
  vector<double> y_mean,x_arr_pos,x_arr_amp0, x_arr_amp1;

  Long64_t nentries = mytree->GetEntries();

    // -------------- 1st loop over events ---------------------
  cout<<"1st loop over events"<<endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  mytree->GetEntry(jentry);
  if (jentry % 10000 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";

  //histograms to determine the cuts to be applied so that we select only the good events
  amp_sipm0_vs_sipm1->Fill(event.ampLo[0], event.ampLo[1]);
  MCP_time_vs_amp->Fill(event.timeMCP,event.ampMCP);
  hprof2D->Fill(event.xtrk, event.ytrk, event.ampLo[1]);
  hBeam->Fill(event.xtrk, event.ytrk);
  amp_vs_x->Fill(event.xtrk, event.ampLo[0]);
  amp_vs_y->Fill(event.ytrk, event.ampLo[0]);
  

  int weight = 0;
  if (event.ampMCP > 100) weight = 1;
  MCP_efficiency->Fill(event.xtrk, event.ytrk, weight);

  //apply cuts -- select only the good events
  bool selections = cuts(event.xtrk, event.ytrk, event.timeMCP);
  if (selections == false) continue;

  MPV_amp0->Fill(event.ampLo[0]);
  MPV_amp1->Fill(event.ampLo[1]);
  amp_sipm0_vs_sipm1Corr->Fill(event.ampLo[0], event.ampLo[1]);
  MCP_amp->Fill(event.ampMCP);
  MCP_hitsProfile->Fill(event.xtrk, event.ytrk, event.ampMCP);


}

  TF1 *l0 = new TF1("l0","landau", MPV_amp0->GetMean() - 2*MPV_amp0->GetRMS() , MPV_amp0->GetMean() + 6*MPV_amp0->GetRMS());
  TF1 *l1 = new TF1("l1","landau", MPV_amp1->GetMean() - 2*MPV_amp1->GetRMS() , MPV_amp1->GetMean() + 6*MPV_amp1->GetRMS());
  MPV_amp0->Fit(l0,"QR");
  MPV_amp1->Fit(l1,"QR");
  TF1 *l2 = new TF1("l2","landau", 0.8*l0->GetParameter(1) , 3*l0->GetParameter(1) );
  TF1 *l3 = new TF1("l3","landau", 0.8*l1->GetParameter(1) , 3*l1->GetParameter(1) );    
  MPV_amp0->Fit(l2,"QR");
  MPV_amp1->Fit(l3,"QR");

  double amp0_ll = round(1/(0.8*l2->GetParameter(1))*100000) / 100000;
  double amp0_ul = round(1/(3*l2->GetParameter(1))*100000) / 100000;
  double amp1_ll = round(1/(0.8*l3->GetParameter(1))*100000) / 100000;
  double amp1_ul = round(1/(3*l3->GetParameter(1))*100000) / 100000;

  cout<<amp0_ll<<" - "<<amp0_ul<<endl;
  cout<<amp1_ll<<" - "<<amp1_ul<<endl;


  TH1D *h1AmpX = amp_vs_x->ProjectionX();
  TH1D *h1AmpY = amp_vs_y->ProjectionX();
  TH1D *h1MCP = MCP_time_vs_amp->ProjectionX();

  fout->cd();
  amp_sipm0_vs_sipm1->Write();
  amp_sipm0_vs_sipm1Corr->Write();
  MCP_time_vs_amp->Write();
  hBeam->Write();
  hprof2D->Write();
  amp_vs_x->Write();
  amp_vs_y->Write();
  h1AmpX->Write();
  h1AmpY->Write();
  h1MCP->Write();
  MPV_amp0->Write();
  MPV_amp1->Write();
  MCP_amp->Write();
  MCP_efficiency->Write();
  MCP_hitsProfile->Write();

  //------ loop over thresholds --------

  for (int ith=40; ith < 131; ith+= 5) {

      Initialize();
      int count = 0;

      //---------- create histogram per channel per threshold ------------

      TH1F *htRes_diff_AmpXposCorrected = new TH1F(Form("htRes_diff_AmpXposCorrected_th%d",ith), "(t0 - t1)/2 after amplitude and position correction; time (ns)", 400, -1, 1);    
      //TH1F *htRes_avg_AmpXposCorrected = new TH1F(Form("htRes_avg_AmpXposCorrected_th%d",ith), "tRes_avg; time (ns)", 200, -2, 2);

      //TH2F *tdiff_vs_x = new TH2F(Form("tdiff_vs_x_th%d",ith), ";x (mm); t0 - t1 (ns)", 50, 22, 29, 100, -1, 1);
      //TH2F *tdiff_vs_x_corrected = new TH2F(Form("tdiff_vs_x_corrected_th%d",ith), ";x (mm); t0 - t1 (ns)", 50, 22, 29, 100, -1, 1);
      TH2F *tdiff_vs_x_ampCorrected = new TH2F(Form("tdiff_vs_x_ampCorrected_th%d",ith), ";x (mm); t0 - t1 (ns)", 40, 20, 30, 100, -1, 1);
      TH2F *tdiff_vs_x_ampXposCorrected = new TH2F(Form("tdiff_vs_x_ampXposCorrected_th%d",ith), ";x (mm); t0 - t1 (ns)", 40, 20, 30, 50, -1, 1);
      //TH2F *tavg_vs_x_ampCorrected = new TH2F(Form("tavg_vs_x_ampCorrected_th%d",ith), ";x (mm); 0.5(t0 + t1) -tMCP (ns)", 50, 22, 29, 100, -0.2, 0.2);
      TH2D *tdiff_vs_ampSIPM0 = new TH2D(Form("tdiff_vs_ampSIPM0_th%d",ith), "; amp (mV); t0 - t1 (ns)", 50, 1/amp0_ll, 1/amp0_ul, 50, -1, 1);
      TH2D *tdiff_vs_ampSIPM1 = new TH2D(Form("tdiff_vs_ampSIPM1_th%d",ith), "; amp (mV); t0 - t1 (ns)", 50, 1/amp1_ll, 1/amp1_ul, 50, -1, 1);
      TH2D *tdiff_vs_MPVampSIPM0 = new TH2D(Form("tdiff_vs_MPVampSIPM0_th%d",ith), "; amp/MPV (mV); t0 - t1 (ns)", 50, 0.8, 3, 50, -1, 1);
      TH2D *tdiff_vs_MPVampSIPM1 = new TH2D(Form("tdiff_vs_MPVampSIPM1_th%d",ith), "; amp/MPV (mV); t0 - t1 (ns)", 50, 0.8, 3, 50, -1, 1);
      TH2D *tdiff_vs_avgMPVamp = new TH2D(Form("tdiff_vs_avgMPVamp_th%d",ith), "; amp/MPV (mV); t0 - t1 (ns)", 50, 0.8, 3, 50, -1, 1);
      //TH2D *t0_minus_tMCP_vs_amp = new TH2D(Form("t0_minus_tMCP_vs_amp_th%d",ith), "; amp (mV); t0 - tMCP (ns)", 100, 1075, 3000, 100, -55.8, -54.6);
      //TH2D *t1_minus_tMCP_vs_amp = new TH2D(Form("t1_minus_tMCP_vs_amp_th%d",ith), "; amp (mV); t1 - tMCP (ns)", 100, 1108, 3000, 100, -55.6, -54.6);
      TH2D *t0_minus_tMCP_vs_amp = new TH2D(Form("t0_minus_tMCP_vs_amp_th%d",ith), "; 1/amp (mV); t0 - tMCP (ns)", 50, amp0_ul, amp0_ll, 260, -58, -47);
      TH2D *t1_minus_tMCP_vs_amp = new TH2D(Form("t1_minus_tMCP_vs_amp_th%d",ith), "; 1/amp (mV); t1 - tMCP (ns)", 50, amp1_ul, amp1_ll, 260, -58, -47);
      TH2D *t0_minus_tMCP_vs_amp_corrected = new TH2D(Form("t0_minus_tMCP_vs_amp_corrected_th%d",ith), "; 1/amp (mV); t0 - tMCP (ns)", 50, amp0_ul, amp0_ll, 100, -0.5, 0.5);
      TH2D *t1_minus_tMCP_vs_amp_corrected = new TH2D(Form("t1_minus_tMCP_vs_amp_corrected_th%d",ith), "; 1/amp (mV); t1 - tMCP (ns)", 50, amp1_ul, amp1_ll, 100, -0.5, 0.5);

      TF1 *fit_xCorrection = new TF1("fit_xCorrection","pol1",26,29);
      //TF1 *fit_xCorrection1 = new TF1("fit_xCorrection1","pol1",26,28.5);
      TF1 *fit_ampwalk0 = new TF1("fit_ampwalk0","pol4",amp0_ul, amp0_ll);
      TF1 *fit_ampwalk1 = new TF1("fit_ampwalk1","pol4",amp1_ul, amp1_ll);

      // -------------- 2nd loop over events ---------------------
      cout<<"2nd loop over events"<<endl;

      for (Long64_t jentry=0; jentry<nentries;jentry++) {
      mytree->GetEntry(jentry);
      if (jentry % 10000 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";

      bool selections = cuts(event.xtrk, event.ytrk, event.timeMCP);
      if (selections == false) continue;
      if(event.ampLo[0] < 0.8*l2->GetParameter(1)) continue;
      if(event.ampLo[0] > 3*l2->GetParameter(1)) continue;
      if(event.ampLo[1] < 0.8*l3->GetParameter(1)) continue;
      if(event.ampLo[1] > 3*l3->GetParameter(1)) continue;
      if (event.ampMCP < 100) continue; 

      t = tot(ith);

      //tRes_avg = abs(0.5*(t[0]+t[1]) - event.timeMCP);  //time resolution with MCP, tMCP is bigger than tbar!
      tdiff[jentry] = (t[count][0] - t[count][1]);  // time difference method (without MCP)
      //tdiff_vs_x->Fill(event.xtrk, tdiff[jentry]);

      t0_minus_tMCP_vs_amp->Fill(1/event.ampLo[0], t[count][0] - event.timeMCP);
      t1_minus_tMCP_vs_amp->Fill(1/event.ampLo[1], t[count][1] - event.timeMCP);

      count += 1;

    }

  if (ith < 55){
    t0_minus_tMCP_vs_amp->GetYaxis()->SetRangeUser(t0_minus_tMCP_vs_amp->GetMean(2) - 2*t0_minus_tMCP_vs_amp->GetRMS(2), t0_minus_tMCP_vs_amp->GetMean(2) + 2*t0_minus_tMCP_vs_amp->GetRMS(2));
    t1_minus_tMCP_vs_amp->GetYaxis()->SetRangeUser(t1_minus_tMCP_vs_amp->GetMean(2) - 2*t1_minus_tMCP_vs_amp->GetRMS(2), t1_minus_tMCP_vs_amp->GetMean(2) + 2*t1_minus_tMCP_vs_amp->GetRMS(2));

  }
  else{
    t0_minus_tMCP_vs_amp->GetYaxis()->SetRangeUser(t0_minus_tMCP_vs_amp->GetMean(2) - 2*t0_minus_tMCP_vs_amp->GetRMS(2), t0_minus_tMCP_vs_amp->GetMean(2) + 2*t0_minus_tMCP_vs_amp->GetRMS(2));
    t1_minus_tMCP_vs_amp->GetYaxis()->SetRangeUser(t1_minus_tMCP_vs_amp->GetMean(2) - 2*t1_minus_tMCP_vs_amp->GetRMS(2), t1_minus_tMCP_vs_amp->GetMean(2) + 2*t1_minus_tMCP_vs_amp->GetRMS(2));
  }
    
    t0_minus_tMCP_vs_amp->Fit(fit_ampwalk0,"QR");
    t1_minus_tMCP_vs_amp->Fit(fit_ampwalk1,"QR");
    //tdiff_vs_x->Fit(fit_xCorrection,"QR");

    count = 0;

    // -------------- 3rd loop over events ---------------------
    cout<< "3rd loop over events"<<endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      mytree->GetEntry(jentry);
      if (jentry % 10000 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";

      //apply cuts -- select only the good events
      bool selections = cuts(event.xtrk, event.ytrk, event.timeMCP);
      if (selections == false) continue;
      if (event.ampLo[0] < 0.8*l2->GetParameter(1)) continue;
      if (event.ampLo[0] > 3*l2->GetParameter(1)) continue;
      if (event.ampLo[1] < 0.8*l3->GetParameter(1)) continue;
      if (event.ampLo[1] > 3*l3->GetParameter(1)) continue;
      if (event.ampMCP < 100) continue; 


     /* y = fit_xCorrection->Eval(event.xtrk);
      tdiffPosCorr = tdiff[jentry] - y;   //x-position correction.  */

      y_amp = fit_ampwalk0->Eval(1/event.ampLo[0]); 
      t0_AmpCorrected[jentry] = (t[count][0] - event.timeMCP) - y_amp;

      y_amp = fit_ampwalk1->Eval(1/event.ampLo[1]); 
      t1_AmpCorrected[jentry] = (t[count][1] - event.timeMCP) - y_amp;

      t0_minus_tMCP_vs_amp_corrected->Fill(1/event.ampLo[0], t0_AmpCorrected[jentry]);
      t1_minus_tMCP_vs_amp_corrected->Fill(1/event.ampLo[1], t1_AmpCorrected[jentry]);

      /*tdiff_vs_ampSIPM0->Fill(event.ampLo[0], t0_AmpCorrected[jentry] - t1_AmpCorrected[jentry]);
      tdiff_vs_ampSIPM1->Fill(event.ampLo[1], t0_AmpCorrected[jentry] - t1_AmpCorrected[jentry]);
      tdiff_vs_MPVampSIPM0->Fill(event.ampLo[0]/l2->GetParameter(1), t0_AmpCorrected[jentry] - t1_AmpCorrected[jentry]);
      tdiff_vs_MPVampSIPM1->Fill(event.ampLo[1]/l3->GetParameter(1), t0_AmpCorrected[jentry] - t1_AmpCorrected[jentry]);
      */
      //tavg_vs_x_ampCorrected->Fill(event.xtrk, 0.5*(t0_AmpCorrected[jentry] + t1_AmpCorrected[jentry]));
      tdiff_vs_x_ampCorrected->Fill(event.xtrk,(t0_AmpCorrected[jentry] - t1_AmpCorrected[jentry]));
      
      //tdiff_vs_x_corrected->Fill(event.xtrk, tdiffPosCorr);
      //htRes_diff_xposCorrected->Fill(tdiffPosCorr/2);

      count+=1;
    }

    tdiff_vs_x_ampCorrected->Fit(fit_xCorrection,"QR");
    //tavg_vs_x_ampCorrected->Fit(fit_xCorrection1,"QR");

    // -------------- 4th loop over events ---------------------
    cout<<"4th loop over events"<<endl;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      mytree->GetEntry(jentry);
      if (jentry % 10000 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";

      //apply cuts -- select only the good events
      bool selections = cuts(event.xtrk, event.ytrk, event.timeMCP);
      if (selections == false) continue;
      if (event.ampLo[0] < 0.8*l2->GetParameter(1)) continue;
      if (event.ampLo[0] > 3*l2->GetParameter(1)) continue;
      if (event.ampLo[1] < 0.8*l3->GetParameter(1)) continue;
      if (event.ampLo[1] > 3*l3->GetParameter(1)) continue;
      if (event.ampMCP < 100) continue; 

      y = fit_xCorrection->Eval(event.xtrk);
      //y1 = fit_xCorrection1->Eval(event.xtrk);
      tdiffAmpPosCorr = (t0_AmpCorrected[jentry] - t1_AmpCorrected[jentry]) - y;   //x-position correction
      //tavgAmpPosCorr = (0.5*(t0_AmpCorrected[jentry] + t1_AmpCorrected[jentry])) - y1;

      tdiff_vs_x_ampXposCorrected->Fill(event.xtrk,tdiffAmpPosCorr);
      htRes_diff_AmpXposCorrected->Fill(tdiffAmpPosCorr/2);

      tdiff_vs_ampSIPM0->Fill(event.ampLo[0], tdiffAmpPosCorr);
      tdiff_vs_ampSIPM1->Fill(event.ampLo[1], tdiffAmpPosCorr);
      tdiff_vs_MPVampSIPM0->Fill(event.ampLo[0]/l2->GetParameter(1), tdiffAmpPosCorr);
      tdiff_vs_MPVampSIPM1->Fill(event.ampLo[1]/l3->GetParameter(1), tdiffAmpPosCorr);
      tdiff_vs_avgMPVamp->Fill(0.5*(event.ampLo[0]/l2->GetParameter(1) + event.ampLo[1]/l3->GetParameter(1)), tdiffAmpPosCorr);
      //htRes_avg_AmpXposCorrected->Fill(tavgAmpPosCorr);

    }

    if(ith == 65){

    TGraphErrors *gr_tRes_vs_amp = new TGraphErrors();
    //TGraphErrors *gr_tdiff_vs_amp = new TGraphErrors();
    vector<double> x_arr_amp = xvalues(0, 0.8, 3);
    //x_arr_amp1 = xvalues(0,1/amp1_ll, 1/amp1_ul);
    vector<double> x_err = xerrors(0.8, 3);
    //vector<double> x_err_1 = xerrors(1/amp1_ll, 1/amp1_ul);

    for(int i = 0; i<xbins_a.size()-1; i++){

        TH1D *t_sipm0 = tdiff_vs_avgMPVamp->ProjectionY("",xbins_a[i],xbins_a[i+1]);
        TF1* f_sipm0 = new TF1("f_sipm0","gaus",t_sipm0->GetMean() - 2.5*t_sipm0->GetRMS(),t_sipm0->GetMean() + 2.5*t_sipm0->GetRMS());
        t_sipm0->Fit(f_sipm0, "QR");
        c->cd();
        t_sipm0->Draw();
        c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46043/c_gaus_test_th%d_%s_%d.png",ith,argv[1],i));
        c->Clear();
        gr_tRes_vs_amp ->SetPoint(gr_tRes_vs_amp->GetN(), x_arr_amp[i], f_sipm0->GetParameter(2)/2);
        gr_tRes_vs_amp ->SetPointError(gr_tRes_vs_amp->GetN()-1, x_err[i], f_sipm0->GetParError(2)/2); 

        //TH1D *t_sipm1 = tdiff_vs_ampSIPM1->ProjectionY("",xbins_a[i],xbins_a[i+1]);
        //TF1* f_sipm1 = new TF1("f_sipm1","gaus",t_sipm1->GetMean() - 2.5*t_sipm1->GetRMS(),t_sipm1->GetMean() + 2.5*t_sipm1->GetRMS());
        //t_sipm1->Fit(f_sipm1, "QR");  
        
        //gr_tdiff_vs_amp1 ->SetPoint(gr_tdiff_vs_amp1->GetN(), x_arr_amp1[i], f_sipm1->GetParameter(2)/2);
        //gr_tdiff_vs_amp1 ->SetPointError(gr_tdiff_vs_amp1->GetN()-1, x_err_1[i], f_sipm1->GetParError(2)/2);  

        delete t_sipm0, f_sipm0;
        //delete t_sipm1, f_sipm1; 

      } 

    //gr_tdiff_vs_amp0->Fit(fit_xCorrection, "QRE");
    c->cd();
    //tdiff_vs_ampSIPM0->Draw("colz");
    gr_tRes_vs_amp->SetMarkerSize(1);
    gr_tRes_vs_amp->SetMarkerStyle(20);
    gr_tRes_vs_amp->SetMarkerColor(kBlack);
    gr_tRes_vs_amp->SetTitle("time resolution vs slices in amp; amp/MPV ; (t_{L} - t_{R})/2 (ns)");
    gr_tRes_vs_amp->Draw("AP same");
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46123/c_tdiff_vs_amp0_sigma_th%d_%s.png",ith,argv[1]));
    c->Clear();

    //c->cd();
    //tdiff_vs_ampSIPM0->Draw("colz");
    //gr_tdiff_vs_amp1->SetMarkerSize(1);
    //gr_tdiff_vs_amp1->SetMarkerStyle(20);
    //gr_tdiff_vs_amp1->SetMarkerColor(kBlack);
    //gr_tdiff_vs_amp1->SetTitle("time resolution vs slices in amp (SiPM1); amp ; (t0 - t1)/2");
    //gr_tdiff_vs_amp1->Draw("AP");
    //c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46123/c_tdiff_vs_amp1_sigma_th%d_%s.png",ith,argv[1]));
    //c->Clear();

    fout->cd();
    gr_tRes_vs_amp->Write("gr_tRes_vs_amp_bestthr");
    //gr_tdiff_vs_amp1->Write("gr_tdiff_vs_amp1_sigma_th40");

    delete gr_tRes_vs_amp;
    //delete gr_tdiff_vs_amp1;

  }


  /* if(ith == 25 || ith == 30){
    
    TGraphErrors *gr_tdiff_vs_x_mean = new TGraphErrors();
    TGraphErrors *gr_tdiff_vs_x_sigma = new TGraphErrors();
    x_arr_pos = xvalues(1);

    for(int i = 0; i<xbins_p.size()-1; i++){

        TH1D *t_ = tdiff_vs_x_ampXposCorrected->ProjectionY("",xbins_p[i],xbins_p[i+1]);
        TF1* f_ = new TF1("f_","gaus",t_->GetMean() - 2.5*t_->GetRMS(),t_->GetMean() + 2.5*t_->GetRMS());
        t_->Fit(f_, "QR");
        c->cd();
        t_->Draw();
        c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46123/c_gaus_pos_test_th%d_%s_%d.png",ith,argv[1],i));
        c->Clear();
        gr_tdiff_vs_x_mean ->SetPoint(gr_tdiff_vs_x_mean->GetN(), x_arr_pos[i], f_->GetParameter(1));
        gr_tdiff_vs_x_mean ->SetPointError(gr_tdiff_vs_x_mean->GetN()-1, 0, f_->GetParError(1)); 
        gr_tdiff_vs_x_sigma ->SetPoint(gr_tdiff_vs_x_sigma->GetN(), x_arr_pos[i], f_->GetParameter(2)/2);
        gr_tdiff_vs_x_sigma ->SetPointError(gr_tdiff_vs_x_sigma->GetN()-1, 0, f_->GetParError(2)/2);   

        delete t_, f_;  

        }    

    gr_tdiff_vs_x_mean->Fit(fit_xCorrection, "QRE");
    c->cd();
    //tdiff_vs_x->Draw("colz");
    gr_tdiff_vs_x_mean->SetMarkerSize(1);
    gr_tdiff_vs_x_mean->SetMarkerStyle(20);
    gr_tdiff_vs_x_mean->SetMarkerColor(kBlack);
    gr_tdiff_vs_x_mean->SetTitle("mean of t0 - t1 vs slices in x position; x(mm) ; t0 - t1");
    gr_tdiff_vs_x_mean->Draw("AP");
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46123/c_tdiff_vs_x_corrected_mean_th%d_%s.png",ith,argv[1]));
    c->Clear(); 

    c->cd();
    //tdiff_vs_x->Draw("colz");
    gr_tdiff_vs_x_sigma->SetMarkerSize(1);
    gr_tdiff_vs_x_sigma->SetMarkerStyle(20);
    gr_tdiff_vs_x_sigma->SetMarkerColor(kBlack);
    gr_tdiff_vs_x_sigma->SetTitle("time resolution vs slices in x position; x(mm) ; (t0 - t1)/2");
    gr_tdiff_vs_x_sigma->Draw("AP");
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46123/c_tdiff_vs_x_corrected_sigma_th%d_%s.png",ith,argv[1]));
    c->Clear(); 

    delete gr_tdiff_vs_x_mean, gr_tdiff_vs_x_sigma;

} */

    TF1 *fit_tRes_diff3 = new TF1("fit_tRes_diff3","gaus",htRes_diff_AmpXposCorrected->GetMean()-3*htRes_diff_AmpXposCorrected->GetRMS(),htRes_diff_AmpXposCorrected->GetMean()+3*htRes_diff_AmpXposCorrected->GetRMS());
    htRes_diff_AmpXposCorrected->Fit(fit_tRes_diff3,"QR");
    TF1 *fit_tRes_diff4 = new TF1("fit_tRes_diff4","gaus",fit_tRes_diff3->GetParameter(1)-4*fit_tRes_diff3->GetParameter(2),fit_tRes_diff3->GetParameter(1)+4*fit_tRes_diff3->GetParameter(2));
    htRes_diff_AmpXposCorrected->Fit(fit_tRes_diff4,"QR");
    TF1 *fit_tRes_diff5 = new TF1("fit_tRes_diff5","gaus",fit_tRes_diff4->GetParameter(1)-3*fit_tRes_diff4->GetParameter(2),fit_tRes_diff4->GetParameter(1)+3*fit_tRes_diff4->GetParameter(2));
    htRes_diff_AmpXposCorrected->Fit(fit_tRes_diff5,"QR");
    
    //tRes_vs_threshold1->SetPoint(tRes_vs_threshold1->GetN(), ith, fit_tRes_diff2->GetParameter(2));
    //tRes_vs_threshold1->SetPointError(tRes_vs_threshold1->GetN()-1, 0, fit_tRes_diff2->GetParError(2));
    tRes_vs_threshold->SetPoint(tRes_vs_threshold->GetN(), ith, fit_tRes_diff5->GetParameter(2));
    tRes_vs_threshold->SetPointError(tRes_vs_threshold->GetN()-1, 0, fit_tRes_diff5->GetParError(2));


   /* c->cd();
    htRes_diff_xposCorrected->Draw("colz");
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46075/c_tRes_diff_xposCorrected_th%d_%s.png",ith,argv[1]));
    c->Clear();
  */



    fout->cd();
    //tdiff_vs_x->Write();
    //tavg_vs_x_ampCorrected->Write();
    //htRes_avg_AmpXposCorrected->Write();
    tdiff_vs_x_ampCorrected->Write();
    tdiff_vs_x_ampXposCorrected->Write();
    t0_minus_tMCP_vs_amp->Write();
    t1_minus_tMCP_vs_amp->Write();
    t0_minus_tMCP_vs_amp_corrected->Write();
    t1_minus_tMCP_vs_amp_corrected->Write();
    tdiff_vs_ampSIPM0->Write();
    tdiff_vs_ampSIPM1->Write();
    tdiff_vs_MPVampSIPM0->Write();
    tdiff_vs_MPVampSIPM1->Write();
    //tdiff_vs_x_corrected->Write();
    //htRes_diff_xposCorrected->Write();
    htRes_diff_AmpXposCorrected->Write();
      
    
    //delete tdiff_vs_x;
    delete tdiff_vs_ampSIPM0;
    delete tdiff_vs_ampSIPM1;
    delete tdiff_vs_MPVampSIPM0;
    delete tdiff_vs_MPVampSIPM1;
    delete tdiff_vs_avgMPVamp;
    delete t0_minus_tMCP_vs_amp;
    delete t1_minus_tMCP_vs_amp;
    delete t1_minus_tMCP_vs_amp_corrected;
    delete t0_minus_tMCP_vs_amp_corrected;
 
    //delete htRes_avg_AmpXposCorrected;
    //delete htRes_diff;
    //delete htRes_diff_xposCorrected;
    //delete fit_tRes_diff;
    delete fit_xCorrection; 
    //delete tdiff_vs_x_corrected;
    delete fit_ampwalk0;
    delete fit_ampwalk1;
    delete tdiff_vs_x_ampCorrected;
    delete htRes_diff_AmpXposCorrected;
    //delete gr_tdiff_vs_x_mean;
    delete fit_tRes_diff3, fit_tRes_diff4, fit_tRes_diff5;
    delete tdiff_vs_x_ampXposCorrected; 
  
}

    tRes_vs_threshold ->SetTitle(" ; threshold (mV/ns); #sigma_{bar} (ns) ");



 
    tRes_vs_threshold->SetMarkerSize(1);
    tRes_vs_threshold->SetMarkerStyle(20);
    tRes_vs_threshold->SetMarkerColor(kBlack);

    tRes_vs_threshold->Draw("ALP");
    tRes_vs_threshold->Write("tRes_vs_threshold");

    delete h1AmpX;
    delete h1AmpY;
    delete h1MCP;
    delete MCP_time_vs_amp;
    delete amp_sipm0_vs_sipm1Corr;
    delete amp_sipm0_vs_sipm1;
    delete hBeam;
    delete hprof2D;
    delete amp_vs_x;
    delete amp_vs_y;
    delete c;
    delete MPV_amp0;
    delete MPV_amp1;
    delete tRes_vs_threshold;
    delete MCP_amp, MCP_efficiency, MCP_hitsProfile;
    delete l0, l1, l2, l3;
    //delete tRes_vs_threshold1;
    //delete gr_tdiff_vs_amp0, gr_tdiff_vs_amp1, t_sipm0, f_sipm0;

    fout->Close();

    return 0;
}


