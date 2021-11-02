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

// Track total memory in use 

#include "sys/types.h"
#include "sys/sysinfo.h"

struct sysinfo memInfo;

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


//TFile *fout = new TFile("./output_timeRes/Run_46076_DItest.root","RECREATE");

double getmin(TGraph *gr) {
  /*double mean = 0;
  for(int i=0; i!=250; ++i)
    mean += gr->GetPointY(i);
  mean /= 250;*/

  double min = 11111111111;
  int p = 0;
  for(int i=250; i!=350; ++i) {
    if(gr->GetPointY(i)<min){
      min = gr->GetPointY(i);
      p = i;
    }
  }
  //return  -( min - mean );
  return p;
}

vector<vector<double>> tot(int th){

//to find the time over threshold for various threshold values via linear interpolation of amplitude values
// just below and above the threshold

      Double_t vbelow0 = -1111, nbelow0 = 200, vabove0 = 1111, nabove0 = 200; 
      Double_t vbelow1 = -1111, nbelow1 = 200, vabove1 = 1111, nabove1 = 200; 
      Double_t m0, cy0, m1, cy1, y0, y1, xb0, xa0, yb0, yb1, ya0, xb1, xa1, ya1, t1, t0; 
      
      //TGraph *grWF0 = new TGraph(1024,event.times,event.channel[0]);
      //TGraph *grWF1 = new TGraph(1024,event.times,event.channel[1]);

      TGraph *gr_t_plus0 = new TGraph();
      TGraph *gr_t_minus0 = new TGraph();
      TGraph *gr_t_plus1 = new TGraph();
      TGraph *gr_t_minus1 = new TGraph();
      TGraph *gr_DI0 = new TGraph();
      TGraph *gr_DI1 = new TGraph();
      //TCanvas *c1 = new TCanvas();
      //TCanvas *c2 = new TCanvas();

      for (int n = 0; n<1024; n++){
       gr_t_plus0->SetPoint(n,event.times[n] + 0.2,event.channel[0][n]);
       gr_t_minus0->SetPoint(n,event.times[n] - 0.2,event.channel[0][n]);
       gr_t_plus1->SetPoint(n,event.times[n] + 0.2,event.channel[1][n]);
       gr_t_minus1->SetPoint(n,event.times[n] - 0.2,event.channel[1][n]);

      }

      for (int j=0; j<1024; j++){
        double yp_0 = gr_t_plus0->Eval(event.times[j]);
        double ym_0 = gr_t_minus0->Eval(event.times[j]);
        double yp_1 = gr_t_plus1->Eval(event.times[j]);
        double ym_1 = gr_t_minus1->Eval(event.times[j]);

        gr_DI0->SetPoint(j,event.times[j],(yp_0 - ym_0)/0.4);
        gr_DI1->SetPoint(j,event.times[j],(yp_1 - ym_1)/0.4);

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
      gr_DI1->SetMarkerSize(0.2);

      c1->Divide(2,2);
      c1->cd(1);
      grWF0->Draw("AP");
      c1->cd(2);
      grWF1->Draw("AP");

      //c2->Divide(2,1);
      c1->cd(3);
      gr_DI0->Draw("AP");
      //gr_DI0->GetXaxis()->SetRangeUser(20,70);
      //gr_DI0->Draw("AP");
      //c1->Update();
      c1->cd(4);
      gr_DI1->Draw("AP");
      //gr_DI1->GetXaxis()->SetRangeUser(20,70);
      //gr_DI1->Draw("AP");
      //c1->Update();
      
     // ifstream file(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46076/c_WFDI_test_th%d.png",th));
     // if(!file){
        //c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46076/c_WF_test_th%d.png",th));
       // c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46076/c_WFDI_test_th%d.png",th));
     // }
      c1->Clear();
      c2->Clear();
      //fout->cd();
      //grWF0->Write("WF0");
      //grWF1->Write("WF1");
      //gr_DI0->Write("DI0");
      //gr_DI1->Write("DI1");
      */

      double point0 = getmin(gr_DI0);
      double point1 = getmin(gr_DI1);

      for (int i = 225; i<point0; i++){
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

      for (int i = 225; i<point1; i++){
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

        return totVals;

        delete gr_t_plus0, gr_t_plus1, gr_t_minus0, gr_t_minus1;

}


vector<int> xbins{18,31,41,51};

/*  vector<double> mean_;
  vector<double> meanerr_;
};

     void slicePositionCorr2D(TH2F *t)
     {

      TCanvas *c1 = new TCanvas();
  
      for(int i = 0; i< xbins.size()-1; i++){

        TH1D *t_ = t->ProjectionY("",xbins[i],xbins[i+1]);
        TF1* f_ = new TF1("f_0","gaus",t_->GetMean() - 3*t_->GetRMS(),t_->GetMean() + 3*t_->GetRMS());
        t_->Fit(f_, "QR");
        g.mean_.push_back(f_->GetParameter(1));
        g.meanerr_.push_back(f_->GetParError(1));

       /* c1->cd();
        t_->Draw();
        c1->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46075/c_test_%d.png",i));
        c1->Clear();


        delete t_, f_;
      }

      delete c1;

      //return g;
     
     }*/


     vector<double> xvalues()
     {
      vector<double> vals;
      double x;
      for(int i = 0; i<xbins.size()-1; i++){
        x = (25+xbins[i]*0.043) + (25+xbins[i+1]*0.043);
        vals.push_back(x/2);
      }
      return vals;

     }


int main(int argc, char* argv[]){ 

  TFile *fout = new TFile(Form("./output_timeRes/Run_%s_DItest1.root",argv[1]),"RECREATE");  
  TFile *infile = TFile::Open(Form("/afs/cern.ch/user/a/amkrishn/private/summer_2021_BTL/BTL_analysis_UVa_Sasha/skims/Run_%s_bar0E00.root",argv[1]), "READ");
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

  TCanvas *c = new TCanvas();
  TH2D *amp_sipm0_vs_sipm1 = new TH2D("amp_sipm0_vs_sipm1", "; amp_SiPM0 (mV); amp_SiPM1 (mV)", 1000, -100, 2000, 1000, -100, 2000);
  TH2D *amp_sipm0_vs_sipm1Corr = new TH2D("amp_sipm0_vs_sipm1Corr", "; amp_SiPM0 (mV); amp_SiPM1 (mV)", 1000, -100, 2000, 1000, -100, 2000);
  TH2D *MCP_time_vs_amp = new TH2D("MCP_time_vs_amp", "; time_MCP; amp_MCP", 1000, 0, 200, 1000, -100, 1000);
  TH2D *amp_vs_x = new TH2D("amp_vs_x", "; x (mm); amp_SiPM0 (mV)", 120, 24, 30, 600, 500, 2000);
  TH2D *amp_vs_y = new TH2D("amp_vs_y", "; y (mm); amp_SiPM0 (mV)", 120, 15, 20, 600, 500, 2000);
  TProfile2D* hprof2D = new TProfile2D("hprof2D", "Profile of Amp vs x and y tracks", 120, 20, 30, 120, 0, 40, -10, 2000);
  TH2D *hBeam = new TH2D("hBeam", "Beam; x (mm); y (mm); events", 120, 10, 40, 120, 0, 40);
  TH1D *MPV_amp0 = new TH1D("MPV_amp0", "SiPM0; amplitude (mV); counts", 400, 200, 2600);
  TH1D *MPV_amp1 = new TH1D("MPV_amp1", "SiPM1; amplitude (mV); counts", 400, 200, 2600);
  TGraphErrors *tRes_vs_threshold1 = new TGraphErrors();
  TGraphErrors *tRes_vs_threshold2 = new TGraphErrors();
   
  Double_t tdiff[6000], tRes_avg, delta_t, y, tdiffPosCorr, y_amp, t0_AmpCorrected[6000], t1_AmpCorrected[6000], tdiffAmpPosCorr, tavgAmpPosCorr, y1; 
  vector<vector<double>> t;
  std::vector<double> y_mean,x_arr_posCorr;

  Long64_t nentries = mytree->GetEntries();
  
  for (int ith=-10; ith > -151; ith+= -5) {
    std::cout << "ith: " << ith << std::endl; 

    Initialize();
    int count = 0;

    sysinfo (&memInfo);

    long long totalPhysMem = memInfo.totalram;
    //Multiply in next statement to avoid int overflow on right hand side...
    totalPhysMem *= memInfo.mem_unit;  

    long long physMemUsed = memInfo.totalram - memInfo.freeram;
    //Multiply in next statement to avoid int overflow on right hand side...
    physMemUsed *= memInfo.mem_unit;

    std::cout << "totalPhysMem: " << totalPhysMem << std::endl;
    std::cout << "physMemUsed: " << physMemUsed << std::endl;

    double pcttotalPhysMemUsed = (float)physMemUsed / (float)totalPhysMem; 

    std::cout << "Percentage of total Available RAM in use: " << pcttotalPhysMemUsed << std::endl; 

    //TH1F *htRes_diff = new TH1F(Form("htRes_diff_th%d",ith), "tRes_diff; time (ns)", 100, -1.5, 1.5);
    TH1F *htRes_diff_xposCorrected = new TH1F(Form("htRes_diff_xposCorrected_th%d",ith), "(t0 - t1)/2 after position correction; time (ns)", 200, -0.2, 0.2);
    TH1F *htRes_diff_AmpXposCorrected = new TH1F(Form("htRes_diff_AmpXposCorrected_th%d",ith), "(t0 - t1)/2 after amplitude and position correction; time (ns)", 200, -0.2, 0.2);    
    TH1F *htRes_avg_AmpXposCorrected = new TH1F(Form("htRes_avg_AmpXposCorrected_th%d",ith), "tRes_avg; time (ns)", 200, -0.4, 0.4);

    TH2F *tdiff_vs_x = new TH2F(Form("tdiff_vs_x_th%d",ith), ";x (mm); t0 - t1 (ns)", 50, 22, 28.23, 50, -0.47, -0.065);
    TH2F *tdiff_vs_x_corrected = new TH2F(Form("tdiff_vs_x_corrected_th%d",ith), ";x (mm); t0 - t1 (ns)", 50, 22, 28.23, 50, -0.2, 0.2);
    TH2F *tdiff_vs_x_ampCorrected = new TH2F(Form("tdiff_vs_x_ampCorrected_th%d",ith), ";x (mm); t0 - t1 (ns)", 50, 22, 28.23, 50, -0.17, 0.19);
    TH2F *tdiff_vs_x_ampXposCorrected = new TH2F(Form("tdiff_vs_x_ampXposCorrected_th%d",ith), ";x (mm); t0 - t1 (ns)", 50, 22, 28.23, 50, -0.2, 0.2);
    TH2F *tavg_vs_x_ampCorrected = new TH2F(Form("tavg_vs_x_ampCorrected_th%d",ith), ";x (mm); 0.5(t0 + t1) -tMCP (ns)", 50, 22, 28.23, 100, -0.2, 0.2);
    TH2D *tdiff_vs_ampSIPM0 = new TH2D(Form("tdiff_vs_ampSIPM0_th%d",ith), "; amp (mV); t0 - t1 (ns)", 400, -10, 2000, 100, -1.5, 1.5);
    TH2D *tdiff_vs_ampSIPM1 = new TH2D(Form("tdiff_vs_ampSIPM1_th%d",ith), "; amp (mV); t0 - t1 (ns)", 400, -10, 2000, 100, -1.5, 1.5);
    TH2D *t0_minus_tMCP_vs_amp = new TH2D(Form("t0_minus_tMCP_vs_amp_th%d",ith), "; amp (mV); t0 - tMCP (ns)", 100, 1200, 3000, 100, -55.8, -54.6);
    TH2D *t1_minus_tMCP_vs_amp = new TH2D(Form("t1_minus_tMCP_vs_amp_th%d",ith), "; amp (mV); t1 - tMCP (ns)", 100, 1200, 3000, 100, -55.6, -54.6);
    TH2D *t0_minus_tMCP_vs_amp_corrected = new TH2D(Form("t0_minus_tMCP_vs_amp_corrected_th%d",ith), "; amp (mV); t0 - tMCP (ns)", 200, 1200, 3000, 100, -0.5, 0.5);
    TH2D *t1_minus_tMCP_vs_amp_corrected = new TH2D(Form("t1_minus_tMCP_vs_amp_corrected_th%d",ith), "; amp (mV); t1 - tMCP (ns)", 200, 1200, 3000, 100, -0.5, 0.5);

    TF1 *fit_xCorrection = new TF1("fit_xCorrection","pol1",26,28);
    TF1 *fit_xCorrection1 = new TF1("fit_xCorrection1","pol1",26,28);
    TF1 *fit_ampwalk0 = new TF1("fit_ampwalk0","pol3",1230,3000);
    TF1 *fit_ampwalk1 = new TF1("fit_ampwalk1","pol3",1290,3000);

    //fit_ampwalk->SetParameter(0,7.7);

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

      mytree->GetEntry(jentry);
      if (jentry % 1000 == 0) cout << "Processing Event " << jentry << " of " << nentries << "\n";

      //histograms to determine the cuts to be applied so that we select only the good events
        amp_sipm0_vs_sipm1->Fill(event.ampLo[0], event.ampLo[1]);
        MCP_time_vs_amp->Fill(event.timeMCP,event.ampMCP);
        hprof2D->Fill(event.xtrk, event.ytrk, event.ampLo[1]);
        hBeam->Fill(event.xtrk, event.ytrk);

      //cuts applied
        if(event.xtrk < 22 ||  event.xtrk > 27.6) continue;   //track cuts
        if(event.ytrk < 16. || event.ytrk > 18.85) continue;
        if(event.timeMCP < 100 ||  event.timeMCP > 118) continue;   //MCP time cuts
        if(event.ampLo[0]< 1230. || event.ampLo[0] > 4605) continue;
        if(event.ampLo[1]< 1290. || event.ampLo[1] > 4740) continue;

        amp_sipm0_vs_sipm1Corr->Fill(event.ampLo[0], event.ampLo[1]);
        amp_vs_x->Fill(event.xtrk, event.ampLo[0]);
        amp_vs_y->Fill(event.ytrk, event.ampLo[0]);
        MPV_amp0->Fill(event.ampLo[0]);
        MPV_amp1->Fill(event.ampLo[1]);

        //cout<<"TEST0"<<endl;

        t = tot(ith);

        //cout<<"TEST"<<endl;
        //tRes_avg = abs(0.5*(t[0]+t[1]) - event.timeMCP);  //time resolution with MCP, tMCP is bigger than tbar!
        tdiff[jentry] = (t[count][0] - t[count][1]);  // time difference method (without MCP)

        //cout<<tdiff[jentry]<<endl;

        //htRes_diff->Fill(tdiff[jentry]);
        tdiff_vs_x->Fill(event.xtrk, tdiff[jentry]);
        tdiff_vs_ampSIPM0->Fill(event.ampLo[0], tdiff[jentry]);
        tdiff_vs_ampSIPM1->Fill(event.ampLo[1], tdiff[jentry]);
        t0_minus_tMCP_vs_amp->Fill(event.ampLo[0], t[count][0] - event.timeMCP);
        t1_minus_tMCP_vs_amp->Fill(event.ampLo[1], t[count][1] - event.timeMCP);

        count += 1;

    }

    TF1 *l0 = new TF1("l0","landau", MPV_amp0->GetMean() - 2*MPV_amp0->GetRMS() , MPV_amp0->GetMean() + 6*MPV_amp0->GetRMS());
    TF1 *l1 = new TF1("l1","landau", MPV_amp1->GetMean() - 2*MPV_amp1->GetRMS() , MPV_amp1->GetMean() + 6*MPV_amp1->GetRMS());
    MPV_amp0->Fit(l0,"QR");
    MPV_amp1->Fit(l1,"QR");
    TF1 *l2 = new TF1("l2","landau", 0.7*l0->GetParameter(1) , 3*l0->GetParameter(1) );
    TF1 *l3 = new TF1("l3","landau", 0.7*l1->GetParameter(1) , 3*l1->GetParameter(1) );    
    MPV_amp0->Fit(l2,"QR");
    MPV_amp1->Fit(l3,"QR");
    //TF1 *l4 = new TF1("l4","landau", l2->GetParameter(1) - 0.7*l2->GetParameter(2) , l2->GetParameter(1) + 3*l3->GetParameter(2));
    //TF1 *l5 = new TF1("l5","landau", l3->GetParameter(1) - 0.7*l3->GetParameter(2) , l3->GetParameter(1) + 3*l3->GetParameter(2));    
    //MPV_amp0->Fit(l4,"QR");
    //MPV_amp1->Fit(l5,"QR");


    TGraphErrors *gr_tdiff_vs_x = new TGraphErrors();
  /*  x_arr_posCorr = xvalues();

   for(int i = 0; i<xbins.size()-1; i++){

        TH1D *t_ = tdiff_vs_x->ProjectionY("",xbins[i],xbins[i+1]);
        TF1* f_ = new TF1("f_","gaus",t_->GetMean() - 3*t_->GetRMS(),t_->GetMean() + 3*t_->GetRMS());
        t_->Fit(f_, "QR");
        c->cd();
        t_->Draw();
        c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46076/c_gaus_test_th%d_%s_%d.png",ith,argv[1],i));
        c->Clear();
        gr_tdiff_vs_x ->SetPoint(gr_tdiff_vs_x->GetN(), x_arr_posCorr[i], f_->GetParameter(1));
        gr_tdiff_vs_x ->SetPointError(gr_tdiff_vs_x->GetN()-1, 0, f_->GetParError(1));   

        }    

    
    gr_tdiff_vs_x->Fit(fit_xCorrection, "QRE");
    /*c->cd();
    tdiff_vs_x->Draw("colz");
    gr_tdiff_vs_x->SetMarkerSize(1);
    gr_tdiff_vs_x->SetMarkerStyle(20);
    gr_tdiff_vs_x->SetMarkerColor(kBlack);
    gr_tdiff_vs_x->Draw("P same");
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46076/c_tdiff_vs_x_with_graph_th%d_%s.png",ith,argv[1]));
    c->Clear(); */

    
    t0_minus_tMCP_vs_amp->Fit(fit_ampwalk0,"QR");
    t1_minus_tMCP_vs_amp->Fit(fit_ampwalk1,"QR");
    tdiff_vs_x->Fit(fit_xCorrection,"QR");
    
    /*c->cd();
    tdiff_vs_x->Draw("colz");
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46076/c_tdiff_vs_x_th%d_%s.png",ith,argv[1]));
    c->Clear();*/

    count = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      mytree->GetEntry(jentry);

      if(event.xtrk < 22 ||  event.xtrk > 27.6) continue;   //track cuts
      if(event.ytrk < 16. || event.ytrk > 18.85) continue;
      if(event.timeMCP < 100 ||  event.timeMCP > 118) continue;  
      if(event.ampLo[0]< 1230. || event.ampLo[0] > 4605) continue;
      if(event.ampLo[1]< 1290. || event.ampLo[1] > 4740) continue;

      //t = tot(ith);


      //y = fit_xCorrection->GetParameter(1)*event.xtrk + fit_xCorrection->GetParameter(0);
      y = fit_xCorrection->Eval(event.xtrk);
      tdiffPosCorr = tdiff[jentry] - y;   //x-position correction

      y_amp = fit_ampwalk0->Eval(event.ampLo[0]); 
      t0_AmpCorrected[jentry] = (t[count][0] - event.timeMCP) - y_amp;

      y_amp = fit_ampwalk1->Eval(event.ampLo[1]); 
      t1_AmpCorrected[jentry] = (t[count][1] - event.timeMCP) - y_amp;

      t0_minus_tMCP_vs_amp_corrected->Fill(event.ampLo[0], t0_AmpCorrected[jentry]);
      t1_minus_tMCP_vs_amp_corrected->Fill(event.ampLo[1], t1_AmpCorrected[jentry]);

      tavg_vs_x_ampCorrected->Fill(event.xtrk, 0.5*(t0_AmpCorrected[jentry] + t1_AmpCorrected[jentry]));
      tdiff_vs_x_ampCorrected->Fill(event.xtrk,(t0_AmpCorrected[jentry] - t1_AmpCorrected[jentry]));
      tdiff_vs_x_corrected->Fill(event.xtrk, tdiffPosCorr);
      htRes_diff_xposCorrected->Fill(tdiffPosCorr/2);

      count+=1;
    }

    tdiff_vs_x_ampCorrected->Fit(fit_xCorrection,"QR");
    tavg_vs_x_ampCorrected->Fit(fit_xCorrection1,"QR");

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      mytree->GetEntry(jentry);

      if(event.xtrk < 22 ||  event.xtrk > 27.6) continue;   //track cuts
      if(event.ytrk < 16.| event.ytrk > 18.85) continue;
      if(event.timeMCP < 100 ||  event.timeMCP > 118) continue; 
      if(event.ampLo[0]< 1230. || event.ampLo[0] > 4605) continue;
      if(event.ampLo[1]< 1290. || event.ampLo[1] > 4740) continue; 

      y = fit_xCorrection->Eval(event.xtrk);
      y1 = fit_xCorrection1->Eval(event.xtrk);
      tdiffAmpPosCorr = (t0_AmpCorrected[jentry] - t1_AmpCorrected[jentry]) - y;   //x-position correction
      tavgAmpPosCorr = (0.5*(t0_AmpCorrected[jentry] + t1_AmpCorrected[jentry])) - y1;

      tdiff_vs_x_ampXposCorrected->Fill(event.xtrk,tdiffAmpPosCorr);
      htRes_diff_AmpXposCorrected->Fill(tdiffAmpPosCorr/2);
      htRes_avg_AmpXposCorrected->Fill(tavgAmpPosCorr);

    }


    TF1 *fit_tRes_diff = new TF1("fit_tRes_diff","gaus",htRes_diff_xposCorrected->GetMean()-3*htRes_diff_xposCorrected->GetRMS(),htRes_diff_xposCorrected->GetMean()+3*htRes_diff_xposCorrected->GetRMS());
    htRes_diff_xposCorrected->Fit(fit_tRes_diff,"QR");
    TF1 *fit_tRes_diff1 = new TF1("fit_tRes_diff1","gaus",fit_tRes_diff->GetParameter(1)-4*fit_tRes_diff->GetParameter(2),fit_tRes_diff->GetParameter(1)+4*fit_tRes_diff->GetParameter(2));
    htRes_diff_xposCorrected->Fit(fit_tRes_diff1,"QR");
    TF1 *fit_tRes_diff2 = new TF1("fit_tRes_diff2","gaus",fit_tRes_diff1->GetParameter(1)-3*fit_tRes_diff1->GetParameter(2),fit_tRes_diff1->GetParameter(1)+3*fit_tRes_diff1->GetParameter(2));
    htRes_diff_xposCorrected->Fit(fit_tRes_diff2,"QR");

    TF1 *fit_tRes_diff3 = new TF1("fit_tRes_diff3","gaus",htRes_diff_AmpXposCorrected->GetMean()-3*htRes_diff_AmpXposCorrected->GetRMS(),htRes_diff_AmpXposCorrected->GetMean()+3*htRes_diff_AmpXposCorrected->GetRMS());
    htRes_diff_AmpXposCorrected->Fit(fit_tRes_diff3,"QR");
    TF1 *fit_tRes_diff4 = new TF1("fit_tRes_diff4","gaus",fit_tRes_diff3->GetParameter(1)-4*fit_tRes_diff3->GetParameter(2),fit_tRes_diff3->GetParameter(1)+4*fit_tRes_diff3->GetParameter(2));
    htRes_diff_AmpXposCorrected->Fit(fit_tRes_diff4,"QR");
    TF1 *fit_tRes_diff5 = new TF1("fit_tRes_diff5","gaus",fit_tRes_diff4->GetParameter(1)-3*fit_tRes_diff4->GetParameter(2),fit_tRes_diff4->GetParameter(1)+3*fit_tRes_diff4->GetParameter(2));
    htRes_diff_AmpXposCorrected->Fit(fit_tRes_diff5,"QR");

    TF1 *fit_tRes_diff6 = new TF1("fit_tRes_diff6","gaus",htRes_avg_AmpXposCorrected->GetMean()-3*htRes_avg_AmpXposCorrected->GetRMS(),htRes_avg_AmpXposCorrected->GetMean()+3*htRes_avg_AmpXposCorrected->GetRMS());
    htRes_avg_AmpXposCorrected->Fit(fit_tRes_diff6,"QR");
    TF1 *fit_tRes_diff7 = new TF1("fit_tRes_diff7","gaus",fit_tRes_diff6->GetParameter(1)-4*fit_tRes_diff6->GetParameter(2),fit_tRes_diff6->GetParameter(1)+4*fit_tRes_diff6->GetParameter(2));
    htRes_avg_AmpXposCorrected->Fit(fit_tRes_diff7,"QR");
    TF1 *fit_tRes_diff8 = new TF1("fit_tRes_diff8","gaus",fit_tRes_diff7->GetParameter(1)-3*fit_tRes_diff7->GetParameter(2),fit_tRes_diff7->GetParameter(1)+3*fit_tRes_diff7->GetParameter(2));
    htRes_avg_AmpXposCorrected->Fit(fit_tRes_diff8,"QR");

    
    tRes_vs_threshold1->SetPoint(tRes_vs_threshold1->GetN(), ith, fit_tRes_diff2->GetParameter(2));
    tRes_vs_threshold1->SetPointError(tRes_vs_threshold1->GetN()-1, 0, fit_tRes_diff2->GetParError(2));
    tRes_vs_threshold2->SetPoint(tRes_vs_threshold2->GetN(), ith, fit_tRes_diff5->GetParameter(2));
    tRes_vs_threshold2->SetPointError(tRes_vs_threshold2->GetN()-1, 0, fit_tRes_diff5->GetParError(2));


   /* c->cd();
    htRes_diff_xposCorrected->Draw("colz");
    c->SaveAs(Form("/eos/user/a/amkrishn/www/MTD_BTL_TB_FNAL_APRIL-MAY2021/46075/c_tRes_diff_xposCorrected_th%d_%s.png",ith,argv[1]));
    c->Clear();
  */

    fout->cd();
    tdiff_vs_x->Write();
    //tdiff_vs_x_1->Write();
    tavg_vs_x_ampCorrected->Write();
    htRes_avg_AmpXposCorrected->Write();
    //htRes_diff->Write();
    tdiff_vs_x_ampCorrected->Write();
    tdiff_vs_x_corrected->Fit(fit_xCorrection,"QR");
    tdiff_vs_x_corrected->Write();
    htRes_diff_xposCorrected->Write();
    htRes_diff_AmpXposCorrected->Write();
    t0_minus_tMCP_vs_amp_corrected->Write();
    t1_minus_tMCP_vs_amp_corrected->Write();
    

    //fit_xCorrection->Write();
    tdiff_vs_ampSIPM0->Write();
    tdiff_vs_ampSIPM1->Write();
    t0_minus_tMCP_vs_amp->Write();
    t1_minus_tMCP_vs_amp->Write();
    tdiff_vs_x_ampXposCorrected->Write();
    
    delete tdiff_vs_x;
    delete tdiff_vs_ampSIPM0;
    delete tdiff_vs_ampSIPM1;
    delete t0_minus_tMCP_vs_amp;
    delete t1_minus_tMCP_vs_amp;
    delete t1_minus_tMCP_vs_amp_corrected;
    delete t0_minus_tMCP_vs_amp_corrected;
 
    delete htRes_avg_AmpXposCorrected;
    //delete htRes_diff;
    delete htRes_diff_xposCorrected;
    //delete fit_tRes_diff;
    delete fit_xCorrection;
    delete tdiff_vs_x_corrected;
    delete fit_ampwalk0;
    delete fit_ampwalk1;
    delete tdiff_vs_x_ampCorrected;
    delete htRes_diff_AmpXposCorrected;
    delete tdiff_vs_x_ampXposCorrected;
  
}
    TH1D *h1AmpX = amp_vs_x->ProjectionX();
    TH1D *h1AmpY = amp_vs_y->ProjectionX();
    TH1D *h1MCP = MCP_time_vs_amp->ProjectionX();

    tRes_vs_threshold1 ->SetTitle(" ; threshold (mV/ns); #sigma_{bar} (ns) ");
    tRes_vs_threshold2 ->SetTitle(" ; threshold (mV/ns); #sigma_{bar} (ns) ");

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

    tRes_vs_threshold1->SetMarkerSize(1);
    tRes_vs_threshold1->SetMarkerStyle(20);
    tRes_vs_threshold1->SetMarkerColor(kBlack);

    tRes_vs_threshold2->SetMarkerSize(1);
    tRes_vs_threshold2->SetMarkerStyle(20);
    tRes_vs_threshold2->SetMarkerColor(kBlack);

    tRes_vs_threshold1->Draw("ALP");
    tRes_vs_threshold2->Draw("ALP");
    tRes_vs_threshold1->Write("tRes_vs_threshold_xCorrOnly");
    tRes_vs_threshold2->Write("tRes_vs_threshold");
  
  
    delete MCP_time_vs_amp;
    delete amp_sipm0_vs_sipm1Corr;
    delete amp_sipm0_vs_sipm1;
    delete hBeam;
    delete hprof2D;
    delete amp_vs_x;
    delete amp_vs_y;
    delete h1AmpX;
    delete h1AmpY;
    delete h1MCP;
    delete c;
    delete MPV_amp0;
    delete MPV_amp1;
    delete tRes_vs_threshold2;
    delete tRes_vs_threshold1;

    fout->Close();



    return 0;
}

