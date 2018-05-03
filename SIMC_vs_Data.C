#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <math.h>

void SIMC_vs_Data() 
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);
  TChain *T = new TChain("T");
  /*
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3892.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3893.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3894.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4073.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074*.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4075*.root");
  */
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074.root");

  T->SetBranchStatus("*",0);
  T->SetBranchStatus("EKL.x_bj",1);
  T->SetBranchStatus("ExtEKL.x_bj",1);
  T->SetBranchStatus("L.tr.tg_y",1);
  T->SetBranchStatus("L.tr.n",1);
  T->SetBranchStatus("L.cer.asum_c",1);
  T->SetBranchStatus("L.prl1.e",1);
  T->SetBranchStatus("L.prl2.e",1);
  T->SetBranchStatus("L.tr.tg_th",1);
  T->SetBranchStatus("L.tr.tg_ph",1);
  T->SetBranchStatus("L.tr.d_th",1);
  T->SetBranchStatus("L.tr.d_ph",1);
  T->SetBranchStatus("L.tr.vz",1);
  T->SetBranchStatus("L.tr.p",1);
  T->SetBranchStatus("L.tr.tg_dp",1);
  T->SetBranchStatus("DBB.evtypebits",1);

  Double_t x_bj=0.,Ext_x_bj=0.,L_tr_tg_y[21],L_tg_th[21],L_tg_ph[21],L_d_th[21],L_d_ph[21],L_cer=0.,L_prl1_e=0.,L_prl2_e=0.,L_tr_p[21],L_tr_tg_dp[21],L_tr_n=0.,L_tr_vz,evtypebits;

  T->SetBranchAddress("EKL.x_bj",&x_bj);
  T->SetBranchAddress("ExtEKL.x_bj",&Ext_x_bj);
  T->SetBranchAddress("L.tr.tg_y",&L_tr_tg_y);
  T->SetBranchAddress("L.tr.n",&L_tr_n);
  T->SetBranchAddress("L.cer.asum_c",&L_cer);
  T->SetBranchAddress("L.prl1.e",&L_prl1_e);
  T->SetBranchAddress("L.prl2.e",&L_prl2_e);
  T->SetBranchAddress("L.tr.tg_th",L_tg_th);
  T->SetBranchAddress("L.tr.tg_ph",L_tg_ph);
  T->SetBranchAddress("L.tr.d_th",L_d_th);
  T->SetBranchAddress("L.tr.d_ph",L_d_ph);
  T->SetBranchAddress("L.tr.vz",L_tr_vz);
  T->SetBranchAddress("L.tr.p",L_tr_p);
  T->SetBranchAddress("L.tr.tg_dp",L_tr_tg_dp);
  T->SetBranchAddress("DBB.evtypebits",&evtypebits);

  TChain *SNT = new TChain("SNT");
  //SNT->Add("/home/skbarcus/Tritium/Simulations/SIMC/simc_gfortran.12/output/current.root");
  //SNT->Add("/home/skbarcus/Tritium/Simulations/SIMC/simc_gfortran.12/output/P0_3.055/Ly50_Lx70_Ld6_wide_2body_1M.root");
  //SNT->Add("/home/skbarcus/Tritium/Simulations/SIMC/simc_gfortran.12/output/P0_3.055_6GeV/Ly50_Lx70_Ld6_wide_2body_no_xp_yp_off_100k.root");
  //SNT->Add("/home/skbarcus/Tritium/Simulations/SIMC/simc_gfortran.12/output/P0_3.055_6GeV/Ly50_Lx70_Ld6_wide_2body_xp_-25_yp_-1_100k.root");
  SNT->Add("/home/skbarcus/Tritium/Simulations/SIMC/simc_gfortran.12/output/P0_3.055_6GeV/Ly50_Lx70_Ld6_wide_2body_xp_-25_yp_-1_10M.root");
  //SNT->Add("/home/skbarcus/Tritium/Simulations/SIMC/simc_gfortran.12/output/current.root");

  Int_t nevts = T->GetEntries();
  cout<<"nevts = "<<nevts<<endl;

  Int_t nbins = 100;

  //Define cuts:
  Double_t vzmin = -0.08, vzmax = 0.08;
  Double_t ymin = -0.028, ymax = 0.028;
  Double_t thmin = -0.07, thmax = 0.07;
  Double_t phmin = -0.05, phmax = 0.05;
  Double_t thdmin = 0.6, thdmax = 1.2;
  Double_t phdmin = -0.06, phdmax = 0.06;
  Double_t dpmin = -0.06, dpmax = 0.06;
  Double_t gcmin = 0.;

  TCut ct_norm = "Weight*Normfac/20000.";
  TCut ct_1tr = "L.tr.n==1";
  TCut ct_trg = "(DBB.evtypebits&1<<3)==1<<3";
  TCut ct_vz = Form("L.tr.vz>%f&&L.tr.vz<%f",vzmin,vzmax);
  TCut ct_y = Form("L.tr.tg_y>%f&&L.tr.tg_y<%f",ymin,ymax);
  TCut ct_y_mc = Form("e_ytar>%f&&e_ytar<%f",ymin*100.,ymax*100.);
  TCut ct_dp = Form("L.tr.tg_dp>%f&&L.tr.tg_dp<%f",dpmin,dpmax);
  TCut ct_dp_mc = Form("e_delta>%f&&e_delta<%f",dpmin*100.,dpmax*100.);
  TCut ct_th = Form("L.tr.tg_th>%f&&L.tr.tg_th<%f",thmin,thmax);
  TCut ct_th_mc = Form("e_xptar>%f&&e_xptar<%f",TMath::ATan(thmin),TMath::ATan(thmax));
  TCut ct_ph = Form("L.tr.tg_ph>%f&&L.tr.tg_ph<%f",phmin,phmax);
  TCut ct_ph_mc = Form("e_yptar>%f&&e_yptar<%f",TMath::ATan(phmin),TMath::ATan(phmax));
  TCut ct_th_d = Form("L.tr.d_th>%f&&L.tr.d_th<%f",thdmin,thdmax);
  TCut ct_ph_d = Form("L.tr.d_ph>%f&&L.tr.d_ph<%f",phdmin,phdmax);
  TCut ct_pr = "L.prl1.e>(-L.prl2.e+2000)";
  TCut ct_gc = Form("L.cer.asum_c>%f",gcmin);

 
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  //c1->SetLogy();

  //Draw x target/theta target.
  T->Draw(Form("TMath::RadToDeg()*TMath::ATan(L.tr.tg_th)>>h1(%d,-5.,5.)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  SNT->Draw(Form("TMath::RadToDeg()*e_xptar>>h2(%d,-5.,5.)",nbins),ct_norm*(ct_y_mc&&ct_th_mc&&ct_ph_mc&&ct_dp_mc),"same");
  h2->SetLineColor(2);
  h1->SetTitle("xtar Data (ATan(target theta)) vs. xtar SIMC");
  h1->GetXaxis()->SetTitle("Theta (degrees)");

  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();

  //Draw y target/phi target.
  T->Draw(Form("21.+TMath::RadToDeg()*TMath::ATan(L.tr.tg_ph)>>h3(%d,18.,24.)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  SNT->Draw(Form("21.+TMath::RadToDeg()*e_yptar>>h4(%d,18.,24.)",nbins),ct_norm*(ct_y_mc&&ct_th_mc&&ct_ph_mc&&ct_dp_mc),"same");
  h4->SetLineColor(2);
  cout<<"# Entries SIMC = "<<h4->GetEntries()<<endl;
  h3->SetTitle("ytar Data (21+ATan(target phi)) vs. ytar SIMC");
  h3->GetXaxis()->SetTitle("Phi (degrees)");

  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();

  //Draw delta p target.
  T->Draw(Form("L.tr.tg_dp>>h5(%d,-0.1,0.1)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  SNT->Draw(Form("e_delta/100.>>h6(%d,-0.1,0.1)",nbins),ct_norm*(ct_y_mc&&ct_th_mc&&ct_ph_mc&&ct_dp_mc),"same");
  h6->SetLineColor(2);
  h5->SetTitle("dp Data vs. e_delta/100 SIMC");
  h5->GetXaxis()->SetTitle("dp=(P-P0)/P0");

  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();

  //Draw xptar data vs xptar SIMC.
  T->Draw(Form("TMath::ATan(L.tr.tg_th)>>h7(%d,-0.1,0.1)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  SNT->Draw(Form("e_xptar>>h8(%d,-0.1,0.1)",nbins),ct_norm*(ct_y_mc&&ct_th_mc&&ct_ph_mc&&ct_dp_mc),"same");
  h8->SetLineColor(2);
  h7->SetTitle("xptar Data (ATan(target theta)) vs. xptar SIMC");
  h7->GetXaxis()->SetTitle("Theta (radians)");

  TCanvas* c5=new TCanvas("c5");
  c5->SetGrid();

  //Draw yptar data vs yptar SIMC.
  T->Draw(Form("TMath::ATan(L.tr.tg_ph)>>h9(%d,-0.1,0.1)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  SNT->Draw(Form("e_yptar>>h10(%d,-0.1,0.1)",nbins),ct_norm*(ct_y_mc&&ct_th_mc&&ct_ph_mc&&ct_dp_mc),"same");
  h10->SetLineColor(2);
  h9->SetTitle("yptar Data (ATan(target phi)) vs. yptar SIMC");
  h9->GetXaxis()->SetTitle("yp target (radians)");

  //T->Draw("EKL.x_bj>>h2(1000,0.5,3.5)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1");

  /*
  //Create and fill histogram with Xbj values with various cuts applied to the data.
  TH1D *h1 = new TH1D("h1","Xbj" , 1000, 0., 4.);
  Int_t evt_136=0,evt_8=0;
  Int_t maxtracks=0;
  for(Int_t i=0;i<nevts;i++) 
    {
      T->GetEntry(i);
      if(L_tr_tg_y[0]>ymin && L_tr_tg_y[0]<ymax && L_tr_n==1 && (evtypebits&1<<3)==1<<3 && L_prl1_e>(-L_prl2_e+2000) && L_tg_ph[0]>phmin && L_tg_ph[0]<phmax && L_tg_th[0]>thmin && L_tg_th[0]<thmax)
	{
	  h1->Fill(Ext_x_bj);
	}
    }
  h1->Draw();
  */

  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
