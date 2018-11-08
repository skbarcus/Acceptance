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

//const Int_t dp_low_bin = 0;
//const Int_t dp_high_bin = 0;

void Final_Cuts() 
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);
  TChain *T = new TChain("T");
  
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3892.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3893.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3894.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4073.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074*.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4075*.root");
  
  //T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074.root");

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

  Double_t f_th_low = -0.042;
  Double_t f_th_high = 0.049;
  Double_t f_ph_low = -0.03;
  Double_t f_ph_high = 0.03;
  Double_t f_dp_low = -0.02;
  Double_t f_dp_high = 0.03;
  //Double_t f_y_low = -0.08;
  //Double_t f_y_high = 0.09;
  Double_t f_y_low = -0.03;//-0.03
  Double_t f_y_high = 0.028;//0.03

  TCut ct_norm = "Weight*Normfac/5500.";
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
  gStyle->SetTitleFontSize(0.08);

  //Draw delta p target.
  T->Draw(Form("L.tr.tg_dp>>h5(%d,-0.1,0.1)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  h5->SetTitle("dP Acceptance");
  h5->GetXaxis()->SetTitle("dp=(P-P0)/P0");
  h5->GetYaxis()->SetTitle("Counts");
  h5->GetYaxis()->CenterTitle(true);
  h5->GetYaxis()->SetLabelSize(0.04);
  h5->GetYaxis()->SetTitleSize(0.06);
  h5->GetYaxis()->SetTitleOffset(0.85);
  h5->GetXaxis()->CenterTitle(true);
  h5->GetXaxis()->SetLabelSize(0.04);
  h5->GetXaxis()->SetTitleSize(0.06);
  h5->GetXaxis()->SetTitleOffset(0.77);

  //cout<<h5->GetBinContent(41)<<endl;
  //cout<<h5->FindBin(f_dp_low)<<endl;

  Int_t dp_low_bin = h5->FindBin(f_dp_low);
  Int_t dp_high_bin = h5->FindBin(f_dp_high);

  //Copy h1 in a clone h1c. Set range and color for h1c
  TH1F *h5c = (TH1F*)h5->Clone();
  h5c->GetXaxis()->SetRange(dp_low_bin,dp_high_bin);
  h5c->SetFillColor(2);
  h5c->Draw("same");

  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();

  //Draw xptar data.
  T->Draw(Form("TMath::ATan(L.tr.tg_th)>>h7(%d,-0.1,0.1)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  h7->SetTitle("#theta Acceptance");
  h7->GetXaxis()->SetTitle("#theta (radians)");
  h7->GetYaxis()->SetTitle("Counts");
  h7->GetYaxis()->CenterTitle(true);
  h7->GetYaxis()->SetLabelSize(0.04);
  h7->GetYaxis()->SetTitleSize(0.06);
  h7->GetYaxis()->SetTitleOffset(0.85);
  h7->GetXaxis()->CenterTitle(true);
  h7->GetXaxis()->SetLabelSize(0.04);
  h7->GetXaxis()->SetTitleSize(0.06);
  h7->GetXaxis()->SetTitleOffset(0.77);

  Int_t th_low_bin = h7->FindBin(f_th_low);
  Int_t th_high_bin = h7->FindBin(f_th_high);

  //Copy h1 in a clone h1c. Set range and color for h1c
  TH1F *h7c = (TH1F*)h7->Clone();
  h7c->GetXaxis()->SetRange(th_low_bin,th_high_bin);
  h7c->SetFillColor(2);
  h7c->Draw("same");

  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();

  //Draw yptar data vs yptar SIMC.
  T->Draw(Form("TMath::ATan(L.tr.tg_ph)>>h9(%d,-0.1,0.1)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  h9->SetTitle("#phi Acceptance");
  //h9->SetTitle("yptar Data (ATan(target phi)) vs. yptar SIMC");
  h9->GetXaxis()->SetTitle("#phi (radians)");
  h9->GetYaxis()->SetTitle("Counts");
  h9->GetYaxis()->CenterTitle(true);
  h9->GetYaxis()->SetLabelSize(0.04);
  h9->GetYaxis()->SetTitleSize(0.06);
  h9->GetYaxis()->SetTitleOffset(0.85);
  h9->GetXaxis()->CenterTitle(true);
  h9->GetXaxis()->SetLabelSize(0.04);
  h9->GetXaxis()->SetTitleSize(0.06);
  h9->GetXaxis()->SetTitleOffset(0.77);

  Int_t ph_low_bin = h9->FindBin(f_ph_low);
  Int_t ph_high_bin = h9->FindBin(f_ph_high);

  //Copy h1 in a clone h1c. Set range and color for h1c
  TH1F *h9c = (TH1F*)h9->Clone();
  h9c->GetXaxis()->SetRange(ph_low_bin,ph_high_bin);
  h9c->SetFillColor(2);
  h9c->Draw("same");

  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();

  //Draw y target data.
  //T->Draw(Form("L.tr.vz>>h1(%d,-0.12,0.12)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_pr&&ct_gc&&ct_dp);
  T->Draw(Form("L.tr.tg_y>>h1(%d,-0.04,0.04)",nbins),ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_pr&&ct_gc&&ct_dp);
  h1->SetTitle("Y Target Acceptance");
  //h9->SetTitle("yptar Data (ATan(target phi)) vs. yptar SIMC");
  h1->GetXaxis()->SetTitle("Y Target (m)");
  h1->GetYaxis()->SetTitle("Counts");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetTitleSize(0.06);
  h1->GetYaxis()->SetTitleOffset(0.85);
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetLabelSize(0.04);
  h1->GetXaxis()->SetTitleSize(0.06);
  h1->GetXaxis()->SetTitleOffset(0.77);

  Int_t y_low_bin = h1->FindBin(f_y_low);
  Int_t y_high_bin = h1->FindBin(f_y_high);

  //Copy h1 in a clone h1c. Set range and color for h1c
  TH1F *h1c = (TH1F*)h1->Clone();
  h1c->GetXaxis()->SetRange(y_low_bin,y_high_bin);
  h1c->SetFillColor(2);
  h1c->Draw("same");


  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
