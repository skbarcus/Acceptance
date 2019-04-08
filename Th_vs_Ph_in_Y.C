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

void Th_vs_Ph_in_Y() 
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);
  TChain *T = new TChain("T");

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

  Int_t nevts = T->GetEntries();
  cout<<"nevts = "<<nevts<<endl;

  Int_t nbins = 200;
  //Define slices on y target for plotting theta vs. phi.
  Int_t nslices = 10;
  Double_t yplot_min=-0.05,yplot_max=0.05;
  Double_t ystep = (yplot_max-yplot_min)/10.;
  Double_t thplot_min=-0.1,thplot_max=0.1;
  Double_t phplot_min=-0.1,phplot_max=0.1;

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
  c1->Divide(3,4);
  //c1->SetLogy();

  //Draw y target.
  c1->cd(1);
  //Draw x target/theta target.
  T->Draw(Form("L.tr.tg_y>>h1(%d,%f,%f)",nbins,yplot_min,yplot_max));//,ct_1tr&&ct_trg&&ct_th&&ct_ph&&ct_y&&ct_pr&&ct_gc&&ct_dp);
  //h1->SetTitle("L.tr.tg_y");
  h1->SetTitle("Y-Target");
  //h1->SetTitleSize(10.06);
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetNdivisions(7);
  h1->GetYaxis()->SetLabelSize(0.12);
  h1->GetYaxis()->SetLabelOffset(-0.02);
  h1->GetYaxis()->SetTitleSize(0.06);
  h1->GetYaxis()->SetTitleOffset(0.75);
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetNdivisions(5);
  h1->GetXaxis()->SetLabelSize(0.12);
  h1->GetXaxis()->SetLabelOffset(0.0);
  h1->GetXaxis()->SetTitleSize(0.06);
  h1->GetXaxis()->SetTitleOffset(0.75);
  //h1->GetXaxis()->SetTitle("L.tr.tg_y");

  //Draw theta vs phi with no cuts.
  c1->cd(2);
  T->Draw(Form("L.tr.tg_th:L.tr.tg_ph>>h2(%d,%f,%f,%d,%f,%f)",nbins,thplot_min,thplot_max,nbins,phplot_min,phplot_max));
  //h2->SetTitle("L.tr.tg_y");
  h2->SetTitle("Total Solid Angle Acceptance");
  h2->GetYaxis()->CenterTitle(true);
  h2->GetYaxis()->SetNdivisions(7);
  h2->GetYaxis()->SetLabelSize(0.12);
  h2->GetYaxis()->SetLabelOffset(-0.02);
  h2->GetYaxis()->SetTitleSize(0.06);
  h2->GetYaxis()->SetTitleOffset(0.75);
  h2->GetXaxis()->CenterTitle(true);
  h2->GetXaxis()->SetNdivisions(5);
  h2->GetXaxis()->SetLabelSize(0.12);
  h2->GetXaxis()->SetLabelOffset(0.0);
  h2->GetXaxis()->SetTitleSize(0.06);
  h2->GetXaxis()->SetTitleOffset(0.75);
  //h2->GetXaxis()->SetTitle("L.tr.tg_ph");
  //h2->GetYaxis()->SetTitle("L.tr.tg_th");

  //Plot theta vs. phi is slices of y target to see how acceptance changes in y.
  for(Int_t i=0;i<nslices;i++)
    {
      c1->cd(3+i);
      T->Draw(Form("L.tr.tg_th:L.tr.tg_ph>>h%d(%d,%f,%f,%d,%f,%f)",3+i,nbins,thplot_min,thplot_max,nbins,phplot_min,phplot_max),Form("L.tr.tg_y>%f&&L.tr.tg_y<%f",yplot_min+i*ystep,yplot_min+(i+1)*ystep));
    }
  //h3->SetTitle("Total Solid Angle Acceptance");

  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
