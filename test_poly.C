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
//#define theta1 21.;
Int_t debug = 0;
Int_t nevts = 10000;
Double_t PI = 3.14159265359;
Double_t deg2rad = PI/180.;
Double_t Y = 0.;
Double_t r = 1.;
Double_t xp = -0.;     //X coordinate of wire1's point.
Double_t yp = 0.;      //Y coordinate of wire1's point.
Double_t xp1 = -37.5;  //Top left hole x (mm).
Double_t yp1 = 75.;    //Top left hole y (mm).
Int_t rows = 7, columns = 7;
const Int_t nholes = rows * columns;
Double_t px[nholes] = {};            //X coordinate of hole centers.
Double_t py[nholes] = {};            //Y coordinate of hole centers.
Double_t qx[nholes] = {};            //X coordinate of random point in hole.
Double_t qy[nholes] = {};            //Y coordinate of random point in hole.
Int_t hole_num = 0;
Double_t horz_sep = 12.5, vert_sep = 25.;
Double_t l = 0., theta = 0.;

Double_t line(Double_t *X, Double_t *par)
{
  Y = m*X[0]+b;
  return Y;
}

Double_t circlep(Double_t *X, Double_t *par)
{
  Y = pow(par[0]*par[0]-(X[0]-par[1])*(X[0]-par[1]),0.5) + par[2];
  return Y;
}

Double_t circlen(Double_t *X, Double_t *par)
{
  Y = -pow(par[0]*par[0]-(X[0]-par[1])*(X[0]-par[1]),0.5) + par[2];
  return Y;
}

void test_poly() 
{
  //Define new random number generator (could use TRandom2,3 instead).
  TRandom *r1=new TRandom();
  r1->SetSeed(0);            //New random seed each time.
  
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Generate the points representing the hole centers.
  for(Int_t i=0;i<rows;i++)
    {
      for(Int_t j=0;j<columns;j++)
	{
	  px[hole_num] = xp1+j*horz_sep;
	  py[hole_num] = yp1-i*vert_sep;
	  hole_num++;
	}
    }

  //Create canvass to draw sieve on.	  
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();

  //TCanvas* c2=new TCanvas("c2");
  //c2->SetGrid();
  Double_t testx[4] = {10,-10,-10,10};
  Double_t testy[4] = {10,10,-10,-10};
  //Define a polygon.
  TGeoPolygon *poly = new TGeoPolygon(4);
  poly->SetXY(testx,testy);
  poly->FinishPolygon();
  //poly->SetLineColor(2);
  //poly->SetMarkerColor(2);
  c1->Update();
  //poly->Paint("same");
  poly->Draw("");
  poly->FinishPolygon();
  //poly->DrawPolygon();
  //poly->SetDrawOption("same");
  c1->Update();
  //poly->Draw("ACp*");
  //poly->Draw("same C");
  cout<<"Poly area = "<<poly->Area()<<endl;
  //cout<<"IsIllegalCheck() = "<<poly->IsIllegalCheck()<<endl;
  cout<<"IsFinished() = "<<poly->IsFinished()<<endl;
  cout<<"GetDrawOption() = "<<poly->GetDrawOption()<<endl;

  TPolyLine *pline = new TPolyLine(4,testx,testy);
  pline->SetLineColor(2);
  //pline->Draw("same");
  //cout<<"GetDrawOption() = "<<pline->GetDrawOption()<<endl;

  st->Stop();
  cout<<"   CPU time = "<<st->CpuTime()<<"   Real time = "<<st->RealTime()<<endl;
}
