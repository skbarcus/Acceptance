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
const Int_t nloop = 100;
Double_t PI = 3.14159265359;
Double_t deg2rad = PI/180.;
Double_t Y = 0.;
Double_t r = 1.;
Double_t xp = -0.;     //X coordinate of wire1's point.
Double_t yp = 0.;      //Y coordinate of wire1's point.
Double_t xp1 = -37.5;  //Top left hole x (mm).
Double_t yp1 = 75.;    //Top left hole y (mm).
Int_t nrows = 7, ncols = 7;
const Int_t nholes = nrows * ncols;
Double_t px[nholes] = {};            //X coordinate of hole centers.
Double_t py[nholes] = {};            //Y coordinate of hole centers.
Double_t qx[nholes] = {};            //X coordinate of random point in hole.
Double_t qy[nholes] = {};            //Y coordinate of random point in hole.
const Int_t nout = 2 * nrows + 2 * ncols - 4;
Double_t outx[nout] = {};            //X coordinate of outer hole point for TGeoPolygon.
Double_t outy[nout] = {};            //Y coordinate of outer hole point for TGeoPolygon.
Double_t linex[nout+1] = {};            //X coordinate of outer hole point for TPolyLine visual.
Double_t liney[nout+1] = {};            //Y coordinate of outer hole point for TPolyLine visual.
Double_t areas[nloop] = {};
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

void Optics_Sieve_Uncertainty_MC() 
{
  //Define new random number generator (could use TRandom2,3 instead).
  TRandom *r1=new TRandom();
  r1->SetSeed(0);            //New random seed each time.
  
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  //Generate the points representing the hole centers.
  for(Int_t i=0;i<nrows;i++)
    {
      for(Int_t j=0;j<ncols;j++)
	{
	  px[hole_num] = xp1+j*horz_sep;
	  py[hole_num] = yp1-i*vert_sep;
	  hole_num++;
	}
    }

  //Create canvass to draw sieve on.	  
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();

  //Plot points for the hole centers.
  //TGraph *gr = new TGraph(2,px,py);
  auto gr = new TGraph(nholes,px,py);
  auto axis = gr->GetXaxis();
  axis->SetLimits(-40.,40.);
  gr->SetTitle("Sieve Plate");
  gr->GetXaxis()->SetTitle("Horizontal Distance from Sieve Center (mm)");
  gr->GetYaxis()->SetTitle("Vertical Distance from Sieve Center (mm)");
  gr->GetHistogram()->SetMaximum(80.);
  gr->GetHistogram()->SetMinimum(-80.);
  gr->Draw("Ap");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.1);

  //Draw radii of circles around the hole centers. 
  //Reset hole_num.
  hole_num = 0;
  TF1 *fc_p[nholes];
  char *circlep_name = new char[100];
  TF1 *fc_n[nholes];
  char *circlen_name = new char[100];
  for(Int_t i=0;i<nrows;i++)
    {
      for(Int_t j=0;j<ncols;j++)
	{
	  sprintf(circlep_name,"fc_p%d",hole_num);
	  sprintf(circlen_name,"fc_n%d",hole_num);
	  fc_p[hole_num] = new TF1(circlep_name, circlep, -r+xp1+j*horz_sep, r+xp1+j*horz_sep, 4);
	  fc_n[hole_num] = new TF1(circlen_name, circlen, -r+xp1+j*horz_sep, r+xp1+j*horz_sep, 4);
	  //cout<<"circlep_name = "<<circlep_name<<endl;
	  //cout<<"circle "<<hole_num<<" -> ("<<xp1+j*horz_sep<<","<<yp1-i*vert_sep<<")"<<endl;
	  
	  fc_p[hole_num]->SetParameter(0,r);
	  fc_p[hole_num]->SetParameter(1,xp1+j*horz_sep);
	  fc_p[hole_num]->SetParameter(2,yp1-i*vert_sep);
	  fc_p[hole_num]->Draw("same");
	  fc_p[hole_num]->SetNpx(1000);
	  fc_p[hole_num]->SetLineColor(1);

	  fc_n[hole_num]->SetParameter(0,r);
	  fc_n[hole_num]->SetParameter(1,xp1+j*horz_sep);
	  fc_n[hole_num]->SetParameter(2,yp1-i*vert_sep);
	  fc_n[hole_num]->Draw("same");
	  fc_n[hole_num]->SetNpx(1000);
	  fc_n[hole_num]->SetLineColor(1);
	  
	  hole_num++;
	}
    }
  
  if(debug==1)
    {
      for(Int_t i=0;i<nholes;i++)
	{
	  cout<<"px["<<i<<"] = "<<px[i]<<"   py["<<i<<"] = "<<py[i]<<endl;
	}
    }

  //Loop over multiple randomly generated areas.
  for(Int_t m=0;m<nloop;m++)
    {
      //Create loop that selects a point randomly in each of the circles.
      for(Int_t i=0;i<nholes;i++)
	{
	  //Generate random length from cricle center inside the circle.
	  l = r1->Uniform(0.,r);
	  //Generate a random angle pivoting about the circle center.
	  theta = r1->Uniform(0.,2*PI);
	  //cout<<"l = "<<l<<"   theta = "<<theta<<endl;
	  
	  //Transform the radial coordinates of the points back into cartesian coordinates.
	  qx[i] = l * TMath::Cos(theta) + px[i];
	  qy[i] = l * TMath::Sin(theta) + py[i];
	  //cout<<"qx["<<i<<"] = "<<qx[i]<<"   qy["<<i<<"] = "<<qy[i]<<endl;
	}
      
      //Plot the random points in the sieve holes.
      if(m<5)  //Only plot the first few random point to save time.
	{
	  auto gr1 = new TGraph(nholes,qx,qy);
	  gr1->Draw("same p");
	  gr1->SetMarkerColor(2);
	  gr1->SetMarkerStyle(20);
	  gr1->SetMarkerSize(0.4);
	}
      //Fill array with the coordinates of the outer edge holes. This array needs to be filled clockwise for the TGeoPolygon to be formed correctly. Start filling from the upper right outer hole.
      
      //Fill top row of sieve's outer holes starting from the top right going towards the top left. 
      for(Int_t i=0;i<ncols;i++)
	{
	  outx[i] = qx[ncols-1-i];
	  outy[i] = qy[ncols-1-i];
	}
      //Fill left column from top to bottom.
      for(Int_t i=0;i<(nrows-2);i++)
	{
	  outx[i+ncols] = qx[(i+1)*ncols];
	  outy[i+ncols] = qy[(i+1)*ncols];
	}
      //Fill bottom row from left to right. 
      for(Int_t i=0;i<(ncols);i++)
	{
	  outx[i+ncols+nrows-2] = qx[ncols*(nrows-1)+i];
	  outy[i+ncols+nrows-2] = qy[ncols*(nrows-1)+i];
	}
      //Fill right column from bottom to top.
      for(Int_t i=0;i<(nrows-2);i++)
	{
	  outx[i+ncols+nrows-2+ncols] = qx[(nrows-1)*ncols-1-(i*ncols)];
	  outy[i+ncols+nrows-2+ncols] = qy[(nrows-1)*ncols-1-(i*ncols)];
	}
      
      //Print coordinates of outer holes.
      if(debug==1)
	{
	  for(Int_t i=0;i<nout;i++)
	    {
	      cout<<"outx["<<i<<"] = "<<outx[i]<<"   outy["<<i<<"] = "<<outy[i]<<endl;
	    }
	}
      
      //Define a polygon.
      TGeoPolygon *poly = new TGeoPolygon(nout);
      poly->SetXY(outx,outy);
      //poly->Draw("same");
      //cout<<"Poly area = "<<poly->Area()<<endl;
      //Fill areas array with the area of this loop's polygon.
      areas[m] = poly->Area();      
      
      //Draw a TPolyLine along the outer holes connecting the points randomly selected in the holes.
      //Need to make a new array for this that adds one more entry where the last entry is the same as the first to connect the final point back to the first visually. 
      if(m<5)   //Only draw the lines for the first few events to save time. 
	{
	  for(Int_t i=0;i<nout;i++)
	    {
	      linex[i] = outx[i];
	      liney[i] = outy[i];
	    }
	  //Fill final line array element with first point from out array to connect the final portion.
	  linex[nout] = outx[0];
	  liney[nout] = outy[0];
	  
	  TPolyLine *pline = new TPolyLine(nout+1,linex,liney);
	  pline->SetLineColor(2);
	  pline->Draw("same");
	}
    }//End loop of randomly gnerated areas.
  
  if(debug==1)
    {
      for(Int_t i=0;i<nloop;i++)
	{
	  cout<<"Areas["<<i<<"] = "<<areas[i]<<endl;
	}
    }

  //Make a histogram of the areas to find the uncertainty in the solid angle due to the sieve.
  TH1* h1 = new TH1D("h1", "Distribution of Areas", 301, 11100., 11400.);
  //TH1* h1 = new TH1D("h1", "Outer Edge Areas for Random Points in Outer Holes", 301, 1700., 2100.);
  //Fill areas histo.
  for(Int_t i=0;i<nloop;i++)
    {
      h1->Fill(areas[i]);
    }

  //Create canvass to draw areas histo on.	  
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  h1->GetXaxis()->SetTitle("Area Inside Connected Points Generated in Each Outer Hole Randomly (mm^{2})");
  h1->GetYaxis()->SetTitle("Counts");
  h1->Draw();

  st->Stop();
  cout<<"   CPU time = "<<st->CpuTime()<<"   Real time = "<<st->RealTime()<<endl;
}
