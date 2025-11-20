#include "padHelper_new.hh"

void test2() 
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  TFile* file;
  TH2Poly* poly;

  file = new TFile("rootfile/allpos1.root", "read");
  poly = (TH2Poly*)file->Get("TPC");


  TCanvas* c2 = new TCanvas("c2","pad", 800,800);
  c2->cd();
  poly->Draw("colz");
}
