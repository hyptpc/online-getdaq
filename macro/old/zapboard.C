void zapboard()
{
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kRainBow);

  TH2D* hist = new TH2D("hist","; Ch; ZAP#", 256,0,256, 31,0,31);
  hist->Fill(3,3);

  TCanvas* c1 = new TCanvas("c1","zap", 1400,800);
  c1->cd()->SetGrid();
  hist->Draw("col");
}
