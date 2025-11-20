void merge() 
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  TFile* file[3];
  TH2D* hist[3];
  TH2D* histall = new TH2D("histall", ";Ch;ZAP#", 260,0,260, 31,1,32);

  for(int i=0; i<3; i++) 
  {
    file[i] = new TFile(Form("rootfile/t%d.root",i+1), "read");
    hist[i] = (TH2D*)file[i]->Get("hist");
  }

  int asad[31] = { 0, 0, 0, 0, 
    		   0, 0, 0, 7,
		   8, 7, 6, 3, 
		   1, 2, 4, 5, 
		   1, 2, 3, 4,
		   5, 6, 7, 8, 
		   1, 2, 3, 4,
		   5, 6, 0};


  for(int k=0; k<3; k++)
    for(int i=0; i<256; i++)
      for(int j=0; j<31; j++)
      {
	if(hist[k]->GetBinContent(i+1, j+1))
	  histall->SetBinContent(i+1, j+1, hist[k]->GetBinContent(i+1, j+1));
      }
  for(int i=0; i<4; i++)
  for(int j=0; j<31; j++)
  {
    histall->SetBinContent(257+i, j+1, asad[j]);
  }

  Double_t levels[8];
  for(int i=0; i<8; i++) levels[i] = i+1;

  TCanvas* c1 = new TCanvas("c1","zap", 1400,800);
  c1->cd()->SetGrid();
  histall->SetContour(sizeof(levels)/sizeof(Double_t), levels);
  histall->SetMaximum(9);
  histall->SetMinimum(1);
  histall->Draw("colz");
}
