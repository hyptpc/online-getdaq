#include "padHelper_new.hh"

#define N_ASAD 30 
void merge_neg() 
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);
  
  int asadNum[N_ASAD];
  int zapNum[N_ASAD];
  for(int i=0; i<N_ASAD; i++) 
  {
    asadNum[i] = 1+i;
    if(i<8) zapNum[i] = 17+i;
    else if( i>=8 && i<24 ) zapNum[i] = i-7;
    else zapNum[i] = 1+i;
    //cout << "i = " << i+1 << ", zap# = " << zapNum[i] << endl;
  }

  TFile* file[4];
  TH2D* hist[4];
  TH2Poly* poly[4];
  TGraph* graph[30][256][4];

  TH2D* histall = new TH2D("histall", ";Ch;ZAP#", 260,0,260, 31,1,32);
  TH2Poly* polyall = new TH2Poly("polyall","TPC",-300,300,-300,300);
  double X[5];
  double Y[5];

  for( int i = 0; i< 32; i++)
  {
    double pLength = padParameter[i][5];
    double st = 180.-(360./padParameter[i][3])*padParameter[i][1]/2.;
    double sTheta  = (-1+st/180.)*TMath::Pi();
    double dTheta  = (360./padParameter[i][3])/180.*TMath::Pi();
    double cRad    = padParameter[i][2];
    int    nPad    = padParameter[i][1];
    for( int j = 0; j< nPad; j++)
    {
      X[1] = (cRad+(pLength/2.))*cos(j*dTheta+sTheta);
      X[2] = (cRad+(pLength/2.))*cos((j+1)*dTheta+sTheta);
      X[3] = (cRad-(pLength/2.))*cos((j+1)*dTheta+sTheta);
      X[4] = (cRad-(pLength/2.))*cos(j*dTheta+sTheta);
      X[0] = X[4];
      Y[1] = (cRad+(pLength/2.))*sin(j*dTheta+sTheta);
      Y[2] = (cRad+(pLength/2.))*sin((j+1)*dTheta+sTheta);
      Y[3] = (cRad-(pLength/2.))*sin((j+1)*dTheta+sTheta);
      Y[4] = (cRad-(pLength/2.))*sin(j*dTheta+sTheta);
      Y[0] = Y[4];
      for(int ii=0; ii<5; ii++) X[ii] -=143;
      polyall->AddBin(5,X,Y);
    }
  }


  for(int i=0; i<4; i++) 
  {
    file[i] = new TFile(Form("rootfile/allneg%d.root",i+1), "read");
    hist[i] = (TH2D*)file[i]->Get("hist");
    poly[i] = (TH2Poly*)file[i]->Get("TPC");
    for(int a=0; a<30; a++)
      for(int c=0; c<256; c++)
      {
	if((TGraph*)file[i]->Get(Form("a%d_c%d", a, c)))
	  graph[a][c][i] = (TGraph*)file[i]->Get(Form("a%d_c%d", a, c));
      }
  }


  for(int k=0; k<4; k++)
    for(int i=0; i<256; i++)
      for(int j=0; j<31; j++)
      {
	//if(hist[k]->GetBinContent(i+1, j+1))
	if(hist[k]->GetBinContent(i+1, j+1)!=5)
	  histall->SetBinContent(i+1, j+1, hist[k]->GetBinContent(i+1, j+1));
      }

  for(int k=0; k<4; k++){
    for(int i=0; i<5768; i++) {
      if(poly[k]->GetBinContent(i+1))
	polyall->SetBinContent(i+1, poly[k]->GetBinContent(i+1));
    }
  }

  Double_t levels[5];
  for(int i=0; i<5; i++) levels[i] = i+0;

  TCanvas* c1 = new TCanvas("c1","zap", 1400,800);
  c1->cd()->SetGrid();
  histall->SetContour(sizeof(levels)/sizeof(Double_t), levels);
  histall->SetMaximum(5);
  histall->SetMinimum(0);
  histall->Draw("colz");
  
  TCanvas* c2 = new TCanvas("c2","pad", 800,800);
  c2->cd();
  polyall->SetContour(sizeof(levels)/sizeof(Double_t), levels);
  polyall->SetMaximum(5);
  polyall->SetMinimum(0);
  polyall->Draw("colz");
  
  TCanvas *c3[4];
  for(int k=0; k<4; k++) {
    c3[k] = new TCanvas(Form("c3_%d",k), Form("c3_%d",k), 1400,600);
    c3[k]->Divide(4,2);
  }
  for(int l=0; l<4; l++) {
    for(int k=0; k<4; k++) {
      for(int i=0; i<8; i++) {
	c3[k]->cd(i+1)->SetGrid();
	c3[k]->cd(i+1)->SetMargin(0.15, 0.1, 0.1, 0.1);
	if(i+k*8>=30) continue;
	for(int j=0; j<256; j++) {
	  if(j==0 && l==0) {
	    graph[i+k*8][j][l]->SetTitle(Form("AsAd#%d + ZAP#%d; Input [V]; ADC Max [ch]", asadNum[i+k*8], zapNum[i+k*8]));
	    graph[i+k*8][j][l]->SetMinimum(0);
	    graph[i+k*8][j][l]->SetMaximum(2500);
	    graph[i+k*8][j][l]->Draw("APL");
	  }
	  else graph[i+k*8][j][l]->Draw("PL");
	}
      }
    }
  }
}
