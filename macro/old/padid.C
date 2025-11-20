#include "padHelper_new.hh"

#define EVENTNUM 1
#define PADNUM 5768

TCanvas *c1;
TH2Poly *poly;

void padid(int padID = 2625)
{
  gStyle->SetOptStat(0);
  poly = new TH2Poly("TPC","TPC",-300,300,-300,300);
  TPolyLine* p[PADNUM];

  //double ll = 586.;
  double ll = 500.;
  double l = (ll/2.)/(1+sqrt(2.));
  Double_t px[9]={-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),
    l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),
    -l*(1+sqrt(2.))};
  Double_t py[9]={l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,
    -l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,
    l};
  TPolyLine* pLine= new TPolyLine(9,px,py);
  pLine->SetLineColor(1);


  int npoint = 5;
  double X[5];
  double Y[5];
  int npl = 0;

  for( int i = 0; i< 32; i++)
  {
    //if(i==10) cout << " npl : " << npl << endl;
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
      poly->AddBin(5,X,Y);
      p[npl] = new TPolyLine(5,X,Y);
      p[npl]->SetUniqueID(npl);
      npl++;
    }
  }
  c1= new TCanvas("c1","",800,800);
  c1->Clear();
  c1->cd();
  /*
  for(int i=1; i<=5768; i++)
  {
    TVector3 pos = getPoint(i);
    double offset=15;
    double x = pos.X();
    double z = pos.Z();
    if ( (x<z+offset && x>z-offset) || (x<-z+offset && x>-z-offset)) poly->SetBinContent(i, 1000);
  }*/
  poly->SetBinContent(padID, 1000);
  poly->Draw("colz");
  pLine->SetLineColor(kGreen);
  pLine->Draw();
  cout << "theta : " << getTheta(padID) << endl;
  //getchar();
  
  TH2D* hist = new TH2D("hist", "GEM Status (Mar 23, 2019)", 100,-300,300,100,-300,300);
  TCanvas* c2= new TCanvas("c2","",800,800);
  c2->cd();
  hist->Draw();
  for(int i=0; i<1344; i++) {
    p[i]->SetLineColor(kSpring+7);
  }
  for(int i=1344; i<5768; i++) {
    p[i]->SetLineColor(kCyan);
  }
  for(int i=0; i<5768; i++) p[i]->Draw();
  pLine->SetLineColor(kRed);
  pLine->Draw();

  TCanvas* c3= new TCanvas("c3","",800,800);
  c3->cd();
  hist->Draw();
  
 

  const double fw = 6;
  const double ew_6div   = 41.3;
  const double ew0_6div  = 30.1;
  const double ew_20div  = 12.3;
  const double ew0_20div = 1.1;
  const double gap = 0.2;
  const double gw = 239.6;
  const double b = 2*fw+gw+94.2;

  double gem_6div_x[6][6]={{0}};
  double gem_6div_y[6][6]={{0}};

  double gem_20div_x[20][6]={{0}};
  double gem_20div_y[20][6]={{0}};

  double gem_mask1_x[35][6]={{0}};
  double gem_mask2_x[14][6]={{0}};
  double gem_mask3_x[21][6]={{0}};
  double gem_mask1_y[35][6]={{0}};
  double gem_mask2_y[14][6]={{0}};
  double gem_mask3_y[21][6]={{0}};

  gem_mask2_y[0][0] = (-43*3-gap/2-gap*3);
  gem_mask2_x[0][0] = gem_mask2_y[0][0]-fw*sqrt(2);
  gem_mask2_y[0][1] = -(gw+fw)/sqrt(2)+fw/sqrt(2);
  gem_mask2_x[0][1] = -(gw+fw)/sqrt(2)-fw/sqrt(2);
  gem_mask2_y[0][2] = gem_mask2_y[0][0];
  gem_mask2_x[0][2] = -gem_mask2_y[0][2]-(fw+gw)*sqrt(2);
  gem_mask2_y[0][3] = gem_mask2_y[0][0];
  gem_mask2_x[0][3] = gem_mask2_x[0][0];
 
  gem_mask2_y[1][0] = (-43*3-gap/2-gap*2);
  gem_mask2_x[1][0] = gem_mask2_y[1][0]-fw*sqrt(2);
  gem_mask2_y[1][1] = (-43*3-gap/2-gap*2);
  gem_mask2_x[1][1] = -gem_mask2_y[1][1]-(fw+gw)*sqrt(2);
  gem_mask2_y[1][2] = gem_mask2_y[0][1]+94.2/sqrt(2);
  gem_mask2_x[1][2] = -gem_mask2_y[1][2]-(fw+gw)*sqrt(2);
  gem_mask2_y[1][3] = (-43*2-gap/2-gap*2);
  gem_mask2_x[1][3] = gem_mask2_x[1][2];
  gem_mask2_y[1][4] = gem_mask2_y[1][3];
  gem_mask2_x[1][4] = gem_mask2_y[1][4]-fw*sqrt(2);
  gem_mask2_y[1][5] = gem_mask2_y[1][0];
  gem_mask2_x[1][5] = gem_mask2_x[1][0];
 
  gem_mask2_y[2][0] = (-43*2-gap/2-gap*1);
  gem_mask2_x[2][0] = gem_mask2_y[2][0]-fw*sqrt(2);
  gem_mask2_y[2][1] = gem_mask2_y[2][0];
  gem_mask2_x[2][1] = gem_mask2_x[1][2];
  gem_mask2_y[2][2] = (-43*1-gap/2-gap*1);
  gem_mask2_x[2][2] = gem_mask2_x[2][1];
  gem_mask2_y[2][3] = (-43*1-gap/2-gap*0);
  gem_mask2_x[2][3] = gem_mask2_y[2][3]-fw*sqrt(2);
  gem_mask2_y[2][4] = gem_mask2_y[2][0];
  gem_mask2_x[2][4] = gem_mask2_x[2][0];

  gem_mask3_y[0][0] = (-43*3-gap/2-gap*3);
  gem_mask3_x[0][0] = gem_mask2_y[0][0]-fw*sqrt(2);
  gem_mask3_y[0][1] = -(gw+fw)/sqrt(2)+fw/sqrt(2);
  gem_mask3_x[0][1] = -(gw+fw)/sqrt(2)-fw/sqrt(2);
  gem_mask3_y[0][2] = gem_mask2_y[0][0];
  gem_mask3_x[0][2] = -gem_mask2_y[0][2]-(fw+gw)*sqrt(2);
  gem_mask3_y[0][3] = gem_mask2_y[0][0];
  gem_mask3_x[0][3] = gem_mask2_x[0][0];
 
  gem_mask3_y[1][0] = (-43*3-gap/2-gap*2);
  gem_mask3_x[1][0] = gem_mask2_y[1][0]-fw*sqrt(2);
  gem_mask3_y[1][1] = (-43*3-gap/2-gap*2);
  gem_mask3_x[1][1] = -gem_mask2_y[1][1]-(fw+gw)*sqrt(2);
  gem_mask3_y[1][2] = gem_mask2_y[0][1]+94.2/sqrt(2);
  gem_mask3_x[1][2] = -gem_mask2_y[1][2]-(fw+gw)*sqrt(2);
  gem_mask3_y[1][3] = (-43*2-gap/2-gap*2);
  gem_mask3_x[1][3] = gem_mask2_x[1][2];
  gem_mask3_y[1][4] = gem_mask2_y[1][3];
  gem_mask3_x[1][4] = gem_mask2_y[1][4]-fw*sqrt(2);
  gem_mask3_y[1][5] = gem_mask2_y[1][0];
  gem_mask3_x[1][5] = gem_mask2_x[1][0];
 
  gem_mask3_y[2][0] = (-43*2-gap/2-gap*1);
  gem_mask3_x[2][0] = gem_mask2_y[2][0]-fw*sqrt(2);
  gem_mask3_y[2][1] = gem_mask2_y[2][0];
  gem_mask3_x[2][1] = gem_mask2_x[1][2];
  gem_mask3_y[2][2] = (-43*1-gap/2-gap*1);
  gem_mask3_x[2][2] = gem_mask2_x[2][1];
  gem_mask3_y[2][3] = (-43*1-gap/2-gap*0);
  gem_mask3_x[2][3] = gem_mask2_y[2][3]-fw*sqrt(2);
  gem_mask3_y[2][4] = gem_mask2_y[2][0];
  gem_mask3_x[2][4] = gem_mask2_x[2][0];

  for(int i=0; i<3; i++) 
    for(int j=0; j<6; j++) 
    {
      gem_mask2_x[7-i][j] = gem_mask2_x[i][j];
      gem_mask2_y[7-i][j] = -gem_mask2_y[i][j];

      gem_mask3_x[20-i][j] = gem_mask3_x[i][j];
      gem_mask3_y[20-i][j] = -gem_mask3_y[i][j];
    }
  
  double gem_sec1_mask_x[35][6]={{0}}; double gem_sec1_mask_y[35][6]={{0}};
  double gem_sec2_mask_x[14][6]={{0}}; double gem_sec2_mask_y[14][6]={{0}};
  double gem_sec3_mask_x[21][6]={{0}}; double gem_sec3_mask_y[21][6]={{0}};
  double gem_sec4_mask_x[14][6]={{0}}; double gem_sec4_mask_y[14][6]={{0}};


  double theta_sec1_mask = 0*TMath::Pi()/180.;
  double theta_sec2_mask = 90*TMath::Pi()/180.;
  double theta_sec3_mask = 180*TMath::Pi()/180.;
  double theta_sec4_mask = 270*TMath::Pi()/180.;
  
  for(int i=0; i<21; i++) 
    for(int j=0; j<6; j++) {
      gem_sec3_mask_x[i][j] = cos(theta_sec3_mask)*gem_mask3_x[i][j] 
			    - sin(theta_sec3_mask)*gem_mask3_y[i][j];
      gem_sec3_mask_y[i][j] = sin(theta_sec3_mask)*gem_mask3_x[i][j] 
	                    + cos(theta_sec3_mask)*gem_mask3_y[i][j];
    }
  for(int i=0; i<14; i++) 
    for(int j=0; j<6; j++) {
      gem_sec4_mask_x[i][j] = cos(theta_sec4_mask)*gem_mask2_x[i][j] 
			    - sin(theta_sec4_mask)*gem_mask2_y[i][j];
      gem_sec4_mask_y[i][j] = sin(theta_sec4_mask)*gem_mask2_x[i][j] 
	                    + cos(theta_sec4_mask)*gem_mask2_y[i][j];
    }
 
  TPolyLine* gem_sec1_mask[35];
  TPolyLine* gem_sec2_mask[14];
  TPolyLine* gem_sec3_mask[21];
  TPolyLine* gem_sec4_mask[14];

  int n_sec3_mask[21] = {4,6,5,5,5,
    		         5,5,5,5,5,
    		         5,5,5,5,5,
    		         5,5,5,5,6,
			 4};
  int n_sec4_mask[14] = {4,6,5,5,5,
    		         5,6,4,5,5,
			 5,5,5,5};
  
  for(int i=0; i<21; i++) 
    gem_sec3_mask[i] = new TPolyLine(n_sec3_mask[i], gem_sec3_mask_x[i], gem_sec3_mask_y[i]);
  for(int i=0; i<14; i++) 
    gem_sec4_mask[i] = new TPolyLine(n_sec4_mask[i], gem_sec4_mask_x[i], gem_sec4_mask_y[i]);
  
  //bool sec1_dead[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0};
  //bool sec2_dead[6]  = {1,0,1,1,1,1};
  //bool sec3_dead[6]  = {0,0,0,0,1,1};
  //bool sec4_dead[6]  = {1,1,0,1,1,1};
//  bool sec1_dead[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0};
//  bool sec2_dead[6]  = {0,0,1,1,1,0};
//  bool sec3_dead[6]  = {0,0,0,0,0,0};
//  bool sec4_dead[6]  = {0,0,0,0,0,0};
//
//  bool sec3_mask_dead[21] = {0,0,0,0,0,
//    			     0,0,0,0,0,
//			     0,0,0,0,0,
//			     0,0,0,0,0,
//			     0};
//  bool sec4_mask_dead[14] = {1,0,1,0,0,
//    			     0,0,0,0,0,
//			     0,0,0,0}; 
  
  bool sec1_dead[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0};
  bool sec2_dead[6]  = {1,0,1,1,1,1};
  bool sec3_dead[6]  = {0,0,0,0,1,1};
  bool sec4_dead[6]  = {1,1,1,1,1,1};

  bool sec3_mask_dead[21] = {1,0,0,0,0,
    			     0,0,0,0,0,
			     0,0,0,0,0,
			     0,0,0,0,0,
			     1};
  bool sec4_mask_dead[14] = {1,0,1,0,0,
    			     0,0,0,0,0,
			     0,0,0,0}; 
  double transparency = 0.3;
  
  for(int i=0; i<21; i++) {
    if(sec3_mask_dead[i]) 
    {
      gem_sec3_mask[i]->SetFillColorAlpha(kBlack, transparency);
      gem_sec3_mask[i]->Draw("f");
    }
    gem_sec3_mask[i]->Draw("");
  }
  for(int i=0; i<14; i++) {
    if(sec4_mask_dead[i]) 
    {
      gem_sec4_mask[i]->SetFillColorAlpha(kBlack, transparency);
      gem_sec4_mask[i]->Draw("f");
    }
    gem_sec4_mask[i]->Draw("");
  }


  for(int i=0; i<6; i++) 
  {
    gem_6div_y[i][0] = fw + (ew_6div+gap)*(5-i);
    gem_6div_y[i][1] = gem_6div_y[i][0];
    
    gem_6div_x[i][0] = -fw;

    if(i!=3) {
      gem_6div_y[i][2] = gem_6div_y[i][0]+ew_6div;
      gem_6div_y[i][3] = gem_6div_y[i][2];
      gem_6div_y[i][4] = gem_6div_y[i][0];
    
      if(i==0) { 
	gem_6div_y[i][2] = gem_6div_y[i][1]+ew0_6div;
	gem_6div_y[i][3] = gem_6div_y[i][2];
	gem_6div_y[i][4] = gem_6div_y[i][0];
      }
      
      gem_6div_x[i][3] = gem_6div_x[i][0];
      gem_6div_x[i][4] = gem_6div_x[i][0];
    } else {
      gem_6div_y[i][2] = fw + 94.2;
      gem_6div_y[i][3] = gem_6div_y[i][0]+ew_6div;
      gem_6div_y[i][4] = gem_6div_y[i][3];
      gem_6div_y[i][5] = gem_6div_y[i][0];
      
      gem_6div_x[i][3] = gem_6div_y[i][3]-b;
      gem_6div_x[i][4] = gem_6div_x[i][0];
      gem_6div_x[i][5] = gem_6div_x[i][0];
    }
    if(i<3) {
      gem_6div_x[i][1] = gem_6div_y[i][1]-b;
      gem_6div_x[i][2] = gem_6div_y[i][2]-b;
    }
    else{
      gem_6div_x[i][1] = -fw-gw;
      gem_6div_x[i][2] = gem_6div_x[i][1];
    }
  }
  for(int i=0; i<20; i++) 
  {
    gem_20div_y[i][0] = fw + (ew_20div+gap)*(19-i);
    gem_20div_y[i][1] = gem_20div_y[i][0];
    
    gem_20div_x[i][0] = -fw;

    if(i!=12) {
      gem_20div_y[i][2] = gem_20div_y[i][0]+ew_20div;
      gem_20div_y[i][3] = gem_20div_y[i][2];
      gem_20div_y[i][4] = gem_20div_y[i][0];
    
      if(i==0) { 
	gem_20div_y[i][2] = gem_20div_y[i][1]+ew0_20div;
	gem_20div_y[i][3] = gem_20div_y[i][2];
	gem_20div_y[i][4] = gem_20div_y[i][0];
      }
      
      gem_20div_x[i][3] = gem_20div_x[i][0];
      gem_20div_x[i][4] = gem_20div_x[i][0];
    } else {
      gem_20div_y[i][2] = fw + 94.2;
      gem_20div_y[i][3] = gem_20div_y[i][0]+ew_20div;
      gem_20div_y[i][4] = gem_20div_y[i][3];
      gem_20div_y[i][5] = gem_20div_y[i][0];
      
      gem_20div_x[i][3] = gem_20div_y[i][3]-b;
      gem_20div_x[i][4] = gem_20div_x[i][0];
      gem_20div_x[i][5] = gem_20div_x[i][0];
    }
    if(i<12) {
      gem_20div_x[i][1] = gem_20div_y[i][1]-b;
      gem_20div_x[i][2] = gem_20div_y[i][2]-b;
    }
    else{
      gem_20div_x[i][1] = -fw-gw;
      gem_20div_x[i][2] = gem_20div_x[i][1];
    }
  }
  double gem_sec1_x[20][6]; double gem_sec1_y[20][6];
  double gem_sec2_x[6][6]; double gem_sec2_y[6][6];
  double gem_sec3_x[6][6]; double gem_sec3_y[6][6];
  double gem_sec4_x[6][6]; double gem_sec4_y[6][6];


  double theta_sec1 = 45*TMath::Pi()/180.;
  double theta_sec2 = 135*TMath::Pi()/180.;
  double theta_sec3 = 225*TMath::Pi()/180.;
  double theta_sec4 = -45*TMath::Pi()/180.;
  
  for(int i=0; i<6; i++) 
    for(int j=0; j<6; j++) {
      gem_sec2_x[i][j] = cos(theta_sec2)*gem_6div_x[i][j] - sin(theta_sec2)*gem_6div_y[i][j];
      gem_sec2_y[i][j] = sin(theta_sec2)*gem_6div_x[i][j] + cos(theta_sec2)*gem_6div_y[i][j];
      
      gem_sec3_x[i][j] = cos(theta_sec3)*gem_6div_x[i][j] - sin(theta_sec3)*gem_6div_y[i][j];
      gem_sec3_y[i][j] = sin(theta_sec3)*gem_6div_x[i][j] + cos(theta_sec3)*gem_6div_y[i][j];
      
      gem_sec4_x[i][j] = cos(theta_sec4)*gem_6div_x[i][j] - sin(theta_sec4)*gem_6div_y[i][j];
      gem_sec4_y[i][j] = sin(theta_sec4)*gem_6div_x[i][j] + cos(theta_sec4)*gem_6div_y[i][j];
    }
  for(int i=0; i<20; i++) 
    for(int j=0; j<6; j++) {
      gem_sec1_x[i][j] = cos(theta_sec1)*gem_20div_x[i][j] - sin(theta_sec1)*gem_20div_y[i][j];
      gem_sec1_y[i][j] = sin(theta_sec1)*gem_20div_x[i][j] + cos(theta_sec1)*gem_20div_y[i][j];
    }
  TPolyLine* gem_sec1[20];
  TPolyLine* gem_sec2[6];
  TPolyLine* gem_sec3[6];
  TPolyLine* gem_sec4[6];
  for(int i=0; i<6; i++) {
    if(i!=3) {
      gem_sec2[i] = new TPolyLine(5, gem_sec2_x[i], gem_sec2_y[i]);
      gem_sec3[i] = new TPolyLine(5, gem_sec3_x[i], gem_sec3_y[i]);
      gem_sec4[i] = new TPolyLine(5, gem_sec4_x[i], gem_sec4_y[i]);
    }
    else {
      gem_sec2[i] = new TPolyLine(6, gem_sec2_x[i], gem_sec2_y[i]);
      gem_sec3[i] = new TPolyLine(6, gem_sec3_x[i], gem_sec3_y[i]);
      gem_sec4[i] = new TPolyLine(6, gem_sec4_x[i], gem_sec4_y[i]);
    }
  }
  for(int i=0; i<20; i++) {
    if(i!=12) {
      gem_sec1[i] = new TPolyLine(5, gem_sec1_x[i], gem_sec1_y[i]);
    }
    else {
      gem_sec1[i] = new TPolyLine(6, gem_sec1_x[i], gem_sec1_y[i]);
    }
  }
  
  double width = 2.0;
  for(int i=0; i<20; i++) {
    if(sec1_dead[i]) 
    {
      //gem_sec1[i]->SetLineColor(kRed);
      //gem_sec1[i]->SetLineWidth(width);
      gem_sec1[i]->SetFillColorAlpha(kBlack, transparency);
      gem_sec1[i]->Draw("f");
    }
    gem_sec1[i]->Draw("");
  }
  for(int i=0; i<6; i++) {
    if(sec2_dead[i]) 
    {
      gem_sec2[i]->SetFillColorAlpha(kBlack, transparency);
      //gem_sec2[i]->SetLineColor(kRed);
      //gem_sec2[i]->SetLineWidth(width);
      gem_sec2[i]->Draw("f");
    }
    if(sec3_dead[i]) 
    {
      gem_sec3[i]->SetFillColorAlpha(kBlack, transparency);
      //gem_sec3[i]->SetLineColor(kRed);
      //gem_sec3[i]->SetLineWidth(width);
      gem_sec3[i]->Draw("f");
    }
    if(sec4_dead[i]) 
    {
      gem_sec4[i]->SetFillColorAlpha(kBlack, transparency);
      //gem_sec4[i]->SetLineColor(kRed);
      //gem_sec4[i]->SetLineWidth(width);
      gem_sec4[i]->Draw("f");
    }
    gem_sec2[i]->Draw("");
    gem_sec3[i]->Draw("");
    gem_sec4[i]->Draw("");
  }

}
