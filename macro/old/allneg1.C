/* Aug 2019
 * Shinhyung */

// AsAd Swapped back to normal
//
#include "padHelper_new.hh"
#include "Unpacker.hh"

#define COBONUM 8
#define ASADNUM 4
#define PADNUM 5768
#define DEBUG 0

#define TBMIN 0
#define TBNUM 512
#define TBCUTMIN 100
#define TBCUTMAX 200

TFile *file[COBONUM];
TTree *tree;
Double_t* adc = new Double_t;
int channelmap[4][31][68];
long eventoffset;
TClonesArray *padArray = nullptr;
TH2Poly *poly;

int run_init=2;

void getPeak(int runNum, int eventNum = 0);

/*
  12 voltages : 1.0, ... , 6.0 V 
  30 AsAd, ZAP boards, 256 ch each
  peak[8][256][6]
*/
#define N_ASAD 30
#define N_CH   256
#define N_VOL  6
double peakid[N_ASAD][N_CH][N_VOL];
double peak[N_ASAD][N_CH][N_VOL];
double asad[N_ASAD][N_CH][N_VOL];

void allneg1()
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);
  
  for(int i=0; i<N_ASAD; i++)
    for(int j=0; j<N_CH; j++)
      for(int k=0; k<N_VOL; k++)
      {
	peak[i][j][k] = -9999;
	asad[i][j][k] = -9999;
      }

  TFile* outFile = new TFile("rootfile/allneg1.root","recreate");
  TH2D* hist = new TH2D("hist", ";Ch;ZAP#", 256,0,256, 31,1,32);
  TH2D* asad_hist = new TH2D("asad_hist", ";Ch;AsAd#", 256,0,256, 31,0,31);
  TGraph *graph[N_ASAD][N_CH];
  
  poly = new TH2Poly("TPC","TPC",-300,300,-300,300);
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
      poly->AddBin(5,X,Y);
    }
  }

  
  //int asadNum[N_ASAD] = {1, 2, 13, 9, 6, 5, 3, 7};
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

  double voltage[N_VOL];
  int runNum[N_VOL];
  for(int i=0; i<N_VOL; i++) {
    voltage[i] = 1.0 + 1.0*i;
    runNum[i]   = run_init+i;
  }

  // CH map ===================================================
  std::ifstream ifs("../HypTPCAna/param/channel_map_20180522.param");
  if(ifs.fail())
  {
    std::cerr<<" --> file open fail ..."<<std::endl;
    exit(-1);
  }
  while(ifs.good())
  {
    std::string buf;
    std::getline(ifs,buf);
    if(buf[0]=='#' || buf.empty())continue;
    std::istringstream is(buf);
    std::istream_iterator<std::string> issBegin(is);
    std::istream_iterator<std::string> issEnd;
    std::vector<std::string> param(issBegin,issEnd);
    if(param.empty()||param[0].empty())continue;
    int aget = atoi(param[0].c_str());
    int asad = atoi(param[1].c_str());
    int ch   = atoi(param[2].c_str());
    int layer  = atoi(param[3].c_str());
    int row  = atoi(param[4].c_str());
    int pid  = atoi(param[5].c_str());
    if(aget!=0)
    {
      channelmap[aget-1][asad-1][ch-1]=pid;
    }

  }
  ifs.close();
  //==================================================


  for(int k=0; k<N_VOL; k++) {
    cout << "k = " << k << endl;
    getPeak(runNum[k]);
  }

  
  // GRAPH ============================
  
  auto GetColor = [](int index)
  {
    Color_t colors[] = {kOrange, kTeal, kViolet, kSpring, kPink, kAzure};
    Color_t color = colors[index%6] + ((index/6)%10);
    return color;
  };

  for(int i=0; i<N_ASAD; i++)
    for(int j=0; j<N_CH; j++)
    {
      graph[i][j] = new TGraph(N_VOL, voltage, peak[i][j]);
      graph[i][j]->SetMarkerStyle(21);
      graph[i][j]->SetMarkerColor(GetColor(j));
      
      if( peak[i][j][N_VOL-1] == -9999 ) {
	hist->Fill(j, zapNum[i], 5);
	poly->SetBinContent(peakid[i][j][N_VOL-1], 5);
      }
      else if( peak[i][j][N_VOL-1] <700 &&  peak[i][j][N_VOL-1]>50 ) {
	hist->Fill(j, zapNum[i], 2);
	poly->SetBinContent(peakid[i][j][N_VOL-1], 2);
      }
      else if( peak[i][j][N_VOL-1] < 50) {
	hist->Fill(j, zapNum[i], 4);
	poly->SetBinContent(peakid[i][j][N_VOL-1], 4);
      }
      else if( peak[i][j][N_VOL-1] > 1700) {
	hist->Fill(j, zapNum[i], 3);
	poly->SetBinContent(peakid[i][j][N_VOL-1], 3);
      }
      else  {
	hist->Fill(j, zapNum[i], 1);
	poly->SetBinContent(peakid[i][j][N_VOL-1], 1);
      }

      if(asad[i][j][N_VOL-1]!=-9999) 
	asad_hist->Fill(j, i, peak[i][j][N_VOL-1]);
    }


  // DRAW ============================
  //TCanvas *c1[N_ASAD];
  TCanvas *c1[4];
  for(int k=0; k<4; k++) {
    c1[k] = new TCanvas(Form("c1_%d",k), Form("c1_%d",k), 1400,600);
    c1[k]->Divide(4,2);
    //for(int i=0; i<N_ASAD; i++) {
    for(int i=0; i<8; i++) {
      c1[k]->cd(i+1)->SetGrid();
      c1[k]->cd(i+1)->SetMargin(0.15, 0.1, 0.1, 0.1);
      for(int j=0; j<N_CH; j++) {
	if(i+k*8>=30) continue;
	if(j==0) {
	  graph[i+k*8][j]->SetTitle(Form("AsAd#%d + ZAP#%d; Input [V]; ADC Max [ch]", asadNum[i+k*8], zapNum[i+k*8]));
	  graph[i+k*8][j]->SetMinimum(0);
	  graph[i+k*8][j]->SetMaximum(2500);
	  graph[i+k*8][j]->Draw("APL");
	}
	else graph[i+k*8][j]->Draw("PL");
      }
    }
  }
  
  TCanvas *c2 = new TCanvas("c2", "zap", 1400,800);
  c2->cd();
  hist->Draw("colz");

  TCanvas *c3 = new TCanvas("c3", "asad", 1400,800);
  c3->cd();
  asad_hist->Draw("colz");

  TCanvas *c4 = new TCanvas("c4", "pad", 800,800);
  c4->cd();
  poly->SetMaximum(5);
  poly->Draw("colz");

  // SAVE ============================
  outFile->cd();
  hist->Write();
  poly->Write();
  asad_hist->Write();
  for(int i=0; i<N_ASAD; i++)
    for(int j=0; j<N_CH; j++)
    {
      graph[i][j]->SetName(Form("a%d_c%d",i,j));
      graph[i][j]->Write();
    }
}

void getPeak(int runNum, int eventNum = 1)
{
  bool flag = true;
  bool eventflag=true;

  int activepad=0;

  for(int coboNum = 0; coboNum<COBONUM ;coboNum++)
  {
    int asadflag[4]={0};
    std::ifstream fin;
    std::stringstream datafile;
    datafile.str("");
      datafile<<"/data"<<coboNum<<"/gemtest/cobo"<<coboNum<<"/run_"<<std::setfill('0')<<std::setw(4)<<runNum<<".dat";

    fin.open(datafile.str(), std::ios::in|std::ios::binary);

    if(!fin.is_open())
    {
      std::cerr<<"#Warning! :"<< datafile.str() << " does not exist !"<<std::endl;
      continue;
    }
    else
    {
      std::cout << "#D " << datafile.str() << " opened " <<std::endl;
    }
    fin.seekg(0, fin.end);
    unsigned long length = fin.tellg();
    fin.seekg(0, fin.beg);
    fin.seekg(12);
    if(!length) continue;
    //std::cout << " start pos : " <<  fin.tellg() << std::endl;
    //std::cout << " length : " <<  length << std::endl;
    //getchar();

    while(!(asadflag[0] &&  asadflag[1] &&  asadflag[2] &&  asadflag[3]) )
    {
      if (coboNum == 7) asadflag[2]++,asadflag[3]++;
      struct GetHeader_t get_header;
      int n_words;
      int padded;
      int asadId;
      while(!fin.eof() && (unsigned long)fin.tellg()<length)
      {
	// Read Header
	fin.read((char*)&get_header,sizeof(get_header));
#if DEBUG
	std::cout << "now : " << std::hex << (unsigned long)fin.tellg() << std::endl;
	std::cout << std::endl;
	//getchar();
#endif
	asadId = (int)get_header.m_AsadIdx;
	n_words = (int)get_header.m_nItems[0]*pow(16,6)
	  +(int)get_header.m_nItems[1]*pow(16,4)
	  +(int)get_header.m_nItems[2]*pow(16,2)
	  +(int)get_header.m_nItems[3]*pow(16,0);

	long eventNumber = (long)((int)get_header.m_EventIdx[0]*pow(16,6)
	  +(int)get_header.m_EventIdx[1]*pow(16,4)
	  +(int)get_header.m_EventIdx[2]*pow(16,2)
	  +(int)get_header.m_EventIdx[3]*pow(16,0));
	unsigned int frameSize = (int)get_header.m_FrameSize[0]*pow(16,4)
	  +(int)get_header.m_FrameSize[1]*pow(16,2)
	  +(int)get_header.m_FrameSize[2]*pow(16,0);

	int headerSize = (int)get_header.m_HeaderSize[0]*pow(16,2)
	  +(int)get_header.m_HeaderSize[1]*pow(16,0);
	int item_size = (int)get_header.m_ItemSize[0]*pow(16,2)
	  +(int)get_header.m_ItemSize[1]*pow(16,0);

	if(eventflag) {
	  eventoffset=eventNumber;
	  //cout << dec<< " eventOffset : " << eventoffset << endl;
	  //getchar();
	  eventflag=false;
	}

	padded = (frameSize-headerSize)*256 - item_size*n_words;
#if DEBUG
//#if 1
	std::cout << std::hex << std::setfill('0') << std::endl;
	std::cout << " ==========================================" << std::endl;
	std::cout << " frameSize  3B : " << std::setw(2) <<
	  (int)get_header.m_FrameSize[0] << " " << std::setw(2) <<
	  (int)get_header.m_FrameSize[1] << " " << std::setw(2) <<
	  (int)get_header.m_FrameSize[2] << std::endl;
	std::cout << " eventIdx   4B : " << std::setw(2) <<
	  (int)get_header.m_EventIdx[0] << " " << std::setw(2) <<
	  (int)get_header.m_EventIdx[1] << " " << std::setw(2) <<
	  (int)get_header.m_EventIdx[2] << " " << std::setw(2) <<
	  (int)get_header.m_EventIdx[3] << " --> " << std::dec << eventNumber << std::endl;
	std::cout << " coboIdx    1B : " << std::setw(2) << (int)get_header.m_CoboIdx << std::endl;
	std::cout << " asadIdx    1B : " << std::setw(2) << (int)get_header.m_AsadIdx << std::endl;
	std::cout << " frameType     : " << (int)get_header.m_FrameType[1] << std::endl;
	std::cout << " itemSize   : " << item_size << std::endl;
	std::cout << " nItems     : " << n_words << std::endl;
	std::cout << " padded     : " << padded << std::endl;
	std::cout << " ==========================================" << std::endl;
	getchar();
#endif


	if(eventNumber-eventoffset == eventNum)
	//if(eventNumber == eventNum)
	{
	  asadflag[asadId]++;
	  cout << "=== Reading CoboNum#" << coboNum << ", Asad#" << asadId << "..."<<endl;
	  //getchar();
	  break;
	}
	else {
	  //std::cout << "before jump : " << std::hex << (unsigned long)fin.tellg() << std::endl;
	  fin.seekg( (unsigned long)fin.tellg() + (unsigned long)(frameSize-headerSize)*256);
	  //if(fin.tellg() == length) {
	    //fin.seekg(12);
	  //}
	  //std::cout << "after jump : " << std::hex << (unsigned long)fin.tellg() << std::endl;
	  //std::cout << std::endl;

	}
      }

      // Read Data
      if( (int)get_header.m_FrameType[1]==2 ) // full
      {
	uint16_t data;
	int coboId = (int)get_header.m_CoboIdx;
	int ntb = n_words/(4*68);

	if(coboId!=coboNum)
	  std::cout<<"#Warning! CoBo Id = "<<coboId <<" and filename "<<coboNum<<" do not match ! " <<std::endl;
	int count[4] = {0,0,0,0};
	double fadc[4][68][512]={{{0}}}; //[aget][ch][tb]
	int hit=0;

	for(int i=0;i<n_words;++i)
	{
	  fin.read((char*)&data,sizeof(data));
	  int agetid   = (data>>6)&0x03;
	  int adc_low  = (data>>8)&0xff;
	  int adc_high = data&0x0f;
	  int adc      = (adc_high<<8)|adc_low;
	  int ch       = count[agetid]%68;
	  int t_bucket = (count[agetid]/68)%ntb;
	  fadc[agetid][ch][t_bucket]=adc;
	  ++count[agetid];
	}
	fin.seekg( (unsigned long)fin.tellg() + padded );
	int padid=0;
	for(int ageti = 0; ageti<4; ageti++)
	{
	  for(int chi = 0; chi<68; chi++)
	  {
	    int chi2;
	    if(chi>=0 && chi<=10) chi2=chi;
	    else if(chi>=12 && chi<=21) chi2 = chi-1;
	    else if(chi>=23 && chi<=44) chi2 = chi-2;
	    else if(chi>=46 && chi<=55) chi2 = chi-3;
	    else if(chi>=57 && chi<=67) chi2 = chi-4;
	    else continue;
	    if(chi==11 || chi==22 ||chi==45||chi==56) padid=-7;
	    else
	    {
	      padid=channelmap[ageti][coboNum*4+asadId][chi2];
	      //cout << "padid : " << padid << endl;
	    }
	    TVector3 pos = getPoint(padid);
	    double  x = pos.Z(), y = pos.X();
	    if(padid>0 && y>=x && y<-x)
	    {
	      double adcSum=0;
	      double adcMax=-999;
	      double ped = 0;
	      for(int tbi = 0; tbi<ntb; tbi++)
	      {
		if(tbi<100) ped+=fadc[ageti][chi][tbi];
		if(tbi >TBCUTMIN && tbi<TBCUTMAX &&adcMax<fadc[ageti][chi][tbi])
		{
		  adcMax=fadc[ageti][chi][tbi];
		}
	      }
	      peak[coboNum*4+asadId][ageti*64+chi2][runNum-run_init] = adcMax-ped/100.;
	      peakid[coboNum*4+asadId][ageti*64+chi2][runNum-run_init] = padid;
	    }
	  }
	}
      }
      else
      {
	std::cerr<<"#E : FrameType error !, FrameType="<<(int)get_header.m_FrameType[1] <<
	  " cobo : " << coboNum << " asad : " << asadId << std::endl;
	exit(-1);
      }

    }
    fin.close();
  }
  cout << "# of active pad : " << activepad << endl;
}
