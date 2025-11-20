/* June 2018
 * Shinhyung */

/* June 2025
 * Modified by Haein */

#ifdef __CINT__
#pragma link C++ enum EMessageTypes;
#pragma link C++ enum EMainCommandIdentifiers;
#endif

#include "padHelper_new.hh"
#include "Unpacker.hh"

#define COBONUM 8
#define ASADNUM 4
#define PADNUM 5768
#define DEBUG 0

#define TBMIN 0
#define TBNUM 170
#define TBCUTMIN 50
#define TBCUTMAX 140

#define TBPEDMIN 140
#define TBPEDMAX 160


bool UserStop=false;
bool endflag=false;

string Data_dir = "../run/";
string Param_dir = "../param/";
string Log_dir = "./debug/log/";

int Header_Size = 256;
int Data_Size = 139264*2;
double threshold = 85;
//double threshold = 0;



class EventDisplay {
 public:
  EventDisplay(const TGWindow *p, UInt_t w, UInt_t h, Int_t runNum = -1, Int_t eventNum = 1);
  virtual ~EventDisplay();
  void LoadEventData(long runNum, long eventNum);
  void GoForward();
  void GoBackward();
  void GoTo();
  void ShowPulse();
  void DoCanvasDraw();
  int Lastentry();

 private:
  TFile *file[COBONUM];
  TTree *tree;
  TGMainFrame *fMain;
  TRootEmbeddedCanvas *fECanvas;
  TGButton *pre_Button, *next_Button, *go_Button;
  TGTextEntry *fEvtHandler;
  TPaveText *fDef;


  TCanvas *fCanvas;
  TH2Poly *poly_time;

  TH1F* tbHist[PADNUM];
  Double_t* adc = new Double_t;
  Double_t padadc[4][4][68][512];
  Int_t fCurrentEvent=0;
  Int_t runNumber=0;
  int channelmap[4][31][68];
  bool asad_status[8][4];
  long eventoffset;
  Int_t fLastEntry = 1000000;



};
TH2Poly *poly_adc;
TH1D *hist_beamprofile;
TH2D *hist_xz;

EventDisplay* gDisplay = nullptr;
void ShowPulse() {
  if (gDisplay) gDisplay->ShowPulse();
}

TClonesArray *padArray = nullptr;

void closeFile(int sig)
{
  std::cout << "#D closeFile" << std::endl;
  UserStop=true;
  gApplication->Terminate(0);
}


EventDisplay::EventDisplay(const TGWindow* p, UInt_t w, UInt_t h, Int_t runNum, Int_t eventNum)
{

  fMain = new TGMainFrame(p,w,h);
  fCurrentEvent= eventNum;

  std::ifstream ifs(Form("%s/channel_map_20180522.param",Param_dir.c_str()));

  if(ifs.fail())
    {
      std::cerr<<" --> channel_map file open fail ..."<<std::endl;
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



  for(int i=0; i<PADNUM; i++)
    {
      tbHist[i]=new TH1F(Form("tbHist%d",i),"",TBNUM,TBMIN,TBNUM-TBMIN);
    }


  signal(SIGINT, closeFile);
  runNumber = runNum;

#if 0
  // ------------
  std::stringstream eventdatafile;
  eventdatafile.str("");
  eventdatafile<<Log_dir<<"eventnum_run"<<std::setfill('0')<<std::setw(4)<<runNum<<".log";
  std::ifstream fevtall(eventdatafile.str());

  if(fevtall.fail())
    {
      std::cerr<<" --> file open fail ..." << eventdatafile.str() <<std::endl;
      sleep(1);
      continue;
    }

  long eventNum =999999999;
  while(fevtall.good())
    {
      std::string buf;
      std::getline(fevtall,buf);
      if(buf[0]=='#' || buf.empty())continue;
      std::istringstream is(buf);
      std::istream_iterator<std::string> issBegin(is);
      std::istream_iterator<std::string> issEnd;
      std::vector<std::string> param(issBegin,issEnd);
      if(param.empty()||param[0].empty())continue;
      int evtmin = atoi(param[0].c_str());
      int evtmax = atoi(param[1].c_str());
      int evtnum = atoi(param[2].c_str());
      if(evtnum < eventNum && evtnum >0) eventNum=evtnum;
    }
  fevtall.close();
  // ------------
  //
  if( eventNum == 999999999 )
    {
      cout << "runNum# " << runNum << " no event.. wait " << endl;
      sleep(1);
      continue;
    }

  while(eventNum < 0)
    {
      std::cerr << " Check eventNum (" << eventNum << ")"<< std::endl;
      //sleep(1);
      continue;
    }
#endif
  std::cout << " * Reading * runNum : " << runNum << " eventNum : " << eventNum << std::endl;
}

EventDisplay::~EventDisplay() {
  fMain->Cleanup();
  delete fMain;
}



void EventDisplay::LoadEventData(long runNum, long eventNum)
{
  bool flag = true;
  bool eventflag=true;

  int min_lastentry = -1;

  for(int coboNum = 0; coboNum<COBONUM ;coboNum++)
    {
      std::ifstream fin;
      std::stringstream datafile;
      datafile.str("");
      datafile<<Data_dir<<"cobo"<<coboNum<<"/run_"<<std::setfill('0')<<std::setw(4)<<runNum<<".dat";

      fin.open(datafile.str(), std::ios::in|std::ios::binary);

      if(!fin.is_open())
	{
	  std::cerr<<"#Warning! :"<< datafile.str() << " does not exist !"<<std::endl;
	  continue;
	}
      else
	{
	  //std::cout << "#D " << datafile.str() << " opened " <<std::endl;
	}

      fin.seekg(0, fin.end);
      unsigned long length = fin.tellg();
      fin.seekg(0, fin.beg);
      fin.seekg(12);

      if(!length) continue;
      min_lastentry = (length - 12) / (Header_Size + Data_Size) / ((coboNum == 7)? ASADNUM -1 : ASADNUM)-1;
      if(min_lastentry < fLastEntry)fLastEntry = min_lastentry;

      for(int asadNum = 0; asadNum<((coboNum == 7)? ASADNUM -1 : ASADNUM) ;asadNum++)
	// for(int asadNum = 1; asadNum<((coboNum == 3)? 2 : 1) ;asadNum++)
	{
	  int asadId;
	  struct GetHeader_t get_header;
	  int n_words;
	  int padded;
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
		eventflag=false;
	      }

	      padded = (frameSize-headerSize)*256 - item_size*n_words;

#if DEBUG
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
	      //getchar();
#endif


	      if(eventNumber-eventoffset == eventNum)
		{
		  //std::cout<<"check"<<std::endl;
		  //cout << "=== Reading CoboNum#" << coboNum << ", Asad#" << asadId << "..."<<endl;
		  //getchar();
		  break;
		}
	      else {
		fin.seekg( (unsigned long)fin.tellg() + (unsigned long)(frameSize-headerSize)*256);
	      }
	    }

	  // Read Data

	  if( (int)get_header.m_FrameType[1]==1 ) //partial
	    {


	      uint32_t data;
	      int asadId = (int)get_header.m_AsadIdx;
	      int coboId = (int)get_header.m_CoboIdx;
	      int ntb = n_words/(4*68);
	      if(coboId!=coboNum)
		std::cout<<"#Warning! CoBo Id = "<<coboId <<" and filename "<<coboNum<<" do not match ! " <<std::endl;
	      double fadc[4][68][TBNUM]={{{0}}}; //[aget][ch][tb]


	      for(int i=0;i<n_words;++i)
		{
		  fin.read((char*)&data,sizeof(data));
		  int agetid  	  = (data>>6)&0x03;
		  int adc_low 	  = (data>>24)&0xff;
		  int adc_high 	  = (data>>16)&0x0f;
		  int adc           = (adc_high<<8)|adc_low;
		  int ch_low 	  = (data>>15)&0x01;
		  int ch_high 	  = data&0x3f;
		  int ch       	  = (ch_high<<1)|ch_low;
		  int t_bucket_low  = (data>>22)&0x03;
		  int t_bucket_high = (data>>8)&0x7f;
		  int t_bucket      = (t_bucket_high<<2)|t_bucket_low;

		  fadc[agetid][ch][t_bucket]=adc;
#if DEBUG
		  std::cout << std::dec << " asad : " << coboNum*4+asadId << std::endl;
		  std::cout << " aget : " << agetid << std::endl;
		  std::cout << " ch : " << ch << std::endl;
		  std::cout << " t_bucket : " << t_bucket << std::endl;
		  std::cout << " adc : " << adc << std::hex<<  std::endl;
		  //getchar();
#endif

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
			}
		      if(padid>0)
			{
			  double adcSum=0;
			  int tbmax=0;
			  double adcMax=-999;
			  std::vector<double> v;
			  for(int tbi = TBMIN; tbi<TBNUM+TBMIN; tbi++)
			    {
			      tbHist[padid-1]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]);
			      if(tbi > TBPEDMIN && tbi < TBPEDMAX)
				adcSum+=fadc[ageti][chi][tbi];
			      v.push_back(fadc[ageti][chi][tbi]);
			      if(adcMax<fadc[ageti][chi][tbi])
				{
				  adcMax=fadc[ageti][chi][tbi];
				  tbmax=tbi;
				}
			    }
			  double stddev = TMath::StdDev(v.size(), v.data());
			  // std::cout << padid << "\t" << stddev << std::endl;
			  // if(stddev == 0 || stddev > 50) continue;

			  double adc_height = adcMax - ((double)adcSum/(TBPEDMAX - TBPEDMIN));

			  //poly_adc->SetBinContent(padid,adc_height);
			  if(adc_height > threshold){
			    int maxBin = tbHist[padid-1]->GetMaximumBin();
			    float maxTime = tbHist[padid-1]->GetBinCenter(maxBin);
			    // if(maxTime > 120 || maxTime < 80) maxTime = 0;
			    //poly_time->SetBinContent(padid,maxTime);
			    TVector3 pad_pos =getPosition(padid);
			    //if(pad_pos.Z()<-143.-40.)
			    if(pad_pos.Z() >-300. && pad_pos.Z() < -200.){
			      hist_beamprofile->Fill(pad_pos.X());

			    }
			    //std::cout<<poly_adc->GetBinContent(padid)<<std::endl;
			    poly_adc->SetBinContent(padid, poly_adc->GetBinContent(padid) + 1);


			    hist_xz->Fill(pad_pos.Z(),pad_pos.X());
			      //hist_beamprofile->Fill(pad_pos.Z());



			  }
			  else{
			    //poly_time->SetBinContent(padid,0);
			  }
			}
		    }
		}
	    }
	  else
	    {
	      std::cerr<<"#E : FrameType error !, FrameType="<<(int)get_header.m_FrameType[1] <<
		" cobo : " << coboNum << " asad : " << asadId << std::endl;
	      //exit(-1);
	      continue;
	    }

	}
      fin.close();
    }

#if 0
  for(int i=0; i<padArray->GetEntries(); i++)
    {
      auto pad = (S2Pad*)padArray->At(i);
      if(pad->GetPadID()==-7) continue;
      //cout << " cobo " << coboNum <<", asad " << asadNum << ", padID : " << pad->GetPadID() << endl;
      //getchar();
      adc = pad->GetADC();
      //if(pad->GetPadID()==3814) cout << " adc[4]" << adc[4] << endl;
      double adcSum=0;
      int tbmax=0;
      double adcMax=-999;
      //for(int tbi = 80; tbi<120; tbi++)
      //cout <<  " padid : " << pad->GetPadID() << endl;
      //getchar();
      for(int tbi = 0; tbi<512; tbi++)
	//for(int tbi = 50; tbi<480; tbi++)
	{
	  //cout << " adc " << adc[tbi] << endl;
	  tbHist[pad->GetPadID()-1]->SetBinContent(tbi+1,adc[tbi]);
	  //cout << "tbi="<< tbi << " adc " << adc[tbi] << endl;
	  //sleep(1);
	  adcSum+=adc[tbi];
	  if(adcMax<adc[tbi])
	    {
	      adcMax=adc[tbi];
	      tbmax=tbi;
	    }
	}
      TVector3 pos = getPoint(pad->GetPadID());
      pos.SetY(tbmax*40*0.001*53-225);
      double y = pos.Y();
      double z = pos.Z();
      if(pad->GetPadID()>0 && adcMax>0)
	{
	  poly_adc->SetBinContent(pad->GetPadID(),adcMax-(adcSum/512.));
	  flag = false;
	}
    }
  padArray->Clear();
#endif

}

void EventDisplay::ShowPulse()
{
  fCanvas->GetPad(1)->cd(1);
  TObject* select = gPad->GetSelected();

  fCanvas->GetPad(2)->cd(1);
  TObject* select_time = gPad->GetSelected();



  if (!select && !select_time){
    return;
  }

  if (!select->InheritsFrom(TH2Poly::Class()) && !select_time->InheritsFrom(TH2Poly::Class())) {
    gPad->SetUniqueID(0);
    return;
  }

  int pyold,px,py, binIdx;
  Float_t upx,upy;
  TH2Poly *pp;
  TH2Poly *pp_time;

  if(select){
    pp = (TH2Poly*)select;
    gPad->GetCanvas()->FeedbackMode(kTRUE);
    pyold = gPad->GetUniqueID();
    px = gPad->GetEventX();
    py = gPad->GetEventY();
    upx = gPad->AbsPixeltoX(px);
    upy = gPad->AbsPixeltoY(py);
    binIdx = pp->FindBin(upx, upy) -1;
  }

  else if(select_time){
    pp_time = (TH2Poly*)select_time;
    gPad->GetCanvas()->FeedbackMode(kTRUE);
    pyold = gPad->GetUniqueID();
    px = gPad->GetEventX();
    py = gPad->GetEventY();
    upx = gPad->AbsPixeltoX(px);
    upy = gPad->AbsPixeltoY(py);
    binIdx = pp_time->FindBin(upx, upy) -1;
  }

  fCanvas->cd(3);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);

  if (binIdx < 0 || binIdx >= PADNUM) return;

  tbHist[binIdx]->SetMaximum(4000);
  tbHist[binIdx]->SetTitle(Form("Pad : %d;Time Bucket;ADC [ch]",binIdx));
  gStyle->SetOptStat("eMR");
  tbHist[binIdx]->Draw("hist");
  fCanvas->Update();
}

void EventDisplay::GoBackward(){
  fCanvas->Clear("D");

  if (fCurrentEvent == 0)
    {
      fCurrentEvent = 0;
      std::cout << "First Event" << std::endl;
    }
  else
    {
      fCurrentEvent--;
    }
  LoadEventData(runNumber, fCurrentEvent);
  DoCanvasDraw();
}

int EventDisplay::Lastentry(){
  return fLastEntry;
}

void EventDisplay::GoForward(){
  fCanvas->Clear("D");

  if(fCurrentEvent == fLastEntry)
    {
      std::cout<< "Last Event" <<std::endl;
    }
  else
    {
      fCurrentEvent++;
    }

  LoadEventData(runNumber, fCurrentEvent);
  DoCanvasDraw();
}

void EventDisplay::GoTo(){
  fCanvas->Clear("D");
  TString eventNumberStr = fEvtHandler ->GetMarkedText();
  Int_t evnum = eventNumberStr.Atoi();

  if(evnum < 0 || evnum > fLastEntry)
    {
      cout<<"No Event" <<endl;
    }

  else
    {
      fCurrentEvent = evnum;
      LoadEventData(runNumber, fCurrentEvent);
      DoCanvasDraw();
    }

}

void EventDisplay::DoCanvasDraw(){

  fCanvas->cd(1);
  fCanvas->GetPad(1)->SetLogz();
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  poly_adc->SetMinimum(0);
  poly_adc->SetStats(false);
  poly_adc->SetTitle(Form("Pulse Height | Run%d Event : %d;Z [mm];X [mm]",runNumber,fCurrentEvent));
  poly_adc->GetZaxis()->SetRangeUser(0, 1000);
  poly_adc->Draw("colz");

  fCanvas->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  poly_time->SetMinimum(0);
  poly_time->SetStats(false);
  poly_time->SetTitle(Form("Pulse Time | Run%d Event : %d;Z [mm];X [mm]",runNumber,fCurrentEvent));
  poly_time->Draw("colz");

  //fCanvas->cd(1);
  fCanvas->Update();
  fCanvas->GetPad(1)->AddExec("dynamic", "ShowPulse()");
  fCanvas->GetPad(2)->AddExec("dynamic", "ShowPulse()");

}

/*
void online_upstream(Int_t runNum = -1, Int_t eventNum = 1)
{
  gStyle->SetOptStat();
  int argc = 0;
  char** argv = nullptr;
  TApplication theApp("App", &argc, argv);

  gDisplay = new EventDisplay(gClient->GetRoot(), 800, 500, runNum, eventNum);
  theApp.Run();
}
*/

void online_upstream(Int_t runNum = -1){

  hist_beamprofile = new TH1D("hist_beamprofile","hist_beamprofile",300,-300,300);
  poly_adc = new TH2Poly("TPCadc","TPCadc",-300,300,-300,300);
  hist_xz = new TH2D("hist_xz","hist_xz",300,-300,300,300,-300,300);


TPolyLine* pl[PADNUM];

  double l = (586./2.)/(1+sqrt(2.));
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
	  poly_adc->AddBin(5,X,Y);
	  //poly_time->AddBin(5,X,Y);
	  pl[npl] = new TPolyLine(5,X,Y);
	  pl[npl]->SetUniqueID(npl);
	  npl++;
	}
    }

  const TGWindow *p = gClient->GetRoot();

  UInt_t w = 1200;
  UInt_t h = 800;
  EventDisplay *ev = new EventDisplay(p, w, h,runNum);

  for(int n=0;n<500;n++){
  //for(int n=0;n<ev->Lastentry();n++){
    cout<<n<<endl;
    //ev->LoadEventData(runNum, n);
    ev->LoadEventData(runNum, n);
  }

  auto c1 =  new TCanvas("c1","c1");
  c1->Divide(2);
  c1->cd(1);
  hist_beamprofile->Draw();
  c1->cd(2);
  //hist_xz->Draw("colz");
  poly_adc->Draw("colz");

}
