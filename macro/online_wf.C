#ifdef __CINT__
#pragma link C++ enum EMessageTypes;
#pragma link C++ enum EMainCommandIdentifiers;
#endif

#include "padHelper_new.hh"
#include "Unpacker.hh"

#define COBONUM 8
#define ASADNUM 4
#define AGETNUM 4
#define TOTALASADNUM 31
#define PADNUM 5768
#define DEBUG 0

#define TBMIN 0
#define TBNUM 170
#define TBCUTMIN 50
#define TBCUTMAX 170


bool UserStop=false;
bool endflag=false;

string Data_dir = "../run/";
string Param_dir = "../param/";
string Log_dir = "./debug/log/";

int Header_Size = 256;
int Data_Size = 139264*2;
double threshold = 50;
double evt_offset = 0;

class EventDisplay {
 public:
  EventDisplay(const TGWindow *p, UInt_t w, UInt_t h, Int_t runNum = -1, Int_t eventNum = 1);
  virtual ~EventDisplay();
  uint32_t FindLastEventId(std::ifstream& fin, bool verbose=false);
  void LoadEventData(long runNum, long eventNum);
  void GoForward();
  void GoBackward();
  void GoLast();
  void GoTo();
  void DoCanvasDraw();

 private:
  TFile *file[COBONUM];
  TTree *tree;
  TGMainFrame *fMain;
  TRootEmbeddedCanvas *fECanvas;
  TGButton *pre_Button, *next_Button, *go_Button, *last_Button;
  TGTextEntry *fEvtHandler;
  TPaveText *fDef;


  TCanvas *fCanvas;
  TH1F* tbHist[PADNUM];
  TH1F* tbHist_aget[TOTALASADNUM][AGETNUM];
  Double_t* adc = new Double_t;
  Double_t padadc[4][4][68][512];
  Int_t fCurrentEvent=0;
  Int_t runNumber=0;
  int channelmap[4][31][68];
  bool asad_status[8][4];
  long eventoffset;
  Int_t fLastEntry = -1;



};

EventDisplay* gDisplay = nullptr;


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
  for(int i=0;i<TOTALASADNUM;i++){
    for(int j=0;j<AGETNUM;j++){
      tbHist_aget[i][j] = new TH1F(Form("tbHist_asad%daget_%d",i,j),"",TBNUM,TBMIN,TBNUM-TBMIN);
      tbHist_aget[i][j]->SetLineColor(j+1);
      tbHist_aget[i][j]->SetTitle(Form("AsAd%d;Time Bucket;ADC [ch]",i+1));

    }
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
      sleep(1);
      continue;
    }
#endif


  std::cout << " * Reading * runNum : " << runNum << " eventNum : " << eventNum << std::endl;

  fECanvas = new TRootEmbeddedCanvas("ECanvas",fMain,1500,1000);
  fMain->AddFrame(fECanvas,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));

  fCanvas = fECanvas->GetCanvas();
  fCanvas->Divide(8,4,0.001,0.001);

  LoadEventData(runNum, eventNum);
  DoCanvasDraw();


  TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);

  pre_Button = new TGTextButton(hframe,"&Previous",1);
  pre_Button->Connect("Clicked()","EventDisplay",this,"GoBackward()");
  hframe->AddFrame(pre_Button, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  next_Button = new TGTextButton(hframe, "&Next",1);
  next_Button->Connect("Clicked()","EventDisplay",this,"GoForward()");
  hframe->AddFrame(next_Button, new TGLayoutHints(kLHintsCenterX ,5,5,3,4));

  last_Button = new TGTextButton(hframe, "&Last",1);
  last_Button->Connect("Clicked()","EventDisplay",this,"GoLast()");
  hframe->AddFrame(last_Button, new TGLayoutHints(kLHintsCenterX ,5,5,3,4));


  fEvtHandler = new TGTextEntry(hframe);
  fEvtHandler ->SetEnabled(true);
  fEvtHandler->SetInsertMode();
  hframe->AddFrame(fEvtHandler, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));

  go_Button = new TGTextButton(hframe, "&Go");
  go_Button->Connect("Clicked()", "EventDisplay", this, "GoTo()");
  hframe->AddFrame(go_Button, new TGLayoutHints(kLHintsCenterX, 5,5,3,4));



  TGTextButton *exit = new TGTextButton(hframe,"&Exit","gApplication->Terminate(0)");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
  fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));


  fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));
  fMain->SetWindowName(Form("run0%d",runNumber));
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->MapWindow();

  std::cout << " * Finished * runNum : " << runNum << " eventNum : " << eventNum << std::endl;

}

EventDisplay::~EventDisplay() {
  fMain->Cleanup();
  delete fMain;
}

struct SkipInfo {
  uint32_t eventId = 0;
  uint32_t dataBytes = 0;
  bool ok = false;
};

static SkipInfo calcSkipFromHeader(const GetHeader_t& h) {
  SkipInfo s;

  const int asadId = (int)h.m_AsadIdx;
  (void)asadId;

  const uint32_t n_words = (uint32_t)h.m_nItems[0] * (uint32_t)std::pow(16, 6)
    + (uint32_t)h.m_nItems[1] * (uint32_t)std::pow(16, 4)
    + (uint32_t)h.m_nItems[2] * (uint32_t)std::pow(16, 2)
    + (uint32_t)h.m_nItems[3] * (uint32_t)std::pow(16, 0);

  const uint32_t eventNumber = (uint32_t)h.m_EventIdx[0] * (uint32_t)std::pow(16, 6)
    + (uint32_t)h.m_EventIdx[1] * (uint32_t)std::pow(16, 4)
    + (uint32_t)h.m_EventIdx[2] * (uint32_t)std::pow(16, 2)
    + (uint32_t)h.m_EventIdx[3] * (uint32_t)std::pow(16, 0);

  const uint32_t frameSize  = (uint32_t)h.m_FrameSize[0] * (uint32_t)std::pow(16, 4)
    + (uint32_t)h.m_FrameSize[1] * (uint32_t)std::pow(16, 2)
    + (uint32_t)h.m_FrameSize[2] * (uint32_t)std::pow(16, 0);

  const uint32_t headerSize = (uint32_t)h.m_HeaderSize[0] * (uint32_t)std::pow(16, 2)
    + (uint32_t)h.m_HeaderSize[1] * (uint32_t)std::pow(16, 0);

  const uint32_t item_size  = (uint32_t)h.m_ItemSize[0] * (uint32_t)std::pow(16, 2)
    + (uint32_t)h.m_ItemSize[1] * (uint32_t)std::pow(16, 0);

  const int64_t padded = (int64_t)( (int64_t)(frameSize) - (int64_t)(headerSize) ) * 256LL
    - (int64_t)item_size * (int64_t)n_words;

  if (frameSize == 0 || item_size == 0) return s;
  if (padded < 0) return s;

  const uint64_t dataBytes = (uint64_t)item_size * (uint64_t)n_words + (uint64_t)padded;

  if (dataBytes > (1ULL<<30)) return s;

  s.eventId = eventNumber;
  s.dataBytes = (uint32_t)dataBytes;
  s.ok = true;
  return s;
}

uint32_t EventDisplay::FindLastEventId(std::ifstream& fin, bool verbose=false) {
  fin.clear();
  fin.seekg(0, std::ios::end);
  const std::streamoff fileSize = fin.tellg();
  fin.clear();
  fin.seekg(12, std::ios::beg);

  GetHeader_t h;
  fin.read(reinterpret_cast<char*>(&h), sizeof(h));

  const uint32_t firsteventNumber = (uint32_t)h.m_EventIdx[0] * (uint32_t)std::pow(16, 6)
    + (uint32_t)h.m_EventIdx[1] * (uint32_t)std::pow(16, 4)
    + (uint32_t)h.m_EventIdx[2] * (uint32_t)std::pow(16, 2)
    + (uint32_t)h.m_EventIdx[3] * (uint32_t)std::pow(16, 0);


  fin.seekg(12, std::ios::beg);

  uint32_t lastEvent = 0;
  size_t nrec = 0;

  while (true) {
    const std::streamoff pos = fin.tellg();
    if (pos < 0) break;
    if (fileSize - pos < (std::streamoff)sizeof(GetHeader_t)) break;

    GetHeader_t h;
    fin.read(reinterpret_cast<char*>(&h), sizeof(h));

    if (!fin) break;

    const auto s = calcSkipFromHeader(h);
    if (!s.ok) {
      if (verbose) {
	std::cerr << "[findLastEventId] header sanity check failed at pos="
                  << (long long)pos << " (0x" << std::hex << (long long)pos << std::dec << ")\n";
      }
      break;
    }

    lastEvent = s.eventId - firsteventNumber;
    evt_offset = firsteventNumber;
    ++nrec;

    fin.seekg((std::streamoff)s.dataBytes, std::ios::cur);
    if (!fin) break;
  }

  if (verbose) {
    std::cerr << "[findLastEventId] scanned records=" << nrec
              << ", lastEventId=" << lastEvent << "\n";
  }

  return lastEvent;
}

bool ReadHeader(std::ifstream& fin, GetHeader_t& h) {
  fin.read(reinterpret_cast<char*>(&h), sizeof(h));
  return (bool)fin;
}

unsigned long GetEventId(const GetHeader_t& h) {
  return (unsigned long)(
			 (int)h.m_EventIdx[0]*pow(16,6) +
			 (int)h.m_EventIdx[1]*pow(16,4) +
			 (int)h.m_EventIdx[2]*pow(16,2) +
			 (int)h.m_EventIdx[3]*pow(16,0)
			 );
}

unsigned long CalcDataBytes(const GetHeader_t& h) {
  int n_words =
    (int)h.m_nItems[0]*pow(16,6) +
    (int)h.m_nItems[1]*pow(16,4) +
    (int)h.m_nItems[2]*pow(16,2) +
    (int)h.m_nItems[3]*pow(16,0);

  unsigned int frameSize =
    (int)h.m_FrameSize[0]*pow(16,4) +
    (int)h.m_FrameSize[1]*pow(16,2) +
    (int)h.m_FrameSize[2]*pow(16,0);

  int headerSize =
    (int)h.m_HeaderSize[0]*pow(16,2) +
    (int)h.m_HeaderSize[1]*pow(16,0);

  int item_size =
    (int)h.m_ItemSize[0]*pow(16,2) +
    (int)h.m_ItemSize[1]*pow(16,0);

  long padded = (long)(frameSize - headerSize) * 256L - (long)item_size * (long)n_words;
  if (padded < 0) return 0;

  return (unsigned long)item_size * (unsigned long)n_words + (unsigned long)padded;
}

void EventDisplay::LoadEventData(long runNum, long eventNum)
{
  bool flag = true;
  for(int i=0;i<TOTALASADNUM;i++){
    for(int j=0;j<4;j++){
      tbHist_aget[i][j]->SetMaximum(0);
    }
  }


  for (int coboNum = 0; coboNum < COBONUM; coboNum++) {

    std::ifstream fin;
    std::stringstream datafile;
    datafile.str("");
    datafile << Data_dir << "cobo" << coboNum
             << "/run_" << std::setfill('0') << std::setw(4) << runNum << ".dat";

    fin.open(datafile.str(), std::ios::in | std::ios::binary);
    if (!fin.is_open()) {
      std::cerr << "#Warning! :" << datafile.str() << " does not exist !" << std::endl;
      continue;
    }

    // file length
    fin.seekg(0, fin.end);
    unsigned long length = (unsigned long)fin.tellg();
    if (!length) { fin.close(); continue; }
    unsigned int lastEventId = FindLastEventId(fin, false);
    fLastEntry = lastEventId;

    fin.clear();
    fin.seekg(12, fin.beg);

    std::cout << "#D " << datafile.str() << " opened " << std::endl;


    GetHeader_t first_header;
    fin.read(reinterpret_cast<char*>(&first_header), sizeof(first_header));
    if (!fin) { fin.close(); continue; }

    long firstEvent =
      (long)((int)first_header.m_EventIdx[0] * pow(16,6) +
             (int)first_header.m_EventIdx[1] * pow(16,4) +
             (int)first_header.m_EventIdx[2] * pow(16,2) +
             (int)first_header.m_EventIdx[3] * pow(16,0));

    long eventoffset = firstEvent;
    long targetEventAbs = eventoffset + eventNum;


    fin.clear();
    fin.seekg(12, fin.beg);


    bool inTarget = false;

    while (!fin.eof() && (unsigned long)fin.tellg() < length) {

      GetHeader_t get_header;
      fin.read(reinterpret_cast<char*>(&get_header), sizeof(get_header));
      if (!fin) break;


      int asadId = (int)get_header.m_AsadIdx;

      int n_words =
        (int)get_header.m_nItems[0] * pow(16,6) +
        (int)get_header.m_nItems[1] * pow(16,4) +
        (int)get_header.m_nItems[2] * pow(16,2) +
        (int)get_header.m_nItems[3] * pow(16,0);

      long eventNumber =
        (long)((int)get_header.m_EventIdx[0] * pow(16,6) +
               (int)get_header.m_EventIdx[1] * pow(16,4) +
               (int)get_header.m_EventIdx[2] * pow(16,2) +
               (int)get_header.m_EventIdx[3] * pow(16,0));

      unsigned int frameSize =
        (int)get_header.m_FrameSize[0] * pow(16,4) +
        (int)get_header.m_FrameSize[1] * pow(16,2) +
        (int)get_header.m_FrameSize[2] * pow(16,0);

      int headerSize =
        (int)get_header.m_HeaderSize[0] * pow(16,2) +
        (int)get_header.m_HeaderSize[1] * pow(16,0);

      int item_size =
        (int)get_header.m_ItemSize[0] * pow(16,2) +
        (int)get_header.m_ItemSize[1] * pow(16,0);

      long padded = (long)(frameSize - headerSize) * 256L - (long)item_size * (long)n_words;


      if (padded < 0) {
	std::cerr << "[ERROR] padded < 0 (sync broken?) cobo=" << coboNum
                  << " pos=" << (unsigned long)fin.tellg() << std::endl;
        break;
      }
      unsigned long dataBytes =
        (unsigned long)item_size * (unsigned long)n_words + (unsigned long)padded;


      if (!inTarget) {
        if (eventNumber < targetEventAbs) {

          fin.seekg((std::streamoff)dataBytes, std::ios::cur);
          continue;
        }
        if (eventNumber == targetEventAbs) {
          inTarget = true;
	  std::cout << "=== Reading CoboNum#" << coboNum << " targetEvent=" << eventNum
                    << " (abs " << targetEventAbs << ") ..." << std::endl;

        } else {
          break;
        }
      } else {

        if (eventNumber != targetEventAbs) break;
      }


      std::cout << "asad : " << asadId << ", eventid : " << (eventNumber - eventoffset) << std::endl;
      if( (int)get_header.m_FrameType[1]==2 ) // full
	{
	  uint16_t data;
	  int coboId = (int)get_header.m_CoboIdx;
	  int ntb = n_words/(4*68);

	  if(coboId!=coboNum)
	    std::cout<<"#Warning! CoBo Id = "<<coboId <<" and filename "<<coboNum<<" do not match ! " <<std::endl;
	  int count[4] = {0,0,0,0};
	  double fadc[4][68][512]={{{0}}}; //[aget][ch][tb]

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
		    }
		  if(padid>0)
		    {
		      double adcSum=0;
		      int tbmax=0;
		      double adcMax=-999;
		      int layerid = getLayerID(padid-1);
		      int rowid = getRowID(padid-1);
		      int asadid = GetASADId(layerid,rowid);
		      for(int tbi = TBMIN; tbi<ntb+TBMIN; tbi++)
			{
			  tbHist[padid-1]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]);

			  adcSum+=fadc[ageti][chi][tbi];
			  if(tbi >TBCUTMIN && tbi<TBCUTMAX &&adcMax<fadc[ageti][chi][tbi])
			    {
			      adcMax=fadc[ageti][chi][tbi];
			      tbmax=tbi;
			    }
			}
		      if((adcMax+700*ageti) > tbHist_aget[asadid][ageti]->GetMaximum()){
			for(int tbi = TBMIN; tbi<TBNUM+TBMIN;tbi++){
			  tbHist_aget[asadid][ageti]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]+700*ageti);
			}
		      }
		      TVector3 position = getPoint(padid);
		      double offset = 15;

		      double adc_height = adcMax - ((double)adcSum/ntb);

		      if(adc_height > threshold){
			int maxBin = tbHist[padid-1]->GetMaximumBin();
			float maxTime = tbHist[padid-1]->GetBinCenter(maxBin);
			// if(maxTime > 120 || maxTime < 80) maxTime = 0;

		      }

		    }
		}
	    }
	}
      else if( (int)get_header.m_FrameType[1]==1 ) //partial
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
#if 0
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
		      int layerid = getLayerID(padid-1);
		      int rowid = getRowID(padid-1);
		      int asadid = GetASADId(layerid,rowid);
		      //int agetid = GetAGETId(asadid, layerid, rowid);
		      for(int tbi = TBMIN; tbi<TBNUM+TBMIN; tbi++)
			{
			  tbHist[padid-1]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]);
			  adcSum+=fadc[ageti][chi][tbi];
			  v.push_back(fadc[ageti][chi][tbi]);
			  if(adcMax<fadc[ageti][chi][tbi])
			    {
			      adcMax=fadc[ageti][chi][tbi];
			      tbmax=tbi;
			    }
			}
		      if((adcMax+700*ageti) > tbHist_aget[asadid][ageti]->GetMaximum()){
			for(int tbi = TBMIN; tbi<TBNUM+TBMIN;tbi++){
			  tbHist_aget[asadid][ageti]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]+700*ageti);
			}
		      }
		      double stddev = TMath::StdDev(v.size(), v.data());
		      // std::cout << padid << "\t" << stddev << std::endl;
		      // if(stddev == 0 || stddev > 50) continue;
		      double adc_height = adcMax - ((double)adcSum/TBNUM);

		      if(adc_height > threshold){
			int maxBin = tbHist[padid-1]->GetMaximumBin();
			float maxTime = tbHist[padid-1]->GetBinCenter(maxBin);
			// if(maxTime > 120 || maxTime < 80) maxTime = 0;

		      }

		    }
		}
	    }
	}
      else
	{
	  std::cerr<<"#E : FrameType error !, FrameType="<<(int)get_header.m_FrameType[1] <<
	    " cobo : " << coboNum << " asad : " << asadId << std::endl;
	  continue;
	}

    }

    fin.close();
  }

  std::cout << " * Finished * runNum : " << runNum << " eventNum : " << eventNum << std::endl;
}


void EventDisplay::GoBackward(){
  fCanvas->Clear();
  fCanvas->Divide(8,4,0.001,0.001);
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

void EventDisplay::GoForward(){
  fCanvas->Clear();
  fCanvas->Divide(8,4,0.001,0.001);

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

void EventDisplay::GoLast(){

  fCanvas->Clear();
  fCanvas->Divide(8,4,0.001,0.001);

  if(fLastEntry == -1){
    cout<<"No Updated Event"<<endl;
  }
  else
    {
      fCurrentEvent = fLastEntry;

      LoadEventData(runNumber, fCurrentEvent);
      DoCanvasDraw();
    }
}

void EventDisplay::GoTo(){

  fCanvas->Clear();
  fCanvas->Divide(8,4,0.001,0.001);
  TString eventNumberStr = fEvtHandler ->GetMarkedText();
  Int_t evnum = eventNumberStr.Atoi();

  /*
  if(evnum < 0 || evnum > fLastEntry)
    {
      cout<<"No Event" <<endl;
    }
  */

  //else
  // {
      fCurrentEvent = evnum;
      LoadEventData(runNumber, fCurrentEvent);
      DoCanvasDraw();
      //  }

}

void EventDisplay::DoCanvasDraw(){

  //waveform
  for(int i=0;i<TOTALASADNUM;i++){
    fCanvas->cd(i+1);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.08);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.02);
    for(int k=0;k<4;k++){
      tbHist_aget[i][k]->SetMinimum(0);
      tbHist_aget[i][k]->SetMaximum(3000);
      tbHist_aget[i][k]->Draw("same");
    }
  }
  fCanvas->Update();
}

void online_wf(Int_t runNum = -1, Int_t eventNum = 1)
{
  gStyle->SetOptStat(0);
  int argc = 0;
  char** argv = nullptr;
  TApplication theApp("App", &argc, argv);
  gDisplay = new EventDisplay(gClient->GetRoot(), 1000, 1000, runNum, eventNum);
  theApp.Run();
}
