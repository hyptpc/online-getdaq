/* June 2018
 * Shinhyung */

#include "padHelper_new.hh"
#include "Unpacker.hh"

#define COBONUM 8
#define ASADNUM 4
#define PADNUM 5768
#define DEBUG 0

#define TBMIN 0
//#define TBNUM 200
//#define TBNUM 512
#define TBNUM 512
#define TBCUTMIN 0
#define TBCUTMAX 512

TFile *outfile;
TFile *file[COBONUM];
TTree *tree;
TCanvas* c1;
TCanvas* c2;
TH2Poly *poly;
TH1F* tbHist[PADNUM];
Double_t* adc = new Double_t;
Double_t padadc[4][4][68][512];
Int_t currentEvent=0;
Int_t runNumber=0;
int channelmap[4][31][68];
long eventoffset;
//bool asad_status[8][4];
TGraph *graph;
TH1D *hist;

bool UserStop=false;
bool endflag=false;

int get_active_pad_num(long runNum, long eventNum);
void next();

TClonesArray *padArray = nullptr;

void closeFile(int sig)
{
  std::cout << "#D closeFile" << std::endl;
  UserStop=true;
}
void event_rate(int runNum)
{
  outfile = new TFile(Form("event_rate/run_%d.root",runNum),"recreate");
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

  for(int i=0; i<PADNUM; i++)
  {
    tbHist[i]=new TH1F(Form("tbHist%d",i),"",TBNUM,TBMIN,TBNUM+TBMIN);
    tbHist[i]->SetMaximum(4100);
    tbHist[i]->SetMinimum(0);
  }

  poly = new TH2Poly("TPC","TPC",-300,300,-300,300);
  TPolyLine* p[PADNUM];

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
      poly->AddBin(5,X,Y);
      p[npl] = new TPolyLine(5,X,Y);
      p[npl]->SetUniqueID(npl);
      npl++;
    }
  }

  c1= new TCanvas("c1","",1200,500);
  c2= new TCanvas("c2","",500,500);
  graph = new TGraph();
  hist = new TH1D("hist","Active Pads", 5000,0,5000);
  signal(SIGINT, closeFile);
  //while(!UserStop)
  signal(SIGINT, closeFile);

  // ------------
  std::stringstream eventdatafile;
  eventdatafile.str("");
  eventdatafile<<"/home/hyptpc/work/shinhyung/debug/log/eventnum_run"<<std::setfill('0')<<std::setw(4)<<runNum<<".log";
  std::ifstream fevtall(eventdatafile.str());

  if(fevtall.fail())
  {
    std::cerr<<" --> file open fail ..." << eventdatafile.str() <<std::endl;
    exit(-1);
  }

  long eventNum =999999999;
  while(fevtall.good())
  {
    std::string buf;
    std::getline(fevtall,buf);
    if(buf[0]=='#' || buf.empty()) continue;
    std::istringstream is(buf);
    std::istream_iterator<std::string> issBegin(is);
    std::istream_iterator<std::string> issEnd;
    std::vector<std::string> param(issBegin,issEnd);
    if(param.empty()||param[0].empty()) continue;
    int evtmin = atoi(param[0].c_str());
    int evtmax = atoi(param[1].c_str());
    int evtnum = atoi(param[2].c_str());
    if(evtnum < eventNum && evtnum >0) eventNum=evtnum;
  }
  fevtall.close();
  // ------------

  while(eventNum < 0) 
  {
    std::cerr << " Check eventNum (" << eventNum << ")"<< std::endl;
    sleep(1);
    continue;
  }
  eventNum = 140;

  int total_active_event=0;
  int total_spark_event=0;
  for( int eventi=0; eventi<eventNum; eventi++)
  {
    int active_pad_num = get_active_pad_num(runNum, eventi);
    graph->SetPoint(eventi, eventi, active_pad_num);
    hist->Fill(active_pad_num);
    if( active_pad_num >= 30 && active_pad_num < 50) total_active_event++;
    if( active_pad_num >= 4500) total_spark_event++;
  }
  cout << " # total event   : " << eventNum << endl;
  cout << " # active events : " << total_active_event << endl;
  //cout << " ratio 	      : " 
  //<< (double)total_active_event/eventNum << endl;
  const double event_rate = 0.433; // Hz 
  double total_time_min = eventNum / event_rate / 60.;
  cout << " active evt rate : " << total_active_event / total_time_min << endl;
  cout << " active spark rate : " << total_spark_event / total_time_min << endl;
  c1->cd();
  graph->SetMarkerColor(kBlue);
  graph->SetMarkerSize(0.3);
  graph->SetMarkerStyle(21);
  graph->Draw("AP");
  c2->cd();
  hist->Draw();

  outfile->cd();
  graph->SetName("graph");
  graph->Write();
  hist->Write();
}

int get_active_pad_num(long runNum, long eventNum)
{
  runNumber = runNum;
  currentEvent = eventNum;
  bool flag = true;
  bool eventflag=true;

  int activepad=0;

//  std::ifstream fasad;
//  std::stringstream asad_file;
//  asad_file.str("");
//  asad_file<<"/home/hyptpc/ganacq_manip/test/acquisition/log/asad_status/run_"
//    <<std::setfill('0')<<std::setw(4)<<runNum<<".log";
//  fasad.open(asad_file.str());
//  if(!fasad.is_open())
//  {
//    std::cerr << asad_file.str() << "  does not exist. " << std::endl;
//    std::exit(-1);
//  }
//  int dum1, dum2;
//  while(fasad >> dum1 >> dum2)
//  {
//    asad_status[dum1/4][dum1%4] = dum2;
//  }

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
      //for(int asadNum = 0; asadNum<ASADNUM ;asadNum++)
    {
      //for(int ai=0; ai<4; ai++)
	//if(!asad_status[coboNum][ai]) asadflag[ai]++;
      if (coboNum == 7) asadflag[3]++, asadflag[2]++;
      struct GetHeader_t get_header;
      int n_words;
      int padded;
      //cout << " padded : " << padded << endl;
      int asadId;
      while(!fin.eof() && (unsigned long)fin.tellg()<length)
      {
	// Read Header
	fin.read((char*)&get_header,sizeof(get_header));
#if DEBUG
	std::cout << "now : " << std::hex << (unsigned long)fin.tellg() << std::endl;
	std::cout << std::endl;
	std::cout << " 5 " << std::endl;
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


	//if(eventNumber-eventoffset == eventNum && asadId == asadNum)  
	if(eventNumber-eventoffset == eventNum)
	  //if(eventNumber == eventNum) 
	{
	  asadflag[asadId]++;
	  //cout << "=== Reading CoboNum#" << coboNum << ", Asad#" << asadId << "..."<<endl;
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
	  //std::cout << " 1 " << std::endl;

	}
	//std::cout << " 2 " << std::endl;
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
	  //int t_bucket = (count[agetid]/68)%512;
	  //int t_bucket = (count[agetid]/68)%200;
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
	    if(padid>0)
	    {
	      double adcSum=0;
	      int tbmax=0;
	      double adcMax=-999;
	      //for(int tbi = 0; tbi<512; tbi++)
	      //for(int tbi = TBMIN; tbi<ntb+TBMIN; tbi++)
	      for(int tbi = 0; tbi<ntb; tbi++)
	      {
		tbHist[padid-1]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]);
		//cout << "tbi="<< tbi << " adc " << fadc[ageti][chi][tbi] << endl;
		//sleep(1);
		adcSum+=fadc[ageti][chi][tbi];
		if(tbi >TBCUTMIN && tbi<TBCUTMAX &&adcMax<fadc[ageti][chi][tbi]) 
		{
		  adcMax=fadc[ageti][chi][tbi];
		  tbmax=tbi;
		}
	      }
	      //poly->SetBinContent(padid,1000);
	      poly->SetBinContent(padid,adcMax-((double)adcSum/ntb));

	      TVector3 position = getPoint(padid);
	      double offset = 16;
	      if( !((position.X() < position.Z()+offset
		      && position.X() > position.Z()-offset)
		    || (position.X() < -position.Z()+offset
		      && position.X() > -position.Z()-offset) )
		  && adcMax - (double)adcSum/ntb > 200 ) activepad++;
	      //poly->SetBinContent(padid,adcMax);
	      //poly->SetBinContent(padid,adcMax);
	      //poly->SetBinContent(padid,adcSum);
	      //cout << dec << " padid : " << padid << endl;
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
	      for(int tbi = TBMIN; tbi<TBNUM+TBMIN; tbi++)
	      {
		tbHist[padid-1]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]);
		//cout << "tbi="<< tbi << " adc " << adc[tbi] << endl;
		//sleep(1);
		adcSum+=fadc[ageti][chi][tbi];
		if(adcMax<fadc[ageti][chi][tbi]) 
		{
		  adcMax=fadc[ageti][chi][tbi];
		  tbmax=tbi;
		}
	      }
	      poly->SetBinContent(padid,adcMax-(adcSum/TBNUM));
	      //cout << " padid = " << padid << " , content : " << adcMax-(adcSum/TBNUM) << endl;
	      if( adcMax-(adcSum/TBNUM) >0 ) activepad++;
	      //cout << " adcMax : " << adcMax << endl;
	    }
	  }
	}
	//poly->SetBinContent(1,1);
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
  cout << "eventNum: " << eventNum << ", # of active pad : " << activepad << endl;
  return activepad;
}
