// -*- C++ -*-

#include <fstream>
#include <iostream>
#include <signal.h>

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2Poly.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>

#include "gem_status.hh"
#include "padHelper_new.hh"
#include "Unpacker.hh"

#define DEBUG 0
#define ASADSTATUS 0

namespace
{
  //const TString dir = "bench_e03_2021feb";
  //const TString dir = "e03_2021feb";
  // const TString dir = "bench_e42_2021may";
  const TString dir = "bench_e72_2024jan";
  const Int_t NumOfCoBo        = 8;
  const Int_t NumOfAsAdPerCoBo = 4;
  const Int_t NumOfAsAd        = 31; // #32 is unused
  const Int_t NumOfPad         = 5768;
  const Int_t NumOfAGETPerAsAd = 4;
  const Int_t NumOfChannelAGET = 68;
  const Int_t MinTimeBucket    = 0;
  //const Int_t NumOfTimeBucket = 200; // himac
  //const Int_t NumOfTimeBucket = 140; // himac
  const Int_t NumOfTimeBucket = 170; // himac
  //const Int_t NumOfTimeBucket = 512; // maximum
  // const Int_t NumOfTimeBucket = 256; // E42
  const Int_t CutMinTimeBucket = 0;
  const Int_t CutMaxTimeBucket = 512;

  Bool_t g_user_stop = false;
  TCanvas*  c1 = nullptr;
  TH2Poly*  poly = nullptr;
  TH2Poly*  poly2 = nullptr;
  std::vector<TH1F*> h_fadc( NumOfPad );
  ULong64_t g_run_number   = 0;
  ULong64_t g_event_number = 0;
  ULong64_t g_total_event_number = 0;
  Int_t     channel_map[NumOfAGETPerAsAd][NumOfAsAd][NumOfChannelAGET];
  Bool_t asad_status[NumOfCoBo][NumOfAsAdPerCoBo];
  Bool_t cobo_status[NumOfCoBo];
  Bool_t endflag = false;
  TClonesArray* padArray = nullptr;
}

//_____________________________________________________________________________
void user_stop( Int_t sig )
{
  std::cout << "Catch signal " << sig << std::endl
	    << "quit" << std::endl;
  // g_user_stop = true;
  gSystem->Exit( EXIT_SUCCESS );
}

//_____________________________________________________________________________
void draw( Int_t run_number, ULong64_t event_number )
{
  std::cout << "draw() run_number : " << run_number
	    << ", event_number : " << event_number << std::endl;
  poly->Reset("ICESM");
  poly2->Reset("ICESM");
  g_run_number = run_number;
  g_event_number = event_number;
  Bool_t flag = true;
  Bool_t eventflag = true;

  Int_t activepad=0;
#if ASADSTATUS
  std::ifstream fasad;
  std::stringstream asad_file;
  asad_file.str("");
  asad_file << "/home/axis/ganacq_manip/"
	    << dir << "/acquisition/log/asad_status/run_"
	    << std::setfill('0') << std::setw(4) << run_number << ".log";
  fasad.open( asad_file.str() );
  if( !fasad.is_open() ){
    std::cerr << asad_file.str() << "  does not exist. " << std::endl;
    std::exit(-1);
  }
  Int_t dum1, dum2;
  while( fasad >> dum1 >> dum2 ){
    asad_status[dum1/4][dum1%4] = dum2;
  }
  for( Int_t i=0; i<NumOfCoBo; ++i ){
    Int_t count = 0;
    for( Int_t j=0; j<NumOfAsAdPerCoBo; ++j ){
      if(asad_status[i][j]) count++;
    }
    if( count ) cobo_status[i] = 1;
    else        cobo_status[i] = 0;
  }
#endif

  for( Int_t icobo=0; icobo<NumOfCoBo; ++icobo ){
#if ASADSTATUS
    if(!cobo_status[icobo]) continue;
#endif
    Int_t asadflag[4] = {0};
    std::ifstream fin;
//     TString raw_data_path = Form( "run_E03Feb/cobo%d/run_%04d.dat",
//     				  icobo, run_number );
//    TString raw_data_path = Form( "bench_e42_2021may/cobo%d/run_%04d.dat",
//				  icobo, run_number );
    // TString raw_data_path = Form( "e42_2021may/cobo%d/run_%04d.dat",
    // 				  icobo, run_number );
    TString raw_data_path = Form( "bench_e72_2024jan/cobo%d/run_%04d.dat",
    				  icobo, run_number );
    fin.open( raw_data_path.Data(), std::ios::in|std::ios::binary );
    if( !fin.is_open() ){
      std::cerr<<"#Warning! :"<< raw_data_path << " does not exist !"
	       << std::endl;
      continue;
    } else {
      std::cout << "#D " << raw_data_path << " opened " << std::endl;
    }
    fin.seekg( 0, fin.end );
    auto length = fin.tellg();
    fin.seekg( 0, fin.beg );
    fin.seekg( 12 );
    if( !length ) continue;
    while( !( asadflag[0] &&  asadflag[1] &&  asadflag[2] &&  asadflag[3]) ){
#if ASADSTATUS
      for( Int_t ai=0; ai<4; ai++ ){
	if( !asad_status[icobo][ai] ) asadflag[ai]++;
      }
#else
      //if (icobo == 7) asadflag[3]++, asadflag[2]++;
      if (icobo == 7) asadflag[3]++;
      //if (icobo == 4) asadflag[3]++;
      //if (icobo == 5) asadflag[1]++;
      //if (icobo == 1) asadflag[3]++;
      if (icobo == 0) asadflag[0]++;
#endif
      GetHeader_t get_header;
      Int_t n_words = -1;
      Int_t padded  = -1;
      Int_t asadId  = -1;
      while(!fin.eof() && fin.tellg()<length){
	// Read Header
	fin.read( (char*)&get_header, sizeof(get_header) );
#if DEBUG
	std::cout << "now : " << std::hex << fin.tellg() << std::endl
		  << std::endl;
	//getchar();
#endif
	asadId = get_header.m_AsadIdx;
	n_words = ( get_header.m_nItems[0] << 24 |
		    get_header.m_nItems[1] << 16 |
		    get_header.m_nItems[2] <<  8 |
		    get_header.m_nItems[3] <<  0 );
	ULong64_t evnum = ( get_header.m_EventIdx[0] << 24 |
			    get_header.m_EventIdx[1] << 16 |
			    get_header.m_EventIdx[2] <<  8 |
			    get_header.m_EventIdx[3] <<  0 );
	UInt_t frameSize = ( get_header.m_FrameSize[0] << 16 |
			     get_header.m_FrameSize[1] <<  8 |
			     get_header.m_FrameSize[2] <<  0 );
	Int_t headerSize = ( get_header.m_HeaderSize[0] << 8 |
			     get_header.m_HeaderSize[1] << 0 );
	Int_t item_size = ( get_header.m_ItemSize[0] << 8 |
			    get_header.m_ItemSize[1] << 0 );
	static ULong64_t event_offset = 0;
	if( eventflag ){
	  event_offset=evnum;
	  eventflag=false;
	}
	padded = (frameSize-headerSize)*256 - item_size*n_words;
#if DEBUG
	std::cout << std::hex << std::setfill('0') << std::endl
		  << TString('=', 80) << std::endl
		  << " frameSize  3B : " << std::setw(2)
		  << get_header.m_FrameSize[0] << " " << std::setw(2)
		  << get_header.m_FrameSize[1] << " " << std::setw(2)
		  << get_header.m_FrameSize[2] << std::endl
		  << " eventIdx   4B : " << std::setw(2)
		  << get_header.m_EventIdx[0] << " " << std::setw(2)
		  << get_header.m_EventIdx[1] << " " << std::setw(2)
		  << get_header.m_EventIdx[2] << " " << std::setw(2)
		  << get_header.m_EventIdx[3] << " --> " << std::dec
		  << evnum << std::endl
		  << " coboIdx    1B : " << std::setw(2)
		  << get_header.m_CoboIdx << std::endl
		  << " asadIdx    1B : " << std::setw(2)
		  << get_header.m_AsadIdx << std::endl
		  << " frameType     : " << get_header.m_FrameType[1]
		  << std::endl
		  << " itemSize   : " << item_size << std::endl
		  << " nItems     : " << n_words << std::endl
		  << " padded     : " << padded << std::endl
		  << TString('=', 80) << std::endl;
	getchar();
#endif
	if( evnum-event_offset == event_number ){
	  asadflag[asadId]++;
	  std::cout << "=== Reading CoboNum#" << icobo << ", Asad#"
		    << asadId << "..." << std::endl;
	  break;
	} else {
	  fin.seekg( (ULong_t)fin.tellg() + (frameSize-headerSize)*256 );
	}
      }
      // Read Data
      if( get_header.m_FrameType[1]==2 ){ // full
	uint16_t data;
	Int_t coboId = get_header.m_CoboIdx;
	Int_t ntb = n_words/(4*NumOfChannelAGET);
	if( coboId != icobo ){
	  std::cout << "#Warning! CoBo Id = " << coboId << " and filename "
		    << icobo << " do not match ! " << std::endl;
	}
	Int_t    count[NumOfAGETPerAsAd] = {};
	Double_t fadc[NumOfAGETPerAsAd][NumOfChannelAGET][NumOfTimeBucket]={{{}}}; //[aget][ch][tb]
	Int_t    hit = 0;
	for( Int_t i=0; i<n_words; ++i ){
	  fin.read( (char*)&data, sizeof(data) );
	  Int_t agetid   = (data>>6)&0x03;
	  Int_t adc_low  = (data>>8)&0xff;
	  Int_t adc_high = data&0x0f;
	  Int_t adc      = (adc_high<<8)|adc_low;
	  Int_t ch       = count[agetid]%NumOfChannelAGET;
	  //Int_t t_bucket = (count[agetid]/NumOfChannelAGET)%NumOfTimeBucket;
	  //Int_t t_bucket = (count[agetid]/NumOfChannelAGET)%200;
	  Int_t t_bucket = (count[agetid]/NumOfChannelAGET)%ntb;
	  fadc[agetid][ch][t_bucket]=adc;
	  ++count[agetid];
	}
	fin.seekg( (ULong_t)fin.tellg() + padded );
	Int_t padid=0;
	for( Int_t ageti = 0; ageti<NumOfAGETPerAsAd; ++ageti ){
	  for( Int_t chi = 0; chi<NumOfChannelAGET; ++chi ){
	    Int_t chi2;
	    if(chi>=0 && chi<=10) chi2=chi;
	    else if(chi>=12 && chi<=21) chi2 = chi-1;
	    else if(chi>=23 && chi<=44) chi2 = chi-2;
	    else if(chi>=46 && chi<=55) chi2 = chi-3;
	    else if(chi>=57 && chi<=67) chi2 = chi-4;
	    else continue;
	    if(chi==11 || chi==22 ||chi==45||chi==56){
	      padid=-7;
	    } else {
	      Int_t id;
	      padid=channel_map[ageti][icobo*4+asadId][chi2];
	    }
	    if( padid>0 ){
	      Double_t adcSum=0;
	      Int_t tbmax=0;
	      Double_t adcMax=-999;
	      std::vector<Double_t> fadc_vec(ntb);
	      for( Int_t tbi = 0; tbi<ntb; ++tbi ){
		h_fadc[padid-1]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]);
		adcSum+=fadc[ageti][chi][tbi];
		fadc_vec[tbi] = fadc[ageti][chi][tbi];
		if( tbi > CutMinTimeBucket &&
		    tbi < CutMaxTimeBucket &&
		    adcMax < fadc[ageti][chi][tbi] ){
		  adcMax = fadc[ageti][chi][tbi];
		  tbmax = tbi;
		}
	      }
	      poly->SetBinContent( padid, adcMax-((Double_t)adcSum/ntb) );
	      TVector3 position = getPoint( padid );
	      Double_t offset = 16;
	      if( !((position.X() < position.Z()+offset
		     && position.X() > position.Z()-offset)
		    || (position.X() < -position.Z()+offset
			&& position.X() > -position.Z()-offset) )
		  && adcMax - (Double_t)adcSum/ntb > 200 ) activepad++;
	    }
	  }
	}
      } else if( get_header.m_FrameType[1]==1 ){ //partial
	uint32_t data;
	Int_t asadId = (Int_t)get_header.m_AsadIdx;
	Int_t coboId = (Int_t)get_header.m_CoboIdx;
	Int_t ntb = n_words/(4*NumOfChannelAGET);
	if( coboId!=icobo ){
	  std::cout << "#Warning! CoBo Id = " << coboId << " and filename "
		    << icobo << " do not match ! " <<std::endl;
	}
	Double_t fadc[NumOfAGETPerAsAd][NumOfChannelAGET][NumOfTimeBucket]={{{0}}}; //[aget][ch][tb]
	for( Int_t i=0; i<n_words;++i ){
	  fin.read((char*)&data,sizeof(data));
	  Int_t agetid        = (data>>6)&0x03;
	  Int_t adc_low       = (data>>24)&0xff;
	  Int_t adc_high      = (data>>16)&0x0f;
	  Int_t adc           = (adc_high<<8)|adc_low;
	  Int_t ch_low 	      = (data>>15)&0x01;
	  Int_t ch_high       = data&0x3f;
	  Int_t ch            = (ch_high<<1)|ch_low;
	  Int_t t_bucket_low  = (data>>22)&0x03;
	  Int_t t_bucket_high = (data>>8)&0x7f;
	  Int_t t_bucket      = (t_bucket_high<<2)|t_bucket_low;
	  // if( fadc[agetid][ch][t_bucket] > 0 ){
	  //   std::cerr << "#W fadc is already filled #" << i << "/" << n_words
	  // 	      << " (aget,ch,tb)=("
	  // 	      << agetid << "," << ch << "," << t_bucket << ") "
	  // 	      << fadc[agetid][ch][t_bucket] << " -> "
	  // 	      << adc << std::endl;
	  // }
	  fadc[agetid][ch][t_bucket]=adc;
#if DEBUG
	  std::cout << std::dec << " asad : " << icobo*4+asadId << std::endl;
	  std::cout << " aget : " << agetid << std::endl;
	  std::cout << " ch : " << ch << std::endl;
	  std::cout << " t_bucket : " << t_bucket << std::endl;
	  std::cout << " adc : " << adc << std::hex<<  std::endl;
	  //getchar();
#endif
	}
	fin.seekg( (unsigned long)fin.tellg() + padded );
	Int_t padid=0;
	for( Int_t ageti=0; ageti<4; ageti++ ){
	  for(Int_t chi = 0; chi<NumOfChannelAGET; chi++){
	    Int_t chi2;
	    if(chi>=0 && chi<=10) chi2=chi;
	    else if(chi>=12 && chi<=21) chi2 = chi-1;
	    else if(chi>=23 && chi<=44) chi2 = chi-2;
	    else if(chi>=46 && chi<=55) chi2 = chi-3;
	    else if(chi>=57 && chi<=67) chi2 = chi-4;
	    else continue;
	    if(chi==11 || chi==22 ||chi==45||chi==56){
	      padid=-7;
	    } else {
	      padid=channel_map[ageti][icobo*4+asadId][chi2];
	    }
	    if( padid > 0 ){
	      if( padid==323 )
		std::cout << "cobo=" << coboId << ", asad=" << asadId
			  << ", aget=" << ageti << ", fech=" << chi
			  << ", fakech=" << chi2 << ", padid=" << padid
			  << std::endl;
	      Double_t adcSum=0;
	      Int_t tbmax=0;
	      Double_t adcMax=-999;
	      for(Int_t tbi=MinTimeBucket; tbi<NumOfTimeBucket+MinTimeBucket; tbi++ ){
		h_fadc[padid-1]->SetBinContent(tbi+1,fadc[ageti][chi][tbi]);
		//cout << "tbi="<< tbi << " adc " << adc[tbi] << endl;
		//sleep(1);
		adcSum+=fadc[ageti][chi][tbi];
		if( adcMax<fadc[ageti][chi][tbi] ){
		  adcMax=fadc[ageti][chi][tbi];
		  tbmax=tbi;
		}
	      }
	      poly->SetBinContent(padid,adcMax-(adcSum/NumOfTimeBucket));
	      //cout << " padid = " << padid << " , content : " << adcMax-(adcSum/NumOfTimeBucket) << endl;
	      if( adcMax-(adcSum/NumOfTimeBucket) >0 ) activepad++;
	      //cout << " adcMax : " << adcMax << endl;
	    }
	  }
	}
	//poly->SetBinContent(1,1);
      } else{
	std::cerr << "#E : FrameType error !, FrameType="
		  << get_header.m_FrameType[1] << " cobo : " << icobo
		  << " asad : " << asadId << std::endl;
	std::exit(-1);
      }
    }
    fin.close();
  }
  std::cout << "# of active pad : " << activepad << std::endl;

  c1->Clear();
  c1->cd()->SetLogz();
  poly->SetMinimum(0);
  poly->SetMaximum(4000);
  //poly->SetMaximum(1000);
  //poly->SetMaximum(500);
  poly->Draw("colz");

  c1->AddExec("dynamic","ShowPulse()");
  //TText *t1 = new TText(-250, 300, Form("run#%ld",run_number));
  //TText *t2 = new TText(-250, 270, Form("evt#%ld",event_number));
  static TText *t1 = new TText;
  static TText *t2 = new TText;
  t1->SetTextAlign(12);
  t2->SetTextAlign(32);
  t1->SetTextSize(0.05);
  t2->SetTextSize(0.05);
  t1->DrawText(-285, 275, Form("Run#%llu", g_run_number));
  t2->DrawText( 285, 275, Form("Event#%llu", g_event_number));
  c1->Update();
}

//_____________________________________________________________________________
void next( void )
{
  g_event_number++;
  draw(g_run_number, g_event_number);
}

//_____________________________________________________________________________
void prev( void )
{
  if( g_event_number > 0 )
    g_event_number--;
  draw( g_run_number, g_event_number );
}

//_____________________________________________________________________________
void current( void )
{
  //  std::ifstream frun("/home/axis/work/shinhyung/switch/runnum.log");
  std::ifstream frun("/home/axis/ganacq_manip/"+dir+"/"+dir+".run_number");
  if(!frun.is_open())
    {
      std::cerr << " test.run_number does not exist. " << std::endl;
      std::exit(-1);
    }
  frun >> g_run_number;
  // ------------
  std::stringstream eventdatafile;
  eventdatafile.str("");
  eventdatafile<<"/home/axis/work/shinhyung/debug/log/eventnum_run"<<std::setfill('0')<<std::setw(4)<<g_run_number<<".log";
  std::ifstream fevtall(eventdatafile.str());

  if(fevtall.fail())
    {
      std::cerr<<" --> file open fail ..." << eventdatafile.str() <<std::endl;
      exit(-1);
    }

  long event_number =999999999;
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
      Int_t evtmin = atoi(param[0].c_str());
      Int_t evtmax = atoi(param[1].c_str());
      Int_t evtnum = atoi(param[2].c_str());
      if(evtnum < event_number && evtnum >0) event_number=evtnum;
    }
  g_total_event_number = event_number-1;
  fevtall.close();
  draw(g_run_number, g_total_event_number);
}

//_____________________________________________________________________________
void current_event( Int_t Runnum )
{
  Int_t run_number = Runnum;
  g_run_number = Runnum;

  // ------------
  std::stringstream eventdatafile;
  eventdatafile.str("");
  eventdatafile<<"/home/axis/work/shinhyung/debug/log/eventnum_run"<<std::setfill('0')<<std::setw(4)<<g_run_number<<".log";
  std::ifstream fevtall(eventdatafile.str());
  std::cout<<eventdatafile.str()<<std::endl;
  if(fevtall.fail())
    {
      std::cerr<<" --> file open fail ..." << eventdatafile.str() <<std::endl;
      exit(-1);
    }

  long event_number =999999999;
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
      Int_t evtmin = atoi(param[0].c_str());
      Int_t evtmax = atoi(param[1].c_str());
      Int_t evtnum = atoi(param[2].c_str());
      if(evtnum < event_number && evtnum >0) event_number=evtnum;
    }
  g_total_event_number = event_number-1;
  fevtall.close();
  draw(run_number, g_total_event_number);
}


void monitor(Int_t Runnum)
{
  while( !gSystem->ProcessEvents() && !g_user_stop ){
    current_event(Runnum);
    for( Int_t i=0; i<5000 && !gSystem->ProcessEvents(); ++i ){
      gSystem->Sleep(1);
    }
  }
}

void monitor( void )
{
  while( !gSystem->ProcessEvents() && !g_user_stop ){
    current();
    for( Int_t i=0; i<5000 && !gSystem->ProcessEvents(); ++i ){
      gSystem->Sleep(1);
    }
  }
}


//_____________________________________________________________________________
void nextrun( void )
{
  g_run_number++;
  draw(g_run_number, 0);
}

//_____________________________________________________________________________
void prevrun( void )
{
  g_run_number--;
  draw(g_run_number, 0);
}

//_____________________________________________________________________________
void print( TString memo )
{
  TString pdffile ="";
  pdffile.Form("pdf/run%llu_ev%llu_", g_run_number, g_event_number);
  pdffile.Append(memo);
  pdffile.Append(".pdf");
  c1->Print(pdffile.Data());
}

//_____________________________________________________________________________
void ShowPulse( void )
{
  TObject* select = gPad->GetSelected();
  if( !select ) return;
  if( !select->InheritsFrom( TH2Poly::Class() ) ){
    gPad->SetUniqueID(0);
    return;
  }
  TH2Poly *pp = dynamic_cast<TH2Poly*>( select );
  gPad->GetCanvas()->FeedbackMode( kTRUE );

  Int_t pyold = gPad->GetUniqueID();
  Int_t px    = gPad->GetEventX();
  Int_t py    = gPad->GetEventY();
  Float_t upy = gPad->AbsPixeltoY(py);
  Float_t upx = gPad->AbsPixeltoX(px);
  TVirtualPad *padsav = gPad;
  TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
  if(c2) delete c2->GetPrimitive("Pulse");
  else c2= new TCanvas("c2","Pulse Canvas", 600,400);
  c2->cd();

  std::cout << "padID : "<< pp->FindBin(upx,upy) << std::endl;
  //h_fadc[pp->FindBin(upx,upy)-1]->GetYaxis()->SetRangeUser(0,4000.);
  //h_fadc[pp->FindBin(upx,upy)-1]->GetYaxis()->SetRangeUser(300,750.);
  h_fadc[pp->FindBin(upx,upy)-1]->Draw();
  Double_t adcSum = 0.;
  for(Int_t i=0; i<NumOfTimeBucket; ++i){
    adcSum+=h_fadc[pp->FindBin(upx,upy)-1]->GetBinContent(i+1);
  }
  Double_t maxadc = h_fadc[pp->FindBin(upx,upy)-1]->GetBinContent(h_fadc[pp->FindBin(upx,upy)-1]->GetMaximumBin());
  std::cout<<"mean:"<<adcSum/(double)NumOfTimeBucket<<", maxadc:"<<maxadc<<std::endl;


  c2->Update();
  padsav->cd();
}

//_____________________________________________________________________________
void DispHypTPC( void )
{
  //gStyle->SetOptStat( 0 );
  // gStyle->SetOptTitle( 0 );
  gStyle->SetPadTopMargin( 0.1 );
  gStyle->SetPadBottomMargin( 0.1 );
  gStyle->SetPadRightMargin( 0.1 );
  gStyle->SetPadLeftMargin( 0.1 );

  {
    std::ifstream ifs("param/channel_map_20180522.param");
    if( !ifs.is_open() ){
      std::cerr<<" --> file open fail ..."<<std::endl;
      std::exit( EXIT_FAILURE );
    }
    while( ifs.good() ){
      std::string buf;
      std::getline(ifs,buf);
      if(buf[0]=='#' || buf.empty())continue;
      std::istringstream is(buf);
      std::istream_iterator<std::string> issBegin(is);
      std::istream_iterator<std::string> issEnd;
      std::vector<std::string> param(issBegin,issEnd);
      if(param.empty()||param[0].empty())continue;
      Int_t aget  = std::atoi( param[0].c_str() );
      Int_t asad  = std::atoi( param[1].c_str() );
      Int_t ch    = std::atoi( param[2].c_str() );
      Int_t layer = std::atoi( param[3].c_str() );
      Int_t row   = std::atoi( param[4].c_str() );
      Int_t pid   = std::atoi( param[5].c_str() );
      if( aget!=0 ){
	channel_map[aget-1][asad-1][ch-1] = pid;
      }
    }
  }

  for( Int_t i=0; i<NumOfPad; ++i ){
    if( h_fadc[i] ) delete h_fadc[i];
    h_fadc[i] =new TH1F( Form("h_fadc%d", i+1),
			 Form("FADC #%d;Time Bucket;ADC", i+1),
			 NumOfTimeBucket, MinTimeBucket,
			 NumOfTimeBucket + MinTimeBucket );
    //h_fadc[i]->SetMaximum( 2048 );
    h_fadc[i]->SetMaximum( 4096 );
    //h_fadc[i]->SetMaximum( 1000 );
    h_fadc[i]->SetMinimum( -10 );
  }

  if( poly ) delete poly;
  poly = new TH2Poly( "TPC",";Z;X", -300, 300, -300, 300 );
  if( poly2 ) delete poly2;
  poly2 = new TH2Poly( "TPC",";Z;X", -300, 300, -300, 300 );

  Double_t X[5];
  Double_t Y[5];
  for( Int_t i=0; i< 32; i++ ){
    Double_t pLength = padParameter[i][5];
    Double_t st      = 180.-(360./padParameter[i][3])*padParameter[i][1]/2.;
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./padParameter[i][3])/180.*TMath::Pi();
    Double_t cRad    = padParameter[i][2];
    Int_t    nPad    = padParameter[i][1];
    for( Int_t j=0; j<nPad; j++ ){
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
      for( Int_t ii=0; ii<5; ii++ ) X[ii] -=143;
      poly->AddBin( 5, X, Y );
      poly2->AddBin( 5, X, Y );
    }
  }

  if( c1 ) delete c1;
  //c1= new TCanvas("c1","",1400,10,500,500);
  c1= new TCanvas( "c1", "c1", 500, 500 );
  //c1= new TCanvas( "c1", "c1", 500, 1000 );
  ::signal( SIGINT, user_stop );
}
