#ifndef PADHELPER_HH
#define PADHELPER_HH
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <TMath.h>
#include <TVector3.h>

const double ZTarget = -143;
//const double ZTarget = 0.;
//#OfPad #division #radius padLength
const Double_t padParameter[32][6]=
 {{0, 48,    14.75, 48, 0,  9.},
  {1, 48,    24.25, 48, 0,  9.},
  {2, 72,    33.75, 72, 0,  9.},
  {3, 96,    43.25, 96, 0,  9.},
  {4, 120,    52.75,120,0,   9.},
  {5, 144,    62.25,144,0,   9.},
  {6, 168,    71.75,168,0,   9.},
  {7, 192,    81.25,192,0,   9.},
  {8, 216,    90.75,216,0,   9.},
  {9, 240,    100.25,240,0,  9.},
  {10,208,    111.5,241, 0,  12.5},
  {11,218,    124.5,271, 0,  12.5},
  {12,230,    137.5,300, 0,  12.5},
  {13,214,    150.5,330, 0,  12.5},
  {14,212,    163.5,360, 0,  12.5},
  {15,214,    176.5,390, 0,  12.5},
  {16,220,    189.5,420, 0,  12.5},
  {17,224,    202.5,449, 0,  12.5},
  {18,232,    215.5,479, 0,  12.5},
  {19,238,    228.5,509, 0,  12.5},
  {20,244,    241.5,539, 0,  12.5},
  {21,232,    254.5,569, 0,  12.5},
  {22,218,    267.5,599, 0,  12.5},
  {23,210,    280.5,628, 0,  12.5},
  {24,206,    293.5,658, 0,  12.5},
  {25,202,    306.5,688, 0,  12.5},
  {26,200,    319.5,718, 0,  12.5},
  {27,196,    332.5,748, 0,  12.5},
  {28,178,    345.5,777, 0,  12.5},
  {29,130,    358.5,807, 0,  12.5},
  {30,108,    371.5,837, 0,  12.5},
  {31,90,     384.5,867, 0, 12.5}};

Int_t getPadID(Int_t layerID, Int_t rowID)
{
  Int_t padID=0;
  for(int layi = 0 ; layi<layerID; layi++) padID += padParameter[layi][1];
  padID+=rowID;
  return padID;

}

Int_t getLayerID(Int_t padID)
{
  padID-=1;
  int layer;
  int sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  return layer;
}

Int_t getRowID(Int_t padID)
{
  padID-=1;
  int layer, row;
  int sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;
  return row;
}

//_____________________________________________________________________________
inline Int_t GetASADId(Int_t layer, Int_t row) //0~30
{
  Int_t flag=layer/4;
  Int_t section;
  if(flag==0) section=0; //layer 0~3
  else if(flag==1) section=1; //layer 4~7
  else if(layer==30||layer==31) section=3; //layer 30~31
  else section=2; //layer 8~29

  Int_t half=padParameter[layer][1]/2;
  Int_t division1=padParameter[layer][1]/6;
  Int_t division2=padParameter[layer][1]*5/6;

  switch(section){
  case 0:
    if(row<half)
      return 0;
    else
      return 1;
  case 1:
    Int_t dummy;
    if(layer%4<2) dummy=2;
    if(2<=layer%4) dummy=5;
    if(row<division1||division2<=row)
      return dummy;
    else if(division1<=row&&row<half)
      return dummy+1;
    else
      return dummy+2;
  case 3:
    return 30;
  default:
    if(row<half)
      return layer-layer%2;
    else
      return layer-layer%2+1;
  }
}
/*
Double_t getTheta(Int_t layerID, Int_t rowID)
{
  Double_t sTheta = 180.-(360./padParameter[layerID][3])*padParameter[layerID][1]/2.;
  Double_t theta = sTheta+(rowID+0.5)*(360.-2*sTheta)/padParameter[layerID][1];
  return theta;
}
*/

Double_t getTheta(Int_t padID)
{
  padID-=1;
  int layer, row;
  int sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;
  Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
  Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
  return theta;
}

Double_t getTheta(Int_t layer, Double_t m_row)
{
  Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
  //Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
  Double_t theta = sTheta+(m_row+0.5)*360./padParameter[layer][3]-180;
  //std::cout<<"theta="<<theta<<std::endl;
  return theta;
}
Double_t getR(Int_t padID)
{
  padID-=1;
  int layer, row;
  int sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;
  Double_t R = padParameter[layer][2];
  return R;
}

TVector3 getPoint(int padID)
{
  padID-=1;
  int layer, row;
  int sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;

  TVector3 result;
  if (row > padParameter[layer][1]){ // out of range
    result.SetX(0);
    result.SetY(-1);
    result.SetZ(0);
  }
  else{
    double x, z;
    Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
    x = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
    z = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) - 143.0;
    result.SetX(x);
    result.SetY(0);
    result.SetZ(z);
  }
  return result;
}

TVector3 getPosition(Int_t padID)
{
  padID-=1;
  Int_t layer, row;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;

  TVector3 result;
  if (row > padParameter[layer][1]){ // out of range
    result.SetX(0);
    result.SetY(-1);
    result.SetZ(0);
  }
  else{
    Double_t x, z;
    //Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;

    //    x = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
    //    z = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) + ZTarget;

    // x = padParameter[layer][2] * sin(getTheta(padID+1)*TMath::Pi()/180.);
    // z = padParameter[layer][2] * cos(getTheta(padID+1)*TMath::Pi()/180.) + ZTarget;
    //std::cout<<"layer="<<layer<<", row"<<row<<std::endl;
    x = padParameter[layer][2] * sin(getTheta(layer,row)*TMath::Pi()/180.);
    z = padParameter[layer][2] * cos(getTheta(layer,row)*TMath::Pi()/180.) + ZTarget;

    // Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
    // Double_t x_ = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
    // Double_t z_ = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) + ZTarget;
    //std::cout<<"x="<<x<<", z="<<z<<", x_="<<x_<<", z_="<<z_<<std::endl;
    result.SetX(x);
    result.SetY(0);
    result.SetZ(z);
  }
  return result;
}



#endif
