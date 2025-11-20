// Author : Hitoshi Sugimura

#ifndef UNPACKER_H
#define UNPACKER_H

struct GetHeader_t{
  uint8_t m_MetaType;         //1
  uint8_t m_FrameSize[3];     //4
  uint8_t m_DataSource;       //5
  uint8_t m_FrameType[2];     //7
  uint8_t m_Revision;         //8
  uint8_t m_HeaderSize[2];    //10
  uint8_t m_ItemSize[2];      //12
  uint8_t m_nItems[4];        //16
  uint8_t m_EventTime[6];     //22
  uint8_t m_EventIdx[4];      //26
  uint8_t m_CoboIdx;          //27
  uint8_t m_AsadIdx;          //28
  uint8_t m_ReadOffset[2];    //30
  uint8_t m_status;           //31
  uint8_t m_HitPatAGet0[9];   //40
  uint8_t m_HitPatAGet1[9];   //49
  uint8_t m_HitPatAGet2[9];   //58
  uint8_t m_HitPatAGet3[9];   //67
  uint8_t m_MultiPlAGet0[2];  //69
  uint8_t m_MultiPlAGet1[2];  //71
  uint8_t m_MultiPlAGet2[2];  //73
  uint8_t m_MultiPlAGet3[2];  //75
  uint8_t m_WindowOut[4];     //79
  uint8_t m_LastCellAGet0[2]; //81
  uint8_t m_LastCellAGet1[2]; //83
  uint8_t m_LastCellAGet2[2]; //85
  uint8_t m_LastCellAGet3[2]; //87
  uint8_t m_reserve1byte;     //88
  uint8_t m_reserve4byte[168];//256             <---(42*4=168)
};
#endif
