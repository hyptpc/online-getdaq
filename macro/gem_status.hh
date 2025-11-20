#ifndef GEM_STATUS_HH
#define GEM_STATUS_HH

bool sec1_dead[20] = {1,0,0,0,0,
  		      0,0,0,0,0,
		      0,0,0,0,0,
		      0,0,0,0,0};
//bool sec2_dead[6]  = {0,0,1,1,1,0};
bool sec2_dead[6]  = {0,0,0,0,0,1};
bool sec3_dead[6]  = {0,0,0,0,0,0};
//bool sec4_dead[6]  = {1,0,0,1,1,1};
bool sec4_dead[6]  = {0};
bool sec3_mask_dead[21] = {0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0,0,
  0};
bool sec2_mask_dead[14] = {0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0};
bool sec4_mask_dead[14] = {0,0,0,0,0,
  0,0,0,0,0,
  0,0,0,0};
/*bool sec4_mask_dead[14] = {1,0,1,0,0,
  0,0,0,0,0,
  0,0,0,0};*/

#endif
