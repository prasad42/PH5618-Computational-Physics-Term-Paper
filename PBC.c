#include"header.h"

float PBC(float x, float L){
  int Lrnd;
  if(x/L<0.5){
    Lrnd = 0;
  }
  if(x/L >= 0.5 && x/L < 1.5){
    Lrnd = 1;
  }
  if(x/L>=1.5){
    Lrnd = 2;
  }
  float xu = x - L*Lrnd;
  return xu;
}
