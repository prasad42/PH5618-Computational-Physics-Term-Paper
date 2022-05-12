#include"header.h"

float LJ(float rij){
  float U;
  if(rij<=2.5){
    U = 4*(1/pow(rij,12)-1/pow(rij,6));
  }
  else{
    U = 0;
  }
  return U;
}
