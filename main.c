#include "header.h"
     
int main(){

  //////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////// 
  //INPUT PARAMETERS
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  
  int N = 300;		//Total Number of Particles
  long double rho = 0.45;	//Total Density
  long double L = pow(N/rho,0.33333);		//Cube Root of Total Volume
  long double Temp = 0.75;	//Temperature
  
  long double l = 1;		//Maximum distance for particle displacement
  long double l1 = 0.01;		//Maximum distance for volume change
  long double V = pow(L, 3);		//Total volume
  
  int MCTot = 100000;		//Total MC Steps
  
  //////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////// 
  //SETUP
  //////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////// 
  long double U = 0;
  ///////////////
  //RANDOM NUMBER	
  ///////////////
  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_mt19937; //MT RAND. NUM. GENERATOR //T =gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  //////////////////////
  //SETUP FOR INITIALIZATION
  //////////////////////
  int n1 = N/2;					// Number of particles in box 1
  long double L1 = pow((pow(L,3)/2),0.33333);	 	// Length of box 1
  int n2 = N-n1;				// Number of particles in box 2
  long double L2 = pow((pow(L,3)-pow(L1,3)),0.33333);	// Length of box 2
  long double V1 = pow(L1, 3);
  long double V2 = pow(L2, 3);
  long double V2o;
  
  ///////
  //FILES
  ///////
  char filename1[50], filename2[50], filename3[50], filename4[50]; //filename5[50];
  sprintf(filename1, "Box1_T=%.2Lf_Rho=%.2Lf.dump", Temp, rho);
  sprintf(filename2, "Box2_T=%.2Lf_Rho=%.2Lf.dump", Temp, rho);
  sprintf(filename3, "Rho1vsMCstep_T=%.2Lf_Rho=%.2Lf.dat", Temp, rho);
  sprintf(filename4, "Rho2vsMCstep_T=%.2Lf_Rho=%.2Lf.dat", Temp, rho);
  //sprintf(filename5, "EvsMCstep_T%.1f_Rho=%.1f.dat", Temp, rho);
  //sprintf(filename5, "exchange.dat");
  
  FILE* fp1 = fopen(filename1, "w");
  FILE* fp2 = fopen(filename2, "w");
  FILE* fp3 = fopen(filename3, "w");
  FILE* fp4 = fopen(filename4, "w");
  //FILE* fp5 = fopen(filename5, "w");
  
  fprintf(fp1, "ITEM: TIMESTEP\n");
  fprintf(fp1, "0\n");
  fprintf(fp1, "ITEM: NUMBER OF ATOMS\n%d\n", n1);
  fprintf(fp1, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(fp1, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L2/2, L1/2, -L2/2, L1/2, -L2/2, L1/2);
  fprintf(fp1, "ITEM: ATOMS id type x y z xu yu zu\n");
  
  fprintf(fp2, "ITEM: TIMESTEP\n");
  fprintf(fp2, "0\n");
  fprintf(fp2, "ITEM: NUMBER OF ATOMS\n%d\n", n2);
  fprintf(fp2, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(fp2, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L2/2, L2/2, -L2/2, L2/2, -L2/2, L2/2);
  fprintf(fp2, "ITEM: ATOMS id type x y z xu yu zu\n");
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //INITIALIZATION
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  long double u;			//Random number variable
  long double x1[N], y1[N], z1[N];	//Box 1 particle arrays
  long double x1u[N], y1u[N], z1u[N];//Box 1 particle arrays(PBC)
  long double x2[N], y2[N], z2[N];	//Box 2 particle arrays
  long double x2u[N], y2u[N], z2u[N];//Box 2 particle arrays(PBC)
  
  ///////
  //BOX 1
  ///////
  for (int i = 0; i < n1; i++){
    //Assign the coordinate
    u = gsl_rng_uniform(r); 	//Random number between 0 and 1
    x1[i] = L1 * (u - 0.5);
    u = gsl_rng_uniform(r);
    y1[i] = L1 * (u - 0.5);
    u = gsl_rng_uniform(r);
    z1[i] = L1 * (u - 0.5);
    //Reject this coordinate if PBC distance between any 2 particles < 1 
    for(int j = 0; j < i && j != i; j++){
      long double dx = x1[i] - x1[j];
      long double dy = y1[i] - y1[j];
      long double dz = z1[i] - z1[j];
      //PBC
      dx = dx - lround(dx/L1) * L1;
      dy = dy - lround(dy/L1) * L1;
      dz = dz - lround(dz/L1) * L1;
      long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
      //printf("rij: %Lf\n",rij2);
      if(rij2 <= 1){
	i--;
	break;
      }
    }
    //PBC coordinates
    x1u[i] = x1[i] - lround(x1[i]/L1) * l;
    y1u[i] = y1[i] - lround(y1[i]/L1) * l;
    z1u[i] = z1[i] - lround(z1[i]/L1) * l;
  }
  //VISUALIZATION OF INITIAL CONFIGURATION
  for (int i = 0; i < n1; i++){ 
    fprintf(fp1, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
    //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
  }
  //Total Energy of initial box 1 
  U = 0;
  for (int i = 0; i < n1; i++){
    for (int j = 0; j < n1; j++){
      long double dx = x1[i] - x1[j];
      long double dy = y1[i] - y1[j];
      long double dz = z1[i] - z1[j];
      //PBC
      dx = dx - lround(dx/L1) * L1;
      dy = dy - lround(dy/L1) * L1;
      dz = dz - lround(dz/L1) * L1;
      long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
      if(i == j){
        U = U;
      }
      else{
        U += LJ(pow(rij2, 0.5));
      }
    }
  }
  printf("Energy for box 1 initial configuration: %Lf\n", U);
  
  ///////
  //BOX 2
  ///////
  for (int i = 0; i < n2; i++){
    //Assign the coordinates
    u = gsl_rng_uniform(r); 	//Random number between 0 and 1
    x2[i] = L2 * (u - 0.5);
    u = gsl_rng_uniform(r);
    y2[i] = L2 * (u - 0.5);
    u = gsl_rng_uniform(r);
    z2[i] = L2 * (u - 0.5);
    //Reject this coordinate if distance between any 2 particles < 1 
    for(int j = 0; j < i && i>0; j++){
      long double dx = x2[i] - x2[j];
      long double dy = y2[i] - y2[j];
      long double dz = z2[i] - z2[j];
      //PBC
      dx = dx - lround(dx/L2) * L2;
      dy = dy - lround(dy/L2) * L2;
      dz = dz - lround(dz/L2) * L2;
      long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
      if(rij2 < 1){
	i--;
	break;
      }
    }
    //PBC coordinates
    x2u[i] = x2[i] - lround(x2[i]/L2) * l;
    y2u[i] = y2[i] - lround(y2[i]/L2) * l;
    z2u[i] = z2[i] - lround(z2[i]/L2) * l;
  }
  //VISUALIZATION OF INITIAL CONFIGURATION
  for (int i = 0; i < n2; i++){ 
    fprintf(fp2, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
    //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);    
  }
  //Total energy of initial box 2 
  U = 0;
  for (int i = 0; i < n2; i++){
    for (int j = 0; j < n2; j++){
      long double dx = x2[i] - x2[j];
      long double dy = y2[i] - y2[j];
      long double dz = z2[i] - z2[j];
      //PBC
      dx = dx - lround(dx/L2) * L2;
      dy = dy - lround(dy/L2) * L2;
      dz = dz - lround(dz/L2) * L2;
      long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
      if(i == j){
        U = U;
      }
      else{
        U += LJ(pow(rij2, 0.5));
      }
    }
  }
  printf("Energy for box 2 initial configuration: %Lf\n", U);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
  //EVOLUTION USING MC
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  int MCStep, iDispTot = 0;
  int ireject = 0;		//Displacement rejection rate
  int iVolTot = 0;		//Total Volume moves
  int irejectVol = 0; 		//Volume rejection rate
  int iExTot = 0;		//Total particle exchange moves
  int irejectEx = 0;		//Particle exchange rejection rate
  for (MCStep = 0; MCStep < MCTot; MCStep++){
    long double choose = gsl_rng_uniform(r);		//Random number between 0 and 1
    //Select with 80-10-10 Probability
    
    //////////////////////////////////////////////////////////
    //Particle Displacement Move
    //////////////////////////////////////////////////////////
    if (choose < 0.8){
      iDispTot++;
      printf("\nDisplacement Step(MC Step: %d)\n", MCStep);
      ///////
      ///////
      ///////
      //Box1
      ///////
      ///////
      ///////
      for (int trial = 0; trial < n1; trial++){
        //pick a particle randomly
        u = gsl_rng_uniform (r);
        int randPart = u * n1;
        //int randPart = 100;
        //printf("\nParticle selection: %d\n", randPart);
        
        //Energy before displacement
        long double eo = 0;
        for (int i = 0; i < n1; i++){
          long double dx = x1[i] - x1[randPart];
          long double dy = y1[i] - y1[randPart];
          long double dz = z1[i] - z1[randPart];
          //PBC
          dx = dx - lround(dx/L1) * L1;
          dy = dy - lround(dy/L1) * L1;
          dz = dz - lround(dz/L1) * L1;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == randPart){
            eo = eo;
          }
          else{
            eo += LJ(pow(rij2, 0.5));
          }
        }
        //printf("The energy of box 1 particle %d before displacement: %Lf\n", randPart, eo);
        //displace the selected particle in a random direction
        iDispTot++;
        long double xold = x1[randPart];
        long double yold = y1[randPart];
        long double zold = z1[randPart];
        u = gsl_rng_uniform (r);
        x1[randPart] = x1[randPart] + l * (u-0.5);
        u = gsl_rng_uniform (r);
        y1[randPart] = y1[randPart] + l * (u-0.5);
        u = gsl_rng_uniform (r);
        z1[randPart] = z1[randPart] + l * (u-0.5);
        //Reject if particles are overlapping
        long double rij2;
        for(int j = 0; j < n1; j++){
          int i = randPart;
          long double dx = x1[i] - x1[j];
          long double dy = y1[i] - y1[j];
          long double dz = z1[i] - z1[j];
          //PBC
          dx = dx - lround(dx/L1) * L1;
          dy = dy - lround(dy/L1) * L1;
          dz = dz - lround(dz/L1) * L1;
          rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if (j == randPart){
            rij2 = 2;
          }
          if(rij2 <= 1){
            x1[randPart] = xold;
            y1[randPart] = yold;
            z1[randPart] = zold;
	    break;
          }
        }
        if (rij2 <= 1){
          ireject++;
          continue;
        }
        //Energy after displacement
        long double en = 0;
        for (int i = 0; i < n1; i++){
          long double dx = x1[i] - x1[randPart];
          long double dy = y1[i] - y1[randPart];
          long double dz = z1[i] - z1[randPart];
          //PBC
          dx = dx - lround(dx/L1) * L1;
          dy = dy - lround(dy/L1) * L1;
          dz = dz - lround(dz/L1) * L1;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == randPart){
            en = en;
          }
          else{
            en += LJ(pow(rij2, 0.5));
          }
        }
        //printf("The energy of box 1 particle %d after displacement: %Lf\n", randPart, en);
        
        //Accept or reject
	long double arg = en - eo;
        long double Bolt = exp(-arg/Temp);
        u = gsl_rng_uniform (r);
        //Reject if acceptance criterion is not satisfied
        if (Bolt < u){
          x1[randPart] = xold; y1[randPart] = yold; z1[randPart] = zold;
          ireject++;
        }
      }
      //VISUALIZATION OF NEW CONFIGURATION
      if(MCStep % 1 == 0){
        fprintf(fp1, "ITEM: TIMESTEP\n");
  	fprintf(fp1, "%d\n", MCStep+1);
  	fprintf(fp1, "ITEM: NUMBER OF ATOMS\n%d\n", n1);
  	fprintf(fp1, "ITEM: BOX BOUNDS pp pp pp\n");
  	fprintf(fp1, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L1/2, L1/2, -L1/2, L1/2, -L1/2, L1/2);
	fprintf(fp1, "ITEM: ATOMS id type x y z xu yu zu\n");
      }
      //PBC coordinates after the trial moves
      for (int i = 0; i < n1; i++){
        x1u[i] = x1[i] - lround(x1[i]/L1) * L1;
        y1u[i] = y1[i] - lround(y1[i]/L1) * L1;
        z1u[i] = z1[i] - lround(z1[i]/L1) * L1;
        if(MCStep % 1 == 0){
          fprintf(fp1, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
          //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
        }
      }
      long double rho1 = n1/pow(L1,3);
      fprintf(fp3,"%d %Lf\n", MCStep, rho1);
      ///////
      ///////
      ///////
      //Box2
      ///////
      ///////
      ///////
      for (int trial = 0; trial < n2; trial++){
        //pick a particle randomly
        u = gsl_rng_uniform (r);
        int randPart = u * n2;
        //int randPart = 100;
        //printf("\nParticle selection: %d\n", randPart);
        
        //Energy before displacement
        long double eo = 0;
        for (int i = 0; i < n2; i++){
          long double dx = x2[i] - x2[randPart];
          long double dy = y2[i] - y2[randPart];
          long double dz = z2[i] - z2[randPart];
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == randPart){
            eo = eo;
          }
          else{
            eo += LJ(pow(rij2, 0.5));
          }
        }
        //printf("The energy of box 2 particle %d before displacement: %Lf\n", randPart, eo);
        //displace the selected particle in a random direction
        iDispTot++;
        long double xold = x2[randPart];
        long double yold = y2[randPart];
        long double zold = z2[randPart];
        u = gsl_rng_uniform (r);
        x2[randPart] = x2[randPart] + l * (u-0.5);
        u = gsl_rng_uniform (r);
        y2[randPart] = y2[randPart] + l * (u-0.5);
        u = gsl_rng_uniform (r);
        z2[randPart] = z2[randPart] + l * (u-0.5);
        //Reject if particles are overlapping
        long double rij2;
        for(int j = 0; j < n2; j++){
          int i = randPart;
          long double dx = x2[i] - x2[j];
          long double dy = y2[i] - y2[j];
          long double dz = z2[i] - z2[j];
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if (j == randPart){
            rij2 = 2;
          }
          if(rij2 <= 1){
            x2[randPart] = xold;
            y2[randPart] = yold;
            z2[randPart] = zold;
	    break;
          }
        }
        if (rij2 <= 1){
          ireject++;
          continue;
        }
        //Energy after displacement
        long double en = 0;
        for (int i = 0; i < n2; i++){
          long double dx = x2[i] - x2[randPart];
          long double dy = y2[i] - y2[randPart];
          long double dz = z2[i] - z2[randPart];
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == randPart){
            en = en;
          }
          else{
            en += LJ(pow(rij2, 0.5));
          }
        }
        //printf("The energy of box 2 particle %d after displacement: %Lf\n", randPart, en);
        
        //Accept or reject
	long double arg = en - eo;
        long double Bolt = exp(-arg/Temp);
        u = gsl_rng_uniform (r);
        //Reject if acceptance criterion is not satisfied
        if (Bolt < u){
          x2[randPart] = xold; y2[randPart] = yold; z2[randPart] = zold;
          ireject++;
        }
      }
      //VISUALIZATION OF NEW CONFIGURATION
      if(MCStep % 1 == 0){
        fprintf(fp2, "ITEM: TIMESTEP\n");
  	fprintf(fp2, "%d\n", MCStep+1);
  	fprintf(fp2, "ITEM: NUMBER OF ATOMS\n%d\n", n2);
  	fprintf(fp2, "ITEM: BOX BOUNDS pp pp pp\n");
  	fprintf(fp2, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L2/2, L2/2, -L2/2, L2/2, -L2/2, L2/2);
	fprintf(fp2, "ITEM: ATOMS id type x y z xu yu zu\n");
      }
      //PBC coordinates after the trial moves
      for (int i = 0; i < n2; i++){
        x2u[i] = x2[i] - lround(x2[i]/L2) * L2;
        y2u[i] = y2[i] - lround(y2[i]/L2) * L2;
        z2u[i] = z2[i] - lround(z2[i]/L2) * L2;
        if(MCStep % 1 == 0){
          fprintf(fp2, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
          //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
        }
      }
      long double rho2 = n2/pow(L2,3);
      fprintf(fp4,"%d %Lf\n", MCStep, rho2);
      //printf("\nRejection ratio in Displacement change (rejection(%d)/Total Displacement Moves(%d)): %Lf\n", ireject, iDispTot, (long double)(ireject)/(long double)(iDispTot));  
    }
    
    
    //////////////////////////////////////////////////////////////////////
    //Volume change move
    /////////////////////////////////////////////////////////////////////
    else if (choose >= 0.8 && choose < 0.9){
      iVolTot++;
      printf("\nVolume Step(MC Step: %d)\n", MCStep + 1);
      printf("\nRejection ratio in volume change (rejection(%d)/Total Volume Moves(%d)): %Lf\n", irejectVol, iVolTot, (long double)(irejectVol)/(long double)(iVolTot));
      long double L1o = L1;
      //Energy of box 1 before volume change
      long double eo = 0;
      for (int i = 0; i < n1; i++){
        for (int j = 0; j < n1; j++){
          long double dx = x1[i] - x1[j];
          long double dy = y1[i] - y1[j];
          long double dz = z1[i] - z1[j];
          //PBC
          dx = dx - lround(dx/L1) * L1;
          dy = dy - lround(dy/L1) * L1;
          dz = dz - lround(dz/L1) * L1;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == j){
            eo = eo;
          }
          else{
            eo += LJ(pow(rij2, 0.5));
          }
        }
      }
      //change the volume randomly
      u = (gsl_rng_uniform(r));
      L1 = L1 + (u - 0.5) * l1;
      V1 = pow(L1, 3);
      //scale the coordinates accordingly
      long double xold[N], yold[N], zold[N];
      for (int i = 0; i < n1; i++){
        //printf("i: %d x1o: %Lf\n", i, x1[i]);
        xold[i] = x1[i];
        yold[i] = y1[i];
        zold[i] = z1[i];
        x1[i] = x1[i]*L1/L1o;
        y1[i] = y1[i]*L1/L1o;
        z1[i] = z1[i]*L1/L1o;
        //printf("i: %d x1n: %Lf\n", i, x1[i]);
      }
      //Reject if particles are overlapping
      long double rij2;
      for (int i = 0; i < n1; i++){
        for(int j = 0; j < n1; j++){
          long double dx = x1[i] - x1[j];
          long double dy = y1[i] - y1[j];
          long double dz = z1[i] - z1[j];
          //PBC
          dx = dx - lround(dx/L1) * L1;
          dy = dy - lround(dy/L1) * L1;
          dz = dz - lround(dz/L1) * L1;
          rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if (j == i){
            rij2 = 2;
          }
          if(rij2 <= 1){
	    break;
          }
        }
        if(rij2 <= 1){
          break;
        }
      }
      if (rij2 <= 1){
        for(int i = 0; i < n1; i++){
          x1[i] = xold[i];
          y1[i] = yold[i];
          z1[i] = zold[i];
        }
        irejectVol++;
        L1 = L1o;
        V1 = pow(L1,3);
        //printf("Box 1 volume move rejected, irejectVol: %d\n", irejectVol);
        continue;
      }
      else{
        //printf("Box 1 volume move accepted, irejectVol: %d\n", irejectVol);
      }
      //Energy of box 1 after volume change
      long double en = 0;
      for (int i = 0; i < n1; i++){
        for (int j = 0; j < n1; j++){
          long double dx = x1[i] - x1[j];
          long double dy = y1[i] - y1[j];
          long double dz = z1[i] - z1[j];
          //PBC
          dx = dx - lround(dx/L1) * L1;
          dy = dy - lround(dy/L1) * L1;
          dz = dz - lround(dz/L1) * L1;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == j){
            en = en;
          }
          else{
            en += LJ(pow(rij2, 0.5));
          }
        }
      }
      long double L2o;
      //Add box 1 Energy before volume change to box 2 before volume change
      for (int i = 0; i < n2; i++){
        for (int j = 0; j < n2; j++){
          long double dx = x2[i] - x2[j];
          long double dy = y2[i] - y2[j];
          long double dz = z2[i] - z2[j];
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == j){
            eo = eo;
          }
          else{
            eo += LJ(pow(rij2, 0.5));
          }
        }
      }
      //change the volume of box 2 such that total volume is the same
      L2o = L2;
      V2o = pow(L2o,3);
      L2 = pow(pow(L,3)-pow(L1,3), 0.33333);
      V2 = pow(L2, 3);
      //printf("V1: %Lf,  V2: %Lf, V = V1 + V2: %Lf, Given V: %Lf\n", V1, V2, V1+V2, V);
      //scale the coordinates accordingly
      long double xold2[N], yold2[N], zold2[N];
      for (int i = 0; i < n2; i++){
        //printf("i: %d x2o: %Lf\n", i, x2[i]);
        xold2[i] = x2[i];
        yold2[i] = y2[i];
        zold2[i] = z2[i];
        x2[i] = x2[i]*L2/L2o;
        y2[i] = y2[i]*L2/L2o;
        z2[i] = z2[i]*L2/L2o;
        //printf("i: %d x2n: %Lf\n", i, x2[i]);
      }
      //Reject if particles are overlapping
      for (int i = 0; i < n2; i++){
        for(int j = 0; j < n2; j++){
          long double dx = x2[i] - x2[j];
          long double dy = y2[i] - y2[j];
          long double dz = z2[i] - z2[j];
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if (j == i){
            rij2 = 2;
          }
          if(rij2 <= 1){
	    break;
          }
        }
        if(rij2 <= 1){
          for (int i = 0; i < n2; i++){ 
            x2[i] = xold2[i];
            y2[i] = yold2[i];
            z2[i] = zold2[i];
          }
          L2 = L2o;
          V2 = pow(L2,3);
	  irejectVol++;
	  //change box 1 coordinates as well
	  L1 = L1o;
	  V1 = pow(L1,3);
	  for(int k = 0; k < n1; k++){
	    x1[k] = xold[k];
            y1[k] = yold[k];
            z1[k] = zold[k];
	  }
          break;
        }
      }
      if (rij2 <= 1){
        //printf("Box 2 volume move rejected, irejectVol: %d\n", irejectVol);
        continue;
      }
      else{
        //printf("No overlap. Box 1 volume move Accepted, irejectVol: %d\n", irejectVol);
      }
      //Total energy of the system after volume change
      //en = 0;
      for (int i = 0; i < n2; i++){
        for (int j = 0; j < n2; j++){
          long double dx = x2[i] - x2[j];
          long double dy = y2[i] - y2[j];
          long double dz = z2[i] - z2[j];
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == j){
            en = en;
          }
          else{
            en += LJ(pow(rij2, 0.5));
          }
        }
      }
      //printf("eo/n2: %Lf en/n2: %Lf\n", eo/2/n2, en/2/n2);
      //////////////////
      //Accept or reject
      //////////////////
      long double arg = en - eo;
      //long double Bolt = (pow(V2,n2)*pow(V1, n1))/(pow(V2o,n2)*pow(V-V2o, n1)) * exp(-arg/Temp);
      long double Bolt = pow((V1/(V-V2o)),(n1+1))*pow((V2/V2o),(n2+1)) * exp(-arg/Temp);
      long double u = gsl_rng_uniform (r);
      //Reject if acceptance criterion is not satisfied
      if (Bolt < u){
        for (int i = 0; i < n2; i++){
          x2[i] = xold2[i]; y2[i] = yold2[i]; z2[i] = zold2[i];
          L2 = L2o;
          V2 = pow(L2,2);
          //change box 1 coordinates as well
	  L1 = L1o;
	  for(int k = 0; k < n1; k++){
	    x1[k] = xold[k];
            y1[k] = yold[k];
            z1[k] = zold[k];
	  }
        }
        irejectVol++;
        //printf("Volume move Bolt rejected, irejectVol: %d\n", irejectVol);
      }
      else{
        //printf("Volume move Bolt accepted, irejectVol: %d\n", irejectVol);
        ////////////////////////////////////
        //VISUALIZATION OF NEW CONFIGURATION
        ////////////////////////////////////
        //Box 1
        if(MCStep % 1 == 0){
          fprintf(fp1, "ITEM: TIMESTEP\n");
  	  fprintf(fp1, "%d\n", MCStep+1);
  	  fprintf(fp1, "ITEM: NUMBER OF ATOMS\n%d\n", n1);
  	  fprintf(fp1, "ITEM: BOX BOUNDS pp pp pp\n");
  	  fprintf(fp1, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L1/2, L1/2, -L1/2, L1/2, -L1/2, L1/2);
	  fprintf(fp1, "ITEM: ATOMS id type x y z xu yu zu\n");
        }
        //PBC coordinates after the trial moves
        for (int i = 0; i < n1; i++){
          x1u[i] = x1[i] - lround(x1[i]/L1) * L1;
          y1u[i] = y1[i] - lround(y1[i]/L1) * L1;
          z1u[i] = z1[i] - lround(z1[i]/L1) * L1;
          if(MCStep % 1 == 0){
            fprintf(fp1, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
            //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
          }
        }
        //Box2
        if(MCStep % 1 == 0){
          fprintf(fp2, "ITEM: TIMESTEP\n");
  	  fprintf(fp2, "%d\n", MCStep+2);
  	  fprintf(fp2, "ITEM: NUMBER OF ATOMS\n%d\n", n2);
  	  fprintf(fp2, "ITEM: BOX BOUNDS pp pp pp\n");
  	  fprintf(fp2, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L2/2, L2/2, -L2/2, L2/2, -L2/2, L2/2);
	  fprintf(fp2, "ITEM: ATOMS id type x y z xu yu zu\n");
        }
        //PBC coordinates after the trial moves
        for (int i = 0; i < n2; i++){
          x2u[i] = x2[i] - lround(x2[i]/L2) * L2;
          y2u[i] = y2[i] - lround(y2[i]/L2) * L2;
          z2u[i] = z2[i] - lround(z2[i]/L2) * L2;
          if(MCStep % 1 == 0){
            fprintf(fp2, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
            //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
          }
        }
        long double rho1 = n1/pow(L1,3);
        fprintf(fp3,"%d %Lf\n", MCStep, rho1);
        long double rho2 = n2/pow(L2,3);
        fprintf(fp4,"%d %Lf\n", MCStep, rho2);
      }
    }
    else if(choose >= 0.9){
      //////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////
      //Particle exchange move
      /////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////
      printf("\nParticle Exchange Move(MC Step: %d)\n", MCStep);
      printf("\nRejection ratio in particle exchange (rejection(%d)/Total Particle Exchange Moves(%d)): %Lf\n", irejectEx, iExTot, (long double)(irejectEx)/(long double)(iExTot));
      iExTot++;
      int randPart;
      //Box 1 Energy before particle exchange
      long double eo = 0;
      for (int i = 0; i < n1; i++){
        for (int j = 0; j < n1; j++){
          long double dx = x1[i] - x1[j];
          long double dy = y1[i] - y1[j];
          long double dz = z1[i] - z1[j];
          //PBC
          dx = dx - lround(dx/L1) * L1;
          dy = dy - lround(dy/L1) * L1;
          dz = dz - lround(dz/L1) * L1;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == j){
            eo = eo;
          }
          else{
            eo += LJ(pow(rij2, 0.5));
          }
        }
      }
      //Add box Energy before particle exchange to box 2 before particle exchange
      for (int i = 0; i < n2; i++){
        for (int j = 0; j < n2; j++){
          long double dx = x2[i] - x2[j];
          long double dy = y2[i] - y2[j];
          long double dz = z2[i] - z2[j];
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(i == j){
            eo = eo;
          }
          else{
            eo += LJ(pow(rij2, 0.5));
          }
        }
      }
      //Pick a particle in system randomly and choose a position for it to be inserted inside other box 
      u = (gsl_rng_uniform(r));
      long double ux, uy, uz;
      long double xold[N], yold[N], zold[N];
      long double xRand, yRand, zRand;
      //Box 1 to Box 2
      if(u < 1){
        //pick a random particle in box 1
        long double u1 = (gsl_rng_uniform(r));
        randPart = u1 * (n1-1);
        //pick a position in the box 2 randomly
        ux = (gsl_rng_uniform(r));
        xRand = L2*(ux-0.5);
        uy = (gsl_rng_uniform(r));
        yRand = L2*(uy-0.5);
        uz = (gsl_rng_uniform(r));
        zRand = L2*(uz-0.5);
        //First check for overlap in box 2
        long double rij2 = 0;
        for(int i = 0; i < n2; i++){
          long double dx = x2[i] - xRand;
          long double dy = y2[i] - yRand;
          long double dz = z2[i] - zRand;
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(rij2 <= 1){
            break;
          }
        }
        if(rij2 <= 1){
          irejectEx++;
          //printf("Particle Exchange move overlap rejected, irejectEx: %d\n", irejectEx);
          continue;
        }
        //If there is no overlap, ASSIGN the coordinates
        //Remove from box 1 and add to box 2
        for(int i = 0; i < n1; i++){
          xold[i] = x1[i];
          yold[i] = y1[i];
          zold[i] = z1[i];
        }
        for(int i = randPart; i < n1; i++){
          x1[i] = x1[i+1];
          y1[i] = y1[i+1];
          z1[i] = z1[i+1];
          //printf("%d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, x1[i], y1[i], z1[i], x1[i+1], y1[i+1], z1[i+1]);
        }
        x2[n2] = xRand;
        y2[n2] = yRand;
        z2[n2] = zRand;
        n1--;
        n2++;
        //Energy after particle exchange
        long double en = 0;
        for (int i = 0; i < n1; i++){
          for (int j = 0; j < n1; j++){
            long double dx = x1[i] - x1[j];
            long double dy = y1[i] - y1[j];
            long double dz = z1[i] - z1[j];
            //PBC
            dx = dx - lround(dx/L2) * L2;
            dy = dy - lround(dy/L2) * L2;
            dz = dz - lround(dz/L2) * L2;
            long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
            if(i == j){
              en = en;
            }
            else{
              en += LJ(pow(rij2, 0.5));
            }
          }
        }
        for (int i = 0; i < n2; i++){
          for (int j = 0; j < n2; j++){
            long double dx = x2[i] - x2[j];
            long double dy = y2[i] - y2[j];
            long double dz = z2[i] - z2[j];
            //PBC
            dx = dx - lround(dx/L2) * L2;
            dy = dy - lround(dy/L2) * L2;
            dz = dz - lround(dz/L2) * L2;
            long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
            long double rij = pow(rij2,0.5);
            if (i==79 && j == n2-1){
              printf("i:%d j:%d rij: %Lf\n", i, j, rij);
            }
            if(i == j){
              en = en;
            }
            else{
              en += LJ(pow(rij2, 0.5));
            }
            //fprintf(fp5 ,"i: %d j: %d en: %Lf\n", i, j, en);
          }
        }
        //Accept or Reject
        long double arg = en - eo;
        long double Bolt = (n1 * V2)/((n2+1)*V1) * exp(-arg/Temp);
        long double u = gsl_rng_uniform (r);
        //Reject if acceptance criterion is not satisfied
        if (Bolt < u){
          for(int i = 0; i < n1; i++){
            x1[i] = xold[i];
            y1[i] = yold[i];
            z1[i] = zold[i];
          }
          n1++;
          n2--;
          irejectEx++;
          //printf("Particle Exchange move Bolt rejected, irejectEx: %d\n", irejectEx);
          continue;
        }
        else{
          //printf("Particle Exchange move Bolt accepted, irejectEx: %d\n", irejectEx);
          ////////////////////////////////////
          //VISUALIZATION OF NEW CONFIGURATION
          ////////////////////////////////////
          //Box 1
          if(MCStep % 1 == 0){
            fprintf(fp1, "ITEM: TIMESTEP\n");
  	    fprintf(fp1, "%d\n", MCStep+1);
  	    fprintf(fp1, "ITEM: NUMBER OF ATOMS\n%d\n", n1);
  	    fprintf(fp1, "ITEM: BOX BOUNDS pp pp pp\n");
  	    fprintf(fp1, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L1/2, L1/2, -L1/2, L1/2, -L1/2, L1/2);
	    fprintf(fp1, "ITEM: ATOMS id type x y z xu yu zu\n");
          }
          //PBC coordinates after the trial moves
          for (int i = 0; i < n1; i++){
            x1u[i] = x1[i] - lround(x1[i]/L1) * L1;
            y1u[i] = y1[i] - lround(y1[i]/L1) * L1;
            z1u[i] = z1[i] - lround(z1[i]/L1) * L1;
            if(MCStep % 1 == 0){
              fprintf(fp1, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
              //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
            }
          }
          //Box2
          if(MCStep % 1 == 0){
            fprintf(fp2, "ITEM: TIMESTEP\n");
  	    fprintf(fp2, "%d\n", MCStep+2);
  	    fprintf(fp2, "ITEM: NUMBER OF ATOMS\n%d\n", n2);
  	    fprintf(fp2, "ITEM: BOX BOUNDS pp pp pp\n");
  	    fprintf(fp2, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L2/2, L2/2, -L2/2, L2/2, -L2/2, L2/2);
	    fprintf(fp2, "ITEM: ATOMS id type x y z xu yu zu\n");
          }
          //PBC coordinates after the trial moves
          for (int i = 0; i < n2; i++){
            x2u[i] = x2[i] - lround(x2[i]/L2) * L2;
            y2u[i] = y2[i] - lround(y2[i]/L2) * L2;
            z2u[i] = z2[i] - lround(z2[i]/L2) * L2;
            if(MCStep % 1 == 0){
              fprintf(fp2, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
              //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
            }
          }
        }
      }
      //Box 2 to Box 1
      else {
        //pick a random particle box 2
        long double u1 = (gsl_rng_uniform(r));
        randPart = u1 * (n2-1);
        //pick a position in the box 1 randomly
        ux = (gsl_rng_uniform(r));
        xRand = L1*(ux-0.5);
        uy = (gsl_rng_uniform(r));
        yRand = L1*(uy-0.5);
        uz = (gsl_rng_uniform(r));
        zRand = L1*(uz-0.5);
        //First check for overlap in box 1
        long double rij2;
        for(int i = 0; i < n1; i++){
          long double dx = x1[i] - xRand;
          long double dy = y1[i] - yRand;
          long double dz = z1[i] - zRand;
          //PBC
          dx = dx - lround(dx/L2) * L2;
          dy = dy - lround(dy/L2) * L2;
          dz = dz - lround(dz/L2) * L2;
          rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
          if(rij2 <= 1){
            break;
          }
        }
        if(rij2 <= 1){
          irejectEx++;
          continue;
        }
        //If there is no overlap, assign the coordinate
        //Remove from box 2 and add to box 1
        for(int i = 0; i < n2; i++){
          xold[i] = x2[i];
          yold[i] = y2[i];
          zold[i] = z2[i];
        }
        for(int i = randPart; i < n2; i++){
          x2[i] = x2[i+1];
          y2[i] = y2[i+1];
          z2[i] = z2[i+1];
        }
        x1[n1] = xRand;
        y1[n1] = yRand;
        z1[n1] = zRand;
        n2--; 
        n1++;
        //Energy after particle exchange
        long double en = 0;
        for (int i = 0; i < n1; i++){
          for (int j = 0; j < n1; j++){
            long double dx = x1[i] - x1[j];
            long double dy = y1[i] - y1[j];
            long double dz = z1[i] - z1[j];
            //PBC
            dx = dx - lround(dx/L2) * L2;
            dy = dy - lround(dy/L2) * L2;
            dz = dz - lround(dz/L2) * L2;
            long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
            if(i == j){
              en = en;
            }
            else{
              en += LJ(pow(rij2, 0.5));
            }
          }
        }
        for (int i = 0; i < n2; i++){
          for (int j = 0; j < n2; j++){
            long double dx = x2[i] - x2[j];
            long double dy = y2[i] - y2[j];
            long double dz = z2[i] - z2[j];
            //PBC
            dx = dx - lround(dx/L2) * L2;
            dy = dy - lround(dy/L2) * L2;
            dz = dz - lround(dz/L2) * L2;
            long double rij2 = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
            if(i == j){
              en = en;
            }
            else{
              en += LJ(pow(rij2, 0.5));
            }
          }
        }
        //Accept or Reject
        long double arg = en - eo;
        long double Bolt = (n2 * V1)/((n1+1)*V2) * exp(-arg/Temp);
        long double u = gsl_rng_uniform (r);
        //Reject if acceptance criterion is not satisfied
        if (Bolt < u){
          for(int i = 0; i < n1; i++){
            x2[i] = xold[i];
            y2[i] = yold[i];
            z2[i] = zold[i];
          }
          n1++;
          n2--;
          irejectEx++;
          //printf("Particle Exchange move Bolt rejected, irejectEx: %d\n", irejectEx);
          continue;
        }
        else{
          //printf("Particle Exchange move Bolt accepted, irejectEx: %d\n", irejectEx);
          ////////////////////////////////////
          //VISUALIZATION OF NEW CONFIGURATION
          ////////////////////////////////////
          //Box 1
          if(MCStep % 1 == 0){
            fprintf(fp1, "ITEM: TIMESTEP\n");
  	    fprintf(fp1, "%d\n", MCStep+1);
  	    fprintf(fp1, "ITEM: NUMBER OF ATOMS\n%d\n", n1);
  	    fprintf(fp1, "ITEM: BOX BOUNDS pp pp pp\n");
  	    fprintf(fp1, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L1/2, L1/2, -L1/2, L1/2, -L1/2, L1/2);
	    fprintf(fp1, "ITEM: ATOMS id type x y z xu yu zu\n");
          }
          //PBC coordinates after the trial moves
          for (int i = 0; i < n1; i++){
            x1u[i] = x1[i] - lround(x1[i]/L1) * L1;
            y1u[i] = y1[i] - lround(y1[i]/L1) * L1;
            z1u[i] = z1[i] - lround(z1[i]/L1) * L1;
            if(MCStep % 1 == 0){
              fprintf(fp1, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
              //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x1[i], y1[i], z1[i], x1u[i], y1u[i], z1u[i]);
            }
          }
          //Box2
          if(MCStep % 1 == 0){
            fprintf(fp2, "ITEM: TIMESTEP\n");
  	    fprintf(fp2, "%d\n", MCStep+2);
  	    fprintf(fp2, "ITEM: NUMBER OF ATOMS\n%d\n", n2);
  	    fprintf(fp2, "ITEM: BOX BOUNDS pp pp pp\n");
  	    fprintf(fp2, "%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n", -L2/2, L2/2, -L2/2, L2/2, -L2/2, L2/2);
	    fprintf(fp2, "ITEM: ATOMS id type x y z xu yu zu\n");
          }
          //PBC coordinates after the trial moves
          for (int i = 0; i < n2; i++){
            x2u[i] = x2[i] - lround(x2[i]/L2) * L2;
            y2u[i] = y2[i] - lround(y2[i]/L2) * L2;
            z2u[i] = z2[i] - lround(z2[i]/L2) * L2;
            if(MCStep % 1 == 0){
              fprintf(fp2, "%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
              //printf("%d %d %Lf %Lf %Lf %Lf %Lf %Lf\n", i, 1, x2[i], y2[i], z2[i], x2u[i], y2u[i], z2u[i]);
            }
          }
          long double rho1 = n1/pow(L1,3);
          fprintf(fp3,"%d %Lf\n", MCStep, rho1);
          long double rho2 = n2/pow(L2,3);
          fprintf(fp4,"%d %Lf\n", MCStep, rho2);
        }
      }
    }
  }
  printf("\nRejection ratio in Displacement change (rejection(%d)/Total Displacement Moves(%d)): %Lf\n", ireject, iDispTot, (long double)(ireject)/(long double)(iDispTot)); 
  printf("\nRejection ratio in Volume change (rejection(%d)/Total Volume Moves(%d)): %Lf\n", irejectVol, iVolTot, (long double)(irejectVol)/(long double)(iVolTot));
  printf("\nRejection ratio in particle exchange (rejection(%d)/Total Particle Exchange Moves(%d)): %Lf\n", irejectEx, iExTot, (long double)(irejectEx)/(long double)(iExTot));
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //END
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  gsl_rng_free (r);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  //fclose(fp5);
  
  return 0;
}
