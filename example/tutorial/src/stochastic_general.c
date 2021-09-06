#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "randomlib.h"

int main(int argc, char *argv[]){
  int length;
  float angle, time, sigma, angle2;
  float *H;
  float *Htemp;
  float *w;
  float deltat;
  float a, b;
  float pi = 3.14159265;
  int i, j, idx, N, nn2;
  FILE *E_FH,*mu_FH;

  length = atoi(argv[1]);
  deltat = atof(argv[2]);
  sigma = atof(argv[3]);
  time = atof(argv[4]);
  angle2 = atof(argv[5]);
  printf("Trajectory length: %i fs\n", length);
  printf("Timestep: %f fs\n", deltat);
  printf("Width: %f cm^-1\n", sigma);
  printf("Coherence time: %f fs\n", time);
  printf("Correlation angle: %f fs\n", angle2);

  nn2 = argc-6;
  N = floor(sqrt(2*nn2));
  printf("NN2: %i, N: %i\n", nn2, N);
  H = (float *)calloc(nn2, sizeof(float));
  Htemp = (float *)calloc(nn2, sizeof(float));
  w = (float *)calloc(N, sizeof(float));

  // Copy the other arguments to a new array
  for (i = 0; i < nn2; i++) {
    H[i] = atof(argv[i + 6]);
  }

  a = exp(-deltat/time);
  b = sqrt(1-a*a);

  // Set initial values of diagonal elements
  RandomInitialise(2511,1974);
  for (i = 0; i < N; i++) {
    w[i] = RandomGaussian(0, sigma);
  }

  E_FH = fopen("Energy.txt","w");
  mu_FH = fopen("Dipole.txt","w");
  
  for (j = 0; j < length; j++) {
    // Reset Htemp
    for (i = 0; i < nn2; i++) {
      Htemp[i] = H[i];
    }
    // Update diagonal elements
    for (i = 0; i < N; i++) {    
      w[i] *= a;
      w[i] += b*RandomGaussian(0, sigma);
    }
    // Assemble temporary Hamiltonian
    for (i = 0; i < N; i++) {
      idx = nn2 - (i+1)*(i+2)/2;
      Htemp[idx] += w[i];
    }
    // Print Hamiltonian to file
    fprintf(E_FH, "%d ", j);
    for (i = 0; i < nn2; i++) {
      fprintf(E_FH, "%f ", Htemp[i]);
    }
    fprintf(E_FH, "\n", 0);

    // Print dipole to file
    fprintf(mu_FH, "%d ", j);
    for (i = 0; i < N; i++) {
      fprintf(mu_FH, "%f %f %f ", 1.0, cos(angle2*pi/180.0), 0.0);
    }
    fprintf(mu_FH, "\n", 0);
  } 

  return 0;
}
