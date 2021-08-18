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
  int i, j, idx, N;
  FILE *E_FH,*mu_FH;

  length = atoi(argv[1]);
  deltat = atof(argv[2]);
  sigma = atof(argv[3]);
  time = atof(argv[4]);
  angle2 = atof(argv[5]);
  
  nn2 = argc-6;
  N = floor(sqrt(2*nn2));
  H = (float *)calloc(nn2, sizeof(float));
  w = (float *)calloc(N, sizeof(float));

  // Copy the other arguments to a new array
  for (i = 0; i < nn2; i++) {
    H[i] = atof(argv[i + 5]);
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
    for (i = 0; i < N; i++) {
      // Update diagonal elements
      w[i] *= a;
      w[i] += b*RandomGaussian(0, sigma);
    }
    // Assemble temporary Hamiltonian
    for (i = 0; i < N; i++) {
      idx = 1 + nn2 + (N-i)*(N-i+1)/2;
      Htemp[idx] += w[i];
    }

    // Print Hamiltonian to file
    fprintf(E_FH, "%d ", j);
    for (i = 0; i < nn2; i++) {
      fprintf(E_FH, "%f ", Htemp[i]);
    }
    fprintf(E_FH, "\n");

    // Print dipole to file
    fprintf(mu_FH, "%d ", j);
    for (i = 0; i < N; i++) {
      fprintf(mu_FH, "%f %f %f ", 1.0, cos(angle2*pi/180.0), 0.0);
    }
    fprintf("\n");

    // Reset Htemp
    for (i = 0; i < nn2; i++) {
      Htemp[i] = H[i];
    }
  } 

  return 0;
}
