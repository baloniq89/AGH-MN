#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "four1.c"
#include "nrutil.c"
#include "nrutil.h"

#define T 1.0
#define T_MAX 3.0*T
#define omega 2.0*M_PI/T
#define sigma T/20.0

float f_rand(float min, float max);
float f0(float t);
float g0(float t);

int main() {
  srand(time(0));
  FILE *fp8 = fopen("file1_k8.dat","w");
  FILE *function = fopen("function.dat", "w");

  //int k = 8;
  //int k = 10;
  int k = 12;
  int N = pow(2,k);
  float dt = T_MAX / N;

  float* f = (float*) calloc(2*N, sizeof(float));
	float* g = (float*) calloc(2*N, sizeof(float));
	float* g1 = (float*) calloc(2*N, sizeof(float));
	float* g2 = (float*) calloc(2*N, sizeof(float));

  for (int i = 1; i < N - 1; ++i) {
			float t = (i-1) * dt;
			f[2*i-1] = f0(t) + f_rand(-0.5, 0.5);
			f[2*i] = 0.0; 

			g1[2*i-1] = g0(t);
			g1[2*i] = 0.0; 

			g2[2*i-1] = g0(t);
			g2[2*i ] = 0.0; 

			fprintf(fp8, "%f\t%f\n", t, f[2*i-1]);
      fprintf(function, "%f\t%f\n",t, f0(t));
		}
    
  four1(f, N, 1);
	four1(g1, N, 1);
	four1(g2, N, -1);


  float a1, a2, b1, b2;
	for (int i = 1; i < N - 1; ++i)
  {
			g[2*i-1] = g1[2*i] + g2[2*i];
			g[2*i] = g1[2*i + 1] + g2[2*i + 1];

			a1 = f[2*i-1];
			b1 = f[2*i];
			a2 = g[2*i-1];
			b2 = g[2*i];

			f[2*i-1] = a1 * a2 - b1 * b2;
			f[2*i] = a1 * b2 + a2 * b1;
	}

  four1(f,N,-1);

  float fmax = -99999.;
	for (int i = 0; i < N - 1; ++i)
			if (fabs(f[2*i]) > fmax) fmax = fabs(f[2*i]);
	
        //printf("f_max(k = %d) = %f\n", k, fmax);
		fprintf(fp8, "\n\n");
		for (int i = 0; i < N - 1; ++i)
			fprintf(fp8, "%f\t%f\n", dt*i, f[2*i]*2.5/fmax);
  printf("%f", fmax);
  free(f);
	free(g);
	free(g1);
	free(g2);
}

float f_rand(float min, float max){
	float r = (float) rand()/RAND_MAX;
	return r*(max-min) + min;
}

float f0(float t) {
	return sin(omega*t) + sin(2*omega*t) + sin(3*omega*t);
}

float g0(float t) {
	float result = 1.0/(sigma*sqrt(2*M_PI));
	result *= exp(-t*t/(2*sigma*sigma));
	return result;
}