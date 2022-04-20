#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nrutil.h"
#include "nrutil.c"
#include "gauher.c"
#include "gammln.c"
#include "gaulag.c"
#include "gauleg.c"




float fun1(float x) {
	return x/(4 * pow(x, 2) + 1);
}

float fun2(float x, int k) {
  return pow(x, k);
}

float fun3_2(float x) {
  return pow(sin(x),2);
}

float fun3_4(float x) {
  return pow(sin(x),4);
}

float factorial(int n)
{
	float result = 1;
	for(float i = 2; i<=n;++i) result*=i;
	
	return result;
}


int main(void) {
  printf("ZAD1\n");
  FILE * fp1 = fopen("roznica.dat","w");
  FILE * fp2 = fopen("roznica2_k5.dat","w");
  FILE * fp3 = fopen("roznica2_k10.dat","w");
  FILE * fp4 = fopen("roznica3.dat","w");

  FILE * fk1 = fopen("sum_kwadratur.dat","w");
  FILE * fk2 = fopen("sum_kwadratur.dat_k5","w");
  FILE * fk3 = fopen("sum_kwadratur_k10.dat","w");
  FILE * fk4 = fopen("sum_kwadratur3.dat","w");
  float c1a = (1./8.)*log(4.*2*2+1.)-(1/8.0)*log(1);
  
	for (int n = 2; n <= 20; ++n) {


	  float* x = vector(1, n);
		float* w = vector(1, n);

		float x1 = 0;
		float x2 = 2;

		gauleg(x1, x2, x-1, w-1, n);
    float suma_k = 0;

		float c1 = 0.0;
		for (int i = 1; i <= n; ++i) {

			c1 += w[i]*fun1(x[i]);
      
		}
    for(int i=0;i<n;i++)
    {
      suma_k +=w[i];
    }

		float err = fabs(c1 - c1a);

		//printf("%d	%f\n", n, err);
    fprintf(fp1,"%d\t%f\n", n, err);
    fprintf(fk1,"%d\t%f\n",n, suma_k);
		free_vector(x, 1, n);
		free_vector(w, 1, n);
  }
  
////////////////////////////////////////////////////////////////////////////////////////
  printf("\n\n");
  printf("ZAD2, k=5");
  printf("\n\n");

  int k = 5;
  float c2a = factorial(k);

	for (int n = 2; n <= 20; ++n) {
		float* x = vector(1, n);
		float* w = vector(1, n);

		gaulag(x-1, w-1, n, 0);
    float suma_k = 0;

		float c2 = 0.0;
		for (int i = 1; i <= n; ++i) {
			c2 += w[i]*fun2(x[i], k);
		}

    for(int i=0;i<n;i++)
    {
      suma_k +=w[i];
    }

		float err = fabs(c2 - c2a);

		//printf("%d	%f\n", n, err);
    fprintf(fp2,"%d\t%f\n", n, err);
    fprintf(fk2,"%d\t%f\n",n, suma_k);
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}


   printf("\n\n");
  printf("ZAD2, k=10");
  printf("\n\n");

  k = 10;
  float c2b = factorial(k);

	for (int n = 2; n <= 20; ++n) {
		float* x = vector(1, n);
		float* w = vector(1, n);

		gaulag(x-1, w-1, n, 0);
    float suma_k = 0;

		float c2 = 0.0;
		for (int i = 1; i <= n; ++i) {
			c2 += w[i]*fun2(x[i], k);
		}

    for(int i=0;i<n;i++)
    {
      suma_k +=w[i];
    }

		float err = fabs(c2 - c2b);

		//printf("%d	%f\n", n, err);
    fprintf(fp3,"%d\t%f\n", n, err);
    fprintf(fk3,"%d\t%f\n",n, suma_k);
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}

//hermit
  float c3a = 0.1919832644;
  for (int n = 2; n <= 15; ++n) {
		float* x = vector(1, n);
		float* w = vector(1, n);

		gauher(x-1, w-1, n);
    float suma_k = 0;
    float s1 = 0, s2 = 0;

		float c3 = 0.0;
		for (int i = 1; i <= n; ++i) 
    {
			s1 += (w[i] * fun3_2(x[i]));
      s2 += (w[i] * fun3_4(x[i]));
		}

    c3 = s1*s2;

    for(int i=0;i<n;i++)
    {
      suma_k +=w[i];
    }

		float err = fabs(c3 - c3a);

		//printf("%d	%f\n", n, err);
    fprintf(fp4,"%d\t%f\n", n, err);
    fprintf(fk4,"%d\t%f\n",n, suma_k);
		free_vector(x, 1, n);
		free_vector(w, 1, n);
	}





  FILE * function1 = fopen("fun1.dat","w");
  //Fukncja podcałkowa 1
  for(float step = 0;step<=2;step +=0.01)
  {
    fprintf(function1,"%f\t%f\n",step, fun1(step));
  }
  

  FILE * function2 = fopen("fun2k5.dat","w");
  //Fukncja podcałkowa 1
  for(float step = 0;step<=3.0;step +=0.01)
  {
    fprintf(function2,"%f\t%f\n",step, fun2(step,5));
  }

  FILE * function3 = fopen("fun2k10.dat","w");
  //Fukncja podcałkowa 1
  for(float step = 0;step<=3.0;step +=0.01)
  {
    fprintf(function3,"%f\t%f\n",step, fun2(step,10));
  }
  



	return 0;

}
