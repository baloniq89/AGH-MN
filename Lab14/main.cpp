#include <iostream>
#include <cmath>
#include <time.h>
#define K 10



double gen_1();
double gen_2();
double gen_3();
int in_range(double p,double u, double pend);

int main() {
  

  int n = 10000;
  //zad1 a
  FILE *fp1 = fopen("gen1_a.dat","w");

  double *u1 = new double[n];

  for(int i=0;i<n;i++)
  {
    u1[i] = gen_1();
  }

  for(int i=0; i<n; i++)
  {
    fprintf(fp1, "%g\t%g\n", u1[i], u1[i+1]);
  }
  //zad1 b

  FILE *fp2 = fopen("gen1_b.dat","w");

  double *u2 = new double[n];

  for(int i=0;i<n;i++)
  {
    u2[i] = gen_2();
  }

  for(int i=0; i<n; i++)
  {
    fprintf(fp2, "%g\t%g\n", u2[i], u2[i+1]);
  }


  //zad2
  n = 1000;
  FILE *fp3 = fopen("gen_troj.dat","w");
  double *u3 = new double[n];
  double step = 6.0/K;

  for(int i=0;i<n;i++)
  {
    u3[i] = gen_3();
  }
  for(int i=0; i<n; i++)
  {
    fprintf(fp3, "%g\t%g\n", u3[i], u3[i+1]);
  }

  for(double p = 1.0; p<=6.5;p+=step)
  {
    int count = 0;
    for(int i=0; i<n; i++)
    {
      count += in_range(p, u3[i], p+step);
    }
    std::cout<<"W przedziale <"<<p<<", "<<p+step<<" >mamy "<<count<<" liczb"<<std::endl;
  }



  delete []u1;
  delete []u2;
  delete []u3;
}


int in_range(double p,double u, double pend){

  

  if(u >= p && u <=pend)
  {
    return 1;
  }
  return 0;
}

double gen_1(){
  srand(time(0));
  static long int x=10;
  int a=123;
  int c=1;
  long int m=pow(2,15);
  x=(a*x+c) % m;
  return x/(m+1.0);
}

double gen_2(){
static long int x=10;
int a=69069;
int c=1;
long int m=pow(2,32);
x=(a*x+c) % m;
return x/(m+1.0);
}

double gen_3(){

  //static long int x = 10;
  double u = 4.0;
  double delta = 3.0;

  double e1 = gen_1();
  double e2 = gen_1();

  double x = u + (e1 + e2 -1)*delta;
  return x;

}