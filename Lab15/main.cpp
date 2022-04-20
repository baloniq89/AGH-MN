#include <iostream>
#include <cmath>
#include <ctime>

double d_rand(double min, double max);
double calc_Z1();
double calc_Z2();
double calculate_V(int n,double x1, double y1, double x2, double y2, double x);

int main() {

  srand(time(NULL));
  FILE *fp = fopen("monteCarlo.dat","w");
  
  double x1=0;
  double y1=0;
  double x2=0;
  double y2=0;
  double V = 0;
  double sigma = 1.0;
  int n = 1000;

  for(double x = 0.0; x<=6.01;x +=0.01)
  {
    x1 = calc_Z1();
    x2 = calc_Z1();
    y1 = calc_Z2();
    y2 = calc_Z2();

    V = calculate_V(n,x1,y1,x2,y2,x);

    fprintf(fp,"%lf\t%lf\n",x,V/n);


  }


  fclose(fp);
  
}
double d_rand(double min, double max){
	double r = (double) rand()/RAND_MAX;
	return r*(max-min) + min;
}


double calc_Z1()
{
  double u1 = d_rand(0,1);
  double u2 = d_rand(0,1);

  return (sqrt((-2) * log(u1)) * cos(2 * M_PI * u2));
}

double calc_Z2()
{
  double u1 = d_rand(0,1);
  double u2 = d_rand(0,1);

  return (sqrt((-2) * log(u1)) * sin(2 * M_PI * u2));
}

double calculate_V(int n,double x1, double y1, double x2, double y2, double x)
{

  double v = 0.0;
  for(int i=0;i<n;i++)
  {
    v +=(1/(sqrt(pow(x1-(x2+x),2) + pow(y1-y2,2))));
  }

  return v;
}