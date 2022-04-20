


#include <iostream>
#include <cmath>
#include <fstream>
#define n 20


double czebyszew(double min, double max, int i);
double fun(double x);
double multiply(double x, int j, double min, double dx);
double W_n(double f[][n+1], double x, double min, double dx);
double W_nCzebyszew(double f[][n+1], double x, double min, double max);
void print_matrix(double F[][n+1]);
double multiply_czebyszew(double x, double j, double min, double max);



int main() {
  double xmin = -5;
  double xmax = 5;
  
  

  FILE *fp = fopen("lab7_wn.txt","w");
  FILE *wezly = fopen("wezly.txt", "w");


  double F[n+1][n+1] = {0};


  double dx = (xmax - xmin)/n;

  //zapis wezlow

  for(int i=0; i<= n; ++i)
  {
    F[i][0] = fun(xmin + i*dx);
    //std::cout<< F[i][0]<< " ,";
    fprintf(wezly,"%lf\t%lf\n",xmin+i*dx, F[i][0]);
  }
  std::cout<<std::endl;
//wartosci ilorazow roznicowych
  for(int j = 1;j <= n; ++j)
  {
    for(int i = j; i <= n; ++i)
    {
      F[i][j] = (F[i][j-1] - F[i-1][j-1]) / ((xmin + i*dx) - (xmin+(i-j)*dx));
    }
  }
//print F
  //print_matrix(F);

//wielomian inerpolacyjny


  double k = xmin;
  
  

  
  while(k<= xmax)
  {
    fprintf(fp, "%lf\t%lf\n",k,W_n(F,k,xmin,dx));
    k += 0.1;
    
  }

  fclose(fp);

  //****************czebyszew**********************************\\
  std::cout<<"\n Czebyszew \n"<<std::endl;

  FILE *wezly2 = fopen("wezly2.txt", "w");
  FILE *fp2 = fopen("lab7_wn_czebyszew.txt", "w");

  for(int i=0; i<= n; ++i)
  {
    F[i][0] = fun(czebyszew(xmin, xmax , i));
    //std::cout<< F[i][0]<< " ,";
    fprintf(wezly2,"%lf\t%lf\n",czebyszew(xmin, xmax , i), fun(czebyszew(xmin, xmax , i)));
  }
  std::cout<<std::endl;



  for(int j = 1;j <= n; ++j)
  {
    for(int i = j; i <= n; ++i)
    {
      F[i][j] = (F[i][j-1] - F[i-1][j-1]) / ((czebyszew(xmin, xmax, i)) - czebyszew(xmin, xmax, i-j));
    }
  }

  print_matrix(F);

  k = xmin;
  while(k<= xmax)
  {
    fprintf(fp2, "%lf\t%lf\n",k, W_nCzebyszew(F, k, xmin, xmax));
    k += 0.1;
    
  }



    

  
    

  return 0;

}

double W_nCzebyszew(double f[][n+1], double x, double min, double max)
{
  double result = 0.0;
  for(int i = 0; i < n+1; i++)
  {
    result += f[i][i] * multiply_czebyszew(x, i, min, max);
  }

  return result;
}

double W_n(double f[][n+1], double x, double min, double dx)
{
  double result = 0.0;
  for(int i = 0; i < n+1; i++)
  {
    result += f[i][i] * multiply(x, i, min, dx);
  }

  return result;
}

double multiply(double x, int j, double min, double dx)
{
  double result = 1.0;

  for(int i=0; i<j; i++)
  {
    result *= (x-(min + i*dx));
  }

  return result;
}

double multiply_czebyszew(double x, double j, double min, double max)
{
  double result = 1.0;

  for(int i=0; i<j; i++)
  {
    result *= (x-czebyszew(min, max, i));
  }

  return result;
}

double fun(double x)
{
    return 1/(1 + x*x);
}


double czebyszew(double min, double max, int i)
{
    return 0.5 * ((min - max) * cos(M_PI * ((2.0 * static_cast<double>(i) + 1.0) / (2.0 * static_cast<double>(n) + 2.0))) + (min + max));
}

void print_matrix(double F[][n+1]){

  for (int i = 0; i <= n; ++i)
  {
      for (int j = 0; j <=n ; ++j)
      {
        printf("%lf    ", F[i][j]);
      }
    std::cout << std::endl;
  }

}