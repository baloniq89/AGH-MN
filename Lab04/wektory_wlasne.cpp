#include <iostream>
#include <cmath>

#define N 50



void normalizacja_wektora(double *wektor);
void multiple(double a[][N], double *c, double *wynik);
void zerowanie_macierzy(double a[][N]);
void initialize( double a[N][N],double L);
void print_matrix(double a[][N]);
void print_tab(double *tab);
double searching_lambda( double left, double right, double h[][N], int wart_wl, int licznik);





int main()
{ 
 
  FILE *fp = fopen("lab4.txt", "w");

  double L = 5.0;
  double H[N][N];
  zerowanie_macierzy(H);
  initialize(H, L);

  //print_matrix(H);
  
  double prawy_koniec = -(H[0][1] + H[2][1]) + H[N-1][N-1];
  double lewy_koniec = H[0][1] + H[2][1] - H[N-1][N-1];

  
  double lambda[6] = {0.0};

  for(int i=1; i< 6; i++)
  {
    lambda[i] = searching_lambda(lewy_koniec, prawy_koniec, H, i,150);
    std::cout<<lambda[i]<<"  ";
  }

  double vec_wlasny [N];



//wykonywalem poniższe operacje 5 razy dla lambda[i] gdzie i <1;5>
  vec_wlasny [0] = 1.0;
  vec_wlasny [1] = lambda[5] - H[0][0] / H[0][1];
  for(int i=2; i < N; i++)
  {
    vec_wlasny[i] = (((lambda[5] - H[i - 1][i - 1]) * vec_wlasny[i - 1] - H[i - 2][i - 1] * vec_wlasny[i - 2])/H[i-1][i]);
  }
  normalizacja_wektora(vec_wlasny);
  for(int i=0; i<N; i++)
  {
    fprintf(fp, "%lf \n", vec_wlasny[i]);
  }

  std::cout<<std::endl;
  
  print_tab(vec_wlasny);






  return 0;
}



void normalizacja_wektora(double *wektor)
{
  double suma = 0.0;

  for(int i=0; i< N; i++)
  {
    suma += pow(wektor[i], 2);
  }
  suma = sqrt(suma);

  for(int i=0; i< N; i++)
  {
    wektor[i] = wektor[i] / suma;
  }
  
}



double searching_lambda( double left, double right, double h[][N], int wart_wl, int dokladnosc )
{
  int sing_change = 0;
  double lambda =  (left + right)/2;

  if(dokladnosc < 1)
  {
    return lambda;
  }
  double w[N] = {0};
  w[0] = 1;
  w[1] = h[0][0] - lambda;
  

  for(int i=2; i < N ; i++)
  {
    w[i] = (h[i - 1][i - 1] - lambda) * w[i - 1] - (h[i - 1][i] * h[i - 1][i]) * w[i - 2];
  }

  for(int i=0; i< N-1; i++)
  {
    if(w[i] * w[i+1] < 0)
      sing_change = sing_change +1;
    
  }
  

  if(sing_change >= wart_wl)
  {
    return searching_lambda(left, lambda, h, wart_wl, dokladnosc - 1);
  }else
  {
    return searching_lambda(lambda, right, h, wart_wl, dokladnosc - 1);
  }

}


void print_matrix( double a[][N]) {

  for (int i = 0; i < N; i++)
  {
      for (int j = 0; j < N; j++)
      {
        printf("%.4f    ", a[i][j]);
      }
    std::cout << std::endl;
  }
}

void print_tab(double *tab) 
{
	
	std::cout<<"[";
    for (int i = 0; i < N; i++) 
	{
        printf("%lf, \n", tab[i]);
  }
  std::cout<<"]"<<std::endl<<std::endl;
	
}




void zerowanie_macierzy(double a[][N])
{
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      a[i][j] = 0;
    }
  }
}


void initialize(double a[N][N], double L){
  
  for(int i=1;i<N;i++)
  {
    a[i][i-1] = a[i-1][i] = (-1)/(2.0*(2*L/N)*(2*L/N));
    
  }
  for(int i=0;i<N;i++)
  {
    a[i][i] = pow(-L + (i+1) * (2.0*(L/N)), 2)/2 + pow(2.0*(L/N), -2);
  }

}