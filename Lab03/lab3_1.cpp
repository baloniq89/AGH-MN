#include <iostream>
#include <cmath>

#define N 1000
#define abs(X) ((X)>0? (X):-(X))


void multiple(double a[][N], double *c, double *wynik);



void zerowanie_macierzy(double a[][N]);

void initialize_a(double a[N][N]);
void initialize_b(double *b);
double iloczyn_skalar_wek(double * wek1, double * wek2);
void wektor_x_skalar (double *wek, double skalar);

double norma_euklidesowa(double * wek);
void najsz_spadek(double b[N], double x[N], double r[N], double a[][N]);




int main()
{ 
 

  double A[N][N];
  double b[N];

  double r[N];
  double x[N];
  for(int i=0;i<N;i++)
  {
    r[i] = 0;
    x[i] = 1;
  }




  initialize_a(A);
  initialize_b(b);

  najsz_spadek(b, x, r, A);


  

  

  
	

    return 0;
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

void initialize_b(double *b){
  for(int i=0;i <N;i++)
  {
    b[i] = i;
  }
}

void initialize_a(double a[N][N]){
  int m = 5;
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      if(abs(i-j) <= m)
      {
        a[i][j] = 1./(1. + abs(i-j));
      }else{
        a[i][j] = 0.0;
      }
    }
  }

}



void najsz_spadek(double b[N], double x[N], double r[N], double a[][N]){

  FILE* fp = fopen("lab3_1.txt", "w");
  double wynik[N], wynik_2[N];
  int k = 0;
  double alfa = 0.0;
  do
  {
    multiple(a,x,wynik);
    for(int i=0;i<N;i++)
    {
      r[i] = b[i] - wynik[i];
    }

    multiple(a, r, wynik);

    alfa = iloczyn_skalar_wek(r,r) / iloczyn_skalar_wek(r, wynik); 

    wektor_x_skalar(r, alfa);

    for (int i = 0; i < N; i++)
		{
			x[i] += r[i];
		}
    
    k++;

    fprintf (fp, " %d     %lf      %lf      %lf\n", k, norma_euklidesowa(r), norma_euklidesowa(x),alfa);

  }while(norma_euklidesowa(r) > 0.0000001);

}

double norma_euklidesowa(double * wek)
{
  return sqrt(iloczyn_skalar_wek(wek, wek));
}



void wektor_x_skalar (double *wek, double skalar){

  for(int i=0;i<N;i++)
  {
    wek[i] = wek[i]*skalar;
  }
}

double iloczyn_skalar_wek(double * wek1, double * wek2)
{
  double suma = 0.0;
  for(int i=0; i<N; i++)
  {
    suma += wek1[i] * wek2[i];
  }
  return suma;
}

void multiple(double a[][N], double *c, double *wynik){

  for(int i = 0; i < N; ++i){
        wynik[i]=0;
  }
 

  for(int i=0;i<N;++i)
  {
    for(int j=0;j<N;++j)
    {
      wynik[i] += a[i][j] * c[j];
    }
  }
  
}
