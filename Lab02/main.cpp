


#include <iostream>
#include <cmath>
#define N 3


void print_matrix(double a[][N]);
void print_tab(double *tab);
void find_y(double *x, double *c, double *y);
void initialize_a(double a[][N], double *c, double *x);;
void zerowanie_macierzy(double a[][N]);
void initialize_LU(double a[][N], double l[][N], double u[][N]);
void forw_subs(double l[][N], double * z, double *y);
void back_subs(double u[][N], double *z, double *x);
double wyznacznik_trojktnej(double a[][N]);
void odwrotna_L(double l[][N],double l_odwr[][N]);
void odwrotna_U(double U_matrix[][N], double result[][N]);
void schemat_hornera(double *w, double *c, double *x);
void mnozenie_macierzy(double a[][N], double b[][N], double c[][N]);
double wsk_uwar(double a[][N], double a_odwr[][N]);




int main()
{ 
	FILE* fp = fopen("lab1_2.txt", "w");
	double a[N][N];
  double l_odwrotna[N][N];
  double u_odwrotna[N][N];
  double a_odwrotna[N][N];

  double x[N] = {2, 5, 7};
  double c[N] = {3.5,  4, -5};
  double y[N];
  double lower[N][N];
  double upper[N][N];
  double z[N];
  double rozw[N];
  double horner_rozw[N];

  zerowanie_macierzy(lower);
  zerowanie_macierzy(upper);


  
  find_y(x,c,y);
  print_tab(y);
  initialize_a(a, c, x);
  std::cout<<"Macierz A"<<std::endl;
  print_matrix(a);

  initialize_LU(a,lower,upper);

  std::cout<<"Macierz L"<<std::endl;
  print_matrix(lower);
  odwrotna_L(lower, l_odwrotna);
  std::cout<<"Macierz odwrotna do L"<<std::endl;
  print_matrix(l_odwrotna); 

  std::cout<<"Macierz U"<<std::endl;
  print_matrix(upper);
  odwrotna_U(upper, u_odwrotna);
  std::cout<<"Macierz odwrotna do U"<<std::endl;
  print_matrix(u_odwrotna);


  //rozwiązywanie układu równan
  forw_subs(lower, z, y);
 // print_tab(z);

  back_subs(upper, z, rozw);
 // print_tab(rozw);
//*********************

  std::cout<<"Wyznacznik macierzy A = "<<wyznacznik_trojktnej(upper)<<std::endl;

  

  mnozenie_macierzy(u_odwrotna, l_odwrotna, a_odwrotna);
  std::cout<<"Macierz odwrotna do macierzy A"<<std::endl;
  print_matrix(a_odwrotna);
  /*
  std::cout<<"Horner"<<std::endl;
  schemat_hornera(horner_rozw, c, x);
  print_tab(horner_rozw);
  */

 std::cout<<"Wskaznik uwarunkowania macierzy wynosi = "<<wsk_uwar(a, a_odwrotna)<<std::endl;

double test[N][N];
mnozenie_macierzy(a, a_odwrotna, test);

print_matrix(test);


  

  
	

    return 0;
}

double wsk_uwar(double a[][N], double a_odwr[][N]) {

  double max_a = 0.0;
  double max_a_odwr = 0.0;
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      if(max_a < a[i][j])
      {
        max_a = a[i][j];
      }
    }
    for(int j=0;j<N;j++)
    {
      if(max_a_odwr < a[i][j])
      {
        max_a_odwr = a[i][j];
      }
    }
  }

  return max_a * max_a_odwr;

}

void schemat_hornera(double *w, double *c, double *x){

  for(int i=0;i<N;i++)
  { 
    w[i] = c[N-1];
    for( int j=N-2;j>=0;j--)
    { 
      w[i] = w[i] * x[i] + c[j];
    }
    
  }
}

void multiple_a_c(double a[][N], double *c){

  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      a[i][j] = a[i][j] * c[j];
    }
  }
}

void find_y(double *x, double *c, double *y)
{
  for(int i=0;i<N;i++)
  {
    
    for(int j=0;j<N;j++)
    {
      y[i] = y[i]+ pow(x[i],j) * c[j];
    }
  }
}

void initialize_a(double a[][N], double *c, double *x){

  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      a[i][j] = c[j] * pow(x[i],j);
    }
  }
}

void initialize_LU(double a[][N], double l[][N], double u[][N])
{
 
  for(int i=0;i<N ;i++)
  {
    for(int j=0;j<N;j++)
    {
      double temp = 0.0;
      for(int k=0;k<i;k++)
      {
        temp = temp + (l[i][k] * u[k][j]);
      }
      u[i][j] = a[i][j] - temp;
    }

    for(int j=0; j<N; j++)
    {
      if(i == j)
      {
        l[i][j] = 1;
      }else {
        double temp = 0.0;
        for(int k=0;k<i;k++)
        {
          temp = temp + l[j][k]*u[k][i];
        }
        l[j][i] = (a[j][i] - temp )/u[i][i];
      }
    }


    
  }
  
}

  
void forw_subs(double l[][N], double * z, double *y){

  z[0] = y[0];
  for(int i=1;i<N;i++)
  {
    double s_temp = 0.0;
    for(int j = 0;j< i-1;j++)
    {
      s_temp += l[i][j] * z[j];
    }
    z[i] = y[i] - s_temp;
  }

}

void back_subs(double u[][N], double *z, double *x){

  x[N-1] = z[N-1]/u[N-1][N-1];

  for(int i= N-2; i >=0;i--)
  {
    double s_temp = 0.0;
    for(int j = i+1;j<N;j++)
    {
      s_temp += u[i][j] * x[j];
    }
    x[i] = (z[i] - s_temp)/u[i][i];
  }

}

double wyznacznik_trojktnej(double a[][N]){

  double det  = a[0][0];
  for(int i=1;i<N;i++)
  {
    det = det *a[i][i];
  }

  return det;
}



void print_matrix(double a[][N]) {
    int i, j; 
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%0.4lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");
}

void print_tab(double *tab) 
{
	
	std::cout<<"[";
    for (int i = 0; i < N; i++) 
	{
        printf("%2.1f, ", tab[i]);
  }
  std::cout<<"]"<<std::endl<<std::endl;
	
}

void zerowanie_macierzy(double a[][N])
{
  for (int i = 0; i < N; i++) 
  {
    for (int j = 0; j < N; j++) 
    {
      a[i][j] = 0;
    }
  }
}

void odwrotna_L(double l[][N],double l_odwr[][N]){

  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      if( i >j)
      {
       double s_temp = 0.0;
       for(int k = j+1;k <= i-1;k++)
       {
         s_temp += l[i][k] * l_odwr[k][j];
       }
       l_odwr[i][j] = -l[i][j] - s_temp;

      }else if(i ==j)
      {
        l_odwr[i][j] = 1;
      }else
      { 
        l_odwr[i][j] = 0;
      }
    }
  }
}


void odwrotna_U(double u[][N], double u_odwr[][N]){

  for(int i=0;i<N;i++)
  {
    u_odwr[i][i] = 1/u[i][i];
    for(int j=i+1; j<N; j++)
    {
      double s_temp = 0.0;
       for(int k =i; k < j; k++)
       {
         s_temp += (u[k][j]*u_odwr[i][k]);
       }
       u_odwr[i][j] = -(s_temp/u[j][j]);
    }
  }
  u_odwr[N-1][N-1] = 1 / u[N-1][N-1];
  
}

void mnozenie_macierzy(double a[][N], double b[][N], double c[][N])
{
  for( int i = 0; i < N; i++ )
    for( int j = 0; j < N; j++ )
    {
      
      for( int k = 0; k < N; k++ )
      {
         c[i][j] += a[i][k] * b[k][j];
      }
      
    }
}