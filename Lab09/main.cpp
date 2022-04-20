


#include <iostream>
#include <cmath>
#define N 2/2
#define M 6/2
#define size M+N


void initialize_C(double c[size + 1]);
int silnia(int x);
double fun(double x);
void printVec(double c[size+1]);
void print_matrix( double a[][M]);
void solve_x_GJ(double* x, double* b, double a[][M]);



int main() {

    FILE *fp = fopen("aproksymacja.txt","w");
    

    double c[size + 1] = {0};
    double A[M][M] = {0};
    double a[N+1] = {0};
    double b[M+1] = {0};
    double x[M] , y[M];

    initialize_C(c);
    std::cout<<"wektor c"<<std::endl;
    printVec(c);

//inicjalizacja macierzy A
    for(int i=0;i<M;i++)
    {
      y[i] = -c[N + 1 + i];

      for(int j=0;j<M;j++)
      {
        A[i][j] = c[N - M + i + j + 1];
      }
      if((N-M)== -2)
      {
        A[0][0] = c[size];
      }
      if((N-M) == -3)
      {
        A[0][0] = c[size-1];
        A[0][1] = c[size];
        A[1][0] = c[size];
      }
    }
 
std::cout<<"Macierz A"<<std::endl;
print_matrix(A);

solve_x_GJ(x,y,A);

b[0] = 1;

for(int i=0;i<M;i++)
{
  b[M-i] = x[i];
}
//print vector b
std::cout<<"wektor b"<<std::endl;
for(int i=0;i<M+1;i++)
{
  std::cout<<b[i]<<std::endl;
}



for(int i=0;i<N+1;i++)
{
  for(int j=0;j<i+1;j++)
  {
    a[i] += c[i-j] * b[j];
  } 
}
std::cout<<"wektor a"<<std::endl;
for(int i=0;i<N+1;i++)
{
  std::cout<<a[i]<<std::endl;
}

double xmin = -5.0;
double xmax = 5.0;
//double P,Q;

double temp = xmin;

while(temp <=xmax)
{
		double P = 0;
    double Q = 0;
		for(int j = 0;j<N+1;j++)
    {
      P += a[j] * pow(temp, 2*j);
    }
    for(int j = 0;j<M+1;j++)
    {
			Q += b[j] * pow(temp, 2*j);
    }
		

    fprintf(fp,"%lf\t%lf\t%lf\n",temp,fun(temp), P/Q);
    temp+=0.01;
}
fclose(fp);


   

return 0;
}

void solve_x_GJ(double* x, double* b, double a[][M]) {
    int i, j, k;
    double w, w2;
    for (i = 0; i < M; i++) {
        w = a[i][i];
        b[i] = b[i] / w;
        for (j = 0; j < M; j++) {
            a[i][j] = a[i][j] / w;
        }
        for (k = 0; k < M; k++) {
            w2 = a[k][i];
            if (k!=i) {
                b[k] = b[k] - b[i] * w2;
                for (j = 0; j < M; j++) {
                    a[k][j] = a[k][j] - w2 * a[i][j];
                }
            }
        }
        
    }
    for (i = 0; i < M; i++) {
        x[i] = b[i];
    }
}

double fun(double x) {
  return exp(-(x*x));
}

int silnia(int x) {

  int silnia = 1;
  for(int i = x;i > 1;i--)
  {
    silnia *=i;
  }
  return silnia;

}

void initialize_C(double c[size + 1])
{
  for(int p = 0; p < size+1 ; p++)
  {
    c[p] = pow(-1,p) / silnia(p);
  }

  
}

void printVec(double c[size+1])
{
  for(int i=0;i<size+1;i++)
  {
    std::cout<<c[i]<<std::endl;
  }
}

void print_matrix( double a[][M]) {

  for (int i = 0; i < M; i++)
  {
      for (int j = 0; j < M; j++)
      {
        printf("%lf    ", a[i][j]);
      }
    std::cout << std::endl;
  }
}