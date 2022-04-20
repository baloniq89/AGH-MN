#include <iostream>
#include <cmath>
#define n 10
;
void solve_x_GJ(double* x, double* b, double a[][n]);
double fun(double x);
double dx_2(double x, double d);
void wyznacz_M(double xw[n], double yw[n],double m[n], double alpha, double beta);
void printVec(double vec[n]);

int main() {
  
  int temp;
  double alpha = 0.0;
  double beta = 0.0;
  double xw[n];  //polozenia wwezlow
  double yw[n];  // wartosci fukncji
  double xmin = -5;
  double xmax = 5;
  double m[n];

  double deltaX = (xmax - xmin) / (n - 1);

  double dx = 0.01;

  //wektor polozenia
  for(int i=0; i<n; i++)
  {
    xw[i] = xmin + deltaX*i;
  }

  //wektor wartosci fukncji

  for(int i=0; i<n; i++)
  {
    yw[i] = fun(xw[i]);
    temp = 1;
  }

  wyznacz_M(xw,yw,m,alpha,beta);

  

  
}

double fun(double x) {
  return 1 /(x*x);
}

double dx_2(double x, double d) {
  return (fun(x - d) - 2*fun(x) + fun(x + d)) / pow(d, 2); 
}

void printVec(double vec[n])
{
  for(int i=0;i<n; i++)
  {
    std::cout<<vec[i]<<"  ";
  }
  std::cout<<std::endl;
}

void wyznacz_M(double xw[n], double yw[n],double m[n], double alpha, double beta) {

  double A[n][n] = {0.0};
  A[0][0] = 1;
  A[n-1][n-1] = 1;
  double vector[n] = {0};
  vector [0] = alpha;
  vector [n-1] = beta;

  double h = 10.0 / (n - 1);
  double lambda = h / (2 * h);
  double my = 1.0 - lambda;
  double v = 0;

  for(int i=1; i < n-1; i++)
  {
    v = (6.0 / (2 * h)) * ((yw[i+1] - yw[i]) / h - (yw[i] - yw[i -1]) / h);
    vector[i] = v;
  }

  printVec(vector);

  for(int i=1; i < n-1; i++)
  {
    A[i][i] = 2.0;
  }
  for(int i=1; i < n-1; i++)
  {
    for(int j=1; j < n-1; j++)
      {
        if(i == j)
        {
          A[i][j-1] = my;
          A[i][j+1] = lambda;
        }
      } 
  }

  solve_x_GJ(m, vector, A);

  printVec(m);
  

}
//funkcja z pierwszych Laboratoriow
void solve_x_GJ(double* x, double* b, double a[][n]) {
    int i, j, k;
    double w, w2;
    for (i = 0; i < n; i++) {
        w = a[i][i];
        b[i] = b[i] / w;
        for (j = 0; j < n; j++) {
            a[i][j] = a[i][j] / w;
        }
        for (k = 0; k < n; k++) {
            w2 = a[k][i];
            if (k!=i) {
                b[k] = b[k] - b[i] * w2;
                for (j = 0; j < n; j++) {
                    a[k][j] = a[k][j] - w2 * a[i][j];
                }
            }
        }
        
    }
    for (i = 0; i < n; i++) {
        x[i] = b[i];
    }
}
