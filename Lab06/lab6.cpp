#include <iostream>
#include <cmath>

#define N 2


void multipleMatrixVector(double matrix[][N],double vec[N], double wynik[N]);
double norm(double wektor[N]);
void initializeVector(double vec[N], double rStart[N]);
void initializeMatrix(double matrix[][N], double vec[N]);
void print_matrix( double a[][N]);
void print_vec(double *tab) ;
//odwracanie macierzy************geekforgeeks
void getCofactor(double A[N][N], double temp[N][N], double p, double q, double n);
double determinant(double A[N][N], double n);
void adjoint(double A[N][N], double adj[N][N]);
void inverse(double A[N][N], double inverse[N][N]);

//odwracanie macierzy ******************
void addVectors(double a[N], double b[N]);


int main() 
{
  double r_o[N]  = {10 ,10};
  double StartMatrix[N][N] ;
  double vecR[N] = {0,0};
  
 

  initializeMatrix(StartMatrix, r_o);
  initializeVector(vecR, r_o);
  
  //print_vec(vecR);
  //print_matrix(StartMatrix);

  

  double deltaR[N]={1.0,1.0};
  double inversed_R[N][N] = {0.0};
  inverse(StartMatrix, inversed_R);
  //multipleMatrixVector(inversed_R, vecR, deltaR);

  //print_matrix(inversed_R);
  //std::cout<<std::endl;
  multipleMatrixVector(inversed_R, vecR, deltaR);




  //print_vec(r_o);
  
  while(norm(deltaR) > pow(10, -6))
  {
    std::cout<<r_o[0]<<"   "<<r_o[1]<<std::endl;



    initializeMatrix(StartMatrix, r_o);
    initializeVector(vecR, r_o);
    inverse(StartMatrix, inversed_R);
    multipleMatrixVector(inversed_R, vecR , deltaR);
    addVectors(r_o, deltaR);
  }
  





  return 0;

}



void addVectors(double a[N], double b[N])
{
  for(int i=0; i<N; i++)
  {
    a[i] += b[i];
  }
}

void print_vec(double *tab) 
{
	
	std::cout<<"[";
    for (int i = 0; i < N; i++) 
	{
        printf("%lf ,", tab[i]);
  }
  std::cout<<"]"<<std::endl<<std::endl;
	
}

double norm(double wektor[N])
{
  double suma = 0.0;

  for(int i=0; i< N; i++)
  {
    suma += pow(wektor[i], 2);
  }
  return sqrt(suma);

  //return suma;
  
}

void multiple_a_c(double a[][N], double c){

  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      a[i][j] = a[i][j] * c;
    }
  }
}

void initializeVector(double vec[N], double rStart[N])
{
  double x = rStart[0];
  double y = rStart[1];

  vec[0] = 2 * x * pow(y,2) - 3 * pow(x,2) * y - 2;
  vec[1] = pow(x,2) * pow(y,3) + 2 * x * y - 12;
}

void initializeMatrix(double matrix[][N], double vec[N])
{

  
  double x = vec[0];
  double y = vec[1];

  matrix[0][0] = -(2 * pow(y, 2) - 6 * x * y);
  matrix[1][0] = -(2 * x * pow(y, 3) + 2 * y);
  matrix[0][1] = -(4 * x * y - 3*pow(x, 2));
  matrix[1][1] = -(3 * pow(y, 2) * pow(x,2) + 2 * x);
  
}

void multipleMatrixVector(double matrix[][N],double vec[N], double wynik[N])
{ 
  wynik[0] = 0;
  wynik[1] = 0;
  for(int i = 0; i < N; i++)
  {
    for(int j=0; j<N; j++)
    {
      wynik[i] += matrix[i][j] * vec[j];
    }
  }
}



void getCofactor(double A[N][N], double temp[N][N], double p, double q, double n) {
  int i = 0, j = 0;

  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {

      if (row != p && col != q) {
        temp[i][j++] = A[row][col];

        if (j == n - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
}

double determinant(double A[N][N], double n) {
  double D = 0; 

  if (n == 1)
    return A[0][0];

  double temp[N][N];
  int sign = 1;

  for (int f = 0; f < n; f++) {
    getCofactor(A, temp, 0, f, n);
    D += sign * A[0][f] * determinant(temp, n - 1);
    sign = -sign;
  }

  return D;
}

void adjoint(double A[N][N], double adj[N][N]) {
  if (N == 1) {
    adj[0][0] = 1;
    return;
  }

  int sign = 1;
  double temp[N][N];

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      getCofactor(A, temp, i, j, N);
      sign = ((i + j) % 2 == 0) ? 1 : -1;
      adj[j][i] = (sign) * (determinant(temp, N - 1));
    }
  }
}

void inverse(double A[N][N], double inverse[N][N]) {

  double det = determinant(A, N);
  if (det == 0) {
    std::cout << "Singular matrix, can't find its inverse";
    
  }

  double adj[N][N];
  adjoint(A, adj);

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      inverse[i][j] = adj[i][j] / double(det);
  
  
}

void print_matrix( double a[][N]) {

  for (int i = 0; i < N; i++)
  {
      for (int j = 0; j < N; j++)
      {
        printf("%lf    ", a[i][j]);
      }
    std::cout << std::endl;
  }
}
