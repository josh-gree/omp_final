#include <stdio.h>
#include <math.h>


#define N 729
#define reps 100
#define N_Threads 6 
#include <omp.h> 

double a[N][N], b[N][N], c[N];
int jmax[N];  


void init1(void);
void init2(void);
void loop1(void);
void loop2(void);
double valid1(void);
double valid2(void);


int main(int argc, char *argv[]) { 

  double start1,start2,end1,end2;
  int r;

  init1(); 

  start1 = omp_get_wtime(); 

  for (r=0; r<reps; r++){ 
    loop1();
  } 

  end1  = omp_get_wtime();  


  double out1 = valid1();
  printf("%f,", (float)(end1-start1)); 


  init2(); 

  start2 = omp_get_wtime(); 

  for (r=0; r<reps; r++){ 
    loop2();
  } 

  end2  = omp_get_wtime(); 


  double out2 = valid2();
  printf("%f,guided,64,%.10f,%.10f\n", (float)(end2-start2), out1, out2); 

} 

void init1(void){
  int i,j; 

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      a[i][j] = 0.0; 
      b[i][j] = 3.142*(i+j); 
    }
  }

}

void init2(void){ 
  int i,j, expr; 

  for (i=0; i<N; i++){ 
    expr =  i%( 3*(i/30) + 1); 
    if ( expr == 0) { 
      jmax[i] = N;
    }
    else {
      jmax[i] = 1; 
    }
    c[i] = 0.0;
  }

  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      b[i][j] = (double) (i*j+1) / (double) (N*N); 
    }
  }
 
} 

void loop1(void) { 
  int i,j; 
#pragma omp parallel for num_threads(N_Threads) schedule(guided,64) default(none) shared(a,b) private(i,j)
  for (i=0; i<N; i++){ 
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    } 
  }

} 



void loop2(void) {
  int i,j,k; 
  double rN2; 

  rN2 = 1.0 / (double) (N*N);  
#pragma omp parallel for num_threads(N_Threads) schedule(guided,64) default(none) shared(b,c,jmax) private(i,j,k) firstprivate(rN2)
  for (i=0; i<N; i++){ 
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){ 
	c[i] += (k+1) * log (b[i][j]) * rN2;
      } 
    }
  }

}

double valid1(void) { 
  int i,j; 
  double suma; 
  
  suma= 0.0; 
  for (i=0; i<N; i++){ 
    for (j=0; j<N; j++){ 
      suma += a[i][j];
    }
  }
  return suma;
} 


double valid2(void) { 
  int i; 
  double sumc; 
  
  sumc= 0.0; 
  for (i=0; i<N; i++){ 
    sumc += c[i];
  }
  return sumc;
} 
 
