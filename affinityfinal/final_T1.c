#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>

#define N_Threads 1
#define N 729
#define Min_Chunk 1

int i;
int j;
double a[N][N], b[N][N], c[N];
int jmax[N];

//Helper function
//----------------------------
int int_ceil_div(int x,int y){
  return (x + y - 1) / y;
} //returns ceil(x/y)
//----------------------------


//Define loops
//---------------------------- 
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

void loop1(int start, int end) { 
  int i,j; 
  for (i=start; i<end; i++){ 
    for (j=N-1; j>i; j--){

      a[i][j] += cos(b[i][j]);

    } 
  }
}

void loop2(int start,int end) {
  int i,j,k; 
  double rN2; 

  rN2 = 1.0 / (double) (N*N);  

  for (i=start; i<end; i++){ 
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
//----------------------------



//Set up Q implementation
//----------------------------

struct Data {
  int start;
  int end;
  int size;
}; //each element of the Q has this struct as its data

struct Node
{
  struct Data* data;
  struct Node* next;
}; //Node in the Q

struct Q {
  struct Node* front;
  struct Node* rear;
}; //Container for the Q -> pointer to front and rear

struct Q* newQ(){
  struct Q* new = (struct Q*)malloc(sizeof(struct Q*));
  new->front = NULL;
  new->rear = NULL;
  return new;
}; //Initialise a new Q

void enq(int xstart, int xend, int xsize,struct Q* q){

  struct Node* temp = (struct Node*)malloc(sizeof(struct Node*));
  struct Data* tempdat = (struct Data*)malloc(sizeof(struct Data*));
  tempdat->start = xstart;
  tempdat->end = xend;
  tempdat->size = xsize;

  temp->data = tempdat;
  temp->next = NULL;

  if(q->front==NULL && q->rear==NULL){
    q->front = temp;
    q->rear = temp;
    return;
  } 
  q->rear->next = temp;
  q->rear = temp;
} //Add element to end of Q

struct Data* deq(struct Q* q) {
  struct Node* temp = q->front;
  if (q->front==NULL) return NULL;
  if (q->front==q->rear){
    struct Data* tempdat = q->front->data;
    q->front = NULL;
    q->rear = NULL;
    return tempdat;
  }
  else {
    struct Data* tempdat = q->front->data;
    q->front = q->front->next;
    return tempdat;
  }
  free(temp);
} //remove element from the front of the Q

struct Data* peek(struct Q* q) {
  if (q->front==NULL) return NULL;
  else {
    struct Data* tempdat = q->front->data;
    return tempdat;
  }
} //view but not remove element from the front of the Q

//!!!Possibly not needed!!!

int myloop(int start,int end){
  int sum = 0;
  for(i=start;i<end;i++){
    sum += i;
  }
  return sum;
} //test loop sum of consecutive integers

void print_list(struct Q* q) {
  struct Node * current = q->front;
  while (current != NULL) {
    int ID = omp_get_thread_num();
    printf("start = %d, end = %d, size = %d ----->%d\n", current->data->start,current->data->end,current->data->size,ID);
    current = current->next;
  }
} //print contents of Q
//----------------------------


//Main Functions
//----------------------------
int * construct_global_partition(){

  static int bounds[N_Threads+1];
  int chunksz = N/N_Threads;
  int rem_its = N % N_Threads;
  int loopcount = 0;

  //Fill bounds array with unadjusted chunksize 

  for(i = 0; i <= N;i+= chunksz){
    bounds[loopcount] = i;
    loopcount++;
  }

  //Add one to the first rem_its chunks 
  for(i = 0;i < rem_its;i++){
    bounds[i+1] += i+1;
  }

  //Adjust leftover bounds 
  for(i = rem_its+1;i < N_Threads+1;i++){
    bounds[i] += rem_its;
  }

  return bounds;
} //returns the upper and lower bounds for each threads local set
  //distributes non divisible iterations one at a time to each thread

void partition_localset(struct Q * list[],int bounds[]){

  int ID = omp_get_thread_num();
  int n = bounds[ID+1] - bounds[ID];
  int sum = bounds[ID];
  int tempn = n;
  
  list[ID] = newQ();
  
  while(sum<bounds[ID+1]){
    int x = int_ceil_div(tempn,N_Threads);
    if(x<Min_Chunk){
      x = Min_Chunk;
    }
  
    if(sum + x > bounds[ID+1]){
      x = bounds[ID+1] - sum;
      enq(sum,sum+x,x,list[ID]);
      break;
    }
    enq(sum,sum+x,x,list[ID]);
    sum += x;
    tempn -= x; 
  }
  
} //forms Q of iterations for each thread
//----------------------------

//Work Functions
//----------------------------
void process_own_Ql1(struct Q * list[]){
  
  int out = 0.0;
  int ID = omp_get_thread_num();
  struct Data* tempdat;
  struct Data* tempdatothers;

    #pragma omp critical
  {
    tempdat = deq(list[ID]);
  }

  while (tempdat != NULL){
    loop1(tempdat->start,tempdat->end);
    #pragma omp critical
    {
      tempdat = deq(list[ID]);
    }
  }
  
} //Process iterations in own Q

void process_own_Ql2(struct Q * list[]){
  
  int out = 0.0;
  int ID = omp_get_thread_num();
  struct Data* tempdat;
  struct Data* tempdatothers;

    #pragma omp critical
  {
    tempdat = deq(list[ID]);
  }

  while (tempdat != NULL){
    loop2(tempdat->start,tempdat->end);
    #pragma omp critical
    {
      tempdat = deq(list[ID]);
    }
  }
  
} //Process iterations in own Q

struct Data* nxt_chunk(struct Q * list[]){
  struct Data* tempdat;
  struct Data* finaldat;
  int large = 0;
  int pos = -1;
  for (i=0;i<N_Threads;i++){
 
 
    tempdat = peek(list[i]);
 
    if (tempdat != NULL){
      if (tempdat->size > large){
        large = tempdat->size;
        pos = i;
      }
    }
  }
  if (pos != -1){
    #pragma omp critical
    {
      finaldat = deq(list[pos]);
    }
    return finaldat;
  }
  else{
    return NULL;
  }
} //Returns the largest available chunk from other Q's


void process_other_Qsl1(struct Q * list[]){
  int out = 0;
  int ID = omp_get_thread_num();
  struct Data* tempdat =  (struct Data*)malloc(sizeof(struct Data*));
 
  
  tempdat = nxt_chunk(list);

  while (tempdat!=NULL){
    loop1(tempdat->start,tempdat->end);

    tempdat = nxt_chunk(list);
  }

} //Process iterations from other Q's until all Q's are empty

void process_other_Qsl2(struct Q * list[]){
  int out = 0;
  int ID = omp_get_thread_num();
  struct Data* tempdat =  (struct Data*)malloc(sizeof(struct Data*));
 
  
  tempdat = nxt_chunk(list);

  while (tempdat!=NULL){
    loop2(tempdat->start,tempdat->end);

    tempdat = nxt_chunk(list);
  }

} //Process iterations from other Q's until all Q's are empty


int main(){
  int * bounds; // Store bounds of threads local set
  struct Q * listl1[N_Threads]; 
  struct Q * listl2[N_Threads]; // Array of Q's one for each thread 
  init1();

  bounds = construct_global_partition();

#pragma omp parallel num_threads(N_Threads) shared(listl1,listl2,bounds)
  {
    partition_localset(listl1,bounds);
    partition_localset(listl2,bounds);
  }
  
  double t0l1= omp_get_wtime();
#pragma omp parallel num_threads(N_Threads) shared(listl1,a,b)
  {
    process_own_Ql1(listl1);
    process_other_Qsl1(listl1); 
  }
  double t1l1= omp_get_wtime();

  init2();

  double t0l2= omp_get_wtime();
#pragma omp parallel num_threads(N_Threads) shared(listl2,a,b,c,jmax)
  {
    process_own_Ql2(listl2);
    process_other_Qsl2(listl2);
  }
  double t1l2= omp_get_wtime();
  
  double out1 = valid1();
  double out2 = valid2();
  printf("%.10f,%.10f,%d,%.10f,%.10f\n",t1l1-t0l1,t1l2-t0l2,N_Threads,out1,out2);
  
  
}
