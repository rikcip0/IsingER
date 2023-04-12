#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXRAND (4294967296ULL)
#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)
#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)
#define sign(a) ( ( (a) < 0 )  ?  -1.   : (((a) > 0 ) ? 1. : 0. ) )

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

int N, M, *graph, **neigh, **neigh_pos, *deg;
double **u, beta;

unsigned rand4init(void) {
  unsigned long long y;
  
  y = myrand * 16807LL;
  myrand = (y & 0x7fffffff) + (y >> 31);
  if (myrand & 0x80000000)
    myrand = (myrand & 0x7fffffff) + 1;
  return myrand;
}

void Init_Random(void) {
  int i;
  
  ip = 128;    
  ip1 = ip - 24;    
  ip2 = ip - 55;    
  ip3 = ip - 61;
  
  for (i = ip3; i < ip; i++)
    ira[i] = rand4init();
}

void error(char* stringa) {
  
  fprintf(stderr, "ERROR: %s\n", stringa);
  exit(1);
}

double gaussRan(double sigma) { //restituisce numeri secondo distribuzione gauss di media nulla e varianza sigma^2
  double fac, rsq, v1, v2;

  do {
    v1 = 2.0 * FRANDOM - 1.0;
    v2 = 2.0 * FRANDOM - 1.0;
    rsq = v1 * v1 + v2 * v2;
  } while (rsq >= 1.0 || rsq == 0.0);
  fac = sqrt(-2.0 * log(rsq) / rsq);
  return v1 * fac * sigma;
}


void allocateMem(void) {

  graph = (int*)calloc(2*M, sizeof(int));
  deg = (int*)calloc(N, sizeof(int));
  neigh = (int**)calloc(N, sizeof(int*));
  neigh_pos = (int**)calloc(N, sizeof(int*));
  u = (double**)calloc(N, sizeof(double*));
}


void makeGraph(void) {
  int i, var1, var2;

  for (i = 0; i < N; i++)
    deg[i] = 0;

  for (i = 0; i < M; i++) {
    var1 = (int)(FRANDOM * N);
    do{
    var2 = (int)(FRANDOM * N);
    }while(var2==var1);
    
    graph[2*i] = var1;
    graph[2*i+1] = var2;
    deg[var1]++;
    deg[var2]++;
  }
  
  for (i = 0; i < N; i++) {
    neigh[i] = (int*)calloc(deg[i], sizeof(int));
    neigh_pos[i] = (int*)calloc(deg[i], sizeof(int));
    u[i] = (double*)calloc(deg[i], sizeof(double));
    deg[i] = 0;
  }
  
  for (i = 0; i < M; i++) {
    var1 = graph[2*i];
    var2 = graph[2*i+1];
    neigh[var1][deg[var1]] = var2;
    neigh_pos[var1][deg[var1]] = deg[var2];

    neigh[var2][deg[var2]] = var1;
    neigh_pos[var2][deg[var2]] = deg[var1];

    deg[var1]++;
    deg[var2]++;
  }
}

void init_u(){

  int i,j;

  for(i=0;i<N;i++){
    for(j=0;j<deg[i];j++){
      u[i][j]=FRANDOM*2.-1.;  //li inizializzo piccoli, tra -1 e 1
    }
  }

}

double f(double h){
  return atanh(tanh(beta)*tanh(beta*h))/beta;
}

void iteration(int t_MAX){
  int t,r,i,j,*n_ch,*pos,n;
  double h,diff=1,eps=1e-30,new_u, average;

  for(t=0;t<t_MAX&&diff>1e-7;t++){

    diff=0.;
    for(r=0;r<N;r++){

      i=(int)(FRANDOM*(double)N);
      n_ch=neigh[i]; pos=neigh_pos[i];

      h=0;

      for(j=0;j<deg[i];j++){  //aggiorno u[i][j]
	h+=u[n_ch[j]][pos[j]];
      }
      
      for(j=0;j<deg[i];j++){
	new_u=f(h-u[n_ch[j]][pos[j]]);
	diff+=fabs((new_u-u[i][j])/(new_u+eps));
	n++;
	u[i][j]=new_u;   
      }

    }
  }

  //mi calcolo la magnetizzazione
  average=0; 
  for(i=0;i<N;i++){
    
    n_ch=neigh[i]; pos=neigh_pos[i];
    
    h=0;

    for(j=0;j<deg[i];j++){  
      h+=u[n_ch[j]][pos[j]];
    }
    average+=tanh(beta*h);
  }

  printf("%g %g %d\n", 1./beta,average/N,t);  
  
}



int main(int argc, char *argv[]) {
  int nIter;
  double c,Tmax,dT,T,Tc;
  
  if (argc != 4&&argc!=5) {
    fprintf(stderr,
	    "usage: %s N c Niter [seed]\n",
	    argv[0]);
    exit(1);
  }
  N = atoi(argv[1]);
  c = atof(argv[2]);
  nIter = atoi(argv[3]);
  
  if(argc==4){  //se non ho messo seed tra gli argomenti lo pesco da un file di numeri random
    FILE *devran = fopen("/dev/urandom","r"); 
    fread(&myrand, 4, 1, devran);
    fclose(devran);
  }
  else myrand=atoi(argv[4]);  // altrimenti uso seed passatomi in input

  if (myrand == 2147483647) {
    fprintf(stderr, "Error: seed must be less than 2147483647\n");
    exit(1);
  }

  M = (int)(0.5 * c * N + 0.5);

  printf("# N = %u   c = %g   M = %d  seed = %u\n", N, c, M, myrand);
  fflush(stdout);
  
  printf("# 1:T 2:m 3:t_convergence\n");

  Init_Random();

  allocateMem();

  makeGraph();

  Tc=1./atanh(1./c);
  dT=1.5*Tc/20;
  Tmax=1.5*Tc; 
  for(T=dT;T<=Tmax;T+=dT){
    beta=1./T;
    init_u();        //inizializzo ogni volta perche' voglio far vedere che sotto Tc BP sceglie in quale stato entrare
    iteration(nIter);
  }
  
  return 0;
}
