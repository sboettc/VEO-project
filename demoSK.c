/* gcc -O3 -lm demoSK.c -o demoSK */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#define min(A,B) ((A) < (B)) ? (A) : (B)
#define max(A,B) ((A) > (B)) ? (A) : (B)
#define MAXVERTEX 1024     /* maximal number of vertices in a graph */


//Type Definitions
typedef int EnergyType;

struct Vertex{
  int index;
  int spin;
  EnergyType fitness;
};

//global variables
/**** Random Number Generator Stuff  ******/
static const int half_int = 1073741824;
static const int rmax = 2147483647;
static const double inv_rmax = .4656612875245797e-9; // =1/rmax 
/**** tau-EO features *****/
float tau;
double tauscale;
short int pdf[MAXVERTEX];
/**** Sherrington-Kirkpatrick Spin Glass ****/
short int magnit,nv;//nv=system size "n"
struct Vertex bestV[MAXVERTEX],V[MAXVERTEX],*p[MAXVERTEX];
short int bonds[MAXVERTEX][MAXVERTEX];
/**** Update Time ****/
int time;


/******************************************************************/
/*                                                                */
/*           Random Number Generator                              */
/*                                                                */
/* This one does the multiplicative LCG of Lewis, Goodman and     */
/* Miller, using the Schrage modulus trick                        */
/*                                                                */
/******************************************************************/
int rng1(int n){
  static const int ia = 16807;
  static const int iq = 127773;
  static const int ir = 2836;
  static int irandom;
  int l;

  if(n > 0) irandom = n;
  l = irandom/iq;
  irandom = (irandom-iq*l)*ia - ir*l;
  if(irandom < 0) irandom = irandom+rmax;
  return(irandom);
}
/************  END Random Number Generator ***********************/

/************  BEGIN of Initialization  ***********************/

/*****************************************************************/
/*   Sets up connectivity matrix with random +/-J bonds          */
/*                                                               */
/*   Although matrix is symmetric, its better to fill matrix.    */
/*   This doubles memory (no bottleneck), but increases          */
/*   speed (big bottleneck).                                     */
/*                                                               */
/*   Returned is the imbalance between +J and -J bonds.          */
/*                                                               */
/*                                                               */
/*****************************************************************/
short int initbonds(void){
  short int i,j,balance=0;
  for(i=0; i<nv; i++){
    for(j=0; j<i; j++)
      balance += bonds[j][i] = bonds[i][j] = (rng1(0) > half_int) ? +1 : -1;
    bonds[i][i] = 0;
  }
  return(balance);
}

/*****************************************************************/
/*   Initiallizes Spin configuration at random.                  */
/*                                                               */
/*   Returned is the magnitization of the configuration.         */
/*                                                               */
/*****************************************************************/
void initspins(void){
  short int i;

  magnit = 0;
  for(i=0; i<nv; i++){
    V[i].spin = (rng1(0) > half_int) ? +1 : -1;
    magnit += V[i].spin;
  }
}


/*****************************************************************/
/*   Initiallizes Fitness for each spin variable                 */
/*                                                               */
/*   "Fitness" here is (twice) the local energy contribution     */
/*   of each spin, ie its good bond weights minus its bad ones.  */
/*                                                               */
/*   Returned is the total energy of the configuration.          */
/*                                                               */
/*****************************************************************/
EnergyType initfitness(void){
  short int i,j;
  EnergyType cost,energy=0;

  for(i=0; i<nv; i++)V[i].fitness = 0;
  for(i=0; i<nv; i++)
    for(j=0; j<i; j++){
      cost = (EnergyType)(V[i].spin == V[j].spin) ? bonds[i][j] : -bonds[i][j];
      V[i].fitness += cost;
      V[j].fitness += cost;
      energy -= cost;
    }
  return(energy);
}


/*****************************************************************/
/*   Initiallizes "ordered" Tree (or "Heap") of Fitnesses        */
/*                                                               */
/*   Initially, the heap is not ordered at all, but just filled  */
/*   sequentially. Order is established DYNAMICALLY with each    */
/*   update.                                                     */
/*                                                               */
/*****************************************************************/

void initlist(void){
  int i;
  for(i=0; i<nv; i++){
    V[i].index = i;//just to make sure!!!
    p[i] = &V[i];
  }
}

/******************************************************************/
/*                                                                */
/*   Initializes tau-probability distribution                     */
/*                                                                */
/*  This provides a discrete approximation for P(k)~k^{-tau},     */
/* 1 <= k <= nv. It maps the distribution onto a large array where*/
/* entries = k occur with frequency in occordance with P(k).      */
/* Shortcoming: if "LENGTH" too small, and/or "tau" too large,    */
/* there will be only values of k<n_0<nv for some cutoff n_0, ie. */
/* certain ranks can't be reached!                                */
/*                                                                */
/*   This is a revised version, where we sort on a heap, mapping  */
/*   P(k) onto Q(l)\sim 2^{-(tau-1)l}, l=0..[ln_2(nv)], where "l" */
/*   refers to a level in the heap, drawn at random from Q(l).    */
/*   For this LENGTH=MAXVERTEX >> ln(MAXVERTEX) suffices.         */
/*                                                                */
/******************************************************************/

void initrank(float tau){
  int i;
  float a,b;
  b = (tau-1)*log(2.0);
  a = (1-exp(-b*(int)(1+log((float)nv)/log(2.0))))/MAXVERTEX;
  for(i = 0; i<MAXVERTEX; i++){
    pdf[i] = (short int)(-log(1-i*a)/b);
  }
}

/************  END of Initialization  ***********************/

/************  BEGIN of Update Procedures  ***********************/

/*****************************************************************/
/*   Picks a spin to update from the heap according to Q(l)      */
/*                                                               */
/*   Somewhat sleazy: recycles old random number to pick level   */
/*   "l" then gets new one to pick a variable within "l".        */
/*   'Hope consecutive random numbers are uncorrelated!!!        */
/*                                                               */
/*   Returned is the index of the vertex (=spin) to be updated   */
/*                                                               */
/*****************************************************************/
short int rank(void){
  static unsigned int oldrandom;
  unsigned int random,base,k;
  
  do{
    random = (unsigned int)((rng1(0)*inv_rmax)*MAXVERTEX);
    base = ldexp(1,pdf[oldrandom])-1; /* =2^level !? */
    oldrandom = random;
    k = base+(random & base);
  }while(k >= nv);//needed if lowest level not complete!!!
  return((short int)k);
}


/*****************************************************************/
/*   Flip selected spin and update magnitization                 */
/*                                                               */
/*****************************************************************/
void flip(struct Vertex *v){
  v->spin = -(v->spin);
  magnit += 2*v->spin;
}

/*****************************************************************/
/*   Updates Fitness for each spin after spin "k" was flipped    */
/*                                                               */
/*   Here we exploit that we have the full bond matrix           */
/*   allowing for an unrestricted loop.                          */
/*                                                               */
/*   Returned is the change to the total energy due to the update*/
/*                                                               */
/*****************************************************************/
EnergyType updatefitness(short int k){
  short int i;

  V[k].fitness = -V[k].fitness;
  for(i=0; i<nv; i++)
    if(V[i].spin == V[k].spin)  
      V[i].fitness += (EnergyType)(bonds[k][i]*2);
    else
      V[i].fitness -= (EnergyType)(bonds[k][i]*2);
  return(-2*V[k].fitness);
}


/*****************************************************************/
/*   Updates the Fitnesses "ordered" on a tree.                  */
/*                                                               */
/*   This is tricky: Each update ALL fitnesses change, but only  */
/*   by O(1), while the flipped spin (and the total energy)      */
/*   can change by O(n). So instead of perfectly ordering the    */
/*   tree (or "heap") each update, we merely parse through all   */
/*   internal nodes and exchange them with one of their two leafs*/
/*   if that leaf has smaller fitness. Bad fitness at the bottom */
/*   trickle to the top after O(log(n)) updates, at least,       */
/*   establishing order DYNAMICALLY. Note that good fitnesses    */
/*   at top move down and get protected IMMEDIATELY, though.     */
/*                                                               */
/*****************************************************************/
void updatelist(void){
  short int n,leaf=0,node=0;
  struct Vertex *help;

  for(;;){
    for(n=0; n<2; n++){
      if(++leaf < nv){
	if(p[node]->fitness > p[leaf]->fitness){
	  help = p[node];
	  p[node] = p[leaf];
	  p[leaf] = help;
	}//end if p
      }else
	return;
    }//end for n
    ++node;
  }//end for()
}

/************ END of Update Procedures  ***********************/

/************ BEGIN of Input/Output Stuff  ***********************/

/*******************************************************************/
/*                                                                 */
/*  This is a bunch of functions making the user interface pretty  */
/*                                                                 */
/*******************************************************************/
void banner(void){
  printf("###############################################################\n");
  printf("#                                                             #\n");
  printf("#                 Extremal Optimization                       #\n");
  printf("#                      for the                                #\n");
  printf("#              Sherrington-Kirkpatrick Model                  #\n");
  printf("#                                                             #\n");
  printf("#                                                             #\n");
  printf("#   This is a Demo for finding Ground-State Energies of the   #\n");
  printf("#   SK-Model with the Extremal Optimization (EO) Heuristic,   #\n");
  printf("#   based on publication Eur. Phys. J. B 46, 501-505 (2005),  #\n");
  printf("#   also at arxiv.org as cond-mat/0407130.                    #\n");
  printf("#                                                             #\n");
  printf("#   You will be able to see EO ordering fitnesses on a tree,  #\n");
  printf("#   print out the current configuration and the state of      #\n");
  printf("#   the optimization so far. You can print out some or all of #\n");
  printf("#   this info each step and hold (slow), or only each time a  #\n");
  printf("#   new minimum E_min was found (fast).                       #\n");
  printf("#                                                             #\n");
  printf("#   You need to specify the size n of an instance, random     #\n");
  printf("#   seeds, a parameter tau for EO, and a number of update     #\n");
  printf("#   steps. Then you can choose how much info you want to      #\n");
  printf("#   print out for that run.                                   #\n");
  printf("#                                                07-10-04     #\n");
  printf("###############################################################\n");
}

int integer_dialog(char request[],int DefaultValue){
  char buffer[120];
  int value;
  do{
    printf("%s (default=%d)",request,DefaultValue) ;
    gets(buffer);
    if(buffer[0] == '\0') return(DefaultValue);
    if(sscanf(buffer,"%d",&value) == 1) return(value);
    printf("I don't understand this input!\n");
  }while (1);
}

float float_dialog(char request[],float DefaultValue){
  char buffer[120];
  float value;
  do{
    printf("%s (default=%5.3f)",request,DefaultValue) ;
    gets(buffer);
    if(buffer[0] == '\0') return(DefaultValue);
    if(sscanf(buffer,"%f",&value) == 1) return(value);
    printf("I don't understand this input!\n");
  }while (1);
}

void printbonds(void){
  int i,j;
  printf("=================Begin Bondmatrix========================================\n"); 
  printf(" |");
  for(j=0; j<nv; j++) printf("%1d",j%10);
  printf(" Spin Fitness");
  printf("\n-+");
  for(j=0; j<nv; j++) printf("-");
  printf("-------------");
  printf("\n");
  for(i=0; i<nv; i++){
    printf("%1d|",i%10);
    for(j=0; j<nv; j++)
      if(bonds[i][j] == 0) printf("0");
      else if(bonds[i][j] > 0) printf("+");
      else printf("-");
    printf(" %4d %7d",bestV[i].spin,bestV[i].fitness);
    printf("\n");
  }
  printf("================= End Bondmatrix ========================================\n");
}

void input(int *maxtime,float *tau,int *printheap,int *printconfig,int *stop_and_go){
  unsigned int seed1,seed2;

  printf("\n########## Input System Settings #############################\n");
  printf("\nInput System Size n(<%d):\n",MAXVERTEX);
  nv = integer_dialog(">> n=",31);
  seed1 = getpid(); // Create a random seed from process-id
  printf("\nInput Random Seed for +/-J Bonds:\n");
  seed1 = integer_dialog(">> bond_seed=",getpid());
  rng1(seed1);
  initbonds();
  printf("\nInput Random Seed for initial Spin Configuration:\n");
  seed2 = integer_dialog(">> spin_seed=",rng1(0)+getpid());
  rng1(seed2);
  initspins();
  printf("\n########## Input EO Settings #############################\n");
  printf("\nInput tau paramter for choosing ranks\n");
  *tau = float_dialog(">> tau=",1.25);
  initrank(*tau);
  initlist();
  printf("\nInput Number of EO update steps (~n^3/10)\n");
  *maxtime = integer_dialog(">> runtime=",pow(nv,3)/10);
  printf("\n########## Output Control #############################\n");
  printf("\nDo you want to stop and print each update?");
  *stop_and_go =  integer_dialog("[Yes=1,no=0]",1);
  if(*stop_and_go){
    printf("\nDo you want to print current configuration each update?");
    *printconfig =  integer_dialog("[Yes=1,no=0]",1);
    printf("\nDo you want to print fitness heap each update?");
  if(nv > 31) printf("(Not pretty for n>31!)");
    *printheap =  integer_dialog("[Yes=1,no=0]",1);
  }
}

void lspace(int blanks){
  int i;
  for(i=0; i<blanks-4; i++) printf(" ");
  if(blanks > 3){
    printf("/");
    for(i=0; i<blanks-1; i++) printf("-");
  }
}

void rspace(int blanks){
  int i;
  if(blanks > 3){
    for(i=0; i<blanks-1; i++) printf("-");
    printf("\\");
  }
  for(i=0; i<blanks-4; i++) printf(" ");
}


void printout(int printconfig,int printheap){
  int i,j,width;
  
  if(printconfig){
    printf("========================Spin Configuration with Fitness============\n");
    for(i=0; i<nv; i++)printf("%3d",V[i].index);
    printf("\n");
    for(i=0; i<nv; i++)printf("%3d",V[i].spin);
    printf("\n");
    for(i=0; i<nv; i++)printf("%3d",V[i].fitness);
    printf("\n");
  }
  
  if(printheap){
    printf("========================Fitness Tree with [vertex#,fitness]========\n");
    printf("\n");
    i = 0;
    j = 2;
    width = 4*ldexp(1,(int)(log((float)nv)/log(2.0)));
    do{
      width /= 2;
      do{
	lspace(width);
	printf("[%2d,%3d]",p[i]->index,(int)(p[i]->fitness));
	rspace(width);
      }while(++i < (min(j-1,nv)));
      printf("\n");
      j *= 2;
    }while(i < nv);
    printf("\n");
  }
}
	 
int main(void){
  int i,k,flag;
  int mintime;
  int maxtime,printheap,printconfig,stop_and_go;
  struct Vertex *worst;
  EnergyType energy,minenergy;

  banner();
  input(&maxtime,&tau,&printheap,&printconfig,&stop_and_go);
  energy = minenergy = initfitness();
  time = mintime = 0;
  for(time=0; time<maxtime; ++time){
    worst = p[k=rank()];
    flip(worst);
    energy += updatefitness(worst->index);
    updatelist();
    if(energy < minenergy){
      minenergy = energy;
      mintime = time;
      printf(">>> New Optimum: n=%4d, tau=%4.3f, E_min=%8d obtained at time=%12d\n",nv,tau,minenergy,mintime);
      for(i=0; i<nv; i++){
	bestV[i].spin = V[i].spin;
	bestV[i].fitness = V[i].fitness;
      }
    }
    if(stop_and_go){
      printf("\nEO:  n=%4d, tau=%4.3f,      E_min=%6d obtained at time=%8d\n",nv,tau,minenergy,mintime);
      printf("           chosen rank=%3d,  E_now=%6d,    current time=%8d\n",k,energy,time);
      printout(printconfig,printheap);
      printf("Continue? [q=Quit]");
      if(getchar() == 'q'){
	getchar();
	goto endgame;
      }
      printf("\n+++++++++++++++++++ NEW UPDATE +++++++++++++++++++++++++++\n");
    }//end if stop
  }//end for time
 endgame:
  printf("\nDo you want to print the best-found solution?");
  if(flag=integer_dialog("[Yes=1,no=0]",1)) printbonds();
  printf("Lowest Energy State found:\n"); 
  if(flag) for(i=0; i<nv; i++)printf(" ");
  printf(" E_min =%7d",minenergy);
  if(flag) printf(" (= -sum of fitnesses/2)");
  printf(", or   E_min/n^1.5 =%8.5f\n",(float)minenergy*pow(nv,-1.5));
  return(0);
}
