//
// cs4900-01
// Project 04
// Brendon Deal
// Due 28 Feb 2020
// System = owens.osc.edu
// Compiler syntax = ./compile.sh proj04
// Job Control File = proj04.pbs
// Results file     = proj04.txt
//


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>


int main(int argc, char** argv) {

  // Standard MPI initialization calls
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int rank, commies, name_len;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commies);
  MPI_Get_processor_name(processor_name, &name_len);

  // Good to know what node and core this rank is running on
  printf("[%d]processor %s, out of %d processors\n", rank, processor_name, commies);
  MPI_Barrier(MPI_COMM_WORLD);

  if (argc <= 4 && rank==0) {
    fprintf(stderr, "%d\n",argc);
    fprintf(stderr, "Usage: This project requires four inputs:\n");
    fprintf(stderr, "The first two variables are vector lengths\n");
    fprintf(stderr, "The second two variables are coefficients\n");
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(1);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Inputs
  int L1 = atoi(argv[1]);
  int L2 = atoi(argv[2]);
  int A = atoi(argv[3]);
  int B = atoi(argv[4]);

  printf("L1: %d \n L2: %d \n A: %d \n B: %d \n", L1, L2, A, B);
  // Constants
  float dx, dy,  pi, norm, norm2;
  pi=acos(-1.0);
  dx=2.0*pi/L1;
  norm=sqrt(2.0/L1);
  dy=2.0*pi/L2;
  norm2=sqrt(2.0/L2);
  int toCalcL1 = 0;
  int toCalcL2 = 0;
  // Calculates which processor has what part of the global elements
  // Block row Partitioned
  float xlast, xlast2, ylast;
  int Nr=(int)(L1/commies);   // All processors have at least this many rows
  int R = (L1%commies);
  if(rank > L1-1 && rank > L2 - 1){
    MPI_Finalize();
  }
  if (rank<R) {
    ++Nr;             // First R processors have one more element
    toCalcL1 = Nr;    // number of places to compute 
    xlast = -dx/2.0+rank*Nr*dx;
  }else{
    toCalcL1 = Nr;
    xlast=-dx/2.0+R*(Nr+1)*dx+(rank-R)*Nr*dx;
  }

  Nr=(int)(L2/commies);   // All processors have at least this many rows
  R = (L2%commies);
  if (rank<R) {
    ++Nr;             // First R processors have one more element
    toCalcL2 = Nr;    // number of places to compute 
    xlast2 = -dx/2.0+rank*Nr*dx;
    ylast = -dy/2.0+rank*Nr*dy;
  }else{
    toCalcL2 = Nr;
    xlast2=-dx/2.0+R*(Nr+1)*dx+(rank-R)*Nr*dx;
    ylast = -dy/2.0+(Nr+1)*dy + (rank-R)*Nr*dy;
  }



  printf("Rank: %d xlast: %f\n", rank, xlast);
  float *vector = (float*)malloc(sizeof(float) * toCalcL1);
  assert(vector != NULL);
  float xi = xlast;
  for(int i = 0; i < toCalcL1; i++){
    xi += dx;
    vector[i] = norm*cos(A*xi);
    //printf("vector1 Rank: %d, %d: %f\n", rank, i, vector[i]);
  }
  xi = xlast2;
  float *vector2 = (float*)malloc(sizeof(float) * toCalcL1);
  assert(vector2 != NULL);
  for(int i = 0; i < toCalcL1; i++){
    xi += dx;
    vector2[i] = norm*sin(B*xi);
    // printf("vector2 %d: %f\n", i, vector2[i]);
  }
  

  //THIS IS NOT WORKING, FIX IT.

  float *vector3 = (float *)malloc(sizeof(float) * toCalcL2);
  assert(vector3 != NULL);
  float yi = ylast; 
  for(int i = 0; i < toCalcL2; i++){
  yi += dy;
  vector3[i] = norm2*cos(B*yi);
  //printf("vector3 %d: %f\n", i, vector3[i]);
 }
  float *v1 = malloc(sizeof(float)*L1);
  printf("working %d\n", rank);
  int prevCalc = 0;
  int totalSend = 0;


  //BIG HAPPY FUN TIME v1 ASSEMBLY AND DISTRIBUTION
  MPI_Barrier(MPI_COMM_WORLD);
  if(R == 0){
    MPI_Allgather(vector, toCalcL1, MPI_FLOAT, v1, toCalcL1, MPI_FLOAT, MPI_COMM_WORLD); //make sure all of the Nodes have vector 1 
  }else{
    if(rank == 0){
      totalSend += toCalcL1;
      MPI_Send(&totalSend, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
      MPI_Send(vector, totalSend, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
    }else{
      MPI_Recv(&totalSend, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(v1, totalSend, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int i = 0; i < toCalcL1; i++){
	printf("adding %f to column %d of v1 in rank %d\n", vector[i], totalSend+i, rank);
	v1[totalSend+i] = vector[i];
      }
      totalSend += toCalcL1;
      if(rank < commies-1){
	MPI_Send(&totalSend, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	MPI_Send(v1, totalSend, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
      }	      
    }
  }
  if(rank == commies-1){
    printf("Rank: %d broadcasting v1 with a total of %d items in it", rank, totalSend);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&totalSend, 1, MPI_INT, commies-1, MPI_COMM_WORLD);
  MPI_Bcast(v1, totalSend, MPI_FLOAT, commies-1, MPI_COMM_WORLD);
 
  //test vectors two and three
  for(int i = 0; i < toCalcL1; i++){
    printf("Rank: %d v2: %d %f\n", rank, i, vector2[i]);
    printf("Rank: %d v3: %d %f\n", rank, i, vector3[i]);
  }


  printf("working2 %d\n", rank);
  /*for(int i = 0; i < L1; i++){
    printf("Rank: %d v: %d %f\n", rank, i, v1[i]);
  }*/ 
  
  float sum = 0.0;
  //make sure the dot product multiplication will be done starting on the right index
  if(rank == 0){
    printf("R = %d", R);
  }
 
  //create an offset and begin adding up the dot multiplication based on the chunks of vectors 2 and 3 known.
 int adder;
  if(rank >= R){
    adder = R*(Nr+1)+(rank-R)*Nr;
    printf("Rank: %d, adder: %d\n", rank, adder);
  }else{
    adder = rank*(Nr+1);
  }
  for (int i = 0; i < toCalcL2; i++) {
    // printf("sum %f is now adding %f on rank %d which is made of %f times %f\n", sum, v1[adder]*vector2[i], rank, v1[adder], vector2[i]);
    sum+= v1[adder]*vector2[i];
    
    adder += 1;
  }
    float *finalsum = malloc(sizeof(float));
    // printf("Sending sum: %f\n", sum);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&sum, finalsum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank == 0){
      printf("final sum: %f\n", finalsum);
    }

 //create the matrix
 float matrix[toCalcL2][L1];



 adder -= toCalcL2;
  
  for(int i = 0; i < toCalcL2; i++){
    for(int j = 0; j < L1; j++){
      matrix[i][j] = vector3[i]*v1[j];
      //printf("Rank: %d at %d %d is %f\n", rank, i, j, matrix[i][j]);
    }
  }
  float diagSum = 0.0;
  float nonDiagSum = 0.0;
  if(rank==0){
   for(int i = 0; i < toCalcL2; i++){
    for(int j = 0; j < L1; j++){
      if(i==j){
	printf("Rank: %d is adding %f to diag sum total\n", rank, matrix[i][j]);
	diagSum+= matrix[i][j];
	printf("Rank: %d new total is %f\n", rank, diagSum);
      }else{
	nonDiagSum += matrix[i][j];
      }
    }
   }
  }else{
    for(int i = 0; i < toCalcL2; i++){
      for(int j = 0; j < L1; j++){
	if(i+adder == j){
	  printf("Rank: %d i is %d adder is %d j is %d\n", rank, i, adder, j);
	  printf("Rank: %d is adding %f to diag sum total\n", rank, matrix[i][j]);
	  diagSum+= matrix[i][j];
	  printf("Rank: %d new total is %f\n", rank, diagSum);
	}else{
	  nonDiagSum += matrix[i][j];
	}
      }
    }
  }

  printf("\n\nRank: %d Final diag sum: %f Final non Diag Sum %f\n", rank, diagSum, nonDiagSum);
  float diagSumFinal = 0.0;
  float nonDiagSumFinal = 0.0;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&diagSum, &diagSumFinal, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&nonDiagSum, &nonDiagSumFinal, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);  
 // if(rank == 0){
  printf("Rank: %d  Final nondiagonal sum: %f\n", rank, &nonDiagSumFinal);
  printf("Rank %d  Final diagonal sum: %f\n", rank, &nonDiagSumFinal);
      // }
//  printf("dot product: %f\n", &finalsum);
  MPI_Finalize();
}

