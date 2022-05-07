# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <stdbool.h>
# include <time.h>
# include <mpi.h>

// ref http://www.tsm-resources.com/alists/mers.html

// Set below to 100 for testing your code
// When ready to run set this to 1 million 1000000
// submit one code called prob1_mpi.c by modifying
//        the code to use MPI and run faster
//        All cores must do work
//        buffer tasks for each core
//        Use 2 nodes and all cores on the node.
// submit one called prob1_omp.c by modifying
//        the code to use openmp and run faster
//        Use 1 node and all cores on the node

#define MAXPRIME 1000000

int main ( int argc, char **argv );
bool is_prime( int n);
void make_prime_vector(int n, int *prime, int *k);
bool quick_is_prime(unsigned long long int j, int *prime, int k);

int main ( int argc, char **argv )
{
  MPI_Init(NULL, NULL);

  int MAXK;
  int MPSize;
  int rank;

  MPI_Comm_size(MPI_COMM_WORLD, &MPSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int PRIMECHUNK = MAXPRIME/MPSize;
  int R = MAXPRIME%MPSize;
  if(rank < R){
    PRIMECHUNK++;
  }
  //printf("rank %d, R: %d\n", rank, R);
  if(rank == MPSize-1){
    MAXK = PRIMECHUNK - 2;
  }else{
    MAXK = PRIMECHUNK;
  }
  int prime[PRIMECHUNK];
  for(int i = 0; i < PRIMECHUNK; i++){
    prime[i] = 0;
  }
  int mersenne[64];
  int k, i;
  // Initial Prime values
  if(rank == 0){
  prime[0]=2;
  prime[1]=3;
  prime[2]=5;
  prime[3]=7;
  prime[4]=11;
  prime[5]=13;
  k = 6;
  }else{
    k = 0;
  } 
  unsigned long long int j;
  int n, m;
  // printf("rank %d, MAXK: %d PRIMECHUNK %d\n", rank, MAXK, PRIMECHUNK);
  clock_t begin = clock();
  if(rank == 0){
    // Create prime vector
    n=17; // starting prime - skip even and factor of 6
    while ( k<MAXK){
      make_prime_vector(n, prime, &k);
      make_prime_vector(n+2, prime, &k);
      n=n+6;
      /* The below many be helpful for debugging only
	 if ((n-17)%10000==0)
	 printf("n=%d k=%d\n",n,k);
      */
    }
    
  }else{
    
    int maxCount = 0;
    if(rank < R){
      maxCount = (PRIMECHUNK+1)*rank;
    }else{
      maxCount = (PRIMECHUNK+1)*R + PRIMECHUNK*(rank-R);
    }

    n = maxCount - ((maxCount - 17)%6);
    n += 6;
    while ( k<MAXK){
      if(is_prime(n)){
	prime[k] = n;
	k++;
      }
      if(is_prime(n+2)){
	prime[k] = n;
	k++
      }
      n=n+6;
      /* The below many be helpful for debugging only
	 if ((n-17)%10000==0)
	 printf("n=%d k=%d\n",n,k);
      */
    }  
   
  }
  
  int isum;
  long int sum=0;
  for (isum=0; isum<k; ++isum){
    if ( sum>1000000000)
      sum=sum-prime[isum];
    else
      sum=sum+prime[isum];

  }
  
  int *tempPrime = (int *)malloc(k*sizeof(int));
  for(int i = 0; i < k; i++){
    tempPrime[i] = prime[i];
  }

  //printf("rank %d TempPrime %d\n", rank, tempPrime[1]);
  int recvCount[MPSize];
  int displace[MPSize];
  displace[0] = 0;
  
  int newPrimeSize = 0;
  //printf("rank: %d, k: %d\n", rank, k);
  MPI_Reduce(&k, &newPrimeSize, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
 
 
  MPI_Allgather(&k, 1, MPI_INT, recvCount, 1, MPI_INT, MPI_COMM_WORLD);

  int *newPrime = (int *)malloc(newPrimeSize*sizeof(int)+2);

  for(int i = 1; i < MPSize; i++){
    displace[i] = recvCount[i];
    displace[i] += displace[i-1];
  }
 

  //if(rank == 0) printf("Working up to gatherv");

  MPI_Gatherv(tempPrime, k, MPI_INT, newPrime, recvCount, displace, MPI_INT, 0, MPI_COMM_WORLD);
  /*if(rank == 0){
    printf("rank %d newPrime %d recvCount %d\n", rank, newPrime[1800], recvCount[1]);
    }*/
  
  //printf("rank: %d k=%d prime[k-1]=%d \n\n",rank,k,prime[k-1]);
  clock_t end = clock();

  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  if(rank == 0){
    printf("time creating prime vector %f \n",time_spent);
  }

 
    
  if(rank == 0){
    n=0;
    //printf("k = %d\n", k);
    for (int i=2; i<64; ++i){
      j=(unsigned long long int)pow(2,i)-1;
      // printf("Working, newPrime[1]: %d\n", newPrime[1]);
      if(quick_is_prime(j, newPrime ,k)){
	mersenne[n]=i;
	++n;
      }
    }
  
  clock_t end2 = clock();
  double time_spent2 = (double)(end2 - end) / CLOCKS_PER_SEC;
  // For mpi only core zero prints the timing
  printf("time creating mersenne primes %f \n",time_spent2);
  }


  //output
  // For mpi only core zero prints this output
   if(rank == 0){
    for (int i=0; i<n; ++i){
      j=(unsigned long long int)pow(2,mersenne[i])-1;
      printf("2^(%d)-1 = %llu \n",mersenne[i],j);
    }
   }
  



  // comment the below for mpi code
  //printf("prime[%d]=%d\n",k-1,prime[k-1]);
  
  /* 
  int isum;
  long int sum=0;
  for (isum=0; isum<k; ++isum){
    if ( sum>1000000000)
      sum=sum-prime[isum];
    else
      sum=sum+prime[isum];
  }
  */

  //comment the below for MPI code
  //printf("sum=%d\n",sum);
  // uncomment the below for mpi code where rank
  // Is the rank for each core
  // All cores must print this.
  
  printf("[%d] prime[%d]=%d\n",rank,k-1,prime[k-1]);
  printf("[%d] sum=%d\n",rank,sum);
  MPI_Barrier(MPI_COMM_WORLD);   
  MPI_Finalize();

}

bool is_prime(int n){
  if (n <= 3)
    return(n > 1);
  if (n%2 == 0 || n%3 == 0)
    return(false);
  int i=5;
  while (i * i <= n){
    if (n%i == 0 || n%(i + 2) == 0)
      return(false);
    i=i+6;
  }
  return (true);
}

void make_prime_vector(int n, int *prime, int *k){

  int i=2;
  while ( i<(*k) ){
    if (n%prime[i]==0)
      return;
    ++i;
  }
  prime[(*k)]=n;
  ++(*k);
  return;
}


bool quick_is_prime(unsigned long long int j, int *prime, int k){
  //printf("j is %lld k is %d prime[1] is %d\n", j, k, prime[1]);
  int i=1;
  while ( i<k ){
    if (j%(unsigned long long int)prime[i]==0)
      return(false);
    //printf("FirstCheckSucceeded\n");
    if ((unsigned long long int)prime[i]*(unsigned long long int)prime[i]>j){
      return(true);
    }
    //printf("SecondCheckSucceeded\n");
    ++i;
  }

  unsigned long long int ii;
  ii=(unsigned long long int)(prime[k-2]+6);
  while (ii * ii <= j){
    if (j%ii == 0 || j%(ii + 2) == 0)
      return(false);
    ii=ii+6;
  }
  return (true);

}
