/*******************************************************************************
 * This file contains a simple test, run and verify for LU Decomposition of 
 * a square Matrix given a block size
 * 
 * Author: Raghav Mohan
 * 
 *******************************************************************************/

#include "matrix.hpp"
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string>
#include <getopt.h>
#include <cassert>
#define NUM_TIMES (1)
float marginErr = 1e-10f;





//global variables

//int numThreads;
int blockSize;
bool verify;
const char *input_file = NULL;


int matrixSize;
Matrix * L, *U, *LU;


void test1();

//method declarations
void getoptions (int argc, char **argv);
void test1();
//void verifyTest1();

int randGen(int tMin,int tMax){
  int f;
  f= 0;
  while(f ==0)
    f =  rand() % (int) tMax;
  return f;
}

//use this if Boost is installed only
/*
   template<typename T>
   T randGen(T tMin,T tMax){
   T f;
   if(boost::is_same<T, double>::value){
   f= 0;
   while(f == 0){
   f = (double)rand() / RAND_MAX;
   f =  tMin + f * (tMax - tMin);
   }
   }
   else{
   f= 0;
   while(f ==0)
   f =  rand() % (int) tMax;
   }
   return f;
   }
   */

using namespace std;
  double 
fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

int intRand(int fMin, int fMax)
{
  int f = (int)rand() % fMax;
  return f;
}


int
create_matrix_from_file(Matrix * LUm, const char* filename, int *size_p){
  int size;
  FILE *fp = NULL;

  fp = fopen(filename, "rb");
  if ( fp == NULL) {
    return -1;
  }

  fscanf(fp, "%d\n", &size);

  if (LUm == NULL) {
    fclose(fp);
    return -1;
  }

  float val;
  for (int i=0; i < size; i++) {
    for (int j=0; j < size; j++) {
      fscanf(fp, "%f ", &val);
      LUm->setElement(i, j, Element(val));
    }
  }


  fclose(fp);


  *size_p = size;
  return 1;
}


void
getoptions (int argc, char **argv) {
  int c;
  while ((c = getopt (argc, argv, "v?hRs:e:A:a:B:b:C:c:D:d:Vo:r:")) != -1){
    switch (c) {
      case 'h':
        std::cout<<"Usage: matrix_ss <matrixsize> <blocksize> <inputfile> [options]\n\toptions\n\t\ti: input file\n\t\tb: blocksize (eg 16 for 4 16x16 blocks for a 256x256 matrix)\n\t\ts: matrix size (eg 16 for a 16x16 matrix\n\t\tv: No verify"<<std::endl;
        exit(1);
        break;
      case 'v':
        verify=false;
        break;
      case 'i':
        input_file=optarg;
        break;
      case 'b':
        blockSize=atoi(optarg);
        break;
      case 's':
        matrixSize=atoi(optarg);
        break;

      default:
        break;
    }
  }
} 




int 
main(int argc, char** argv,char ** envp){
  //total time
  struct timeval startTotal, endTotal; 
  gettimeofday(&startTotal, NULL);

  //set in case not defined
  blockSize=2;
  matrixSize =8;
  verify = true;

  getoptions(argc, argv);
  test1();
  //total timing
  gettimeofday(&endTotal, NULL);

  long mtime, seconds, useconds;
  seconds  = endTotal.tv_sec  - startTotal.tv_sec;
  useconds = endTotal.tv_usec - startTotal.tv_usec;
  mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
  cout << "\nTotal running time : "<< mtime<< " ms."<<endl; 
  return 0;
}


//Test Block LU
void 
test1(){
  std::cout << "Running Test Block LU"<<endl;
  std::cout << "Matrix Size: " << matrixSize <<endl;
  std::cout << "Block Size: "<< blockSize << endl;


  //declarations
  Matrix * Lm = new Matrix(matrixSize, matrixSize);
  Matrix * Um = new Matrix(matrixSize, matrixSize);
  Matrix * LUm;
  LUm =  new Matrix(matrixSize, matrixSize);

  long double diffVal = 0;

  if(input_file == NULL){
    for(int i= 0 ; i < matrixSize; ++i){
      for(int j= 0 ; j < matrixSize; ++j){
        if(i == j){
          Lm->setElement(i, j,Element( (double) 1.0 ));
          Um->setElement(i, j,Element( (double) randGen(-20, 20)));
        }
        else if(i < j){
          Lm->setElement(i, j,Element( (double) 0. ));
          Um->setElement(i, j,Element( (double) randGen(-20, 20)));
        }
        else{
          Lm->setElement(i, j,Element( (double) randGen(-20, 20)));
          Um->setElement(i, j,Element( (double) 0. ));
        }

        LUm->setElement(i, j,Element( (double) randGen(-20, 20)));
        Lm->setElement(i, j,Element( (double) randGen(-20, 20)));
      }
    }
  }
  else{
    printf("Reading matrix from file %s\n", input_file);
    int test = 0;
    int ret = create_matrix_from_file(LUm, input_file,&test);
    if (ret != 1) {
      LUm = NULL;
      fprintf(stderr, "error create matrix from file %s\n", input_file);
      exit(0);
    }
  }



  Matrix * Lresult = new Matrix(matrixSize, matrixSize);
  Matrix * Uresult = new Matrix(matrixSize, matrixSize);
  Matrix * LUresult =  new Matrix(matrixSize, matrixSize);

  LUm->partition();
  struct timeval start, end; 
  gettimeofday(&start, NULL);
  LUm->blockLU();
  gettimeofday(&end, NULL);

  long mtime, seconds, useconds;
  seconds  = end.tv_sec  - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;
  mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

  LUm->getBlockResults(Lresult, Uresult);

  cout << "\nTook "<< mtime<< " milliseconds to compute."<<endl; 

  if(verify){
    cout << "\nStarting verification ... "<<endl; 

    (*LUresult) = (*Lresult) * (*Uresult);

    LUm->freeBlocks();

    for(int i= 0 ; i < matrixSize ; ++i){
      for(int j= 0 ; j < matrixSize; ++j){
        if(!(LUresult->getElement(i,j) == LUm->getElement(i,j))){
          cout << "assertion failed at LU("<<i<<","<<j<<") == LUresult("<<i<<","<<j<<") " << endl;
          std::cout.unsetf(ios::floatfield); 
          std::cout.precision(20);
          diffVal =  LUresult->getElement(i,j).getValue() -  LUm->getElement(i,j).getValue();
          diffVal = fabs(diffVal);

          if(diffVal > marginErr )
            std::cout << "diffVal: "<< diffVal<<endl;     
        }
        assert( (!(LUresult->getElement(i,j).getValue() == 0)) &&(diffVal < marginErr) );
      }
    }
  }
  cout << "Verification Complete!"<<endl; 
  cout << "Test OK \t (Block LU)"<<endl; 

  delete LUm, delete Lm; delete Um;delete LUresult, delete Lresult; delete Uresult;  
}
