#include "matrix2.hpp"
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string>


#define NUM_TIMES (1)
//float marginErr = 1e-2f;
float marginErr = 1e-10f;





//global variables

//int numThreads;
int blockSize;

int matrixSize;
Matrix * L, *U, *LU;


void test1();

//method declarations
void test1();
//void verifyTest1();


//generate random num
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


using namespace std;
double fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

int intRand(int fMin, int fMax)
{
  int f = (int)rand() % fMax;
  return f;
}



int main(int argc, char** argv,char ** envp){
  //total time
  struct timeval startTotal, endTotal; 
  gettimeofday(&startTotal, NULL);

  //set in case not defined
  blockSize=2;
  matrixSize =8;
  
  if(argc > 1){
    matrixSize = atoi(argv[1]);
  }
  
  //get num threads via cpp ...
  //char * threadString;
  //threadString = getenv ("PROMETHEUS_THREADS");
  //if (threadString != NULL)
  //  numThreads = atoi(threadString);

  if(argc > 2){
    //----------------------------------------//
    blockSize=atoi(argv[2]);
    //----------------------------------------//
  }

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
//TODO: Add Timing
void test1(){
  
  //declarations
  Matrix * Lm = new Matrix(matrixSize, matrixSize);
  Matrix * Um = new Matrix(matrixSize, matrixSize);
  Matrix * LUm =  new Matrix(matrixSize, matrixSize);
  
  Matrix * Lresult = new Matrix(matrixSize, matrixSize);
  Matrix * Uresult = new Matrix(matrixSize, matrixSize);
  Matrix * LUresult =  new Matrix(matrixSize, matrixSize);
  

  /*
  (*LUm) = (*LU) ;
  (*Lm) = (*L) ;
  (*Um) = (*U) ;
  */

  long double diffVal = 0;
   
  for(int i= 0 ; i < matrixSize; ++i){
    for(int j= 0 ; j < matrixSize; ++j){
      if(i == j){
	Lm->setElement(i, j,Element( (double) 1.0 ));
	Um->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
      }
      else if(i < j){
	Lm->setElement(i, j,Element( (double) 0. ));
	Um->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
      }
      else{
	Lm->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
	Um->setElement(i, j,Element( (double) 0. ));
      }
      
      LUm->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
      Lm->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
      
      //Um->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
    }
  }
  //(*LUm) = (*Lm) * (*Um);


  /*
  for(int i= 0 ; i < matrixSize; ++i){
    for(int j= 0 ; j < matrixSize; ++j){
      //Lm->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
      //Um->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
      LUm->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
      Lm->setElement(i, j,Element( (double) randGen<int>(-20, 20)));
      //Lm->setElement(i, j,Element( (double) Larr[i][j]));
      //Um->setElement(i, j,Element( (double) Uarr[i][j]));
    }
  }
  */	
  
  //(*LUm) = (*Lm);
 

  
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
  cout << "\nStarting verification ... "<<endl; 
  
  (*LUresult) = (*Lresult) * (*Uresult);
  
  LUm->freeBlocks();
  
  /*
 
  cout << "LU"<<endl;
  LUm->print();
  cout << "LUresult"<<endl;
  LUresult->print();
  
  cout << "L"<<endl;
  Lm->print();
  cout << "U"<<endl;
  Um->print();
  cout << "Lresult"<<endl;
  Lresult->print();
  cout << "Uresult"<<endl;
  Uresult->print();
  cout << "origLU"<<endl;
  LUm->print();
  
  cout << "LUresult"<<endl;
  LUresult->print();
  */
  
   for(int i= 0 ; i < matrixSize ; ++i){
    for(int j= 0 ; j < matrixSize; ++j){
      if(!(LUresult->getElement(i,j) == LUm->getElement(i,j))){
	cout << "assertion failed at LU("<<i<<","<<j<<") == LUresult("<<i<<","<<j<<") " << endl;
	   std::cout.unsetf(ios::floatfield); 
	   std::cout.precision(20);
	   diffVal =  LUresult->getElement(i,j).getValue() -  LUm->getElement(i,j).getValue();
	   diffVal = fabs(diffVal);
	   //std::cout << "\t"<< (float)this->array[i][j].getValue();
	   //std::cout << "\t"<<  fabs(diffVal); 
	   //std::cout << "diffVal: "<< diffVal;     
	   
	   if(diffVal > marginErr )
	     std::cout << "diffVal: "<< diffVal<<endl;     
	   
	   
	
	
	   //std::cout << "LU("<<i<<","<<j<<"): "<< LUm->getElement(i,j).getValue() << endl;
	   //std::cout << "LUresult("<<i<<","<<j<<"): "<< LUresult->getElement(i,j).getValue() << endl;
      }
      //assert( (!(LUresult->getElement(i,j).getValue() == 0)) &&(LUresult->getElement(i,j) == LUm->getElement(i,j)) );
      assert( (!(LUresult->getElement(i,j).getValue() == 0)) &&(diffVal < marginErr) );
    }
  }
 
  
   
   //cout << "\nTook "<< mtime<< " milliseconds."<<endl; 
   cout << "Verification Complete!"<<endl; 
   cout << "Test OK \t (Block LU)"<<endl; 
  
  //LUm->print();
  delete LUm, delete Lm; delete Um;delete LUresult, delete Lresult; delete Uresult;
  
  
}
