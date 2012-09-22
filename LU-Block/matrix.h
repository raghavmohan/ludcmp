/*******************************************************************************
 * This file contains class definitions for class
 * Matrix - implemented as a 2d array of Element Objects
 * 
 *******************************************************************************/
#include "element.hpp"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>
#include <vector>
#include <set>

#ifndef SS_CPP_ONLY
#include "parakram.hpp"
#endif

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

//#define DEBUG (0)
#define DEBUG (2)

#define TINY 1e-10f
#define CALIGN 64 

class Matrix;

#ifndef SS_CPP_ONLY
typedef set<Matrix*> obj_set_t;
#endif




int (Matrix::*fcnptr)(void *) = NULL;


typedef struct MGetSetArgs{
  Matrix * MGet, *MSet;
  int rowG, columnG, rowS, columnS;
} MGetSetArgs;

typedef struct VerifyTestArgs{
  Matrix * LUm, *Ltest, *Utest;
} VerifyTestArgs;






typedef struct BlockLUArgs{
  Matrix * LUm;
  Matrix * Lm;
  Matrix* Um;
}BlockLUArgs;


//struct for inverse function args
typedef struct MInverseArgs{
  Matrix * passM;
  Matrix * resultM;
  int flag;
}MInverseArgs;


//struct for inverse function args
typedef struct MOpArgs{
  Matrix * result;
  Matrix * M1;
  Matrix * M2;
  ElementOp::OpType type;
  
}MOpArgs;


//struct for a general function
typedef struct MFcnArgs{
  Matrix * LUm;
  void * fcn;
  void * fcnArgs;
}MFcnArgs;


//variable declarations
extern int blockSize;



typedef set < MOpArgs * > MOpArgsSet;
MOpArgsSet * MOpArgsSetToFree ;

typedef set < BlockLUArgs * > MLUArgsSet;
MLUArgsSet * MLUArgsSetToFree ;

typedef set < MInverseArgs * > MInvArgsSet;
MInvArgsSet * MInvArgsSetToFree;

typedef set < MGetSetArgs * > MGetSetArgsSet;
MGetSetArgsSet *  MGetSetArgsSetToFree;

typedef set< Matrix * > MatrixSet;
MatrixSet * matricesToFree;

#ifndef SS_CPP_ONLY
typedef set< prometheus::obj_set_t * > ObjSet;
ObjSet * objsToFree;
#endif

class Matrix
#ifndef SS_CPP_ONLY
: public prometheus::static_xact_t
#endif
{
private:
  Element ** array;
  int numRows;
  int numCols;
  
  //for blockLU
  int tinyRowSize, tinyColSize;
  int bigRowSize, bigColSize;


public:
  //for block LU
  vector< vector <Matrix*> > Ls;
  vector< vector<Matrix*> > Us;
  vector< vector< Matrix *> > matrixBlocks;

  
  //constructor that constructs a Dense Square 2-D array of Matrix Elements
  //Matrix needs to be square for LU- thus a NxM matrix is made into a 
  //N x N matrix where N > M, the values are zero filled.
  //Matrix(int numRows=0, int numCols=0);
  Matrix(int numRows, int numCols);
 
  //copy constructor
  Matrix(const Matrix &M); 
  
  //destructor
  ~Matrix();
  
  
  //constructor that takes a 2d double array - makes main testing easier
  Matrix(int numRows, int numCols , double** initArray);
  

  //Matrix Operators
  Matrix& operator= (const Matrix &M);
  Matrix& operator+=(const Matrix &M);
  Matrix& operator-=(const Matrix &M);
  Matrix& operator*=(const Matrix &M);
  const Matrix operator+(const Matrix &other) const;
  const Matrix operator-(const Matrix &other) const;
  const Matrix operator*(const Matrix &other) const;
  
  /*
  void * operator new(size_t size);
  void operator delete(void * ptr);
  */

  //accessors
  int getNumRows()const;
  int getNumCols()const;
  void getElement(int row, int column, Element * element);
  Element * getElementPtr(int i, int j);
  Element& getElement(int row, int column)const ;
  
  //setters
  bool setElement(int row, int column, Element element);
  void setElementPtr(int row, int column, Element* element);
  
  void getSetElement(MGetSetArgs getSetArgs);
  
 
  

  //current LU Function. Only solves Positive definite square matrices without pivot.
  //most benchmark specs use this kind of matrix, and then perform Cholesky on it
  //for additional performance measures
  //void Doolittle(int d, Matrix* L, Matrix * U);
  void Doolittle(Matrix* L, Matrix * U);

  //block LU - blocked LU, blocks are chosen by cmd line args
  void blockLU();

  //vector dot product - can be made a static function
  Element vectorDot(Element el1[], Element el2[], int numEls ) const;
  
  // Routine that swaps a[i][*] with a[j][*]
  void swapRows(int i, int j);
  
  //transpose the current matrix
  void transposeMatrix();
 
  //functions that compute what they say but need for df_exec calls
  void multiply(const Matrix* M1, const Matrix* M2);
  void subtract(const Matrix* M1, const Matrix* M2);
  

  //partition the matrix for blockLU 
  void partition();
  //get the results from blockLU factorization
  void getBlockResults(Matrix * Lresult, Matrix * Uresult);
  //cleanup from blockLU
  void freeBlocks();
  
  //Debug functions
  void print();
  
  //print wolfram alpha style - makes it easy to debug wih wolfram
  void printWolfram();

}; //end Matrix class


//other static functions
void getSetElement(MGetSetArgs* getSetArgs);

#ifndef SS_CPP_ONLY
void execLU(prometheus::obj_set_t *wr_set, prometheus::obj_set_t *rd_set, void * arg); 
void execMInverse(prometheus::obj_set_t *wr_set, prometheus::obj_set_t *rd_set, void * arg);
void execMOp(prometheus::obj_set_t *wr_set, prometheus::obj_set_t *rd_set, void * arg); 
void execMGetSet(prometheus::obj_set_t *wr_set, prometheus::obj_set_t *rd_set, void * arg);
#endif

void getPassArray(Element ** original, Element ** arrayToPass, int expand,int expand2, int n);
Element computeDeterminant(Element ** a, int n, int flag);
void transposeMatrix(Matrix* M);
void transposeArray(Element ** array, int numRows, int numCols);
void getTriangularInverse(Matrix * originalM, Matrix * resultM, int flag);
void CoFactor(Element **a,int n,Element **b);
