/*******************************************************************************
 * This file contains class definitions for class
 * Matrix - implemented as a 2d array of Element Objects
 * 
 *******************************************************************************/
#include "matrix2.h"
#include <unistd.h>

//-------------------------Start Matrix Class Functions-----------------------//


//constructor that constructs a Dense Square 2-D array of Matrix Elements
//Matrix needs to be square for LU- thus a NxM matrix is made into a 
//N x N matrix where N > M, the values are zero filled.
Matrix::Matrix(int numRows=0, int numCols=0)
{
  this->numRows =0;this->numCols=0;
  if(numRows < 0  || numCols < 0)
    return;
  //check if rows and cols are equal, if not zero fill
  if(numRows != numCols){
    this->numRows = max(numRows, numCols);
    this->numCols = this->numRows;
	
  }
  else{
    this->numRows = numRows; 
    this->numCols= numCols;
  }

  this->array = new Element*[this->numRows];
  for(int i = 0; i <this->numCols; ++i)
    this->array[i] = new Element[this->numCols];
}

//copy constructor
Matrix::Matrix(const Matrix &M){
    
  this->numRows = M.numRows; 
  this->numCols= M.numCols;
  this->array = new Element*[numRows];
  for(int i = 0; i < numCols; ++i)
    this->array[i] = new Element[numCols];
      

  if( (M.array != NULL)){
    for(int i = 0; i < numRows; ++i){
      for(int j = 0; j < numCols; ++j)
	this->setElement(i, j, M.array[i][j]);
    }
  }
}
  
//destructor
Matrix::~Matrix(){
  //need to delete the 2d array
  for(int i = 0; i < numRows; ++i){
    delete[] this->array[i];
  }
  delete[] this->array;
 
}
  
//constructor that takes a 2d double array - makes main testing easier
Matrix::Matrix(int numRows, int numCols , double** initArray)
{
  //TODO: Check if numRows == numCols, if not zero Fill
  if(numRows != numCols){
    this->numRows = max(numRows, numCols);
    this->numCols = this->numRows;
	
  }
  else{
    this->numRows = numRows; 
    this->numCols= numCols;
  }
      
  this->array = new Element*[this->numRows];
  for(int i = 0; i < this->numCols; ++i)
    this->array[i] = new Element[this->numCols];
      

  if( (initArray != NULL)){
    for(int i = 0; i < this->numRows; ++i){
      for(int j = 0; j < this->numCols; ++j)
	this->setElement(i, j, (initArray[i][j]));
    }
  }
}

//Matrix Operators
Matrix& Matrix::operator= (const Matrix &M)
{
  if (this == &M)
    return *this;
   
  if(numRows != M.getNumRows() || numCols != M.getNumCols())
    return *this;

      
  for(int i = 0; i < M.getNumRows(); ++i){
    for(int j = 0; j < M.getNumCols(); ++j)
      array[i][j] =(M.getElement(i,j));	  
  }
  return *this;
}

Matrix& Matrix::operator+=(const Matrix &M) {
  for(int i = 0; i < M.getNumRows(); ++i){
    for(int j = 0; j < M.getNumCols(); ++j){
      this->array[i][j] += M.getElement(i,j);
    }
  }
  return *this;
}

Matrix& Matrix::operator-=(const Matrix &M) {
  for(int i = 0; i < M.getNumRows(); ++i){
    for(int j = 0; j < M.getNumCols(); ++j){
      this->array[i][j] -= M.getElement(i,j);
    }
  }
  return *this;
}

Matrix& Matrix::operator*=(const Matrix &M) {
  Matrix* newMatrix = new Matrix(numRows, numCols);
  *newMatrix = *this;
  for(int i = 0; i < M.getNumRows(); ++i){
    for(int j = 0; j < M.getNumCols(); ++j){
      //Element * sum  = new Element(0.0);
      Element sum (0.);
      for(int k = 0; k < numCols; ++k){
	sum += newMatrix->array[i][k] * M.array[k][j];
	//*sum += newMatrix->array[i][k] * M.array[k][j];
      }
      this->array[i][j] = sum;
      //this->array[i][j] = *sum;
      // delete sum;
    }
  }

  delete newMatrix;
  return *this;
}

const Matrix Matrix::operator+(const Matrix &other) const {
  Matrix result = *this;
  result += other;
  return result;
}

const Matrix Matrix::operator-(const Matrix &other) const {
  Matrix result = *this;
  result -= other;
  return result;
}

const Matrix Matrix::operator*(const Matrix &other) const {
  Matrix result = *this;
  result *= other;
  return result;
}

//Cache align operators. 
//TODO:
//      Hardcoded CALIGN for now, figure out a way to get
//      this information pragmatically
/*
void * Matrix::operator new(size_t size){
  void *mem = malloc(size+CALIGN+sizeof(void*));
  void **ptr = (void**)((long)(mem+CALIGN+sizeof(void*)) & ~(CALIGN-1));
  ptr[-1] = mem;
  return ptr;
}

void  Matrix::operator delete(void *ptr){
  free(((void**)ptr)[-1]);
}
*/

  
Element Matrix::vectorDot(Element el1[], Element el2[], int numEls ) const{
  Element el(0.0);
  for(int i =0; i < numEls; i++){
    el += el1[i] * el2[i];
  }
  return el;
}
  
Element * Matrix::getElementPtr(int i, int j){
  if(i >= 0 && i < numRows && j >= 0 && j < numCols)
    return &(array[i][j]);
  else
    return NULL;
}


bool Matrix::setElement(int row, int column, Element element){
  if(row >= 0 && row < numRows && column >= 0 && column < numCols){
    array[row][column] = element;
    return true;  
  }
  else
    return false;
}

void Matrix::setElementPtr(int row, int column, Element* element){
  if(row >= 0 && row < numRows && column >= 0 && column < numCols){
    array[row][column] = *element;
  }
}


Element& Matrix::getElement(int row, int column)const {
  //if(row >= 0 && row < numRows && column >= 0 && column < numCols)
  return array[row][column]; 
  //else
  //return NULL;
  //return Element(0.);
}

void Matrix::getElement(int row, int column, Element * element){
  if( (row < numRows) && (row >= 0) && (column < numCols) && (column >= 0) )
    (*element) =  array[row][column]; 
}


void Matrix::partition(){
  MOpArgsSetToFree = new MOpArgsSet ;
  MLUArgsSetToFree = new MLUArgsSet;
  MInvArgsSetToFree = new MInvArgsSet;
  MGetSetArgsSetToFree = new  MGetSetArgsSet;
  matricesToFree = new MatrixSet;
  

  //Since always a square matrix...
  tinyRowSize = ( (blockSize < this->numRows) ? blockSize : this->numRows );
  tinyColSize = ( (blockSize < this->numRows) ? blockSize : this->numRows );
  
  bigRowSize = (int)ceil(( (double) this->numRows)/(tinyRowSize)) ;
  bigColSize = (int)ceil(( (double) this->numRows)/(tinyColSize));
  
  if(DEBUG == 2){
    std::cout << "this->numRows "<<this->numRows<< endl;
    std::cout << "this->numCols "<<this->numCols<< endl;
    std::cout << "tinyRowSize "<< tinyRowSize<< endl;
    std::cout << "tinyColSize "<< tinyColSize<< endl;
    std::cout << "bigRowSize "<< bigRowSize<< endl;
    std::cout << "bigColSize "<< bigColSize<< endl;
    std::cout << "Block Size "<< blockSize<< endl;
  }



  //allocate space for block matrices
  matrixBlocks.resize(bigRowSize);
  Ls.resize(bigRowSize);
  Us.resize(bigRowSize);
  for(int i = 0; i < bigRowSize; ++i){
    Ls[i].resize(bigColSize);
    Us[i].resize(bigColSize);
    matrixBlocks[i].resize(bigColSize);
  }

  //initialize each element in the block matrix
  for(int i = 0; i < bigRowSize; ++i){
    for(int j = 0; j < bigColSize; ++j){
      Ls[i][j] = new Matrix(tinyRowSize, tinyColSize);
      Us[i][j] = new Matrix(tinyRowSize, tinyColSize);
      matrixBlocks[i][j] = new Matrix(tinyRowSize, tinyColSize);
 
      for(int entI = 0; entI < tinyRowSize; ++entI){
	for(int entJ = 0; entJ < tinyColSize; ++entJ){
	  int arrayI, arrayJ;
	  arrayI =  ((i* tinyRowSize) + entI) % (this->numRows);
	  arrayJ = ((j* tinyColSize) + entJ) % (this->numCols) ;
	  (matrixBlocks[i][j])->setElement(entI, entJ, array[arrayI][arrayJ]);
	}
      }
    }
  }
}


void Matrix::blockLU(){
  //these two matrices will correspond to inv(L00) and inv(U00)
  
  Matrix * L00inv = new Matrix(tinyRowSize, tinyColSize);
  Matrix * U00inv = new Matrix(tinyRowSize, tinyColSize);

   
  matricesToFree->insert(U00inv);
  matricesToFree->insert(L00inv);

  //start computation after setup
  (matrixBlocks[0][0])->Doolittle(Ls[0][0],Us[0][0]);

  getTriangularInverse(Ls[0][0], L00inv, 1);
  getTriangularInverse(Us[0][0], U00inv, 2); 
   
  for(int i = 1; i < bigRowSize; ++i){  
    //compute L[i][0] && U[0][i]   
    Ls[i][0]->multiply(matrixBlocks[i][0], U00inv);
    Us[0][i]->multiply(L00inv,matrixBlocks[0][i]);
  }
  

  for(int i = 1; i < bigRowSize; ++i){
    for(int j = 1; j < bigColSize; ++j){
     	
      Matrix * temp1 = new Matrix(tinyRowSize, tinyColSize);
      matricesToFree->insert(temp1);
      
      Matrix * temp2 = new Matrix(tinyRowSize, tinyColSize);
      matricesToFree->insert(temp2);
      
      //write explanation for this later
      Matrix * junk = new Matrix(tinyRowSize, tinyColSize);
      matricesToFree->insert(junk);

      (*temp1) =  *(matrixBlocks[i][j]);

      //go through all cols of Ls, rows of Us
      for(int k = 0; k < bigRowSize; ++k){
	//end case 
	if(( (i == bigRowSize && j == bigColSize) &&(k != (bigRowSize - 1) ) ) || (!(i == k && j == k)) ) {
	//if(( (i == bigRowSize && j == bigColSize) &&(k != (bigRowSize - 1) ) )) {
	  temp2->multiply(Ls[i][k], Us[k][j]);
	  temp1->subtract(temp1, temp2);
	  //(*temp2) = ( *(Ls[i][k])) * (*(Us[k][j]));
	  //(*temp1) -= (*temp2) ; 
	}
      }

      if(i == j){
	temp1->Doolittle( Ls[i][j], Us[i][j]);
      }

      
      else if(i > j){
	getTriangularInverse(Us[j][j], junk, 2);
	Ls[i][j]->multiply(temp1, junk);
	// *(Ls[i][j]) = (*temp1) *(*(junk));
      }


      else{
	getTriangularInverse(Ls[i][i], junk, 1);
	(Us[i][j])->multiply(junk, temp1);
      }

    }
  }
}

void Matrix::getBlockResults(Matrix *Lresult, Matrix * Uresult){
  //Since always a square matrix...
  for(int i = 0; i < bigRowSize; ++i){
    for(int j = 0; j < bigColSize; ++j){

      for(int entI = 0; entI < tinyRowSize; ++entI){
	for(int entJ = 0; entJ < tinyColSize; ++entJ){
	  int arrayI, arrayJ;
	  arrayI =  ((i*tinyRowSize) + entI) % (this->numRows);
	  arrayJ = ((j*tinyColSize) + entJ) % (this->numCols) ;

	  Lresult->setElement(arrayI, arrayJ, Ls[i][j]->getElement(entI,entJ));
	  Uresult->setElement(arrayI, arrayJ, Us[i][j]->getElement(entI,entJ));	  
	}
      }
    }
  }
}

int Matrix::getNumRows()const {
  return this->numRows;
}
int Matrix::getNumCols()const {
  return this->numCols;
}


void Matrix::freeBlocks(){
  for(int i = 0; i <(int) Ls.size() ; ++i){
    for(int j = 0; j <(int) Ls[i].size(); ++j){
      delete Ls[i][j];
      delete Us[i][j];
      delete matrixBlocks[i][j];
    }
  }
  //cleanup
  set<Matrix *>::iterator matricesIt;
  for(matricesIt = matricesToFree->begin(); matricesIt != matricesToFree->end(); ++matricesIt){
    delete (*matricesIt);
  }
  
  MOpArgsSet::iterator mOpIt;
  for(mOpIt = MOpArgsSetToFree->begin(); mOpIt != MOpArgsSetToFree->end(); ++mOpIt){
    delete (*mOpIt);
  }

  MLUArgsSet::iterator mLUIt;
  for(mLUIt = MLUArgsSetToFree->begin(); mLUIt != MLUArgsSetToFree->end(); ++mLUIt){
    delete (*mLUIt);
  }
  MInvArgsSet::iterator mInvIt;
  for(mInvIt = MInvArgsSetToFree->begin(); mInvIt != MInvArgsSetToFree->end(); ++mInvIt){
    delete (*mInvIt);
  }

  MGetSetArgsSet::iterator mGSIt;
  for(mGSIt = MGetSetArgsSetToFree->begin(); mGSIt != MGetSetArgsSetToFree->end(); ++mGSIt){
    delete (*mGSIt);
  }
  
  delete MOpArgsSetToFree;delete MInvArgsSetToFree; delete MLUArgsSetToFree;delete MGetSetArgsSetToFree;
  delete matricesToFree;
}


//current LU Function. Only solves Positive definite square matrices without pivot.
//most benchmark specs use this kind of matrix, and then perform Cholesky on it
//for additional performance measures
void Matrix::Doolittle(Matrix* L, Matrix * U)
{

  int d = this->numRows;
  //TODO:
  // - Implement Cholesky, Crout
  Matrix * D = new Matrix(d,d);          
  set<Element *> elementSet;
  set<ElementOp *> elementOpSet;
    
  for(int k=0;k<d;++k){
    for(int j=k;j<d;++j){
      Element* sum = new Element(0.);
      elementSet.insert(sum);
	
      for(int p=0;p<k;++p){

	Element * temp1 = new Element(0.);
	ElementOp *elOpMul = new ElementOp(ElementOp::MUL,temp1, D->getElementPtr(k,p), D->getElementPtr(p, j));
	ElementOp * sumOp = new ElementOp(ElementOp::ADD,sum, sum, temp1);
	  
	elementSet.insert(temp1);elementOpSet.insert(elOpMul);elementOpSet.insert(sumOp);
	  
	elOpMul->execute();
	sumOp->execute();
      }


      ElementOp *elOpSub = new ElementOp(ElementOp::SUB,D->getElementPtr(k, j), &(array[k][j]) , sum);
	
      elementOpSet.insert(elOpSub);
	
      elOpSub->execute();
	
    }


    for(int i=k+1;i<d;++i){
      Element * sum2 = new Element(0.0);
	
      elementSet.insert(sum2);
	
      for(int p=0;p<k;++p){
	Element * temp2 = new Element(0.);
	ElementOp *elOpMul2 = new ElementOp(ElementOp::MUL,temp2, D->getElementPtr(i,p), D->getElementPtr(p, k));
	ElementOp * sumOp2 = new ElementOp(ElementOp::ADD,sum2, sum2, temp2);
	  
	elementSet.insert(temp2);elementOpSet.insert(elOpMul2);elementOpSet.insert(sumOp2);
	  
	elOpMul2->execute();
	sumOp2->execute();
      }


      Element * r1 = new Element (0.);
      ElementOp *elOpSub2 = new ElementOp(ElementOp::SUB,r1, &(array[i][k]), sum2);
      ElementOp *elOpDiv1 = new ElementOp(ElementOp::DIV,D->getElementPtr(i,k), r1, D->getElementPtr(k,k));
	
      elementSet.insert(r1);elementOpSet.insert(elOpSub2);elementOpSet.insert(elOpDiv1);
	
      elOpSub2->execute();
      elOpDiv1->execute();
	
    }
  }

    
  for(int i = 0; i < numRows; ++i){
    for(int j =0; j < numCols; ++j){	
      if(j<i){
	L->setElement(i, j, D->getElement(i,j));
	U->setElement(i, j, Element(0.0));
      }

      else if(i==j){
	L->setElement(i, j, Element(1.0));
	U->setElement(i, j, D->getElement(i,j));
      }
      else{
	L->setElement(i, j, Element(0.0));
	U->setElement(i, j, D->getElement(i,j));
      }
    }
  }    

  delete D;
    
  set<Element * >::iterator elIt;
  set<ElementOp* >::iterator elOpIt;
    
  for(elIt = elementSet.begin(); elIt != elementSet.end(); ++elIt)
    delete *elIt;
    
  for(elOpIt = elementOpSet.begin(); elOpIt != elementOpSet.end(); ++elOpIt)
    delete *elOpIt;
    
}

// Routine that swaps a[i][*] with a[j][*]
void Matrix::swapRows(int i, int j){
  if(i < 0 || i > numRows ||j < 0 || j > numRows  )
    return;

  Element * temp = array[i];
  array[i]= array[j];
  array[j] = temp;
}

//transpose the current matrix
void Matrix::transposeMatrix(){
  for(int i = 0; i < numRows ; ++i){
    for(int j = 0; j < i; ++j){
      if(i != j){
	Element temp(0.);
	temp = array[i][j];
	array[i][j] = array[j][i];
	array[j][i] = temp;
      }
    }
  }
}


void Matrix::multiply(const Matrix* M1, const Matrix* M2){
  (*this) = (*M1) *(*M2);
}

void Matrix::subtract(const Matrix* M1, const Matrix* M2){
  (*this) = (*M1) - (*M2);
}

 


//Debug functions
void Matrix::print(){
  for(int i = 0; i < numRows; i++){
    for(int j = 0; j < numCols; j++){
      long double diffVal = 0;
      std::cout.unsetf(ios::floatfield); 
      std::cout.precision(20);
      diffVal = (float)this->array[i][j].getValue() - 0;
      //std::cout << "\t"<< (float)this->array[i][j].getValue();
      //std::cout << "\t"<<  fabs(diffVal); 
      std::cout << "\t"<< diffVal;     
    }
    std::cout << std::endl;
  }
}
  
//print wolfram alpha style - makes it easy to debug with wolfram
void Matrix::printWolfram(){
  printf("\n M : {");
  for(int i = 0; i < numRows; ++i){
    printf("{");
    for(int j =0; j < numCols; ++j){
      printf("%g", array[i][j].getValue());      
      if(j != numCols-1)
	printf(",");
    }
    
    printf("}");
    if(i != numRows-1){
      //if(flag)
      //  printf(",");
      //else
      printf(",");
    }
  }
  printf("}\n");
}

//-------------------------End Matrix Class Functions-------------------------//

using namespace std;
void getSetElement(MGetSetArgs* getSetArgs){
  //if( (getSetArgs->rowS < numRows) && (getSetArgs->rowS >= 0) && (getSetArgs->columnS < numCols) && (getSetArgs->columnS >= 0) ){
    
  getSetArgs->MSet->setElement(getSetArgs->rowS, getSetArgs->columnS, getSetArgs->MGet->getElement(getSetArgs->rowG, getSetArgs->columnG));
  //array[getSetArgs->rowS][getSetArgs->columnS] = getSetArgs->MGet->getElement(getSetArgs->rowG, getSetArgs->columnG);
  //}
}




//square matrix triangular inverse function 
//NOTE: This uses more memory than it technically needs, but it makes it easier
//       because we do not have to add it to the write set for df_execute
//         flag => 1 lower, 2=> upper
void getTriangularInverse(Matrix * passM, Matrix * resultM, int flag){
  
  //base cases
  if(passM->getNumRows() == 1 && passM->getNumCols() == 1 ){
    long double passVal =passM->getElement(0,0).getValue() ;
    if( passVal == 0)
      return;
    else{
      resultM->setElement(0, 0, Element((long double) (1/passVal)) );
    } 
  }


  //check for flag - only want to invert triangular matrices 
  if(flag != 1 && flag != 2) return;
  
  int n = passM->getNumRows();

  Matrix * originalM = new Matrix(n, n); 
  (*originalM) = (*passM);

  


  //upper triangular, so transpose it
  if(flag==2){
    transposeMatrix(originalM);
  }
  

  //Invert the diagonal elements of the upper triangular matrix U.
  for(int i = 0; i < n; ++i){
    if(originalM->getElement(i,i).getValue() == 0){
      delete originalM;
      return;
    }
    else{
      //shouldn't need this
      for(int j = 0; j < n; ++j){
	Element temp(0.);
	resultM->setElement(i,j,temp);
      }
      Element temp2(0.);
      temp2.setValue( (1.0)/ (originalM->getElement(i,i).getValue()));
      resultM->setElement(i,i, temp2);
    }
  }
 

  for(int j = 0; j < n; ++j){
    for(int i = j+1; i < n; ++i){
      
      Element temp(0.);
      for(int k = j; k < i; ++k){
	temp += originalM->getElement(i,k) * resultM->getElement(k,j);
      }
      
      Element temp2 (0.);
      temp2.setValue( ((-1.0 / originalM->getElement(i,i).getValue())) * temp.getValue());
      resultM->setElement(i, j, temp2);  
    }
  }

  //transpose the result back if upper triangular
  if(flag == 2){
    transposeMatrix(resultM);    
  }
  delete originalM;

}



void getPassArray(Element ** original, Element ** arrayToPass, int expand,int expand2, int n){
  for(int iOrig = 0,iPass = 0; iOrig < n; ++iOrig){
    if( (iOrig == expand)){
      continue;
    }
    for(int jOrig = 0,jPass = 0; jOrig < n; ++jOrig){
      if( (jOrig== expand2)){
	continue;
      }
      
      arrayToPass[iPass][jPass] = original[iOrig][jOrig];
      ++jPass;
    }
    ++iPass;
  }
}


//TODO: Add a flag for lower, upper or regular matrix
//  Find the cofactor matrix of a square matrix
//flag : 0 => regular, 1 => lower, 2 => upper
void CoFactor(Element **a,int n,Element **b)
{
  //regular
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      // if(a[i][j].getValue() != 0){
      Element ** arrayToPass;

      arrayToPass = new Element*[n-1];
      for(int i1 = 0; i1 < n-1; ++i1)
	arrayToPass[i1] = new Element[n-1];
      
      for(int i2 = 0; i2 < n-1; ++i2){
	for(int j2 = 0; j2 < n-1; ++j2)
	  arrayToPass[i2][j2] =  Element(0.);
      }
      
      
      getPassArray(a, arrayToPass, i, j, n);
    
      Element temp(0.);
      cout << "a["<<i<<"]["<<j<<"]: "<< a[i][j]<<endl; 
      temp = computeDeterminant(  arrayToPass, n-1, 0) ;
      
      b[i][j] =temp ;
      if( (i+j) % 2 != 0)
	b[i][j] = b[i][j] * Element(-1.0);
      
      delete2dArray<Element>(arrayToPass, n-1);
    }
  }
  transposeArray(b, n , n);
}

void transposeMatrix(Matrix * M){
  for(int i = 0; i < M->getNumRows() ; ++i){
    for(int j = 0; j < i; ++j){
      if(i!= j){
	Element temp(0.);
	temp = M->getElement(i, j);
	M->setElement(i, j, M->getElement(j, i));
	M->setElement(j, i, temp);
      }
    }
  }
}


void transposeArray(Element ** array, int numRows, int numCols){
  for(int i = 0; i < numRows; ++i){
    for(int j = 0; j < i; ++j){
      
      if(i!= j){
	Element temp(0.);
	temp = array[i][j];
	array[i][j] = array[j][i];
	array[j][i] = temp;
      }
    }
  }
}


//TODO: Auto detect triangular matrices
//Function that recursively computes the determinant
//a => 2d Element array, n=> size of array, flag => 0 regular, 1 triangular 
Element computeDeterminant(Element ** a, int n, int flag){
  //base case
  if( n == 2){
    return ( (a[0][0] * a[1][1]) - (a[1][0] * a[0][1]) );   
  }
  
  Element result = Element(0.);
  
  //base case for triangular
  if(flag){
    result = Element(1.0);
    for(int i = 0;i < n; ++i){
      result *= a[i][i];
    }
    return result;
  }

  //TODO this will change depending upon the flag
  int expand = 0;

  //'n' is numCols, could have written this the other way  
  for(int i = 0; i < n; ++i){
    if(a[expand][i].getValue() != 0){

      Element ** arrayToPass;
      arrayToPass = new Element*[n-1];
      for(int i1 = 0; i1 < n-1; ++i1)
	arrayToPass[i1] = new Element[n-1];
    
      for(int i2 = 0; i2 < n-1; ++i2){
	for(int j2 = 0; j2 < n-1; ++j2)
	  arrayToPass[i2][j2] =  Element(0.);
      }
      getPassArray(a, arrayToPass, expand,i, n);
    
      if(expand + i % 2 == 0 )
	result += a[expand][i] * computeDeterminant(arrayToPass, n-1, 0);
      else 
	result -= a[expand][i] * computeDeterminant(arrayToPass, n-1, 0);
  
      delete2dArray<Element>(arrayToPass, n-1);
    }
  }
    
  return result;
}

