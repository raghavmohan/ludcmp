#include <algorithm> // for std::swap
#include <cstddef>
#include <cassert>
#include <iostream>
#include <vector> 

using namespace std;
template <class T>
class Matrix
{
public:

  Matrix(int rows, int columns){
    this->rows = rows;
    this->columns = columns;  
    
    //allocate memory
    this->matrixArray.resize(rows);
    for (int i = 0; i < rows; ++i)
      matrixArray[i].resize(columns);

    //TODO: more efficient memory management than previous impl
  } 
  
  //TODO: make this constructor? seems like a bit of work
  //Matrix(int rows, int columns, void * values){
  //}
 
  ~Matrix() { ~matrixArray(); }
  

  int numRows(){ return this->rows;}
  int numColumns(){ return this->columns;}


  bool addValue(int row, int column, const T& element){
    // if(row > this->rows || column > this->columns )
    //return false;
    if(!checkIndex(row, column))
      return false;
    
    matrixArray[row][column] = element;
    return true;
    
  }

  T getValue(int row, int column){
    // if(row > this->rows || column > this->columns )
    if(!checkIndex(row, column))  
      return NULL;
    return this->matrixArray[row][column];


  }

  // ~Matrix() { delete [][] matrixArray; }
  /*
  int push(const T&); 
  int pop(T&) ;  // pop an element off the stack
  int isEmpty()const { return top == -1 ; } 
  int isFull() const { return top == size - 1 ; } 
  */
private:
  int rows ;  // Number of rows in Matrix
  int columns ;  // Number of columns in Matrix
  //T** matrixArray ;//generic array used to implement the Matrix  
  //T* matrixArray ;//generic array used to implement the Matrix
  vector< vector <T > > matrixArray;
  
  bool checkIndex(int row, int column){
    if(row < 0 || column < 0)
      return false;

    if(row > this->rows || column > this->columns )
       return false;
     return true;
  }
  
};
/*
//constructor with the default size 10
template <class T>
Stack<T>::Stack(int s)
{
  size = s > 0 && s < 1000 ? s : 10 ;  
  top = -1 ;  // initialize stack
  stackPtr = new T[size] ; 
}
// push an element onto the Stack 
template <class T>
int Stack<T>::push(const T& item)
{
  if (!isFull())
    {
      stackPtr[++top] = item ;
      return 1 ;  // push successful
    }
  return 0 ;  // push unsuccessful
}

// pop an element off the Stack
template <class T> 
int Stack<T>::pop(T& popValue) 
{
  if (!isEmpty())
    {
      popValue = stackPtr[top--] ;
      return 1 ;  // pop successful
    }
  return 0 ;  // pop unsuccessful
}
*/
int main(int argc, char* argv[]){
  typedef Matrix<int> MatrixInt;
  typedef Matrix<double> MatrixDouble;
  

  MatrixDouble * md = new MatrixDouble(3,3);
  cout << "Rows: " << md->numRows();
  cout << "Cols: " << md->numColumns();
  
  if(md->addValue(0, 1, 5.0))
    cout << "Value at 0, 1: "<< md->getValue(0, 1)<<endl;
  else
    cout << "Erorr adding Value";
  
  cout << "Value at 0, 0: "<< md->getValue(0, 0)<<endl;
  return 0;
}
