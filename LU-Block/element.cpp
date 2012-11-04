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
  } 
  

  ~Matrix() { ~matrixArray(); }
  

  int numRows(){ return this->rows;}
  int numColumns(){ return this->columns;}


  bool addValue(int row, int column, const T& element){
    if(!checkIndex(row, column))
      return false;
    
    matrixArray[row][column] = element;
    return true;
    
  }

  T getValue(int row, int column){
    if(!checkIndex(row, column))  
      return NULL;
    return this->matrixArray[row][column];


  }
private:
  int rows ;  // Number of rows in Matrix
  int columns ;  // Number of columns in Matrix
  vector< vector <T > > matrixArray;
  
  bool checkIndex(int row, int column){
    if(row < 0 || column < 0)
      return false;

    if(row > this->rows || column > this->columns )
       return false;
     return true;
  }
  
};
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
