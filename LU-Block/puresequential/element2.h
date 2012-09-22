/*******************************************************************************
 * This file contains class definitions for 
 * Element
 * ElementOp (Element Operations)
 * 
 *******************************************************************************/

#ifndef SS_CPP_ONLY
//#include "parakram.hpp"
#endif

#include <climits>

//Need these for various statistics
int shelved_xact = 0;
int tokens_passed = 0;
int tokens_to_acquire = 0;
int total_delegates = 0;
int total_xact_executed = 0;
int total_xact_pushed = 0;
int total_xact_unshelved = 0;
int maxx_list_size = 0;
float avg_list_size = 0;
int counter_lock = 0;
int total_tokens_requested = 0;
int hit_under_miss = 0;



//simple routine that allocates a 2d array
template <class T> 
void allocate2dArray(T** arg1, int rows, int cols ){
  arg1 = new T*[rows];
  for(int i = 0; i < rows; ++i){
    arg1[i] = new T[cols];
  } 
}


//simple routine that deallocates a 2d array
template <class T> 
void delete2dArray(T** arg1, int cols ){
  for(int i = 0; i < cols; ++i){
    delete[] arg1[i];
  }
  delete[] arg1;
}

#include <iostream>
using namespace std;

/******************************************************************************/
/*
 * Simple wrapper class Element. Wrapper for Double, but can easily make 
 * into a wrapper for T
 */
/******************************************************************************/
class Element
#ifndef SS_CPP_ONLY
//: public prometheus::static_xact_t
#endif
{
 public:
  Element(){
    this->value = new double (0.);
  }
 
  //Element(double valueInit = 0.0) :
 Element(double valueInit) :
  value(new double(valueInit))
    {}

  // Copy constructor
  Element(const Element &el)
    {
      value = new double (*el.value);
    }
  
  
  Element& operator= (const Element &el)
    {
      // check for self-assignment by comparing the address of the
      // implicit object and the parameter
      if (this == &el)
        return *this;

      *value = *(el.value);
      return *this;
    }
  
  Element& operator= (const double val)
    {  
      *value = val;
      return *this;
    }

  
  //Define Operators. Can one call df_execute on an operator?
  Element& operator+=(const Element &elrhs) {
    *(this->value) += *(elrhs.value); 
    return *this;
  }

  Element& operator-=(const Element &elrhs) {
    *(this->value) -= *(elrhs.value); 
    return *this;
  }
  Element& operator*=(const Element &elrhs) {
    *(this->value) *= *(elrhs.value); 
    return *this;
  }
  

  Element& operator/=(const Element &elrhs) {
    if(*(elrhs.value) != 0){
      *(this->value) /= *(elrhs.value); 
    }
    return *this;
  }
  

  const Element operator+(const Element &other) const {
    Element result = *this;
    result += other;
    return result;
  }

  const Element operator-(const Element &other) const {
    Element result = *this;
    result -= other;
    return result;
  }

  const Element operator*(const Element &other) const {
    Element result = *this;
    result *= other;
    return result;
  }
  
  const Element operator/(const Element &other) const {
    Element result = *this;
    result /= other;
    return result;
  }
  
  
  bool operator==(const Element &el) const {
    return (  (float) ( *(this->value)) ==  (float)(*(el.value))); 
  }
  
  bool operator!=(const Element &el) const {
    return (*(this->value) != *(el.value)); 
  }
  
  bool operator>(const Element &el) const {
    return (*(this->value) > *(el.value)); 
  }

  bool operator<(const Element &el) const {
    return (*(this->value) < *(el.value)); 
  }
  
  operator double() const {return (*value);  }

  operator int() const {return (int) (*value);  }

  friend ostream &operator << ( ostream &out, const Element &el){
    out << *(el.value);
    return out;
  }
  
  //destructor
  ~Element()
    {
      delete value;
    }

  double getValue(){
    return *(this->value);
  }
 private:
  double* value;
};


/******************************************************************************/
/*
 * Simple class of ElementOp (Element Operations). Need this because if one
 * calls df_execute on an Element, then a thread needs to have an operation
 * to call df_execute on where it can execute the operation and get the
 * the result later
 */
/******************************************************************************/
class ElementOp
#ifndef SS_CPP_ONLY
//: public prometheus::static_xact_t
#endif
{
 public:
  typedef enum OpType{
    ADD,
    ADDEQ,
    SUB,
    MUL,
    DIV,
    SPEC1, // el1 += el2 * el3
    SET 
  }OpType;


 private:
  Element* el1;
  Element* el2;
  Element* el3;
  
  OpType opType;

 public:  
  ElementOp(OpType opType, Element *el1, Element *el2, Element *el3){
    this->el1 = el1; 
    this->el2 = el2;
    this->el3 = el3;

    this->opType = opType;
  }
  
 
  //routine that executes the Operation
  void execute(){
    // *(result) = 0.0;
    /*****************************************************/
    //spin loop for measuring performance differences
    //int a;
    // for(int i = 0; i < INT_MAX; ++i) ++a ;
    // for(int i = 0; i < 1000000000; ++i) ;
    /*****************************************************/
    
    switch(this->opType){
    case ADD:
      *(this->el1) = *(this->el2) + *(this->el3);
      break;
    case SUB:
      *(this->el1) = *(this->el2) - *(this->el3);
      break;
    case MUL:
      *(this->el1) = (*(this->el2)) * (*(this->el3));
      break;
    case DIV:
      *(this->el1) = (*(this->el2)) / (*(this->el3)) ;
      break;
      //   case SPEC1:
      //this->result = this->el1 + this->el2 * this->el3;
      // break;
      //case ADDEQ:
      //this->el1 += this->el2 + this->el3;
      //break;
      /*TODO
	case ADDEQ:
	*(this->result) += *
	break;
      */
    default:
      break;
    }
    
  }
  
  
  void getResult(Element * el){
    *(el) = *(this->el1);
  }
  
  /*
    Element * getResultPtr(){
    return  this->result;
    }
  
  
    Element * getElement1(){
    return  this->el1;
    }
    Element * getElement2(){
    return  this->el2;
    }  
  */


  Element getResultVal(){
    return *(this->el1);
  }

};


