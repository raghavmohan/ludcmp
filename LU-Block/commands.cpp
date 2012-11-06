/* system example : DIR */
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <map>

#define backDir ".."
#define path "src"
#define DS '/'

//to turn debugging on, set to 1 or 2
#define DEBUG 0



using namespace std;
void constCharToChar(const char * a, char * b);
string getStdoutFromCommand(std::string cmd);
string intToString(int i);
int main(int argc, char *argv[])
{

  if (argc != 4){
    cout << "Usage: commands <size of matrix> <blocksize> <numthreads>" << endl;
    return 1;
  }

  int matrixSize = atoi(argv[1]);
  int numThreads = atoi(argv[3]);
  int numRuns = 3;
  int blockSize = atoi(argv[2]);
  int numBlocks = 1;


  int midRange = ((int) ceil(sqrt(matrixSize)));
  int startRange = (int) midRange/2;
  int endRange = (int) midRange*2;



  if(DEBUG > 0){  
    cout << "matrixSize: " << matrixSize <<  endl;
    cout << "numthreads: " << numThreads <<  endl;
    cout << "startRange: " << startRange <<  endl;
    cout << "midRange: " << midRange <<  endl;
    cout << "endRange: " << endRange <<  endl;
  }

  char cwd[MAXPATHLEN]; 
  getcwd(cwd, MAXPATHLEN);

  if (!system(NULL)) 
    exit (1);


  system("rm -f results.txt");
  system("touch results.txt");
  ostringstream stringStream;
  string cmd, output;

  //Declare and initialize vector
  //TODO:
  //there should be a better way to implement this (aka map)
  // values[numthread] = array
  //                [blocksize] = array
  //                        [runnumber] = timevalue
  //vector< vector < vector <int> > > values; 
  map<int, map<int, map<int, int> > > values;

  for(int i = 0; i <numThreads; ++i){
    map < int, map<int, int> > tempMap1; 
    values[i] =  tempMap1; 
    for(int j = 0; j < numBlocks; ++j){
      map< int , int> tempMap2; 
      values[i][j] =  tempMap2; 
      for(int k = 0; k < numRuns; ++k){
        values[i][j][k] =  0;
      }
    }
  }



  for(int i = 0; i <numThreads; ++i){
    for(int j = 0; j < numBlocks; ++j){
      for(int k = 0; k < numRuns; ++k){
        stringStream.str("");

        //threads can't be 0
        cmd = stringStream.str();
        if(DEBUG == 2)
          cout<<cmd<<endl;  

        char * writable = new char[cmd.size() + 1];
        std::copy(cmd.begin(), cmd.end(), writable);
        writable[cmd.size()] = '\0'; // don't forget the terminating 0
        putenv(writable);

        system("rm -f out.txt");
        stringStream.str("");

        stringStream<<"matrix_ss ";
        stringStream<< matrixSize <<" "<< blockSize<<" > out.txt";
        cmd = stringStream.str();
        if(DEBUG == 2)
          cout << cmd << endl;

        system(cmd.c_str());

        stringStream.str("");
        string temp1, temp2;

        temp1 = getStdoutFromCommand("cat out.txt | grep \"num_threads\" | awk -F'=' '{print $2}'");
        temp2 = getStdoutFromCommand("cat out.txt | grep \"Took\" | awk '{print $2}'");


        int a; 
        a = atoi(temp1.c_str());
        int b = atoi(temp2.c_str());

        values[i][j][k] = b;
        delete[] writable;
      }
    }
  } 





  if(DEBUG > 0){
    for(int i = 0; i <numThreads; ++i){
      for(int j = 0; j < numBlocks; ++j){
        for(int k = 0; k < numRuns; ++k){
          cout << "values["<<i<<"]["<< j <<"]["<<k<<"]: " << values[i][j][k] << endl;
        }
      }
    } 
  }

  map < int, map <int , double> > averages;
  for(int i =0; i < numThreads; ++i){
    map< int, double > tempMap;
    averages[i] = tempMap;
    for(int j=0; j < numBlocks ; ++j){
      averages[i][j] = 0.0;
    }
  }

  for(int i =0; i < numThreads; ++i){
    for(int j=0; j < numBlocks ; ++j){
      double sum = 0.0;
      for(int k = 0; k < numRuns; ++k){
        sum +=  values[i][j][k];
      }
      averages[i][j] = ((long double)sum / numRuns);
    }
  }

  for(int i = 0; i <numThreads; ++i){
    for(int j = 0; j < numBlocks; ++j){
      cout << "averages["<<i+1<<"]["<< j <<"]: " << averages[i][j] << endl;
    }
  }
  system("rm -f results.txt out.txt"); 
  return 1;

}

std::string getStdoutFromCommand(std::string cmd) {
  std::string data;
  FILE *stream;
  int MAX_BUFFER = 256;
  char buffer[MAX_BUFFER];
  cmd.append(" 2>&1");
  stream = popen(cmd.c_str(), "r");
  if (!stream){
    std::cout<<"popen error"<<endl;
    exit(1);
  }
  while (!feof(stream)){
    if (fgets(buffer, MAX_BUFFER, stream) != NULL){
      data.append(buffer);
    }
  }
  pclose(stream);
  return data;
}

std::string intToString(int i){
  std::string s;
  std::stringstream out;
  out << i;
  s = out.str();
  return s;
}

void constCharToChar(const char * a, char * b){
  while(*a != '\0'){
    *b = *a;
    b++;
    a++;
  }
}

