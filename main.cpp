/* 
 * File:   main.cpp
 * Author: Jonathan
 *
 * Created on 08 April 2014, 01:45
 * 
 * An unsigned long int should be big enough to store any lattice within interest
 * i.e q = 4 n = 4.
 */

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <vector>
#include <string>
#include <time.h>
#include "groupCounter.h"
#include <fstream>
#include <list>
#include "matConvert.h"
#include <ctime>
#include <armadillo>
//DEBUG FLAG, verbose output, -D command line.
bool debug = false;


bool checkSymmetry(matConvert mat);
void writeMatGroup(matConvert &mat, unsigned long int group);
void buildFileNamePrefix();
bool parityCheck(std::vector<int> qList, matConvert mat, matConvert comp, matConvert ref, int cur);
bool parWrap(std::vector<int> qList, matConvert tmat, matConvert comp, matConvert ref);
std::vector<int> calcCat(groupCounter &paraCount, matConvert &pLat, int* catHash[]);
int boundary(int tmp, int Ntemp);
void outputGod(std::vector<std::vector<std::vector< int > > > &God);
void newtonRhapsonSolver(std::vector<std::vector<std::vector< int > > > &God);
using namespace std;
using namespace arma;
/*
 * 
 */
//Global Variables, sue me.
bool Z3 = true, Z2 = true, ROT = true;
int N = 0, Q = 0, B = 2;
string filePrefix;
double tol = 1e-2;
std::vector<matConvert > symmetryList;

int main(int argc, char** argv) {
    char *outputFile;
    clock_t t;
    t = clock();
    srand(time(0));
    //Command line interface.
    for (int i = 0; i < argc; i++) {
        if (std::strcmp(argv[i], "-o") == 0) {
            outputFile = argv[i + 1];
            filePrefix.append(outputFile);
            i += 1;
        } else if (std::strcmp(argv[i], "-N") == 0) {
            N = atoi(argv[i + 1]);
            i += 1;
        } else if (std::strcmp(argv[i], "-Q") == 0) {
            Q = atoi(argv[i + 1]);
            i += 1;
        } else if (std::strcmp(argv[i], "-Z3") == 0) {
            Z3 = false;
        } else if (std::strcmp(argv[i], "-Z2") == 0) {
            Z2 = false;
        } else if (std::strcmp(argv[i], "-R") == 0) {
            ROT = false;
        } else if (std::strcmp(argv[i], "-D") == 0) {
            debug = true;
        } else if (std::strcmp(argv[i], "-B") == 0) {
            B = atoi(argv[i + 1]);
            i += 1;
        } else if (std::strcmp(argv[i], "-t") == 0) {
	  tol = atof(argv[i + 1]);
	  i += 1;
	}
    }

    buildFileNamePrefix();
    cout << filePrefix << endl;
    //Output Command line values.
    if (debug) {
        cout << "N = " << N << endl;
        cout << "Q = " << Q << endl;
        cout << "Z2 = " << Z2 << endl;
        cout << "Z3 = " << Z3 << endl;
        cout << "Rot = " << ROT << endl;
        //cout << "Output File = " << outputFile << endl;
        // cout << "WARNING: Debug mode ENABLED, verbose output" << endl;
    }

    //
    groupCounter paraCount = groupCounter();
    matConvert mat = matConvert(1, Q, N / B);
    matConvert mat2 = matConvert(1, Q, N / B);
    int catHash[(int) (pow(Q, N * N / (B * B)))][2];

    bool found;
    //Find all N/b lattices 
    for (unsigned long int i = 1; i <= pow(Q, N * N / pow(B, 2.0)); i++) {
        if (i % 1000 == 0) {
            cout << 100 * i / pow(Q, N * N / pow(B, 2.0)) << endl;
        }
        found = false;
        mat.convertVal(i);
        checkSymmetry(mat);
        for (unsigned long int j = 0; j < paraCount.currentGroupCount(); j++) {
            mat2.convertVal(paraCount.matValAt(j));

            for (int k = 0; k < symmetryList.size(); k++) {
                if (symmetryList[k].equalTo(mat2)) {
                    //writeMatGroup(mat, j);
                    catHash[i % (unsigned long int) (pow(Q, N * N / pow(B, 2.0)))][0] = j;
                    catHash[i % (unsigned long int) (pow(Q, N * N / pow(B, 2.0)))][1] = 0;
                    found = true;
                    break;
                }
            }

        }
        if (found == false) {
	  mat.printMat();
            paraCount.newGroup(i);
            catHash[i % (unsigned long int) (pow(Q, N * N / pow(B, 2.0)))][0] = paraCount.currentGroupCount() - 1;
            catHash[i % (unsigned long int) (pow(Q, N * N / pow(B, 2.0)))][1] = 1;
            //writeMatGroup(mat,paraCount.currentGroupCount() - 1);

        }

    }

    std::vector<int> renormHash;
    //create renormalisation hash table (BxB lattices)
    matConvert renormMatHash = matConvert(1, Q, B);
    renormMatHash.convertVal(0);
    renormHash.push_back(renormMatHash.findRN());
    for (int i = 1; i <= pow(Q, B * B); i++) {
        renormMatHash.convertVal(i);
        renormHash.push_back(renormMatHash.findRN());
    }
    matConvert parent = matConvert(1, Q, N);
    int renorm[N / B][N / B];
    int patch[N / B][N / B];
    unsigned long int hash = 0;
    int* pcatHash[(int) (pow(Q, N * N / (B * B)))];
    *pcatHash = catHash[0];
    int x, y;
    int group;
    int curCat;
    int parentCounter[paraCount.currentGroupCount()];
    //int God[paraCount.currentGroupCount()][(unsigned long int)pow(Q,N*N*(1 - 1/(B*B)))][paraCount.currentGroupCount()];
    std::vector<std::vector<std::vector< int > > > God;
    std::vector<int> tmp;
    God.resize(paraCount.currentGroupCount());
    tmp.resize(paraCount.currentGroupCount());
    cout << "God has been initialised" << endl;
    for (unsigned long int n = 1; n <= pow(Q, N * N); n++) {
        //parent loop


      parent.convertVal(n);
        for (int i = 0; i < N / B; i++) {
            for (int j = 0; j < N / B; j++) {
                hash = 0;
                for (int k = 0; k < B; k++) {
                    for (int l = 0; l < B; l++) {
                        hash += parent.lat[B * i + k][B * j + l] * pow(Q, k*B + l);
                    }
                }
                renorm[i][j] = renormHash[hash];

                //if is paradigm, do

            }
        }
        hash = 0;
        for (int m = 0; m < N / B; m++) {
            for (int p = 0; p < N / B; p++) {
                hash += renorm[m][p] * pow(Q, m*N/B + p);
            }
        }
        //if renomarlised is paradigm


        if (catHash[hash][1]) {
            curCat = catHash[hash][0];
            for (int b = 0; b < paraCount.currentGroupCount(); b++) {
                tmp[b] = 0;
            }

            for (int m = 0; m < N; m++) {
                for (int p = 0; p < N; p++) {
                    hash = 0;
                    x = m;
                    y = p;
                    for (int k = 0; k < N / B; k++) {
                        for (int l = 0; l < N / B; l++) {
                            hash += pow(Q, k*N/B + l) * parent.lat[(x + k) % N][(y + l) % N];
                        }
                    }

                    group = catHash[hash][0];
                    tmp[group]++;
                }
            }
            God[curCat].push_back(tmp);
        }
        if (n % 3000000 == 0) {
            cout << n / pow(Q, N * N) << endl;
        }
    }
    cout << "God matrix sizing: " << endl;
    for(int i = 0; i < God.size(); i++){
      cout << God[i].size() << endl;
    }
    cout << endl;
    cout << "Outputting God...";
    outputGod(God);
    cout << "Done" << endl;
    newtonRhapsonSolver(God);
    return 0;
}

void outputGod(std::vector<std::vector<std::vector<int > > > &God) {
    ostringstream groupss;
    groupss << ".txt";
    stringstream ss;
    string tempN;
    for (int i = 0; i < God.size(); i++) {
        ss.clear();
        ss.str(std::string());
        ss << i;


        tempN.clear();
        tempN = ss.str();
        tempN.append(groupss.str());
        std::fstream out;
        out.open(tempN.c_str(), std::ios::out | std::ios::app);
        for (int j = 0; j < God[i].size(); j++) {
            for (int k = 0; k < God[i][j].size(); k++) {
                out << God[i][j][k] << "\t";
            }
            out << endl;
        }
        out.close();
        //output parent x cat matrix in separate file for each
    }
}

void newtonRhapsonSolver(std::vector<std::vector<std::vector< int  > > > &God){
  
  double coupling[God.size()];
  double couplingDeriv[God.size()];
  double Jacobian[God.size()][God.size()];
  double JacobianInverse[God.size()][God.size()];
  double bSum[God.size()];
  double bSumDeriv[God.size()][God.size()];
  double identity[God.size()][God.size()];
  double premult = (double)B/(double)N;
  double con =  2;
  double mult = 2;
  srand(time(0));
  for(int i = 0; i < God.size(); i++){
    //coupling[i] = 0;
    coupling[i] = mult*rand()/(double)RAND_MAX - con;
    for(int j = 0; j < God.size(); j++){
      if(i == j){
	identity[i][j] = 1;
      }
      else{
	identity[i][j] = 0;
      }
    }

  }
  //coupling[0] = -log(2);
  bool fin = false;
  //loop
  //coupling[0] = -0.843;
  //coupling[1] = -0.178;
  //coupling[2] = -0.453;
  //coupling[3] = 0.774;
  while(!fin){



  for(int i = 0; i < God.size(); i++){

  }



  for(int i = 0; i < God.size(); i++){
    bSum[i] = 0;
    for(int j = 0; j < God.size(); j++){
      bSumDeriv[i][j] = 0;

    }
 }
  double temp = 0;

  for(int i = 0; i < God.size(); i++){
    for(unsigned long int j = 0; j < God[i].size(); j++){
      temp = (double)N*N*coupling[0];

      for(int k = 1; k < God.size(); k++){
	temp += God[i][j][k]*coupling[k];
      }
      bSum[i] += exp(temp);
    }
  }

  for(int i = 0; i < God.size(); i++){

  }

  for(int i = 0; i < God.size(); i++){
    for(int m = 0; m < God.size(); m++){
      for(unsigned long int j = 0; j < God[i].size(); j++){
	temp = (double)N*N*coupling[0];
	for(int k = 1; k < God.size(); k++){
	  temp += God[i][j][k]*coupling[k];
	}
	bSumDeriv[i][m] += God[i][j][m]*exp(temp);
      }
    }
  }
  couplingDeriv[0] = pow(premult,2.0)*log(bSum[0]);

  for(int i = 1; i < God.size(); i++){
    couplingDeriv[i] = pow(premult,2.0)*log(bSum[i]/bSum[0]);

  }


  fin = true;
  for(int i = 0; i < God.size(); i++){
    if(abs(couplingDeriv[i] - coupling[i]) > tol){
      fin = false;
    }
  }
    
      
 
    

  double sum = 0;
  bool skip = false;
  if(fin){
    for(int i = 0; i < God.size(); i++){
      sum += abs(coupling[i]);
    }
    if(abs(sum - log((double)Q)) < 1e-1){

      for(int i = 0; i < God.size(); i++){
	coupling[i] = mult*rand()/(double)RAND_MAX - con;
      }
      fin = false;
      skip = true;
    }
    
  }
  

  if(fin){
    break;
  }
  if(!skip){
  for(int i = 0; i < God.size(); i++){
    Jacobian[0][i] = pow(premult,2.0)*bSumDeriv[0][i]/bSum[0];
  }
  
  for(int i = 1; i < God.size(); i++){
    for(int j = 0; j < God.size(); j++){
      Jacobian[i][j] = pow(premult,2.0)*(bSumDeriv[i][j]/bSum[i]) - Jacobian[0][j];
    }
  }
  
  for(int i = 0; i < God.size(); i++){
    for(int j = 0; j < God.size(); j++){
      Jacobian[i][j] = Jacobian[i][j] - identity[i][j];
    }
  }


  mat coup2;
  coup2.set_size(God.size(),1);
  mat Jac2;
  Jac2.set_size(God.size(),God.size());

  for(int i = 0; i < God.size(); i++){

    coup2(i,0) = couplingDeriv[i] - coupling[i];
    for(int j = 0; j < God.size(); j++){

      Jac2(i,j) = Jacobian[i][j];
    }
    
  }

  Jac2 = Jac2.i();
  coup2 = Jac2*coup2;

  double temp2 = 0;

  for(int i = 0; i < God.size(); i++){
    coupling[i] -= coup2(i,0);
  }

  for(int i = 0; i < God.size(); i ++){
    if(abs(coupling[i]) > 3){
      for(int j = 0; j < God.size(); j++){
	coupling[j] = mult*rand()/(double)RAND_MAX - con;
	
      }
      break;
    }
  }
  }
  skip = false;

  }
  
  for(int i = 0; i < God.size(); i++){
    Jacobian[i][i]++;
  }
    
  mat JacT;
  JacT.set_size(God.size(), God.size());
  for(int i = 0; i < God.size(); i++){
    for(int j = 0; j < God.size(); j++){
      JacT(i,j) = Jacobian[i][j];
      cout << Jacobian[i][j] << "\t";
    }
    cout << endl;
  }
  cout << endl;
  //calculate eigen vectors
  cx_vec eigval;
  cx_mat eigvec;
  eig_gen(eigval, eigvec, JacT);
  for(int i = 0; i < eigval.n_elem; i++){
    cout << eigval(i) << " ";
    //cout << log((double)B)/log(eigval(i)) << " ";
  }
}

std::vector<int> calcCat(groupCounter &paraCount, matConvert &pLat, int* catHash[]) {
    std::vector<int> paraMemberCounter;
    paraMemberCounter.resize(paraCount.currentGroupCount());
    matConvert paraLatTemp = matConvert(1, Q, N / B);
    //Initialise
    for (int i = 0; i < paraMemberCounter.size(); i++) {
        paraMemberCounter[i] = 0;
    }
    matConvert curLat = matConvert(1, Q, N / B);
    int x, y;
    int group;
    for (int i = 0; i < pLat.N; i++) {
        for (int j = 0; j < pLat.N; j++) {
            x = i;
            y = j;
            for (int l = 0; l < N / B; l++) {
                for (int p = 0; p < N / B; p++) {
                    curLat.lat[l][p] = pLat.lat[(x + l) % N][(y + p) % N];

                }

            }
            group = (*catHash)[curLat.returnValue()*2];
            paraMemberCounter[group]++;

        }
    }
    return paraMemberCounter;
}

int boundary(int tmp, int Ntemp) {
    if (tmp == Ntemp - 1) {
        return 0;
    } else {
        return tmp++;
    }
}

void writeMatGroup(matConvert &mat, unsigned long int group) {
    ostringstream groupss;
    groupss << group;
    groupss << ".txt";
    string tempN;
    tempN = filePrefix;
    tempN.append(groupss.str());
    std::fstream out;

    out.open(tempN.c_str(), std::ios::out | std::ios::app);

    for (int i = 0; i < N / B; i++) {
        for (int j = 0; j < N / B; j++) {
            out << mat.lat[i][j];
        }
    }

    out << "\n";
    out.close();

}

bool checkSymmetry(matConvert mat) {
    symmetryList.clear();
    matConvert temp = mat;
    for (int rotNum = 0; rotNum < 4; rotNum++) {
        for (int parNum = 0; parNum < Q; parNum++) {

            for (int bndNum = 0; bndNum < N; bndNum++) {
                for (int w = 0; w < N; w++) {
                    for (int oneTwo = 0; oneTwo < 2; oneTwo++) {

                        temp = mat;
                        symmetryList.push_back(temp);
                        mat.oneTwoSwitch();
                    }
                    mat.nextBoundaryY();
                }

                mat.nextBoundaryX();
            }

            mat.nextParity();

        }

        mat.rotate();
    }

    return false;
}

bool parWrap(std::vector<int> qList, matConvert tmat, matConvert comp, matConvert ref) {
    for (int i = 0; i < Q; i++) {
        if (parityCheck(qList, tmat, comp, ref, i)) {
            return true;
        }
    }
    return false;
}

bool parityCheck(std::vector<int> qList, matConvert tmat, matConvert comp, matConvert ref, int cur) {
    matConvert mat = tmat;
    std::vector<int> tList = qList;
    for (int i = 0; i < qList.size(); i++) {
        tList = qList;
        mat = tmat;
        mat.exchangeSpin(qList[i], cur, ref);
        tList.erase(tList.begin() + i);
        cout << qList[i] << endl;
        //need a way to check that the number of spin exchanges is "closed"
        if (mat.countQ(ref) && mat.equalTo(comp)) {
            return true;
        } else if (tList.size() == 0) {
            return false;
        } else if (parityCheck(tList, mat, comp, ref, qList[i])) {
            return true;
        }
    }
    return false;
}

void buildFileNamePrefix() {
    ostringstream temp;
    filePrefix.append("_N");
    temp << N;
    filePrefix.append(temp.str());
    filePrefix.append("_Q");
    temp.clear();
    temp.str("");
    temp << Q;
    filePrefix.append(temp.str());
    temp.clear();
    temp.str("");
    temp << Z3;
    filePrefix.append("_Z3");
    filePrefix.append(temp.str());
    temp.clear();
    temp.str("");
    temp << Z2;
    filePrefix.append("_Z2");
    filePrefix.append(temp.str());
    temp.clear();
    temp.str("");
    temp << ROT;
    filePrefix.append("_ROT");
    filePrefix.append(temp.str());
    temp.clear();
    temp.str("_GROUP_");
    filePrefix.append(temp.str());
}

