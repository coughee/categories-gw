/* 
 * File:   matConvert.cpp
 * Author: Jonathan
 * 
 * Created on 09 April 2014, 00:08
 */

#include "matConvert.h"
#include "Random.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

matConvert::matConvert(unsigned long int value, int Q, int Nval) {
    this->base = Q;
    this->val = value;
    this->N = Nval;
    this->curBoundaryX = 0;
    this->curBoundaryY = 0;
    parList.resize(0);
    boundaryList.resize(0);
    convertVal(value);
    makeBoundaryList();
    makeParityList();
    srand(time(0));
}

int matConvert::returnValue(){
    int returnVal = 0;
    for(int i = 0; i < this->N; i++){
        for(int j = 0; j < this->N; j++){
            returnVal += lat[i][j]*pow(this->base,i*N + j);
        }
    }
    return returnVal;
}

void matConvert::convertVal(unsigned long int value){
    this->curBoundaryX = 0;
    this->curBoundaryY = 0;
    int temp [N*N];
    this->val = value;
    int rem;
    int count = 0;
    lat.clear();
    lat.resize(N);
    for(int i = 0; i < N; i++){
        lat[i].resize(N);
    }
    
    for (int i = 0; i < N * N; i++) {
        temp[i] = 0;
    }

    while (value) {
        rem = value % this->base;
        temp[count] = rem;
        value /= this->base;
        count++;
    }

    //Integer values have now been converted and saved to N*N 1d Array. Either we perform
    //Symmetry checks on those or we now convert them to 2D arrays;

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            lat[i][j] = temp[i*N + j];
        }
    }
}

matConvert matConvert::findRenormalised(int B,std::vector<int> &renormHash){
    matConvert returnMat = matConvert(1,this->base,N/B);
    matConvert tempMat = matConvert(1,this->base,B);
    for(int i = 0; i < N/B; i++){
        for(int j = 0; j < N/B; j++){
            for(int k = 0; k < B; k++){
                for(int l = 0; l < B; l++){
                    tempMat.lat[k][l] = lat[i*B + k][j*B + l];
                }
            }
            returnMat.lat[i][j] = renormHash[tempMat.returnValue()];
        }
    }
    return returnMat;
}
int matConvert::findRN(){
    int qVals[base];
    for(int i = 0; i < base; i++){
        qVals[i] = 0;
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            qVals[lat[i][j]]++;
        }
    }
    int max = qVals[0];
    int qMaxVal = 0;
    for(int i = 1; i < base; i++){
        if(qVals[i] > max){
            max = qVals[i];
            qMaxVal = i;
        }
    }
    std::vector< int > Degens;
    for(int i = 0; i < base; i++){
        if(qVals[i] == max){
            Degens.push_back(i);
        }
    }
    
    if (Degens.size() > 1) {
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                for(int k = 0; k < Degens.size(); k++){
                    if(lat[i][j] == Degens[k]){
                        return Degens[k];
                    }
                }
            }
        }
    }
    
    return qMaxVal;
}


bool matConvert::equalTo(matConvert &mat){
    if (mat.N != N){
        std::cout << "Error: Matrices must be same size for comparison" << std::endl;
    }
    for(int i = 0; i < this->N; i++){
        for(int j = 0; j < this->N; j++){
            if(this->lat[i][j] != mat.lat[i][j]){
                return false;
            }
        }
    }
    return true;
}

void matConvert::exchangeSpin(int s1, int s2, matConvert ref){
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(lat[i][j] == s1 && ref.lat[i][j] == s1){
                lat[i][j] = s2;
            }
        }
    }
}
bool matConvert::countQ(matConvert ref){
    std::vector<int> QcounterM;
    std::vector<int> QcounterR;
    QcounterM.resize(base);
    QcounterR.resize(base);
    for(int i = 0; i < base; i++){
        QcounterM[i] = 0;
        QcounterR[i] = 0;
    }
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                QcounterM[lat[i][j]]++;
                QcounterR[ref.lat[i][j]]++;
            }
        }
    bool Qcor[base];
    for(int i = 0; i < base; i++){
        Qcor[i] = false;
    }
    for(int i = 0; i < base; i++){
        for(int j = 0; j < QcounterR.size(); j++){
            if(QcounterM[i] == QcounterR[j]){
                Qcor[i] = true;
                QcounterR.erase(QcounterR.begin() + j);
                j--;
            }
        }
    }
    for(int i = 0; i < base; i++){
        if(!Qcor[i]){
            return false;
        }
    }
    
    return true;
}
void matConvert::nextBoundaryX(){
    matConvert temp = matConvert(1,base,N);
    temp.lat = this->lat;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if(j > 0){
                lat[i][j] = temp.lat[i][j - 1];
            }
            else{
                lat[i][j] = temp.lat[i][N - 1];
            }

        }
    }
    
        
}

void matConvert::nextBoundaryY(){
    matConvert temp = matConvert(1,base,N);
    temp.lat = this->lat;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if(i > 0){
                lat[i][j] = temp.lat[i - 1][j];
            }
            else{
                lat[i][j] = temp.lat[N - 1][j];
            }

        }
    }
    
}
void matConvert::nextParity(){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            lat[i][j]++;
            if(lat[i][j] > base - 1){
                lat[i][j] = 0;
            }
        }
    }
}

void matConvert::oneTwoSwitch(){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(lat[i][j] == 0){
                lat[i][j] = 1;
            }
            else if (lat[i][j] == 1){
                lat[i][j] = 0;
            }
        }
    }
}

void matConvert::rotate(){
    matConvert temp = matConvert(1,base,N);
    temp.lat = this->lat;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            this->lat[i][j] = temp.lat[N - j - 1][i];
        }
    }
}

void matConvert::makeBoundaryList(){
    bndTot = 1;
}

void matConvert::makeParityList(){
    parTot = base;
}

matConvert::matConvert(const matConvert& orig) {

    this->val = orig.val;
    this->base = orig.base;
    this->N = orig.N;
    this->parTot = orig.N;
    this->bndTot = orig.bndTot;
    this->curBoundaryX = orig.curBoundaryX;
    this->curBoundaryY = orig.curBoundaryY;
    this->boundaryList = orig.boundaryList;
    this->lat = orig.lat;
    this->parList = orig.parList;
}


matConvert::~matConvert() {
}

void matConvert::printMat(){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            std::cout << lat[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

