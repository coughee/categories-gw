/* 
 * File:   groupCounter.cpp
 * Author: Jonathan
 * 
 * Created on 08 April 2014, 06:59
 */

#include "groupCounter.h"

groupCounter::groupCounter() {
}

groupCounter::groupCounter(const groupCounter& orig) {
}

groupCounter::~groupCounter() {
}

unsigned long int groupCounter::currentGroupCount(){
    return matVals.size();
}

unsigned long int groupCounter::matValAt(int id){
    return matVals[id];
}

void groupCounter::newGroup(unsigned long int matVal){
    matVals.push_back(matVal);
}
