/* 
 * File:   Random.cpp
 * Author: Jonathan
 * 
 * Created on 09 April 2013, 03:13
 */

#include "Random.h"
#include "math.h"
#include <stdlib.h>
#include <time.h>

Random::Random() {
    srand(time(0));
}

Random::Random(const Random& orig) {
}

Random::~Random() {
}

Random::Random(int max_, int _min){
    this->max = max_;
    this->min = _min;
    srand(time(0));
}

int Random::Next(){
    return rand() % (max - min + 1) + min;
}

double Random::Rand(){
    return (double)rand()/(double)RAND_MAX;
}

