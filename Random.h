/* 
 * File:   Random.h
 * Author: Jonathan
 *
 * Created on 09 April 2013, 03:13
 */

#ifndef RANDOM_H
#define	RANDOM_H

class Random {
public:
    int min, max;
    
    Random(int min, int max);
    Random();
    Random(const Random& orig);
    virtual ~Random();
    int Next();
    double Rand();
    
private:

};

#endif	/* RANDOM_H */

