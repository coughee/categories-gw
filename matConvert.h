/* 
 * File:   matConvert.h
 * Author: Jonathan
 *
 * Created on 09 April 2014, 00:08
 */

#ifndef MATCONVERT_H
#define	MATCONVERT_H
#include <vector>
class matConvert {
public:
    unsigned long int val;
    int base;
    int N;
    int parTot;
    int bndTot;
    int curBoundaryX;
    int curBoundaryY;
    std::vector<std::vector<int > > boundaryList;
    std::vector<std::vector<int > > lat;
    std::vector<std::vector<int > > parList;
    std::vector<matConvert > symmetryTable;
    matConvert(unsigned long int value, int Q, int Nval);
    matConvert(const matConvert& orig);
    virtual ~matConvert();
    void convertVal(unsigned long int value);
    bool equalTo(matConvert &mat);
    void nextBoundaryX();
    void nextBoundaryY();
    void nextParity();
    void nextParity(int Qmax);
    void rotate();
    void makeBoundaryList();
    void makeParityList();
    int findRN();
    void exchangeSpin(int s1, int s2, matConvert ref);
    bool countQ(matConvert ref);
    void printMat();
    void oneTwoSwitch();
    matConvert findRenormalised(int B, std::vector<int> &renormHash);
    int returnValue();
private:

};

#endif	/* MATCONVERT_H */

