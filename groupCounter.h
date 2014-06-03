/* 
 * File:   groupCounter.h
 * Author: Jonathan
 *
 * Created on 08 April 2014, 06:59
 */

#ifndef GROUPCOUNTER_H
#define	GROUPCOUNTER_H
#include <vector>
class groupCounter {
public:
    std::vector<unsigned long int> matVals;
    groupCounter();
    groupCounter(const groupCounter& orig);
    virtual ~groupCounter();
    unsigned long int matValAt(int id);
    unsigned long int currentGroupCount();
    void newGroup(unsigned long int matVal);
private:

};

#endif	/* GROUPCOUNTER_H */

