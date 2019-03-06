#ifndef CCREC_H
#define CCREC_H
#include <iostream>
#include <vector>
#include <string>
//#include <fftw3.h>
#include "InfoLists.h"



struct CC_todo
{
    StaInfo station1, station2;
    std::string monthdir, chan1, chan2;
    std::vector <int> daylst;
    
    //std::vector < std::string >  infname1, infname2;
    /* ------------------------------ con/destructors and operators ------------------------------ */
    /* constructors */
    CC_todo();
    CC_todo(const StaInfo & sta1, const StaInfo & sta2, const std::string & MDIR, const std::vector <int> & DLst);
    /* operators */

    /* destructor */
    ~CC_todo();
};










struct hole_rec
{
    int rec_b, rec_e;
    hole_rec();
    ~hole_rec();
    hole_rec(int BEG, int END);
};


#endif
