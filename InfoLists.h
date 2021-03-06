#ifndef INFOLISTS_H
#define INFOLISTS_H

#include "SysTools.h"
#include <iostream>
#include <vector>
#include <string>


/* seed record wrapper consists of filename, year, month, and day */
struct SeedInfo
{
    std::string name;
    int year, month, day;
    SeedInfo() : year(0.), month(0.), day(0.) {}
    SeedInfo(const char* inname, const int& yy, const int& mm, const int& dd)
    {
        if(inname)name.assign(inname);
        year    = yy;
        month   = mm;
        day     = dd;
    }
    friend std::ostream& operator<< (std::ostream& o, const SeedInfo& sr)
    {
        o<<"( "<<sr.name<<" "<<sr.year<<" "<<sr.month<<" "<<sr.day<<" )";
        return o;
    }
};


/* station wrapper consists of stationname, network code, longitude, and latitude */
struct StaInfo
{
    std::string name, netcode;
    float lon, lat;
    int CCflag;
    StaInfo() : lon(0.), lat(0.) {}
    
    StaInfo(const char* inname, const char* innetcode, const float& lonin, const float& latin)
    {
        if(inname)name.assign(inname);
        if(innetcode)netcode.assign(innetcode);
        lon = lonin;
        lat = latin;
    }
    
    friend bool operator== (StaInfo& a, StaInfo& b)
    {
        return ( a.name.compare(b.name)==0 && a.lon==b.lon && a.lat==b.lat && a.netcode.compare(b.netcode)==0 );
    }
    
    friend std::ostream& operator<< (std::ostream& o, const StaInfo& sr)
    {
        o<<"( "<<sr.name<<"."<<sr.netcode<<" "<<sr.lon<<" "<<sr.lat<<" "<<" )";
        return o;
    }
    
    bool checkdoCC(const StaInfo & STA )
    {
        if  (CCflag==STA.CCflag && CCflag==0)
            return false;
        else if(CCflag != STA.CCflag && CCflag != 0 && STA.CCflag != 0)
            return false;
        else if (CCflag< 0 && name == STA.name && netcode == STA.netcode)
        {
            std::cout << "GROUP: "<<CCflag<<std::endl;
            std::cout<<"DO NOT DO AUTO"<<std::endl;
            return false;
        }
        else
            return true;
    }
    
};


/* Daily Info contains StaInfo and SeedInfo */
struct DailyInfoData
{
    SeedInfo seed;
    StaInfo sta;
    std::string rdsexe, evrexe;
    std::string chname;
    std::string osac_outname, fsac_outname;
    std::string rec_outname, resp_outname;
    std::string monthdir; 
    std::string outdir;
    float sps, perl, perh, t1, tlen;
    int tnorm_flag;
    float Eperl, Eperh, timehlen;
    float frechlen, memomax;

};

struct DailyInfo : public DailyInfoData
{
    //const std::string MonthName[13] {
    const std::vector<std::string> MonthName
    {
        "INVALID",
        "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
        "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
    };

    std::string& seedname   = seed.name;
    std::string& staname    = sta.name;
    std::string& netname    = sta.netcode; // added on 2019-03-04
    
    int& year   = seed.year;
    int& month  = seed.month;
    int& day    = seed.day;

    DailyInfo() {}
    DailyInfo( const SeedInfo sei, const StaInfo& sti, const std::string& chi )
    {
        Update(sei, sti, chi);
    }
    DailyInfo( const DailyInfo& di2 ) : DailyInfoData(di2) {}
    DailyInfo& operator= ( const DailyInfo& di2 )
    {
        *( reinterpret_cast<DailyInfoData*>(this) ) = di2;
        return *this;
    }

    void Update( const SeedInfo sei, const StaInfo& sti, const std::string& chi )
    {
        seed    = sei;
        sta     = sti;
        chname  = chi;
        std::string mdir    = std::to_string(year) + "." + MonthName[month];
        std::string ddir    = mdir + "." + std::to_string(day);
        outdir              = mdir + "/" + ddir;
        MKDirs( outdir.c_str() );
        std::string outname = ddir + "." + netname+"."+staname + "." + chname + ".SAC";
        osac_outname        = mdir + "/" + ddir + "/" + outname;
        fsac_outname        = mdir + "/" + ddir + "/ft_" + outname;
        rec_outname         = fsac_outname + "_rec"; 
        resp_outname        = mdir + "/" + ddir + "/" + "RESP." + netname+"."+staname + ".." + chname; // changed on 2019-03-04
        monthdir            = mdir;
    }

    friend std::ostream& operator<< (std::ostream& o, DailyInfo& di)
    {
        o<<di.seed<<"   "<<di.sta<<"   ("<<di.chname<<")";
        return o;
    }

};



/* ------------------------------ ListBase ------------------------------ */
template< typename T >
class ListBase
{
public:

    ListBase( std::ostream& reportin = std::cerr )
        : report( &reportin ), icurrent( list.end() ) {}
    /* call Load if input_info is provided */
    ListBase( const std::string& input_info, std::ostream& reportin = std::cerr )
        : report( &reportin )
    {
        Load(input_info);
    }
    /* load records either from an input file or directly from the input_info string */
    virtual void Load( const std::string& ) = 0;
    /* check if icurrent is meaningful */
    bool IsEnded()
    {
        return icurrent>=list.end() || icurrent<list.begin();
    }
    bool NotEnded()
    {
        return icurrent<list.end() && icurrent>=list.begin();
    }
    /* return an iterator to the current item in the list */
    typename std::vector<T>::iterator GetRec()
    {
        return icurrent;
    }
    /* rewind */
    void Rewind()
    {
        icurrent = list.begin();
    }
    /* get to the next record and return true on success */
    bool NextRec()
    {
        icurrent++;
        if( icurrent<list.begin() || icurrent>=list.end() ) return false;
        return true;
    }

protected:
    std::vector<T> list;
    typename std::vector<T>::iterator icurrent;
    std::ostream* report = &(std::cerr);
};


/* ------------------------------ Channellist ------------------------------ */
class Channellist : public ListBase<std::string>
{
public:
    /* load channel list from an input string */
    void Load( const std::string& chinfo );
};


/* ------------------------------ Seedlist ------------------------------ */
/* seed list with the 'NextRec' and the 'ReLocate' operation */
class Seedlist : public ListBase<SeedInfo>
{
public:
    /* load station records from an input file */
    void Load( const std::string& );
    /* search for the first match of the input Seed date and return true on success
       icurrent will be moved to the match on succed and to the first rec after on failure */
    bool ReLocate( int year, int month, int day );
    //~Seedlist() { if( seedrec ) delete seedrec; }
};


/* ------------------------------ Stationlist ------------------------------ */
/* station list with the 'NextRec' and the 'ReLocate' operation */
class Stationlist : public ListBase<StaInfo>
{
public:
    /* load station records from an input file */
    void Load( const std::string& );
    /* search for the first match of the input staname and return true on success
       icurrent=.end() if no such match is found */
    bool ReLocate( const std::string& staname, const std::string& netname );
};

#endif
