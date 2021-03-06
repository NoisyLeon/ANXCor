#ifndef SACREC_H
#define SACREC_H

#include "mysac64.h"
#include "MyOMP.h"
//#include <cstddef>
#include <iostream>
#include <string>
#include <memory>
#include <limits>

#ifndef FuncName
#define FuncName __FUNCTION__
#endif
/*
// a workaround for nullptr when compiled by old gcc compiler
const class {				// this is a const object...
public:
   template<class T>			// convertible to any type
   operator T*() const { return 0 }	// of null non-member pointer...
   template<class C, class T>		// or any type of null
   operator T C::*() const { return 0; }// member pointer...
private:
  void operator&() const;		// whose address can't be taken
} nullptr = {};				// and whose name is nullptr
*/


class SacRec
{
public:
    std::string fname;			// input file name
    SAC_HD shd;				// sac header
    std::unique_ptr<float[]> sig;	// pointer to the signal
    //std::auto_ptr<float> sig;
public:
    /* ------------------------------ con/destructors and operators ------------------------------ */
    /* constructors */
    SacRec( std::ostream& reportin = std::cerr );	// default 1
    SacRec( const std::string& fnamein, std::ostream& reportin = std::cerr );	// default 2
    SacRec( const SacRec& recin );		// copy
    SacRec( SacRec&& recin );			// move
    /* operators */
    SacRec &operator= ( const SacRec& recin );	// assignment
    SacRec &operator= ( SacRec&& recin );	// move
    /* destructor */
    ~SacRec();

    /* assign header and allocate memory */
    void MutateAs ( const SacRec& recin );

    /* ------------------------------ sac file read/write ------------------------------ */
    /* load sac header from file 'fname' */
    void LoadHD( const std::string& fnamein )
    {
        fname = fnamein;
        LoadHD();
    }
    void LoadHD();
    /* read sac header+signal from file 'fname', memory is allocated on heap */
    void Load( const std::string& fnamein )
    {
        fname = fnamein;
        Load();
    }
    void Load();
    /* write to file '*fname' */
    void WriteHD( const std::string& fname );
    void Write( const std::string& fname );
    /* dump signal to txt */
    void Dump( const std::string fname );
    /* clear sac and release memory */
    void clear()
    {
        sig.reset();
        shd = sac_null;
        fname.clear();
    }

    /* ------------------------------ header operations ------------------------------ */
    void ChHdr(const std::string& field, const std::string& value);

    /* ------------------------------ header/signal information ------------------------------ */
    /* compute the absolute time in sec relative to 1900.01.00 */
    double AbsTime ();
    /* update/reformat header time if shd.nzmsec is modified and is out of the range [0,1000) */
    void UpdateTime();
    /* search for min&max signal positions and amplitudes */
    void MinMax (int& imin, int& imax, float tbegin, float tend);
    void MinMax ( float tbegin, float tend, float& tmin, float& min, float& tmax, float& max );
    /* compute the root-mean-square average in a given window */
    void RMSAvg ( float tbegin, float tend, float& rms )
    {
        RMSAvg( tbegin, tend, 1, rms);
    }
    void RMSAvg ( float tbegin, float tend, int step, float& rms );
    bool MeanStd ( float tbegin, float tend, float& mean, float& std )
    {
        return MeanStd(tbegin, tend, 1, mean, std);
    }
    bool MeanStd ( float tbegin, float tend, int step, float& mean, float& std );
    //Check precossor
    bool CheckPrecNoise();
    /* ------------------------------ single-sac operations ------------------------------ */
    void Mul( const float mul );
    void Addf( const SacRec& sac2 );
    void Divf( const SacRec& sac2 );
    void PullUpTo( const SacRec& sac2 );
    void ToAm()
    {
        ToAm(*this);
    }
    void ToAm( SacRec& sac_am )
    {
        SacRec sac_ph;
        ToAmPh( sac_am, sac_ph );
    }
    void ToAmPh( SacRec& sac_am, SacRec& sac_ph );	// in series when sig is large
    void ToAmPh_p( SacRec& sac_am, SacRec& sac_ph );	// always parallel
    void FromAmPh( SacRec& sac_am, SacRec& sac_ph, const short outtype = 0 );		// in series when sig is large
    void FromAmPh_p( SacRec& sac_am, SacRec& sac_ph, const short outtype = 0 );	//	always parallel
    /* filters */
    void LowpassFilt( double fh1, double fh2 )
    {
        LowpassFilt(fh1, fh2, *this);
    }
    void LowpassFilt( double fh1, double fh2, SacRec& srout )
    {
        Filter(-1., -1., fh1, fh2, srout);
    }
    void BandpassFilt( double f1, double f2, double f3, double f4 )
    {
        BandpassFilt(f1, f2, f3, f4, *this);
    }
    void BandpassFilt( double f1, double f2, double f3, double f4, SacRec& srout )
    {
        Filter(f1, f2, f3, f4, srout);
    }
    void GaussianFilt( double fc, double fhlen )
    {
        GaussianFilt(fc, fhlen, *this);
    }
    void GaussianFilt( double fc, double fhlen, SacRec& srout )
    {
        Filter(-1., fc, fhlen, -1., srout);
    }
    //else if( f1==-1. && f4==-1. ) pimpl->TaperGaussian( f2, f3, dom, nk, sf );
    /* method that performs 3 different types of filters
     * lowpass when ( (f1==-1. || f2==-1.) && (f3>0. && f4>0.) )
     * bandpass when ( f1>=0. && f2>0. && f3>0. && f4>0. )
     * gaussian when ( f1==-1. && f4==-1. ) where f2 = center freqency and f3 = frequency half length */
    void Filter ( double f1, double f2, double f3, double f4 )
    {
        Filter(f1, f2, f3, f4, *this);    // in-place (in series when sig is large)
    }
    void Filter ( double f1, double f2, double f3, double f4, SacRec& srout );				// out-of-place (in series when sig is large)
    void Filter_p ( double f1, double f2, double f3, double f4, SacRec& srout );			// always parallel
    /* cosine tapers */
    void cosTaperL( const float fl, const float fh );
    void cosTaperR( const float fl, const float fh );
    /* remove mean and trend */
    void RTrend();
    /* remove response and apply filter */
    void RmRESP( const std::string& fresp, float perl, float perh )
    {
        std::string evrexe;
        RmRESP( fresp, perl, perh, evrexe );
    }
    void RmRESP( const std::string& fresp, float perl, float perh, const std::string& evrexe );
    /* resample (with anti-aliasing filter) the signal to given sps */
    void Resample( float sps );
    /* smoothing ( running average ) */
    void Smooth( float timehlen, SacRec& sacout ) const;
    void Hilbert()
    {
        Hilbert(*this);
    }
    void Hilbert( SacRec& sacout );
    void Envelope()
    {
        Envelope(*this);
    }
    void Envelope( SacRec& sacout );

    void cut( float tb, float te )
    {
        cut(tb, te, *this);
    }
    void cut( float tb, float te, SacRec& );
    /* ------------------------------ inter-sac operations ------------------------------ */
    /* merge a second sacrec to the current */
    void Merge( SacRec sacrec2 )
    {
        merge( sacrec2 );
        const char *RECNAME = nullptr;
        arrange(RECNAME, 1000);
    }
    void merge( SacRec sacrec2 );
    int arrange( const char *recname , int FileNo );
    //int arrange( const char *recname = nullptr);
    void Correlate( SacRec& sac2 )
    {
        Correlate(sac2, *this);
    }
    /* Correlate with another sac record
    	ctype=0: Correlate (default)
    	ctype=1: deconvolve (sac.am/sac2.am)
    	ctype=2: deconvolve (sac2.am/sac.am) */
    void Correlate( SacRec& sac2, SacRec& sacout, int ctype = 0 );

    /* ------------------------------- cut by event ---------------------------------- */
    void ZoomToEvent( const std::string etime, float evlon, float evlat, float tb, float tlen, std::string ename = "" );
    void ZoomToEvent( const SAC_HD& eshd, float evlon, float evlat, float tb, float tlen, std::string ename );

    /* ------------------------------- temporal normalizations ------------------------------- */
    void OneBit();
    void RunAvg( float timehlen, float Eperl, float Eperh );

    /* ------------------------------- memory consumed ------------------------------- */
    float MemConsumed() const;
    void AlwaysParallel()
    {
        maxnpts4parallel = std::numeric_limits<int>::max();
    }
    // to run the fftw, 16 times the original npts is required ( in&out complex double array with size doubled for specturm ). 20 is used to be safe
    void SetMaxMemForParallel( float MemInMb )
    {
        maxnpts4parallel = (MemInMb * 1024. * 1024. - 1000.) / (4. * 20.);
    }

protected:
    static constexpr float NaN = -12345.;
    int maxnpts4parallel = 1e6;

private:
    /* impl pointer */
    struct SRimpl;
    std::unique_ptr<SRimpl> pimpl;
    /* reporting stream */
    std::ostream* report = &(std::cerr);

};

/*
class AMPH {
private:
   class APimpl;
   std::unique_ptr<APimpl> pimpl;
   std::unique_ptr<float[]> am, ph;
};
*/

namespace ErrorSR
{

class Base : public std::runtime_error
{
public:
    Base( const std::string funcname, const std::string message )
        : runtime_error(funcname + ": " + message)
    {
        //PrintStacktrace();
    }
};

class BadFile : public Base
{
public:
    BadFile(const std::string funcname, const std::string info = "")
        : Base(funcname, "Invalid or non-accessable file ("+info+").") {}
};

class BadParam : public Base
{
public:
    BadParam(const std::string funcname, const std::string info = "")
        : Base(funcname, "Bad parameters ("+info+").") {}
};

class EmptySig : public Base
{
public:
    EmptySig(const std::string funcname, const std::string info = "")
        : Base(funcname, "No sac signal loaded in the memory ("+info+").") {}
};

class InsufData : public Base
{
public:
    InsufData(const std::string funcname, const std::string info = "")
        : Base(funcname, "Insufficient data points ("+info+").") {}
};

class SizeMismatch : public Base
{
public:
    SizeMismatch(const std::string funcname, const std::string info = "")
        : Base(funcname, "Size mismatch ("+info+").") {}
};

class ExternalError : public Base
{
public:
    ExternalError(const std::string funcname, const std::string info = "")
        : Base(funcname, "External error ("+info+").") {}
};

class MemError : public Base
{
public:
    MemError(const std::string funcname, const std::string info = "")
        : Base(funcname, "Memory error ("+info+").") {}
};
};

#endif
