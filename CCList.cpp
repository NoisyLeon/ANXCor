#include "SacRec.h"
#include "CCDatabase.h"
#include "CCRec.h"
#include "MyOMP.h"
#include <iostream>
#include <deque>
#include <unistd.h>
#include "MyLogger.h"

bool FileExists(const char* filename);
bool FileExists(const std::string& filename);
bool doCor(SacRec & out_CC, const SacRec & sac_am1, const SacRec & sac_ph1, const SacRec & sac_am2, const SacRec & sac_ph2,
           const int & lagtime, const int & fs, const int & mintlen, const float & tlen, const int & ftlen, std::string & fname1, std::string & fname2)   ;


void GCCtodoList(std::vector < StaInfo > & StaList, std::string & MDIR,
        std::vector <CC_todo> & CC_List, std::vector < std::string > & CHAN, int skipccflag)
{
    const int CHsize    = CHAN.size();
    const int CCsize    = CHsize*CHsize;
    for (int s1=0; s1<StaList.size(); s1++)
        for (int s2=0; s2<StaList.size(); s2++)
        {
            if (StaList[s1].netcode+"."+StaList[s1].name > StaList[s2].netcode+"."+StaList[s2].name
                    || !StaList[s1].checkdoCC(StaList[s2]) )
                continue;
            bool skip_mcc   = true;
            // check the existence of monthly stacked cc file
            for (int f=0; f<CCsize; f++)
            {
                int ch1 = (int)(f/CHsize);
                int ch2 = f%CHsize;
                std::string month_fout  = MDIR+"/COR/"+StaList[s1].netcode+"."+StaList[s1].name+"/COR_" +
                                            StaList[s1].netcode+"."+StaList[s1].name + "_" +CHAN[ch1] + "_" +
                                                StaList[s2].netcode+"."+StaList[s2].name + "_" + CHAN[ch2] +".SAC";
                if (!FileExists(month_fout))
                {
                    skip_mcc    = false;
                    break;
                }
            }
            if(skip_mcc && skipccflag == 1)
            {
                //std::cout<<"CC File: "<<MDIR<<":"<<StaList[s1].netcode+"."+StaList[s1].name
                //        <<"_"<<StaList[s2].netcode+"."+StaList[s2].name<<" Exist! "<<std::endl;
                continue;
            }
            std::vector <int> daylst;
            for (int d=1; d<=31; d++)
            {
                bool all_existflag  = true;
                std::string infname1;
                std::string infname2;
                for (int ch=0; ch<CHsize; ch++)
                {
                    std::string daydir  = MDIR+"."+std::to_string(d);
                    infname1            = MDIR+"/"+daydir+"/ft_"+daydir+"."+
                                            StaList[s1].netcode+"."+StaList[s1].name+"."+CHAN[ch]+".SAC";
                    infname2            = MDIR+"/"+daydir+"/ft_"+daydir+"."+
                                            StaList[s2].netcode+"."+StaList[s2].name+"."+CHAN[ch]+".SAC";
                    if ( !FileExists(infname1+".am") || ! FileExists(infname1+".ph") ||
                            !FileExists(infname2+".am") || !FileExists(infname2+".ph"))
                    {
                        all_existflag   = false;
                        break;
                    }
                }
                if ( !all_existflag )
                    continue;
                daylst.push_back(d);
            }
            if (daylst.size() > 0)
                CC_List.push_back( CC_todo(StaList[s1], StaList[s2], MDIR, daylst) );
        }
}


void CCList2CC( std::vector<CC_todo> & CC_List, const std::vector < std::string > & CHAN, const CCPARAM & cdbParams)
{
    const int CHsize    = CHAN.size();
    const int CCsize    = CHsize*CHsize;
    int fskipcc         = cdbParams.fskipcrco;
    int ftlen           = cdbParams.ftlen;
    float tlen          = cdbParams.tlen;
    int mintlen         = cdbParams.mintlen;
    int SIZE            = CC_List.size();
    
    if (fskipcc != 2 && SIZE != 0)
    {
        int i               = -1;
        int si              = -1;
        const int fs        = cdbParams.sps;
        const int LagTime   = cdbParams.lagtime;
        const int Coutflag  = cdbParams.CorOutflag;
        const int checkprec = cdbParams.fprcs;
        
        #pragma omp parallel private ( i ) shared( CC_List, si)
        {
            while (1)
            {
                #pragma omp critical
                {
                    si++;
                    i   = si;
                }
                if(i >= SIZE)
                    break;
                //#pragma omp critical 
                //{
                //    std::cout<<i<<std::endl;
                //}
                std::string staid1              = CC_List[i].station1.netcode+"."+CC_List[i].station1.name;
                std::string staid2              = CC_List[i].station2.netcode+"."+CC_List[i].station2.name;
                std::string cor_dir             = CC_List[i].monthdir + "/COR/" + staid1;
                const std::vector <int>  daylst = CC_List[i].daylst;
                
                std::vector < SacRec > m_sacV(CCsize);          // monthly sac vector
                bool init_month_sac_flag        = false;        // if the monthly sac vector has been initialized or not
                std::vector <std::string> m_foutname(CCsize);   // monthly output file name list
                // set montly output file name
                for (int j=0; j<CCsize; j++)
                {
                    int ch1         = (int)(j/CHsize); 
                    int ch2         = j%CHsize;
                    m_foutname[j]   = CC_List[i].monthdir+ "/COR/" + staid1 + "/COR_" +
                                                staid1 + "_" +CHAN[ch1] + "_" + staid2 + "_" + CHAN[ch2] +".SAC";
                }
                // loop over days
                for( int iday=0; iday < daylst.size(); iday++ )
                {
                    std::vector < SacRec > d_sacV(CCsize);          // vector includes daily cross-correlation data
                    std::vector < std::string> foutnameD(CCsize);   // daily output file name list
                    std::vector < SacRec > sac_amV1(CHsize), sac_phV1(CHsize), sac_amV2(CHsize), sac_phV2(CHsize);
                    
                    const int day           = daylst[iday];
                    const std::string daydir= CC_List[i].monthdir+"."+std::to_string(day); 
                    bool skipccflag         = false;            // skip this day or not
                    // Read amp/ph files
                    for (int ch=0; ch<CHsize; ch++)
                    {
                        const std::string infname1  = CC_List[i].monthdir+"/"+daydir+
                                                        "/ft_"+daydir+"."+staid1+"."+CHAN[ch]+".SAC";
                        const std::string infname2  = CC_List[i].monthdir+"/"+daydir+
                                                        "/ft_"+daydir+"."+staid2+"."+CHAN[ch]+".SAC";
                        sac_amV1[ch].Load(infname1+".am");
                        sac_phV1[ch].Load(infname1+".ph");
                        sac_amV2[ch].Load(infname2+".am");
                        sac_phV2[ch].Load(infname2+".ph");
                    }
                    // Loop over cross-correlation channels
                    for (int j=0; j<CCsize; j++)
                    {
                        int ch1     = (int)(j/CHsize); 
                        int ch2     = j%CHsize;
                        foutnameD[j]= CC_List[i].monthdir+ "/COR_D/" + staid1 + "/COR_" +
                                      staid1 + "_" +CHAN[ch1] + "_" + staid2 + "_" + CHAN[ch2]  +
                                      "_" + std::to_string(day) + ".SAC";
                        //auto & d_sac                = d_sacV[j];
                        auto & foutname             = foutnameD[j];
                        const std::string infname1  = CC_List[i].monthdir+"/"+daydir+
                                                        "/ft_"+daydir+"."+staid1+"."+CHAN[ch1]+".SAC";
                        const std::string infname2  = CC_List[i].monthdir+"/"+daydir+
                                                        "/ft_"+daydir+"."+staid2+"."+CHAN[ch2]+".SAC";
                        // Start to do CC                        
                        //d_sacV[j].LoadHD(infname1+".am");
                        d_sacV[j].MutateAs(sac_amV1[ch1]);
                        std::string frec1           = infname1 + "_rec";
                        std::string frec2           = infname2 + "_rec";
                        
                        if  ( !doCor(d_sacV[j], sac_amV1[ch1], sac_phV1[ch1], sac_amV2[ch2], sac_phV2[ch2],
                                     LagTime, fs, mintlen, tlen, ftlen, frec1, frec2) )
                        {
                            
                            std::string mesg    = "Error am&ph File CC: " + staid1 + "_" + CHAN[ch1]  +
                                                " with " + staid2 + "_" + CHAN[ch2] + " Date: " +
                                                    daydir + " ---- List No. " + std::to_string(i+1);
                            #pragma omp critical 
                            {
                                std::cout<<mesg<<std::endl;
                            }
                            skipccflag  = true;
                            break;
                        }
                        else
                        {
                            std::string mesg    = " Do CC: " + staid1 + "_" + CHAN[ch1]  +
                                  " with " + staid2 + "_" + CHAN[ch2] + " Date: " +
                                  daydir + " ---- List No. " + std::to_string(i+1);
                            //logger.Hold( INFO, mesg, FuncName );
                            #pragma omp critical 
                            {
                                std::cout<<mesg<<std::endl;
                            }
                        }
                        
                        if (checkprec==1) // NEED TEST!
                        {
                            if (d_sacV[j].CheckPrecNoise())
                                skipccflag  = true;
                            break;
                        }
                        // save daily cross-correlation results
                        if (Coutflag != 0 )
                        {
                            MKDirs((CC_List[i].monthdir + "/COR_D/" + staid1).c_str());
                            d_sacV[j].Write(foutname); // Save daily CC data
                        }
                    }
                    ////if (skipccflag)
                    ////{
                    ////    for (int j=0; j<CCsize; j++)
                    ////    {
                    ////        if (d_sacV[j].sig)
                    ////            d_sacV[j].sig.reset();
                    ////    }
                    ////}
                    // append monthly cross-correlation results
                    if( Coutflag != 1 && !skipccflag )
                    {
                        if (init_month_sac_flag)
                        {
                            for (int j=0; j<CCsize; j++)
                            {
                                m_sacV[j].shd.user0 += 1;
                                for (int k=0; k < d_sacV[j].shd.npts; k++)
                                    m_sacV[j].sig[k]+= d_sacV[j].sig[k];
                            }
                        }
                        else
                        {
                            for (int j=0; j<CCsize; j++)
                                 m_sacV[j]      = std::move(d_sacV[j]);
                            init_month_sac_flag = true;
                        }
                    }
                    // clear memory
                    for (int ch=0; ch<CHsize; ch++)
                    {
                        sac_amV1[ch].clear();
                        sac_phV1[ch].clear();
                        sac_amV2[ch].clear();
                        sac_phV2[ch].clear();
                    }
                    sac_amV1.clear();
                    sac_phV1.clear();
                    sac_amV2.clear();
                    sac_phV2.clear();
                    
                    for (int j=0; j<CCsize; j++)
                    {
                        d_sacV[j].clear();
                    }
                    d_sacV.clear();
                }
                // save monthly cross-correlation results
                if( Coutflag != 1 && init_month_sac_flag)
                {
                    MKDirs((CC_List[i].monthdir + "/COR/" + staid1).c_str());
                    for (int j=0; j<CCsize; j++)
                    {
                        auto & m_sac    = m_sacV[j];
                        m_sac.Write(m_foutname[j]);
                    }
                }
                #pragma omp critical 
                {
                    std::string stapair_str = staid1 + "_"+staid2;
                    std::cout<<stapair_str<<" ";
                }
                // clear memory
                for (int j=0; j<CCsize; j++)
                {
                    m_sacV[j].clear();
                }
            }
        }
        std::cout<<std::endl;
    }
}

