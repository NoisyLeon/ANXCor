#define DEBUG
#define MAIN
#include "InfoLists.h"
#include "SeedRec.h"
#include "SacRec.h"
#include "CCDatabase.h"
#include "MyLogger.h"
#include "MyOMP.h"
#include "SysTools.h"
#include "CCRec.h"
#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <unistd.h>

/* normalize all sac files in sacV (simultaneously if SyncNorm==true) by Ye Tian*/
void TNormAll( std::deque<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm);
/* normalize and apply taper*/
void FNormAll( std::deque<SacRec>& sacV, const std::vector<DailyInfo>& dinfoV, bool SyncNorm);
bool FileExists(const char* filename);
bool FileExists(const std::string & filename);

void GCCtodoList(std::vector < StaInfo > & StaList, std::string & MDIR,
                 std::vector<CC_todo> & CC_List, std::vector < std::string > & CHAN, int skipccflag);
/* convert amp & ph files to xcorr */
void CCList2CC( std::vector<CC_todo> & CC_List, const std::vector < std::string > & CHAN, const CCPARAM & cdbParams);
/* daily component correction*/
//void Rotation_daily(std::deque < SacRec > & SACV, std::vector < std::string > Rout);

extern MyLogger logger;
MyLogger logger;
extern MEMO memo;
MEMO memo;



int main(int argc, char *argv[])
{
    if(argc!=2 && argc!=3)
    {
        std::cerr<<"Usage: "<<argv[0]<<" [Parameters_file] [monthdir](optional, e.g. 2019.JAN)"<<std::endl;
        return -1;
    }
    bool SyncNorm = true;
    if (argc == 2)
        logger.Rename("logs_Seed2Cor/Test");
    else
    {
        std::string name_logger(argv[2]);
        name_logger = "logs_Seed2Cor/Test_" + name_logger;
        logger.Rename(name_logger);
    }
    try
    {
        /* Initialize the CC Database with the input parameter file */
        CCDatabase cdb(argv[1]);
        const CCPARAM& cdbParams = cdb.GetParams();
        /* check total memory available */
        float MemTotal = memo.MemTotal();
        logger.Hold( INFO, "Estimated total memory = "+std::to_string(MemTotal)+" Mb", FuncName );
        logger.flush();
        /* iterate through the database and handle all possible events */
        #pragma omp parallel
        {
            // parallel region S
            while(1)   // main loop
            {
                //int ithread = omp_get_thread_num();
                /* dynamically assign events to threads, one at a time */
                bool got;
                std::vector<DailyInfo> dinfoV;
                #pragma omp critical(cdb)
                {
                    // critical S
                    got     = cdb.GetRec_AllCH(dinfoV);
                    cdb.NextEvent();
                }   // critical E
                if( !got ) break;
                if( cdbParams.fskipesac==2 && cdbParams.fskipresp==2 && cdbParams.fskipamph==2 )
                    break;
                try   // handle current event
                {
                    std::deque<SacRec> sacV;
                    std::vector < std::string > sacfout_R;
                    std::vector < std::string > sacfout_O;
                    std::vector<std::stringstream> reportV( dinfoV.size() );
                    
                    /*----  STEP 1: seed to SAC file(removed response) ----*/
                    for( int ich=0; ich < dinfoV.size(); ich++ )
                    {
                        auto& dinfo     = dinfoV[ich];
                        bool extract_flag=false;
                        /* daily info from the database */
                        //logger.Hold( INFO, dinfo.seedname + " " + dinfo.staname + " " + dinfo.chname, FuncName );
                        /* stringstream for reporting */
                        auto& report    = reportV[ich];
                        
                        /* extract the original sac from seed */
                        SacRec sac( report );
                        if ( cdbParams.fskipesac == 0 || ( ! FileExists(dinfo.osac_outname) && cdbParams.fskipesac==1 ) ) 
                        {
                            float gapfrac;
                            sac.SetMaxMemForParallel( MemTotal * dinfo.memomax * 0.8 / omp_get_num_threads() );
                            SeedRec seedcur( dinfo.seedname, dinfo.rdsexe, report );
                            if( seedcur.ExtractSac( dinfo.staname, dinfo.netname, dinfo.chname, dinfo.sps, dinfo.rec_outname,
                                                    dinfo.resp_outname, gapfrac, sac ) )
                            {
                                extract_flag    = true;
                                sac.Write( dinfo.osac_outname );
                                ////logger.Hold( INFO,"Extracting SAC file: "+dinfo.osac_outname, FuncName );
                                sacfout_O.push_back(dinfo.osac_outname);
                            }
                        }

                        /* remove response and cut */
                        // if do not do removing response step
                        if (  cdbParams.fskipresp==2 || (FileExists(dinfo.fsac_outname) && cdbParams.fskipresp==1)  ) 
                        {
                            
                            if (FileExists(dinfo.fsac_outname))
                            {
                                sac.Load(dinfo.fsac_outname);
                                sacfout_R.push_back(dinfo.fsac_outname);
                            }
                            sacfout_O.push_back(dinfo.osac_outname);
                            sacV.push_back( std::move(sac) );
                            continue;
                        }
                        
                        //if extraction step is skipped
                        if (extract_flag==false) 
                        {
                            // Get the path to resp file
                            std::vector<std::string> resp_list;
                            if (List( (dinfo.outdir).c_str(), ("RESP*."+dinfo.netname+"."+dinfo.staname+".*."+dinfo.chname).c_str(), 2, resp_list))
                                dinfo.resp_outname  = resp_list[0];
                            else
                            {
                                // if resp file NOT found, check if ft_*SAC file exists, continue then
                                if (FileExists(dinfo.fsac_outname))
                                {
                                    sac.Load(dinfo.fsac_outname);
                                    sacfout_R.push_back(dinfo.fsac_outname);
                                }
                                sacfout_O.push_back(dinfo.osac_outname);
                                sacV.push_back( std::move(sac) );
                                continue;
                            }
                            
                            if (FileExists(dinfo.osac_outname))
                            {
                                ////logger.Hold( INFO,"Reading existing SAC file: "+dinfo.osac_outname, FuncName );
                                sac.Load(dinfo.osac_outname);
                                sacfout_O.push_back(dinfo.osac_outname);
                            }
                            else
                            {
                                // if original SAC file NOT found, check if ft_*SAC file exists, continue then
                                if (FileExists(dinfo.fsac_outname))
                                {
                                    sac.Load(dinfo.fsac_outname);
                                    sacfout_R.push_back(dinfo.fsac_outname);
                                }
                                sacfout_O.push_back(dinfo.osac_outname);
                                sacV.push_back( std::move(sac) );
                                continue;
                            }
                        }
                        
                        if ( sac.shd.npts <=1 ) { 
                                std::cout<<"ATTENTION: Skip Preprocessing for: "<<dinfo.fsac_outname<<std::endl;
                                continue;
                        }
                        
                        //  remove response
                        sac.RmRESP( dinfo.resp_outname, dinfo.perl*0.8, dinfo.perh*1.3, dinfo.evrexe );
                        char evtime[15];
                        sprintf( evtime, "%04d%02d%02d000000\0", dinfo.year, dinfo.month, dinfo.day );
                        sac.ZoomToEvent( evtime, -12345., -12345., dinfo.t1, dinfo.tlen );
                        sac.Write( dinfo.fsac_outname );
                        sacV.push_back( std::move(sac) );
                        sacfout_R.push_back(dinfo.fsac_outname);
                    }
                    
                    /*--- Component correction ---*/
                    /*float del=3.0;
                    if (sacV.size()==dinfoV.size() && sacV.size()> 1 &&sacV.size()<4)
                    {
                        // Do not take do esac, not do rsac, do amp into consideration!!!
                        if ( abs(sacV[0].shd.cmpaz-90.0 ) > del && sacV[0].sig && sacV[1].sig)
                            Rotation_daily(sacV, sacfout_R);
                    } */
                    
                    /*----  STEP 2: convert SAC to amp&ph file ----*/
                    if ( cdbParams.fskipamph == 2)
                        continue;
                    if (cdbParams.fskipamph == 1)
                    {
                        bool skip_amph_flag  = false;
                        for ( int ich=0; ich < dinfoV.size(); ich++ )
                        {
                            if ( FileExists(dinfoV[ich].fsac_outname+".am") == false)
                            {
                                skip_amph_flag  = true;
                                break;
                            }
                            if ( FileExists(dinfoV[ich].fsac_outname+".ph") == false)
                            {
                                skip_amph_flag  = true;
                                break;
                            }
                        }
                        if (skip_amph_flag)
                        continue;
                    }
                    
                    /* time-domain normalization */
                    TNormAll( sacV, dinfoV, SyncNorm );
                    /* fre-domain normalization */
                    // convert sacs to am&ph and store amp in sacV
                    for( int isac=0; isac<sacV.size(); isac++ )
                    {
                        auto& dinfo = dinfoV[isac];
                        auto& sac   = sacV[isac];
                        if( ! sac.sig ) continue;
                        SacRec sac_am, sac_ph;
                        sac.ToAmPh( sac_am, sac_ph );
                        sac         = std::move(sac_am); 
                        sac_ph.Write( dinfo.fsac_outname + ".ph" );
                    }
                    // normalize
                    FNormAll( sacV, dinfoV, SyncNorm );
                    // write am
                    for( int isac=0; isac<sacV.size(); isac++ )
                    {
                        auto& dinfo = dinfoV[isac];
                        auto& sac = sacV[isac];
                        if( ! sac.sig ) continue;
                        sac.Write( dinfo.fsac_outname + ".am" );
                    }
                    
                    // Remove raw or removed-resp SAC files
                    if (cdbParams.fdelosac==1)
                        for (int isac=0; isac<sacfout_O.size(); isac++ )
                            if (FileExists(sacfout_O[isac]))
                                fRemove(sacfout_O[isac].c_str());
                    if (cdbParams.fdelosac==2)
                        for (int isac=0; isac<sacfout_R.size(); isac++ )
                            if (FileExists(sacfout_R[isac]))
                                fRemove(sacfout_R[isac].c_str());
                    if (cdbParams.fdelosac==3)
                    {
                        for (int isac=0; isac<sacfout_R.size(); isac++ )
                            if (FileExists(sacfout_R[isac]))
                                fRemove(sacfout_R[isac].c_str()); 
                        for (int isac=0; isac<sacfout_O.size(); isac++ )
                            if (FileExists(sacfout_O[isac]))
                                fRemove(sacfout_O[isac].c_str());
                    }
                    /* log if any warning */
                    for( const auto& report : reportV )
                    {
                        std::string warning = report.str();
                        if( ! warning.empty() )
                            logger.Hold( WARNING, "\n" + warning, FuncName );
                        logger.flush();
                    } // for dinfo
                
                std::string mesg    = "Preprocess done ::: "+ std::to_string(dinfoV[0].year)+"-"+std::to_string(dinfoV[0].month)+"-"+std::to_string(dinfoV[0].day)
                                +" "+ dinfoV[0].netname +"."+dinfoV[0].staname;
                logger.Hold( INFO, mesg, FuncName );
                }
                catch ( std::exception& e )
                {
                    logger.Hold( ERROR, e.what(), FuncName );
                } // current event done
            } // main while loop
        } // parallel region E
        
        /*---- CC code start here ----*/
        clock_t time_before;
        time_before = clock();
        if(cdbParams.fskipcrco == 3) return 0;
        int fskipcc = cdbParams.fskipcrco;
        //------------------------
        // Get channel list
        //------------------------
        std::vector < std::string > channel;
        Channellist temp_chlst  = cdb.GetchList();
        for( temp_chlst.Rewind(); !temp_chlst.IsEnded(); temp_chlst.NextRec() )
            channel.push_back(*(temp_chlst.GetRec()));
        std::string CH_pre; // channel prefix
        CH_pre.append(channel[0].begin(),channel[0].end()-1);
        //------------------------
        // Get station list
        //------------------------
        std::vector < StaInfo > stationlist;
        Stationlist temp_stalst = cdb.GetstaList();
        for( temp_stalst.Rewind(); !temp_stalst.IsEnded(); temp_stalst.NextRec() )
            stationlist.push_back( *(temp_stalst.GetRec()) );
        ////for( int ista=0; ista < stationlist.size(); ista++ )
        ////    std::cout << stationlist[ista].netcode<<"."<<stationlist[ista].name<<std::endl;
        //------------------------
        // Get month list
        //------------------------
        std::vector < std::string > monthdir;
        const std::vector<std::string> month_dict
        {
            "INVALID",
            "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
            "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
        };
        Seedlist temp_seedlst   = cdb.GetSeedList();
        bool new_month  = true;
        int temp_month  = -1;
        int temp_year   = -1;
        for( temp_seedlst.Rewind(); !temp_seedlst.IsEnded(); temp_seedlst.NextRec() )
        {
            SeedInfo temp_seedinfo  = *(temp_seedlst.GetRec());
            if (temp_year != temp_seedinfo.year || temp_month != temp_seedinfo.month)
                new_month   = true;
            if (new_month)
            {
                monthdir.push_back( std::to_string(temp_seedinfo.year)+"."+month_dict[temp_seedinfo.month]);
                new_month   = false;
                temp_year   = temp_seedinfo.year;
                temp_month  = temp_seedinfo.month;
            }
        }
        bool in_monlst      = false;
        if (argc == 3)
        {
            for( int imon=0; imon < monthdir.size(); imon++ )
                if (monthdir[imon] == argv[2])
                    in_monlst   = true;
            if (! in_monlst)
                std::cout << argv[2]<<" NOT in the month list"<<std::endl;
        }
        //------------------------------------------------
        // perform cross-correlation month by month
        //----------------------------------------------
        for( int imon=0; imon < monthdir.size(); imon++ )
        {
            if (in_monlst && monthdir[imon] != argv[2])
                continue;
            std::vector< CC_todo > CC_todolist;
            GCCtodoList( stationlist, monthdir[imon], CC_todolist, channel, fskipcc );
            std::cout <<"===== "<<monthdir[imon]<<":"<<CC_todolist.size()<<" pairs"<<std::endl;
            //for ( int icc=0; icc < CC_todolist.size(); icc++ )
            //    std::cout << CC_todolist[icc].station1.netcode+"."+CC_todolist[icc].station1.name+"_"
            //            + CC_todolist[icc].station2.netcode+"."+CC_todolist[icc].station2.name+" ";
            //std::cout<<std::endl;
            CCList2CC( CC_todolist, channel, cdbParams );
        }
        /*---- CC code end here ----*/
        logger.Hold( INFO, "All threads finished.", FuncName );
        std::cout<<"Elapsed Time: "<<(float(clock()-time_before))/CLOCKS_PER_SEC<<" secs"<<std::endl;
    }
    catch ( std::exception& e )
    {
        logger.Hold(FATAL, e.what(), FuncName);
        return -2;
    }
    catch (...)
    {
        logger.Hold(FATAL, "unknown exception", FuncName);
        return -2;
    }
    return 0;
}
