import os
import obspy
import glob
datadir     = '/work2/leon/temp_working_2004_2008'
start_date  ='20040401'
end_date    ='20080531'
stime       = obspy.UTCDateTime(start_date)
etime       = obspy.UTCDateTime(end_date)
monthdict   = {1: 'JAN', 2: 'FEB', 3: 'MAR', 4: 'APR', 5: 'MAY', 6: 'JUN', 7: 'JUL', 8: 'AUG', 9: 'SEP', 10: 'OCT', 11: 'NOV', 12: 'DEC'}
month       = stime.month - 1
while(stime < etime):
    # print stime
    if month < stime.month:
        month   += 1
        print 'Deleting data: '+str(stime.year)+'.'+monthdict[stime.month]
    dailydir    = datadir+'/'+str(stime.year)+'.'+monthdict[stime.month]+'/'+str(stime.year)+'.'+monthdict[stime.month]+'.'+str(stime.day)
    if os.path.isdir(dailydir):
        ftamsaclst  = glob.glob(dailydir+'/ft*SAC.am')
        ftphsaclst  = glob.glob(dailydir+'/ft*SAC.ph')
        ftrecsaclst = glob.glob(dailydir+'/ft*SAC_rec')
        for phfname in ftphsaclst:
            os.remove(phfname)
        for amfname in ftamsaclst:
            os.remove(amfname)
        for recfname in ftrecsaclst:
            os.remove(recfname)
    else:
        print dailydir+' NOT exists!'
    if stime.month == 12 and stime.day == 31:
        stime   = obspy.UTCDateTime(str(stime.year + 1)+'0101')
        month   = 0
    else:
        stime.julday += 1
