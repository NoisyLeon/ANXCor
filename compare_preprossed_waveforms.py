import obspy
import os
import glob
import numpy as np

datadir1    = '/work2/leon/benchmark_old'
datadir2    = '/work2/leon/benchmark_new'

subdir_lst  = glob.glob(datadir1+'/20*/*')

for subdir1 in subdir_lst:
    subdir2     = datadir2+'/'+subdir1.split('/')[-2]+'/'+subdir1.split('/')[-1]
    if not os.path.isdir(subdir2):
        print subdir2 + ' NOT exist!'
        continue
    
    ftsaclst2   = glob.glob(subdir2+'/'+'ft_*SAC')
    
    for ftsac2 in ftsaclst2:
        ftrec2  = ftsac2+'_rec'
        tr2     = obspy.read(ftsac2)[0]
        network = tr2.stats.network
        try:
            ftsac1  = glob.glob(subdir1+'/ft_*'+tr2.stats.station+'.'+tr2.stats.channel+'.SAC')[0]
            # ftsac1  = glob.glob(subdir1+'/ft_*'+network+'.'+tr2.stats.station+'.'+tr2.stats.channel+'.SAC')[0]
            tr1     = obspy.read(ftsac1)[0]
        except IndexError:
            print subdir1+'/ft_*'+tr2.stats.station+'.'+tr2.stats.channel+'.SAC'+' NOT exist!'
            continue
        if not np.sum(abs(tr1.data - tr2.data)) < tr1.data.max()/100.:
            # try:
            #     if network == 'CN':
            #         ftsac1  = glob.glob(subdir1+'/ft_*XL.'+tr2.stats.station+'.'+tr2.stats.channel+'.SAC')[0]
            #     else:
            #         ftsac1  = glob.glob(subdir1+'/ft_*CN.'+tr2.stats.station+'.'+tr2.stats.channel+'.SAC')[0]
            #     tr1     = obspy.read(ftsac1)[0]
            #     if not np.sum(abs(tr1.data - tr2.data)) < tr1.data.max()/100.:
            #         print ftsac1, ftsac2, ' tr NOT SAME!'
            # except:
            print ftsac1, ftsac2, ' tr NOT SAME!'
        # compare rec files
        if os.path.isfile(ftrec2):
            ftrec1  = ftsac1+'_rec'
            
            if not os.path.isfile(ftrec1):
                print ftrec1+ ' NOT exist!'
            else:
                inarr1  = np.loadtxt(ftrec1)
                inarr2  = np.loadtxt(ftrec2)
                if not np.allclose(inarr1, inarr2):
                    print ftrec1, ftrec2, ' NOT SAME!'
            # print ftrec2, ftrec1
        else:
            try:
                tram1   = obspy.read(ftsac1 + '.am')[0]
            except:
                print ftsac1 + '.am'+' NOT exist!'
                
            try:
                trph1   = obspy.read(ftsac1 + '.ph')[0]
            except:
                print ftsac1 + '.ph'+' NOT exist!'
                
            try:
                tram2   = obspy.read(ftsac2 + '.am')[0]
            except:
                print ftsac2 + '.am'+' NOT exist!'
                
            try:
                trph2   = obspy.read(ftsac2 + '.ph')[0]
            except:
                print ftsac2 + '.ph'+' NOT exist!'
                
            try:
                if not np.allclose(tram1.data, tram2.data):
                    print ftsac1, ftsac2, ' amp NOT SAME!'
            except:
                pass
            
            try:
                if not np.allclose(trph1.data, trph2.data):
                    print ftsac1, ftsac2, ' ph NOT SAME!'
            except:
                pass
        