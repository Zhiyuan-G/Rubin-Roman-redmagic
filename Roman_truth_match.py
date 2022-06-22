# Code from Troxel with slight changes.

import os
from astropy.io import fits
import h5py
import numpy as np
import gzip
import glob
import fitsio as fio
from scipy.spatial import cKDTree
Roman_output_dir = '/global/cscratch1/sd/zg64/Rubin-Roman-Redmagic/dc2_sim_output/'

def get_match(m,t):
    tree = cKDTree(np.c_[t['ra'].ravel(), t['dec'].ravel()])
    dd,ii = tree.query(np.c_[m['alphawin_j2000'].ravel(), m['deltawin_j2000'].ravel()], k=3)
    m1=dd[:,0]*60.*60.<1
    m2=dd[:,1]*60.*60.<1
    m3=dd[:,2]*60.*60.<1
    dm1=(rtd['mag_F184'][ii[:,0]]-rd['mag_auto_F184'])
    dm2=(rtd['mag_F184'][ii[:,1]]-rd['mag_auto_F184'])
    dm3=(rtd['mag_F184'][ii[:,2]]-rd['mag_auto_F184'])
    mask = np.ones(len(rd))*-1
    mask[m1] = ii[m1,0]
    mask[m2&(np.abs(dm2)<np.abs(dm1))]= ii[m2&(np.abs(dm2)<np.abs(dm1)),1]
    mask[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2))]= ii[m3&(np.abs(dm3)<np.abs(dm1))&(np.abs(dm3)<np.abs(dm2)),2]
    return m[mask>=0],t[mask[mask>=0].astype(int)]


mask = []
start  = 0
rtd=None
for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz'))[5:-4]):
    print(i)
    try:
        tmp = fio.FITS(f)[-1].read()
        if rtd is None:
            print('test')
            rtd = np.zeros(100000000,dtype=tmp.dtype)
        for col in rtd.dtype.names:
            rtd[col][start:start+len(tmp)] = tmp[col]
        mask.append(i)
        start+=len(tmp)
    except:
        pass

rtd=rtd[rtd['x']>0]

start  = 0
rd = None
for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz'))[5:-4]):
    if i in mask:
        tmp = fio.FITS(f)[-1].read()
        if rd is None:
            rd = np.zeros(100000000,dtype=tmp.dtype)
        for col in rd.dtype.names:
            rd[col][start:start+len(tmp)] = tmp[col]
        mask.append(i)
        start+=len(tmp)

rd=rd[rd['number']>0]

rd=rd[(rd['flags']==0)&(rd['flux_auto']/rd['fluxerr_auto']>5)]
rd_match,rtd_match = get_match(rd,rtd)
fr = ['Y106','J129','H158','F184']
for i in range(4):
    rd_match['mag_auto_'+fr[i]] += rtd_match['dered_'+fr[i]]
    rtd_match['mag_'+fr[i]] += rtd_match['dered_'+fr[i]]
    mask = (rd_match['mag_auto_'+fr[i]]>17)&(rd_match['mag_auto_'+fr[i]]<20)&(rtd_match['gal_star']==1)
    rd_match['mag_auto_'+fr[i]] -= np.mean(rd_match['mag_auto_'+fr[i]][mask]-rtd_match['mag_'+fr[i]][mask])
    
# now write some data
fits = fio.FITS('m_t_match_rtd.fits','rw')
# create a new table extension and write the data
fits.write(rtd_match)
fits.close()

# now write some data
fits = fio.FITS('m_t_match_rt.fits','rw')
# create a new table extension and write the data
fits.write(rd_match)
fits.close()