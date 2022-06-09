import h5py
import numpy as np
import gzip
import matplotlib.pyplot as plt
from lsst.sims.photUtils import BandpassDict, Sed, Bandpass
from lsst.utils import getPackageDir
from lsst.sims.photUtils import cache_LSST_seds, getImsimFluxNorm
import os

class Finder():
    
    def __init__(self):
        
        
        self.lsst_bplist = np.array(['u','g','r','i','z','y'])
        self.file_s_e = np.load('/global/cscratch1/sd/zg64/Rubin-Roman-Redmagic/lookup_file_list_range.npy')
        self.sed_look_dir = "/global/projecta/projectdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sedLookup/"
        self._galaxy_sed_dir = os.path.join(getPackageDir('sims_sed_library'))
        self.lookupfile = os.listdir(self.sed_look_dir)
        self.lsst_bp_dict,self.dummy_bp_dict = BandpassDict.loadBandpassesFromFiles()
        
    def find_h5_file (self, gal_id):

        for i,f in enumerate(self.file_s_e):
            if gal_id > f[0] and gal_id < f[1]:
                if i > 57:
                    which_h5file = i + 1
                else:
                    which_h5file = i
                print('SED LookupFile: ', self.lookupfile[which_h5file])
        return(which_h5file)
    
    
    
    def find_target_sed(self,gal_id, plot = True, calc = True):
    
        h5f_id = self.find_h5_file(gal_id=gal_id)
        look_file = self.sed_look_dir + self.lookupfile[h5f_id]
        hf = h5py.File(look_file, 'r')
        galaxy_id = hf['galaxy_id'][:]
        bulge_sed = hf['bulge_sed'][:]
        disk_sed = hf['disk_sed'][:]
        sed_names = hf['sed_names'][:]

        idx = np.where(galaxy_id == gal_id)[0][0]

        target_bulge = sed_names[bulge_sed[idx]]
        target_bulge_sed = target_bulge.decode("utf-8") 

        target_disk = sed_names[disk_sed[idx]]
        target_disk_sed = target_disk.decode("utf-8") 
        bulge_spec = Sed()
        bulge_spec.readSED_flambda(os.path.join(self._galaxy_sed_dir, target_bulge_sed))
        disk_spec = Sed()
        disk_spec.readSED_flambda(os.path.join(self._galaxy_sed_dir, target_disk_sed))
        print('SED template (bulge): ',target_bulge_sed)
        print('SED template (disk): ',target_disk_sed)
        if plot == True:
            fig,ax = plt.subplots(1,1,figsize = (8,6))
            ax2 = ax.twinx()


            m = bulge_spec.wavelen < 2000
            ax.plot(bulge_spec.wavelen[m],bulge_spec.flambda[m],color = 'c', label = 'bulge')
            ax.legend(loc = 'upper right',fontsize = 14)
            ax.set_ylabel('bulge',fontsize = 16)

            ax2.plot(disk_spec.wavelen[m],disk_spec.flambda[m],color = 'm',label = 'disk') 
            ax2.legend(loc = 'upper left',fontsize = 14)
            ax2.set_ylabel('disk',fontsize = 16)
            ax.set_xlabel('Wavelength [nm]',fontsize = 16)
            
        if calc == True:
            return(hf,idx,bulge_spec,disk_spec)
        else:
            return(hf,idx)


    def calc_color (self,gal_id,colors):
        lsst_bplist = np.array(['u','g','r','i','z','y'])
        hfile,ID,bulge_spec,disk_spec = self.find_target_sed(gal_id=gal_id, plot = False)
        color_comb = list()
        for c in colors:
            color_comb.append( np.where(lsst_bplist == c)[0][0])

        mags = list()
        for c_idx in color_comb:

            bulge_s = Sed(wavelen =  bulge_spec.wavelen , flambda = bulge_spec.flambda)
            disk_s = Sed(wavelen =  disk_spec.wavelen , flambda = disk_spec.flambda)
            disk_magnorm = hfile['disk_magnorm'][c_idx][ID]
            bulge_magnorm = hfile['bulge_magnorm'][c_idx][ID]

            # disk
            fnorm_disk = getImsimFluxNorm(disk_s, disk_magnorm)
            disk_s.multiplyFluxNorm(fnorm_disk)
            ax, bx = disk_s.setupCCM_ab()
            disk_s.addDust(ax, bx, A_v=hfile['disk_av'][ID], R_v=hfile['disk_rv'][ID])
            # bulge
            fnorm_bulge = getImsimFluxNorm(bulge_s, bulge_magnorm)
            bulge_s.multiplyFluxNorm(fnorm_bulge)
            ax, bx = bulge_s.setupCCM_ab()
            bulge_s.addDust(ax, bx, A_v=hfile['bulge_av'][ID], R_v=hfile['bulge_rv'][ID])


            # Run for various redshift and calculate colors
            spec_f = Sed(wavelen =  bulge_s.wavelen , flambda = bulge_s.flambda+disk_s.flambda)
            spec_f.redshiftSED(hfile['redshift'][ID],dimming = True)

            inband_mag = spec_f.calcMag(self.lsst_bp_dict[lsst_bplist[c_idx]])
            mags.append(inband_mag)
        return(mags)

    def calc_colors (self,gal_id,file_save = True):
        hfile,ID,bulge_s,disk_s = self.find_target_sed(gal_id=gal_id, plot = False)
        disk_magnorm = hfile['disk_magnorm'][0][ID]
        bulge_magnorm = hfile['bulge_magnorm'][0][ID]
        redshift = hfile['redshift'][ID]

        fnorm = getImsimFluxNorm(disk_s, disk_magnorm)
        disk_s.multiplyFluxNorm(fnorm)
        ax, bx = disk_s.setupCCM_ab()
        disk_s.addDust(ax, bx, A_v=hfile['disk_av'][ID], R_v=hfile['disk_rv'][ID])
        disk_s.redshiftSED(hfile['redshift'][ID], dimming=True)

        fnorm = getImsimFluxNorm(bulge_s, bulge_magnorm)
        bulge_s.multiplyFluxNorm(fnorm)
        ax, bx = bulge_s.setupCCM_ab()
        bulge_s.addDust(ax, bx, A_v=hfile['bulge_av'][ID], R_v=hfile['bulge_rv'][ID])


        # Run for various redshift and calculate colors
        spec_f = Sed(wavelen =  bulge_s.wavelen , flambda = bulge_s.flambda+disk_s.flambda)
        redshift =  np.arange(0,1.2,0.01)
        color = list()
        mags = list()
        for z in redshift:
            spec_f = Sed(wavelen =  bulge_s.wavelen , flambda = bulge_s.flambda+disk_s.flambda)
            spec_f.redshiftSED(z,dimming = True)
            maglist = lsst_bp_dict.magListForSed(spec_f)
            color.append(maglist[:-1]-maglist[1:])
            mags.append(maglist)
        if file_save == True:
            np.save('data/single_gal_cali_u_%i'%gal_id,np.array(color))
        else:   
            return(np.array(mags),bulge_s.wavelen,bulge_s.flambda+disk_s.flambda)