import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import h5py
import multiprocessing
import time
import get_template

Finder = get_template.Finder()

Roman_joint_hf = h5py.File('data/Rubin_Roman_gal_cross_match.h5', 'r')
joint_gal_id = Roman_joint_hf['gal_id'][:]
Roman_gal_idx = Roman_joint_hf['roman_id_idx'][:]

def worker(return_dict,gal_id,gal_idx,tag):
    Roman_mag_list = []
    Rubin_mag_list = []
    name_list = []
    cross_idx = []
    for i in range(len(gal_id)):
        Roman_mags = Finder.calc_roman_color(gal_id=gal_id[i], roman_gal_idx=gal_idx[i],check = False)
        Rubin_mags = Finder.calc_color(gal_id=gal_id[i],colors = ['u','g','r','i','z','y'])
        Roman_mag_list.append(Roman_mags)
        Rubin_mag_list.append(Rubin_mags)
        name_list.append(gal_id)
        cross_idx.append(gal_idx)
        
        
    return_dict[tag] = (Roman_mag_list,Rubin_mag_list,name_list)
    
    
if __name__ == "__main__":
    start_time = time.perf_counter()
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    for i_start in range(0,10,10):
        s = slice(i_start,i_start + 10)
    #for i in range(64*8):
        p = multiprocessing.Process(target=worker, 
                                    args=(return_dict,joint_gal_id[s],
                                        Roman_gal_idx[s],i_start))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()
    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time-start_time} seconds")
    
    Roman_mags = list()
    Rubin_mags = list()
    Gal_id = list()
    cross_id = list()
    for i_start in return_dict.keys():
        Roman_mags.extend(return_dict[i_start][0])
        Rubin_mags.extend(return_dict[i_start][1])
        Gal_id.extend(return_dict[i_start][2][0])
        cross_id.extend(return_dict[i_start][3][0])
        
    np.save('data/gal_id', np.array(Gal_id))
    np.save('data/Rubin_mags', np.array(Rubin_mags))
    np.save('data/roman_mags', np.array(Roman_mags))
    np.save('data/cross_id', np.array(cross_id))