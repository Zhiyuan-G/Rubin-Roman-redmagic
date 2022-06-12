import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from tqdm import tqdm
os.chdir('gcr-catalogs-master/')
import GCRCatalogs

print('Loading cosmoDC2_v1.1.4_redmagic_v0.8.1_highdens')
# Load CosmoDC2 Redmagic catalog
gc = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_redmagic_v0.8.1_highdens')
data = gc.get_quantities([ 'id','ra', 'dec','redshift','zspec'])


# Apply Rubin-Roman joint imaging footprint mask to CosmoDC2 redmagic catalog
joint_area_mask = (data['ra'] > 51) & (data['ra'] < 56) & (data['dec'] > -42) & (data['dec'] < -38)

print('Loading Roman truth fits file.')
# Load Roman sim truth catalog
Roman_truth_file = 'dc2_truth_gal_icrs.fits'
hdu = fits.open(Roman_truth_file)
roman_data = hdu[1].data
gind = np.array(roman_data['gind'])
# Narrow gal_id range for shorter computing time
Roman_gal_id = gind

print('Start storing joint galaxy info.')
# Make joint galaxy list: Store galaxy_id and indx in Roman truth catalog
Joint_gal_id = list()
Roman_joint_gal_idx = list()

# Loop over CosmoDC2 Redmagic
for ids in tqdm(data['id'][joint_area_mask]):
    # Cross match check
    if ids in gind:
        # Store info
        Joint_gal_id.append(ids)
        Roman_joint_gal_idx.append(np.where(Roman_gal_id == ids)[0][0])

print('Saving file..')
hf = h5py.File('../data/Rubin_Roman_gal_cross_match.h5', 'w')
hf.create_dataset('gal_id', data=np.array(Joint_gal_id))
hf.create_dataset('roman_id_idx', data=np.array(Roman_joint_gal_idx))
hf.close()