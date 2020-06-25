import numpy as np
import matplotlib as mpl
from astropy.table import Table, Column, join 
#mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
import sys
import math
import os
import glob
import sys
from sortedcontainers import SortedDict
import datetime as dt
import imageio
import os
from PIL import Image
from matplotlib.colors import LogNorm
from astropy.nddata.utils import Cutout2D
from astropy import units as u
import datetime as dt 
import astropy.units as u
from astroML.crossmatch import crossmatch_angular 
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import glob
import csv

field = 'Dusty10'
year = '2018'
month = '06'
ran = '500'

def getstars(filename):
    ra = [] ; dec = [] ; par = []; epar = []
    bprp = [] ; gmag = [] ; gaiaid = []; b = []
    df = open(filename,"r")
    for line in df:
        if len(line)>0:
            line_split = line.split()
            if len(line_split)>=25:
                ra.append(float(line_split[0]))
                #print(line_split[4])
                dec.append(float(line_split[2]))
                #gaiaid.append(int(line_split[2]))
                        #par.append(float(line_split[5]))
                        #epar.append(float(line_split[6]))
                        #gmag.append(float(line_split[14]))
                        #bprp.append(float(line_split[24]))
                        #gaiaid.append(line_split[4])
                
    ra = np.array(ra)
    dec = np.array(dec)
    #gaiaid = np.array(gaiaid,dtype=str)

    return ra,dec


filename1 = '/home/swebb/oz100/NOAO_archive/archive_NOAO_data/scripts/create_lc/total_field_sources/fields/'+ year + '_'+month+'_'+ field + '_full_field_sources.ascii'
filename2 = '/home/swebb/oz100/NOAO_archive/archive_NOAO_data/scripts/GAIA_volume_searchs/'+ran+'pc/'+field+'field.'+ran+'pc.sample'

ra1, dec1 = np.loadtxt(filename1, skiprows=1, unpack = True)
ra2, dec2 = getstars(filename2)
#print(ra)
#ra1 = []
#for i in ra:
	#ra1.append(i/15)
first_cat = np.empty((len(ra1), 2), dtype = np.float64)
first_cat[:, 0] = ra1
first_cat[:, 1] = dec1

new_cat = np.empty((int(len(ra2)), 2), dtype = np.float64)
new_cat[:,0] = ra2
new_cat[:,1] = dec2

print(print(new_cat))
print(print(first_cat))
max_radius = 1/3600 #.5 arc second

dist_between, ind_row = crossmatch_angular(first_cat, new_cat, max_radius)
match =~np.isinf(dist_between)

match_table = Table()
match_table['matched_true_false'] = match 
match_table['match_ID'] = ind_row
match_table['matched_firstcat_RA'] = first_cat[:,0]
match_table['matched_firstcat_DEC'] = first_cat[:,1] 



new_cat_line_ID = range(len(match))
#print('LENGTH OF MATCH: ' + str(len(match)))
match_table['match_table_false_id'] = new_cat_line_ID 

new_cat_match_true = []
new_cat_match_false = []
new_cat_ID_true = []
new_cat_ID_false = []
first_cat_RA_true = []
first_cat_DEC_true = []


first_cat_RA_false= []
first_cat_DEC_false = []

for row in match_table: 
	if row['matched_true_false'] == True:
 		new_cat_match_true.append(row['matched_true_false'])
 		new_cat_ID_true.append(row['match_ID'])
 		first_cat_RA_true.append(row['matched_firstcat_RA'])
 		first_cat_DEC_true.append(row['matched_firstcat_DEC'])
				
for row in match_table: 
	if row['matched_true_false'] == False:
		new_cat_match_false.append(row['matched_true_false'])
		new_cat_ID_false.append(row['match_table_false_id'])
		first_cat_RA_false.append(row['matched_firstcat_RA'])
		first_cat_DEC_false.append(row['matched_firstcat_DEC'])


#print('new_cat_match_false:' + str(len(new_cat_match_false))) 

new_cat_RA_true  = []
new_cat_DEC_true  = []
new_cat_gmag_auto_true  = []
new_cat_gmag_auto_err_true  = []
new_cat_gmag_aper_true  = []
new_cat_gmag_aper_err_true  = []
new_cat_mjd_true = []

for j in new_cat_ID_true: 
	new_cat_RA_true.append(new_cat[j,0])
	new_cat_DEC_true.append(new_cat[j,1])



new_cat_RA_false = []
new_cat_DEC_false  = []
new_cat_gmag_auto_false  = []
new_cat_gmag_auto_err_false  = []
new_cat_gmag_aper_false  = []
new_cat_gmag_aper_err_false  = []
new_cat_mjd_false = []

for k in new_cat_ID_false: 
	new_cat_RA_false.append(first_cat[k,0])
	new_cat_DEC_false.append(first_cat[k,1])

match_true_table = Table () 
match_true_table['RA'] = new_cat_RA_true
match_true_table['DEC'] = new_cat_DEC_true

print(match_true_table)
match_false_table = Table()
match_false_table['RA'] = new_cat_RA_false
match_false_table['DEC'] = new_cat_DEC_false

output = year + '_' + month+ '_'+field+'_'+ran+'_gaia_matched.ascii'

match_true_table.write(output, format='ascii', overwrite=True)
#print('TRUE: ' + str(len(match_true_table)))
#print('FALSE: ' + str(len(match_false_table)))
