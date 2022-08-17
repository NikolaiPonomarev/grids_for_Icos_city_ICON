import netCDF4
from netCDF4 import Dataset
import numpy as np
from math import fsum
import time
start_time = time.time()
from pylab import rcParams
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
font = {'family' : 'normal',
        'size'   : 6}

plt.rc('font', **font)
widthP=int(12/2.54)
heightP=int(10/2.54)
rcParams['figure.figsize'] = (widthP, heightP)
rcParams.update({'figure.autolayout': True})
rcParams.update({'figure.autolayout': True})

#Read file with borders
with open ('/home/yury/SILAM/Apps571/input/CAMS-REG-AP_v4_2a-MSK_CB5V2/london_city') as f:
      coords = [(float(z[1:-1].split(',')[0]), float(z[1:-1].split(',')[1])) for z in next(f).split(' ')]
lclon = (-0.1138211-0.0727493)/2       
lclat = (51.5068696+51.5233122)/2

poly = Polygon(coords)

print('Read polyshape: ', np.shape(coords))

dicts=[]
#n species = 28
sp = ('NOX_CAMS', 'CO_CAMS','PM2_5_10_CAMS', 'PM2_5_CAMS', 'SO2_CAMS', 'ALD2_CAMS', 'ALDX_CAMS', 'BENZENE_CAMS', 'CH3Cl_CAMS',  'EC_COARSE_CAMS', 'EC_FINE_CAMS', 'EC_FINEff_CAMS', 'EC_FINEwb_CAMS', 'ETHA_CAMS', 'ETH_CAMS', 'ETOH_CAMS', 'HCHO_CAMS', 'IOLE_CAMS', 'MINERAL_COARSE_CAMS', 'MINERAL_FINE_CAMS', 'NH3_CAMS', 'NMVOC_CAMS', 'OC_COARSE_CAMS', 'OC_FINE_CAMS', 'OLE5_CAMS', 'PAR5_CAMS', 'TOLUENE_CAMS', 'XYLENE_CAMS')
#sp = ('ALD2_CAMS', 'ALDX_CAMS', 'BENZENE_CAMS', 'CH3Cl_CAMS', 'CO_CAMS', 'EC_COARSE_CAMS', 'EC_FINE_CAMS', 'EC_FINEff_CAMS')
for i, s in enumerate(sp):
 #Paths to the original emission files
 #"/home/yury/SILAM/Apps/input/emis/CAMS-Lockdown/CO_CAMS-REG-AP_v4_2.sa2.hdr"
 NCfileEms = '/home/yury/SILAM/Apps/input/emis/CAMS-Lockdown/'+ str(sp[i]) + '-REG-AP_v4_2.sa2.nc4'
 HDRfile = '/home/yury/SILAM/Apps/input/emis/CAMS-Lockdown/'+ str(sp[i]) + '-REG-AP_v4_2.sa2.hdr'
 HDRname = str(sp[i]) + '-REG-AP_v4_2.sa2.hdr'
 #Read the HDR and create a dictionary (with the space as a separator for lines.split())
 with open(HDRfile) as HDR:
   lineslist=[]
   for l in HDR:
    lineslist.append(l.split())
 dictHDR={}   
 for i in range (0, len(lineslist)):
   if 'AREA_SOURCE_2' in lineslist[i]:
    k=0
    while 'END_AREA_SOURCE_2' not in lineslist[i+k]:
     k=k+1
 #    print('k', k, 'i+k', i+k, 'len', len(lineslist))
     if len(lineslist[i+k])!= 0:
      if lineslist[i+k][0] not in dictHDR:
       dictHDR[lineslist[i+k][0]]=list()
      dictHDR[lineslist[i+k][0]].extend(lineslist[i+k][2:])
# print('header', dictHDR)
# print('Lengths', len(dictHDR['source_name']), len(dictHDR['hour_in_day_index']), len(dictHDR['hour_in_day_index'])/len(dictHDR['source_name']))
 

 fh = Dataset(NCfileEms, mode='r')
 lonsE = fh.variables['lon'][:]
 latsE = fh.variables['lat'][:]
 emis = fh.variables['value'][:]
 src_size = fh.variables['src_size'][:]
 source_name = fh.variables['source_name'][:]
 src_start_index = fh.variables['src_start_index'][:]
 fh.close()
 em_sc={}

 
 for isr, sr in enumerate(dictHDR['source_name']):
     em=0
     
     if '_GBR' in sr:
         #print('Full name ', sr, 'Short name ', sr[2:])
         sind = src_start_index[int(dictHDR['netcdf_data'][3*isr])]
         ssize = src_size[int(dictHDR['netcdf_data'][3*isr])]
         for j in range (sind, sind+ssize):
             #if Point((lonsE[j], latsE[j])).within(poly):
             if (abs(lonsE[j]-lclon)<0.05) and (abs(latsE[j]-lclat)<0.05):
                 em=em+emis[j]
         if sr[2:] in em_sc:       
             em_sc[sr[2:]]+=em
         else:
             em_sc[sr[2:]]=em            
         #print('name, number, start, size, i ', sr, int(dictHDR['netcdf_data'][3*isr]), src_start_index[int(dictHDR['netcdf_data'][3*isr])], src_size[int(dictHDR['netcdf_data'][3*isr])], 3*isr)
 #print(em_sc, sum(em_sc.values()))
 sum_em = sum(em_sc.values())
 if sum_em>0:
     em_sc.update({n: 100*em_sc[n]/sum_em for n in em_sc.keys()})
 dicts.append(em_sc)

colors = {
          'A' : 'red',
          'B' : 'blue',
          'C' : 'green',
          'D' : 'yellow',
          'E' : 'orange',
          'F' : 'black',
          #'F2' : 'olive',
          #'F3' : 'darkolivegreen',
          #'F4' : 'dimgrey',
          'G' : 'violet',
          'H' : 'grey',
          'I' : 'brown',
          'J' : 'pink',
          'K' : 'aqua',
          'L' : 'lightgoldenrodyellow',
          }

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), fontsize=4, loc='upper center', bbox_to_anchor=(0.5, -0.07), fancybox=True, shadow=True, ncol=5).set_zorder(102)


for t in range (0, len(dicts), 4):
 fig, ax = plt.subplots()
 spnames=[x[:-5] for x in sp[t:t+4]]
 #print(spnames)
 for i in range (0, 4):
     valsum=0
     vtr=0
     for k, v in dicts[t+i].items():
         if 'exhaust' not in k: 
            #print(k, dicts, len(k.split('_')), k.split('_')) 
            ax.bar(i+1, v, bottom=valsum, width = 0.5, color=colors[k.split('_')[1]], label = str(k.split('_')[0]), zorder=10)#label = str(k.split('_')[0])
            valsum+=v 
         else:
            vtr+=v        
            #ax.bar(i+1, vtr, bottom=valsum, width = 0.5, color=colors[k.split('_')[-3]], label = str(k.split('_')[0])+str(k.split('_')[2]), zorder=10)#label = str(k.split('_')[0]) 
             
     ax.bar(i+1, vtr, bottom=valsum, width = 0.5, color=colors['F'], label = 'RoadTransport', zorder=10)  
 ax.set_xlim([0, 5])
 ax.set_xticks(range(1, 5)) 
 print(np.shape(range(1, 5)), np.shape(spnames), spnames) 
 ax.set_xticklabels(spnames)
 #ax.legend(loc='upper right', fontsize=8)
 legend_without_duplicate_labels(ax)
 ax.grid(axis='y', linestyle='--', zorder=1)
 ax.set_title('City of London')
 ax.set_ylabel('Emission sector distribution, %')
 fig.tight_layout()
 name = 'CityLondon_Annual_Emission_sectors_distribution_'+str("_".join(spnames))+'.png'
 fig.savefig(name, dpi=500) 
 plt.close(fig)     
           