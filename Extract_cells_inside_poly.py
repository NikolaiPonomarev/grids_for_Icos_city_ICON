from copyreg import clear_extension_cache
from netCDF4 import Dataset
from netCDF4 import Dataset
import numpy as np
from math import fsum
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt

#Read file with borders
with open ('/users/nponomar/grids/Zurich_borders.txt') as f:
      coords = [(float(z[1:-1].split(',')[0]), float(z[1:-1].split(',')[1])) for z in next(f).split(' ')]

poly = Polygon(coords)

print('Read polyshape: ', np.shape(coords))

#Read nc grid file
fileN='/users/nponomar/grids/icon_Zurich_1_DOM02.nc'
dataset = Dataset(fileN, maskandscale=True, mmap=False)

lats = dataset.variables['clat'][:]
lons = dataset.variables['clon'][:]

clon = np.rad2deg( lons )
clat = np.rad2deg( lats )
print('Centers shape, lon, lat', np.shape(clon), np.shape(clat), clon, clat)

#Search for indecies inside poly

#clon_poly = np.zeros(np.shape(clon), dtype=np.float64)
#clat_poly = np.zeros(np.shape(clat), dtype=np.float64)
indcs_poly=[]

for ix in range(0, len(clon)):
    if Point((clon[ix], clat[ix])).intersects(poly):
            indcs_poly.append(ix)

print('Centers inside poly shape, indecies', np.shape(indcs_poly), indcs_poly)

#Extract cells 

clon_poly = clon[indcs_poly]
clat_poly = clat[indcs_poly]

figure, ax  = plt.subplots()
ax.plot(clon, clat, 'o', color='black')
ax.plot(np.asarray(coords)[:, 0], np.asarray(coords)[:, 1], 'o', color='black')
plt.savefig('check.png')
