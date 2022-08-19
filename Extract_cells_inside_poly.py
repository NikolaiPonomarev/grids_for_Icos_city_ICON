from copyreg import clear_extension_cache
from netCDF4 import Dataset
from netCDF4 import Dataset
import numpy as np
from math import fsum
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import xarray as xr
from cdo import *

def plot_icongrid(ds, ax):
      n_cells = len(dataset.variables["clon"])
      corners = np.zeros((n_cells, 3, 2))
      corners[:, :, 0] = dataset.variables["vlon"][dataset.variables["vertex_of_cell"] - 1].T
      corners[:, :, 1] = dataset.variables["vlat"][dataset.variables["vertex_of_cell"] - 1].T
      

      poly_coll = PolyCollection(
                  np.rad2deg(corners),
                  #cmap=cmap,
                  edgecolors="black",
                  linewidth=0.04,
                  antialiased=True, # AA will help show the line when set to true
                  alpha=0.6,
            )
      # Add the collection to the ax
      ax.add_collection(poly_coll)


#Read file with borders
with open ('/users/nponomar/grids/Zurich_borders.txt') as f:
      coords = [(float(z[1:-1].split(',')[0]), float(z[1:-1].split(',')[1])) for z in next(f).split(' ')]

poly = Polygon(coords)

print('Read polyshape: ', np.shape(coords))

#Read nc grid file
fileN='/users/nponomar/grids/icon_Zurich_1_DOM02.nc'
dataset = Dataset(fileN, maskandscale=True, mmap=False)
dataset = xr.load_dataset(fileN)

lats = dataset['clat'][:]
lons = dataset['clon'][:]

clon = np.rad2deg( lons )
clat = np.rad2deg( lats )
#print('Centers shape, lon, lat', np.shape(clon), np.shape(clat), clon, clat)

#print(dataset['vlat'] < dataset['vlat'].mean())
#dropped = dataset.where(dataset['vlat'] < dataset['vlat'].mean(), drop=True)
#Search for indecies inside poly
#dropped.to_netcdf("/users/nponomar/grids/icon_Zurich_1_DOM02_droppedcells_python.nc" )

#clon_poly = np.zeros(np.shape(clon), dtype=np.float64)
#clat_poly = np.zeros(np.shape(clat), dtype=np.float64)

indcs_poly=[]
x = np.unique(clon)
epsilon = (x[1] - x[0])
print('Unique Longitudes, epsilon', x, epsilon)
for ix in range(0, len(clon)):
    if Point((clon[ix], clat[ix])).distance(poly)<epsilon:
            indcs_poly.append(ix)

print('Centers inside poly shape, indecies', np.shape(indcs_poly), indcs_poly)

#Extract cells 
""""

ifile = fileN
ofile = "/users/nponomar/grids/icon_Zurich_1_DOM02_extractedcells_python.nc" 

cdo = Cdo()
inds = ','.join(str(e+1) for e in indcs_poly)
cdo.selgridcell(inds, input=ifile, output=ofile)

cdo.selgridcell()
"""
dataset_cells = dataset.isel(cell=indcs_poly) #extract vars which have cell dimension

inds_verticies = np.array([[x for x in y] for y in dataset.vertex_of_cell[:, indcs_poly].values-1]).flatten()
print('Extracted inds, inds_v', indcs_poly, inds_verticies)
inds_edges = np.array([[x for x in y] for y in dataset.edge_of_cell[:, indcs_poly].values-1]).flatten()

dataset_cells_verts = dataset_cells.isel(vertex=inds_verticies)
dataset_cells_verts_edges = dataset_cells_verts.isel(edge=inds_edges)

dataset_cells_verts_edges.to_netcdf("/users/nponomar/grids/icon_Zurich_1_DOM02_extracted_vars.nc" )

#Plot the grid and city's shape file
clon_poly = clon[indcs_poly]
clat_poly = clat[indcs_poly]

figure, ax  = plt.subplots()
ax.plot(clon, clat, 'o', color='black', markersize=2.5)
ax.plot(np.asarray(coords)[:, 0], np.asarray(coords)[:, 1], 'o', color='black', markersize=2.5)
plt.savefig('check.png')


#Read newly created nc grid file
fileN='/users/nponomar/grids/icon_Zurich_1_DOM02_extracted_vars.nc'
dataset = Dataset(fileN, maskandscale=True, mmap=False)
dataset = xr.load_dataset(fileN)
print(dataset)
clats = dataset.variables['clat'][:]
clons = dataset.variables['clon'][:]

vlats = dataset.variables['vlat'][:]
vlons = dataset.variables['vlon'][:]

clon1 = np.rad2deg( clons )
clat1 = np.rad2deg( clats )

vlon1 = np.rad2deg( vlons )
vlat1 = np.rad2deg( vlats )

clon_poly = clon[indcs_poly]
clat_poly = clat[indcs_poly]
#Plot the grid and city's shape file
figure, ax  = plt.subplots()
ax.plot(clon, clat, 'o', color='black', markersize=2.5)
ax.plot(clon_poly, clat_poly, 'o', color='green', markersize=2.5)
ax.plot(np.asarray(coords)[:, 0], np.asarray(coords)[:, 1], 'o', color='black', markersize=2.5)
ax.plot(clon1, clat1, 'o', color='red', markersize=2.5)
ax.plot(vlon1, vlat1, 'o', color='gray', markersize=2.5)
#plot_icongrid(dataset, ax)
plt.savefig('Extracted_points_with_xr_isle.png')