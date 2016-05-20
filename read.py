#!/usr/bin/python
import os, sys
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
# ================Read NC file ====================#
#   Dimension:
#       File 1: GHCND 1949 ~ 2011
#           (12, 63, 73, 96)
#           (month, year, lat, lon)
#       File 2:  EX2  1901 ~ 2010
#           (12, 110, 73, 96)
#           (month, year, lat, lon)
#
#		time = [19490100, ....]
#
# ==================================================#
nc1 = '../data/HadGHCND_1949-2011_TXx.nc'
nc2 = '../data/H2_TXx_1901-2010_RegularGrid_global_3.75x2.5deg_LSmask.nc'
nc3 = '../data/HadGHCND_1949-2011_FD.nc'
nc4 = '../data/HadGHCND_1949-2011_ID.nc'
root1 = Dataset(nc1, mode='r')
root2 = Dataset(nc2, mode='r')
root3 = Dataset(nc3, mode='r')
root4 = Dataset(nc4, mode='r')
GHCND = np.array(root1.variables['Annual'][2:62][:][:], dtype=np.float64)
EX2 = np.array(root2.variables['Ann'][50:][:][:], dtype=np.float64)
FD = np.array(root3.variables['Annual'][2:62][:][:], dtype=np.float64)
ID = np.array(root4.variables['Annual'][2:62][:][:], dtype=np.float64)

#quit()
# ================ Region Coordinates ====================#
#		Greenland 60~80N, 60~20W		
#		China	  105~120E, 25~30N
#		GHCND	  (60, 73, 96) 
#		EX2		  (60, 73, 96)
# ========================================================#

# test[0]: greenland, test[1]: china
# test[0][0]: lat start, test[0][2]: lon start
lat1 = np.array(root1.variables['lat'])
lon1 = np.array(root1.variables['lon'])
lat2 = np.array(root2.variables['lat'])
lon2 = np.array(root2.variables['lon'])
lat3 = np.array(root3.variables['lat'])
lon3 = np.array(root3.variables['lon'])
lat4 = np.array(root4.variables['lat'])
lon4 = np.array(root4.variables['lon'])

test = np.zeros((2, 4))

for i in range(0, 96):
    if lon2[i] >= 180:  
        lon2[i] -= 360
#print 'lat'
#print lat1
#print lat3

#print 'lon'
#print lon1
#print lon3

for i in range(0, 73):
	# Greenland
	if lat1[i] == 60:
		test[0][0] = i
	if lat1[i] == 80:
		test[0][1] = i
# China
	if lat1[i] == 25:
		test[1][0] = i
	if lat1[i] == 30:
		test[1][1] = i

for i in range(0, 96):
	if lon1[i] == 105:
		test[1][2] = i
	if lon1[i] == 120:
		test[1][3] = i
	if lon1[i] == -60:
		test[0][2] = i
	if lon1[i] == -20:
		test[0][3] = i
# GHCND

time = np.linspace(1951, 2010, 60)


#GHCND = np.array([np.array(root1.variables['January'][2:62][:][:], dtype=np.float32),  
#                np.array(roo2.variables['February'][2:62][:][:], dtype=np.float32), 
#                np.array(roo2.variables['March'][2:62][:][:], dtype=np.float32), 
#                np.array(roo2.variables['April'][2:62][:][:], dtype=np.float32),
#                np.array(roo2.variables['May'][2:62][:][:], dtype=np.float32),
#                np.array(roo2.variables['June'][2:62][:][:], dtype=np.float32),
#                np.array(roo2.variables['July'][2:62][:][:], dtype=np.float32),
#                np.array(roo2.variables['August'][2:62][:][:], dtype=np.float32),
#                np.array(roo2.variables['September'][2:62][:][:], dtype=np.float32),
#                np.array(roo2.variables['October'][2:62][:][:], dtype=np.float32),
#                np.array(roo2.variables['November'][2:62][:][:], dtype=np.float32),
#                np.array(roo2.variables['December'][2:62][:][:], dtype=np.float32)])


#EX2 = np.array([np.array(root2.variables['Jan'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Feb'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Mar'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Apr'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['May'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Jun'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Jul'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Aug'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Sep'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Oct'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Nov'][:][:][50:], dtype=np.float32),
#                np.array(root2.variables['Dec'][:][:][50:], dtype=np.float32)])
#quit()
#for i in range(4, 13):
#    for j in range(32, 43):
#        g = np.zeros(60)
#        for k in range(0, 60):
#            if GHCND[k][i][j] > -998:
#                g[k] = GHCND[k][i][j]
#        plt.plot(time, g)
#        plt.savefig('./test/'+str(i)+str(j)+'.png')
#        plt.clf()
