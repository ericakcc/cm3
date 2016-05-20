#!/usr/bin/python
import read
import model3
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

print 'plot.py'

lons1 = read.root1.variables['lon'][:]
lats1 = read.root1.variables['lat'][:]
lats1, lons1 = np.meshgrid(lats1, lons1)

print lons1
lons2 = read.root2.variables['lon'][:]
lats2 = read.root2.variables['lat'][:]
etrend, lons2 = shiftgrid(180, model3.ex2_trend, lons2, start=False)
lats2, lons2 = np.meshgrid(read.lat2, lons2)
print lons2
gtrend = model3.ghcnd_trend
#etrend = model3.ex2_trend
fdtrend = model3.fd_trend
idtrend = model3.id_trend
ghcnd_extend = model3.ghcnd_extend
# masked array
z = ma.masked_invalid(gtrend)
zz = ma.masked_invalid(etrend)
zzz = ma.masked_invalid(fdtrend)
Z = ma.masked_invalid(idtrend)
s = ma.masked_invalid(ghcnd_extend)

cmap = plt.get_cmap('coolwarm', 10)

fig = plt.figure(3)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
#m = Basemap(lon_0=0, lat_0=0, llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='mill')
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
m.plot(model3.ghcnd_sig_test_lon, model3.ghcnd_sig_test_lat, '+')
im1 = m.pcolormesh(lons1, lats1, z.T, shading='flat', cmap=cmap, vmin=-1.5, vmax=1.5, latlon=True)
cbar = m.colorbar(im1, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadGHCND 1951-2010')
plt.savefig('./img/ghcnd60.png')
#quit()
fig = plt.figure(4)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
m.plot(model3.ex2_sig_test_lon, model3.ex2_sig_test_lat, '+')
im2 = m.pcolormesh(lons2, lats2, zz.T, shading='flat', cmap=cmap, vmin=-1.5, vmax=1.5, latlon=True)

cbar = m.colorbar(im2, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadEX2 1951-2010')
plt.savefig('./img/hadex60.png')

fig = plt.figure(5)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
im3 = m.pcolormesh(lons1, lats1, s[0].T, shading='flat', cmap=cmap, vmin = -1, vmax = 1, latlon=True)

cbar = m.colorbar(im3, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadGHCND 1951-1960')
plt.savefig('./img/ghcnd5160.png')

fig = plt.figure(6)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
im3 = m.pcolormesh(lons1, lats1, s[1].T, shading='flat', cmap=cmap, vmin = -1, vmax = 1, latlon=True)

cbar = m.colorbar(im3, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadGHCND 1961-1970')
plt.savefig('./img/ghcnd6170.png')

fig = plt.figure(7)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
im3 = m.pcolormesh(lons1, lats1, s[2].T, shading='flat', cmap=cmap, vmin = -1, vmax = 1, latlon=True)

cbar = m.colorbar(im3, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadGHCND 1971-1980')
plt.savefig('./img/ghcnd7180.png')

fig = plt.figure(8)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
im3 = m.pcolormesh(lons1, lats1, s[3].T, shading='flat', cmap=cmap, vmin = -1, vmax = 1, latlon=True)

cbar = m.colorbar(im3, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadGHCND 1981-1990')
plt.savefig('./img/ghcnd8190.png')

fig = plt.figure(9)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
im3 = m.pcolormesh(lons1, lats1, s[4].T, shading='flat', cmap=cmap, vmin = -1, vmax = 1, latlon=True)

cbar = m.colorbar(im3, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadGHCND 1991-2000')
plt.savefig('./img/ghcnd9100.png')

fig = plt.figure(10)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
im3 = m.pcolormesh(lons1, lats1, s[5].T, shading='flat', cmap=cmap, vmin = -1, vmax = 1, latlon=True)

cbar = m.colorbar(im3, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadGHCND 2001-2010')
plt.savefig('./img/ghcnd0110.png')


fig = plt.figure(11)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
#m = Basemap(lon_0=0, lat_0=0, llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='mill')
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawcoastlines()
m.drawcountries()
im3 = m.pcolormesh(lons1, lats1, s[5].T, shading='flat', cmap=cmap, vmin = -1, vmax = 1, latlon=True)

cbar = m.colorbar(im3, location='bottom', pad="10%")
cbar.set_label('Trend($^\circ$C/decade)')
ax.set_title('HadGHCND 2001-2010')
plt.savefig('./img/ghcnd0110.png')


fig = plt.figure(12)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
#m = Basemap(lon_0=0, lat_0=0, llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,projection='mill')
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
im1 = m.pcolormesh(lons1, lats1, zzz.T, shading='flat', cmap=cmap, vmin=-1.5, vmax=1.5, latlon=True)

cbar = m.colorbar(im1, location='bottom', pad="10%")
cbar.set_label('Trend(days)')
ax.set_title('HadGHCND 1951-2010')
plt.savefig('./img/ghcndforst60.png')

fig = plt.figure(13)
ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
m = Basemap(projection='cyl', lon_0=0, lat_0=0)
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
m.drawmapboundary(fill_color='white')
m.drawcoastlines()
m.drawcountries()
im1 = m.pcolormesh(lons1, lats1, Z.T, shading='flat', cmap=cmap, vmin=-1.5, vmax=1.5, latlon=True)

cbar = m.colorbar(im1, location='bottom', pad="10%")
cbar.set_label('Trend(days)')
ax.set_title('HadGHCND 1951-2010')
plt.savefig('./img/ghcndice60.png')

plt.show()
