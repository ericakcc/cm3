#!/usr/bin/python
import read
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import math

# =====================================
#   From read routine:
#      GHCND_dim :  
#           (12, 60, 73, 96)
#      EX2_dim   : 
#           (12, 60, 73, 96)
#     Total days = 52560
# =====================================

# =====================================
# Theils : return y = ax + b
#          a = res[0], b = res[1]
# Conditions:
#          1. grid has at least 40 yrs
#          2. end no earlier than 2003
# =====================================

def global_mean(array, invalid, lat, cond1):
    mean = 0
    W = 0
    for i in range(array.shape[0]):
        s = 0
        w = 0
        for j in range(array.shape[1]):
            if cond1[i][j] >= 1  and array[i][j] > invalid:
                s = s + array[i][j]*math.cos(math.pi*lat[i]/180)
                w = w + math.cos(math.pi*lat[i]/180)
        mean += s
        W += w
    mean = mean / W
    return mean

print "model3.py"

ghcnd = read.GHCND
ex2 = read.EX2
fd = read.FD
idd = read.ID

ghcnd_trend = np.zeros((73, 96))
ex2_trend = np.zeros((73, 96))
fd_trend = np.zeros((73, 96))
id_trend = np.zeros((73, 96))

ghcnd_sig_test_lat = list()
ghcnd_sig_test_lon = list()
ex2_sig_test_lat = list()
ex2_sig_test_lon = list()
fd_sig_test_lat = list()
fd_sig_test_lon = list()
id_sig_test_lat = list()
id_sig_test_lon = list()

ghcnd_cond1 = np.zeros((73, 96))
ex2_cond1 = np.zeros((73, 96))
fd_cond1 = np.zeros((73, 96))
id_cond1 = np.zeros((73, 96))

ghcnd_cond2 = np.zeros((73, 96))
ex2_cond2 = np.zeros((73, 96))
fd_cond2 = np.zeros((73, 96))
id_cond2 = np.zeros((73, 96))

ghcnd_cond3 = np.zeros((73, 96))
ex2_cond3 = np.zeros((73, 96))
fd_cond3 = np.zeros((73, 96))
id_cond3 = np.zeros((73, 96))

ghcnd_ana = np.zeros(60)
ex2_ana = np.zeros(60)
fd_ana = np.zeros(60)
id_ana = np.zeros(60)

fd_ave = np.zeros(60)
id_ave = np.zeros(60)
time = np.linspace(1951, 2010, 6)

# ----------------
# Set Condition

for i in range(0, 73):
    for j in range(0, 96):
        for k in range(0, 60):
            if ghcnd[k][i][j] > -998:
                ghcnd_cond1[i][j] = ghcnd_cond1[i][j] + 1
            if ex2[k][i][j] > -90:
                ex2_cond1[i][j] = ex2_cond1[i][j] + 1
            if fd[k][i][j] > -998:
                fd_cond1[i][j] = fd_cond1[i][j] + 1
            if idd[k][i][j] > -998:
                id_cond1[i][j] = id_cond1[i][j] + 1

for i in range(0, 73):
    for j in range(0, 96):
        for k in range(52, 60):
            if ghcnd[k][i][j] > -998:
                ghcnd_cond2[i][j] = 1
            if ex2[k][i][j] > -90:
                ex2_cond2[i][j] = 1
            if fd[k][i][j] > -998:
                fd_cond2[i][j] = 1
            if idd[k][i][j] > -998:
                id_cond2[i][j] = 1
print 'cond3'
for i in range(0, 73):
    for j in range(0, 96):
        if ghcnd_cond1[i][j] >= 54 and ex2_cond1[73-i-1][(j+48)%96] >= 54:
            ghcnd_cond3[i][j] = 1
            ex2_cond3[73-i-1][(j+48)%96] = 1

# ------------------------
# Trend Calculation

for i in range(0, 73):
    for j in range(0, 96):
        if ghcnd_cond1[i][j] >= 40 and ghcnd_cond2[i][j] == 1:
            temp = list()
            tt = list()
            for k in range(0, 60):
                if ghcnd[k][i][j] > -998:
                    temp.append(ghcnd[k][i][j])
                    tt.append(1951+k)
            # tmp is a time series of a grid
            gres = stats.theilslopes(temp, tt, 0.90)
            ghcnd_trend[i][j] = gres[0]*10
            # significant test
            if 0 < gres[2] or 0 > gres[3]:
                ghcnd_sig_test_lat.append(read.lat1[i])
                ghcnd_sig_test_lon.append(read.lon1[j])
        else:
            ghcnd_trend[i][j] = np.nan

for i in range(0, 73):
    for j in range(0, 96):
        if ex2_cond1[i][j] >= 40 and ex2_cond2[i][j] == 1:
            temp = list()
            tt = list()
            for k in range(0, 60):
                if ex2[k][i][j] > -90:
                    temp.append(ex2[k][i][j])
                    tt.append(1951+k)
            # tmp is a time series of a grid
            exres = stats.theilslopes(temp, tt, 0.90)
            ex2_trend[i][j] = exres[0]*10
            # significant test
            if 0 < exres[2] or 0 > exres[3]:
                ex2_sig_test_lat.append(read.lat2[i])
                ex2_sig_test_lon.append(read.lon2[j])
        else:
            ex2_trend[i][j] = np.nan

for i in range(0, 73):
    for j in range(0, 96):
        if fd_cond1[i][j] >= 40 and fd_cond2[i][j] == 1:
            fday = list()
            tt = list()
            for k in range(0, 60):
                if fd[k][i][j] > -998:
                    fday.append(fd[k][i][j])
                    tt.append(1951+k)
            # tmp is a time series of a grid
            gres = stats.theilslopes(fday, tt, 0.90)
            fd_trend[i][j] = gres[0]
            # significant test
            if 0 < gres[2] or 0 > gres[3]:
                fd_sig_test_lat.append(read.lat3[i])
                fd_sig_test_lon.append(read.lon3[j])
        else:
            fd_trend[i][j] = np.nan

for i in range(0, 73):
    for j in range(0, 96):
        if id_cond1[i][j] >= 40 and id_cond2[i][j] == 1:
            iday = list()
            tt = list()
            for k in range(0, 60):
                if idd[k][i][j] > -998:
                    iday.append(idd[k][i][j])
                    tt.append(1951+k)
            # tmp is a time series of a grid
            gres = stats.theilslopes(iday, tt, 0.90)
            id_trend[i][j] = gres[0]
            # significant test
            if 0 < gres[2] or 0 > gres[3]:
                id_sig_test_lat.append(read.lat4[i])
                id_sig_test_lon.append(read.lon4[j])
        else:
            id_trend[i][j] = np.nan

for k in range(0, 60):
    ghcnd_ana[k] = global_mean(ghcnd[k][:][:], -998, read.lat1, ghcnd_cond3) 
    fd_ana[k] = global_mean(fd[k][:][:], -998, read.lat3, fd_cond1) 
    id_ana[k] = global_mean(idd[k][:][:], -998, read.lat4, id_cond1) 
    ex2_ana[k] = global_mean(ex2[k][:][:], -98, read.lat2, ex2_cond3)

id_ave = id_ana
fd_ave = fd_ana

plt.figure(1)
plt.plot(read.time, ghcnd_ana, 'b', label='HadGHCND')
plt.plot(read.time, ex2_ana, 'r', label='HadEX2')
plt.legend(loc='upper left')
plt.title('Global Land Mean')
plt.ylabel('Temperature ($^\circ$C)')
plt.savefig('./img/mean.png')

plt.figure(2)
#plt.plot(read.time, fd_ana, 'b', label='Anomaly')
plt.plot(read.time, id_ave, 'r', label='Mean')
plt.legend(loc='upper left')
plt.title('Global Land Mean')
plt.ylabel('Ice Days')
plt.savefig('./img/icedaysaveana.png')

gmean = np.mean(ghcnd_ana[10:39])
emean = np.mean(ex2_ana[10:39])
fmean = np.mean(fd_ana[10:39])
imean = np.mean(id_ana[10:39])

for k in range(0, 60):
    ex2_ana[k] = ex2_ana[k] - emean
    ghcnd_ana[k] = ghcnd_ana[k] - gmean
    fd_ana[k] = fd_ana[k] - fmean
    id_ana[k] = id_ana[k] - imean


plt.figure(14)
#plt.plot(read.time, fd_ana, 'b', label='Anomaly')
plt.plot(read.time, fd_ave, 'r', label='Mean')
plt.legend(loc='upper left')
plt.title('Global Land Mean')
plt.ylabel('Frost Days')
plt.savefig('./img/frostdaysaveana.png')

plt.figure(15)
plt.plot(read.time, id_ave, 'r', label='Mean')
plt.legend(loc='upper left')
plt.title('Global Land Mean')
plt.ylabel('Ice Days')
plt.savefig('./img/icedaysanaave.png')
# extension
ghcnd_extend = np.zeros((6, 73, 96), dtype=float)
ex2_extend = np.zeros((6, 73, 96), dtype=float)
threshold = 2
for i in range(0, 73):
    for j in range(0, 96):
        if ghcnd_cond1[i][j] >= 40 and ghcnd_cond2[i][j] == 1:
            temp = list()
            tt = list()
            for k in range(0, 60):
                if k % 10 <= 9:
                    if ghcnd[k][i][j] > -998:
                        temp.append(ghcnd[k][i][j])
                        tt.append(1951+k)
                        #print temp
                        #print tt
                if k % 10 == 9 and len(temp) > threshold and len(tt) > threshold:
                    # tmp is a time series of a grid
                    gres = stats.theilslopes(temp, tt, 0.90)
                    ghcnd_extend[k/10][i][j] = gres[0]
                    temp = list()
                    tt = list()
                elif k % 10 == 9 and len(temp) <= threshold and len(tt) <= threshold:
                    ghcnd_extend[k/10][i][j] = np.nan
        else:
            for m in range(0, 6):
                ghcnd_extend[m][i][j] = np.nan

# correlation
ghcnd_corrf = np.zeros(1096)
ghcnd_corri = np.zeros(1143)
fd_corr = np.zeros(1096)
id_corr = np.zeros(1143)

n = 0
k = 0
for i in range(0, 73):
    for j in range(0, 96):
        if ~np.isnan(ghcnd_trend[i][j]) and ~np.isnan(fd_trend[i][j]):
            ghcnd_corrf[n] = ghcnd_trend[i][j]
            fd_corr[n] = fd_trend[i][j]
            n = n + 1
        if ~np.isnan(ghcnd_trend[i][j]) and ~np.isnan(id_trend[i][j]):
            ghcnd_corri[k] = ghcnd_trend[i][j]
            id_corr[k] = id_trend[i][j]
            k = k + 1

plt.figure(16)
res = stats.theilslopes(fd_corr, ghcnd_corrf, 0.90)
lsq_res = stats.linregress(ghcnd_corrf, fd_corr)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ghcnd_corrf, fd_corr, 'b.')
ax.plot(ghcnd_corrf, res[1] + res[0] * ghcnd_corrf, 'r-')
ax.plot(ghcnd_corrf, lsq_res[1] + lsq_res[0] * ghcnd_corrf, 'g-')

plt.figure(17)
res = stats.theilslopes(id_corr, ghcnd_corri, 0.90)
lsq_res = stats.linregress(ghcnd_corri, id_corr)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ghcnd_corri, id_corr, 'b.')
ax.plot(ghcnd_corri, res[1] + res[0] * ghcnd_corri, 'r-')
ax.plot(ghcnd_corri, lsq_res[1] + lsq_res[0] * ghcnd_corri, 'g-')

print np.corrcoef(ghcnd_corrf, fd_corr)
print np.corrcoef(ghcnd_corri, id_corr)
plt.show() 
