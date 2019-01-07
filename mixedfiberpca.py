import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize':'xx-large',
         'axes.labelsize':'xx-large',
         'axes.titlesize':'xx-large',
         'ytick.labelsize':'xx-large',
         'xtick.labelsize':'xx-large'}#,
         #'axes.facecolor':'white',
         #'figure.facecolor':'white'}
pylab.rcParams.update(params)
#plt.style.use(['dark_background','seaborn-pastel'])

import glob
import os

from astropy.io import fits

def gauss(a,b,c,x):
    return a*np.exp(-(x-b)**2/c)**2

pattern = '/media/maja/Elements/rebinned20180822/rebinned20180822v02*_exp03_multi*.fits.fits'
pattern2 = '/media/maja/Elements/rebinned20180822/rebinned*_exp0?_multi*_042_*LL.fits.fits'
#pattern = '/media/maja/Elements/work/03946/hetdex/maverick/red1/reductions/201801*/virus/*/exp0?/virus/multi*_042_*LL.fits'
ff = glob.glob(pattern)
ff = np.array(ff)
ff2 = glob.glob(pattern2)
#print(len(ff))

fiber = []
here = [x[-20:-18]!='09' for x in ff]

#print(here)
for f in ff[np.where(here)]:
    hdu = fits.open(f)
    fiber.append(hdu['sky_spectrum'].data)
    hdu.close()
fiber = np.array(fiber)
fiber.shape

fiberarray = np.concatenate([x[:,4:] for x in fiber])
fiberarray[np.where(np.isnan(fiberarray))] = 0

i = 5000
fiberarray[i] = fiberarray[i]+gauss(10,600,8, np.arange(0,fiberarray.shape[-1],1))

fiberarray.shape

meanspec = np.nanmean(fiberarray[:], axis=0)
stdspec = np.nanstd(fiberarray[:],axis=0)

fiberarray_0 = (fiberarray-meanspec)/stdspec

#stdspec.shape
plt.figure(figsize=(20,4))
plt.plot(meanspec)
plt.title('mean spectrum')
plt.show()

cov = np.cov(fiberarray_0.T)

eigenvals, eigenvecs = np.linalg.eig(cov)
eigenvals = np.real(eigenvals)
eigenvecs = np.real(eigenvecs)

eigenpairs = [(np.abs(eigenvals[i]), eigenvecs[:,i]) for i in np.argsort(abs(eigenvals))[::-1]]

"""plt.figure(figsize=(20,4))
plt.plot([x[0] for x in eigenpairs])
plt.axvline(100)
plt.axvline(150)
plt.yscale('log')"""

ncomp = 150

imp = np.array([eigenpairs[i][1] for i in range(ncomp)])

fiberpca = np.dot(fiberarray_0, imp.T)

new = np.dot(fiberpca, imp)

newspec = new*stdspec+meanspec

plt.figure(figsize=(20,4))
plt.plot(fiberarray[i,580:620]-newspec[i,580:620])
plt.plot(gauss(10,600,8, np.arange(0,fiber.shape[-1],1))[580:620])
plt.show()

plt.figure(figsize=(20,4))
plt.plot(fiberarray[i,:])
plt.plot(newspec[i,:])
plt.show()

re = (fiberarray-newspec)/fiberarray

stddevi = np.nanstd(re[i])
print('\nstd of re[i]: ',stddevi)

plt.figure(figsize=(20,4))
plt.plot(re[i])
plt.axhline(stddevi)
plt.axhline(-stddevi)
plt.title('std {:.3f}'.format(stddevi))
plt.show()


"""plt.figure(figsize=(25, 17))
plt.subplot(341)
plt.scatter(fiberpca[:,0], fiberpca[:,1], alpha=0.1)
plt.xlabel('PC 1')
plt.ylabel('PC 2')
plt.subplot(342)
plt.scatter(fiberpca[:,1], fiberpca[:,2], alpha=0.1)
plt.xlabel('PC 2')
plt.ylabel('PC 3')
plt.subplot(343)
plt.scatter(fiberpca[:,2], fiberpca[:,3], alpha=0.1)
plt.xlabel('PC 3')
plt.ylabel('PC 4')
plt.subplot(344)
plt.scatter(fiberpca[:,3], fiberpca[:,4], alpha=0.1)
plt.xlabel('PC 4')
plt.ylabel('PC 5')
plt.subplot(345)
plt.scatter(fiberpca[:,4], fiberpca[:,5], alpha=0.1)
plt.xlabel('PC 5')
plt.ylabel('PC 6')
plt.subplot(346)
plt.scatter(fiberpca[:,5], fiberpca[:,6], alpha=0.1)
plt.xlabel('PC 6')
plt.ylabel('PC 7')
plt.subplot(347)
plt.scatter(fiberpca[:,6], fiberpca[:,7], alpha=0.1)
plt.xlabel('PC 7')
plt.ylabel('PC 8')
plt.subplot(348)
plt.scatter(fiberpca[:,7], fiberpca[:,8], alpha=0.1)
plt.xlabel('PC 8')
plt.ylabel('PC 9')
plt.subplot(349)
plt.scatter(fiberpca[:,8], fiberpca[:,9], alpha=0.1)
plt.xlabel('PC 9')
plt.ylabel('PC 10')
plt.subplot(3,4,10)
plt.scatter(fiberpca[:,9], fiberpca[:,10], alpha=0.1)
plt.xlabel('PC 10')
plt.ylabel('PC 11')
plt.subplot(3,4,11)
plt.scatter(fiberpca[:,10], fiberpca[:,11], alpha=0.1)
plt.xlabel('PC 11')
plt.ylabel('PC 12')
plt.subplot(3,4,12)
plt.scatter(fiberpca[:,11], fiberpca[:,12], alpha=0.1)
plt.xlabel('PC 12')
plt.ylabel('PC 13')
#plt.savefig('fiberpcas.png', bbox_inches='tight')

for j in range(0,10):
    plt.figure(figsize=(20,4))
    plt.plot(imp[j])
"""
