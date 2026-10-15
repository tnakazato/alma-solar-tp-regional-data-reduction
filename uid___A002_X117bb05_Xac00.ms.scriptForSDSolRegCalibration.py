import os
import glob
import shutil
import numpy as np
import analysisUtils as aU
from scipy import interpolate
from datetime import datetime, timezone
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

#==============Parameter Setup============================

#-- Epehemris Table
# Name of Calibrated MS
asdm ='uid___A002_X117bb05_Xac00'

#--------- Measurement Set Name
mso = asdm + '.ms'
vis = mso + '.split.cal'
ephem = vis+'/FIELD/EPHEM2_Sun_Regional_60454.6.tab'

#-- Calibrated Full Sun image for Flux Calibration
FSfile1= 'member.uid___A001_X3666_X2b0.fullSun-num0_sci.spw2_104021MHz.sd.TP.cont.I.manual.sd.rescl'

#-- directories of the regional maps
#"member."+MOUS ID+'regSun.scn*.image'
orgfile = glob.glob('member.uid___A001_X3666_X2b0.regSun.scn*.image')

#-- Size of FoV (Diameter of double-circle pattern)
dcDia = 5.*60.

## Don't change w/o the comments from specialits.
#-- Threshold for seraching the Refernce point on the regional map
thTb = 5000.
#-- size of seraching area comparing with FoV
safeF = 0.2

#============= Search regional maps
orgfile.sort()

#============= Determining the reference point automatically

#----get dimensions and create narray
forShape = imhead(orgfile[0])
imgsX = forShape['shape'][0]
imgsY = forShape['shape'][1]
strImgsX = str(imgsX-1)
strImgsY = str(imgsY-1)
imgNum = len(orgfile)
imgCube = np.zeros([imgsX, imgsY, imgNum])

#---- read images
for imgCnt in range(0,len(orgfile)):
  inpfile = orgfile[imgCnt]
  imgCASA = imval(inpfile, box = '0, 0, '+strImgsX+', '+strImgsY)
  imgCube[0:,0:,imgCnt] = imgCASA['data']

#---- calcuration for determining the refarance point
mskImg = np.average(imgCube, axis=2)
mskImg[mskImg < thTb] = 0
mskImg[mskImg >= thTb] = 1
stdImg = np.std(imgCube, axis=2)

radFcen = np.zeros([imgsX, imgsY])

for i in range(0, imgsX):
    for j in range(0,imgsY):
        radFcen[i,j] = np.sqrt(((i - imgsX/2.)*forShape['incr'][0]*3600.*180./np.pi)**2.
                             + ((j - imgsY/2.)*forShape['incr'][1]*3600.*180./np.pi)**2.)

radFcen[radFcen < dcDia/2. * safeF] = 1
radFcen[radFcen >= dcDia/2. * safeF] = 0

filStdImg = stdImg*mskImg*imgCASA['mask']*radFcen
filStdImg[filStdImg == 0 ] = 1e8
FCRefP = np.unravel_index(np.argmin(filStdImg),filStdImg.shape)

print('####################################################')
print('Reference Point Coordinate [Pixels in Image] X:'+str(FCRefP[0])+'  Y:'+str(FCRefP[1]))
print('####################################################')

#============= re-scaling
#----- Read the date-obs of the Full-Sun images
fSTime1 = imhead(FSfile1, mode='get', hdkey='date-obs')[:-7]
#fSTime2 = imhead(FSfile2, mode='get', hdkey='date-obs')[:-7]
fSTims1 = qa.convert(qa.quantity(fSTime1), outunit='s')['value']
#fSTims2 = qa.convert(qa.quantity(fSTime2), outunit='s')['value']

#----- Read the ephemeris of the target region
tb.open(ephem, nomodify=True)
epheMJ = tb.getcol("MJD")
epheRA = tb.getcol("RA")
epheDE = tb.getcol("DEC")
tb.close()

epheTI = qa.convert(qa.quantity(epheMJ, 'd'), outunit='s')['value']

fRA = interpolate.interp1d(epheTI, epheRA, kind='cubic')
fDE = interpolate.interp1d(epheTI, epheDE, kind='cubic')

#------ Read the parameters of the reginal maps
rCdel1 = imhead(orgfile[0], mode='get', hdkey='cdelt1')
rCdel2 = imhead(orgfile[0], mode='get', hdkey='cdelt2')
rCrpi1 = imhead(orgfile[0], mode='get', hdkey='crpix1')
rCrpi2 = imhead(orgfile[0], mode='get', hdkey='crpix2')
#rCval1 = imhead(orgfile[0], mode='get', hdkey='crval1')
#rCval2 = imhead(orgfile[0], mode='get', hdkey='crval2')
rBmaj1 = imhead(orgfile[0], mode='get', hdkey='bmaj')

#----- Calculate the RADEC coordinate of the reference point at the full-sun images
RefFSRA1 = fRA(fSTims1)*np.pi/180. +  (FCRefP[0] - rCrpi1) * rCdel1['value']
RefFSDe1 = fDE(fSTims1)*np.pi/180. +  (FCRefP[1] - rCrpi2) * rCdel2['value']
#RefFSRA1 = rCval1['value']+(FCRefP[0] - rCrpi1) * rCdel1['value']
#RefFSDe1 = rCval2['value']+(FCRefP[1] - rCrpi2) * rCdel2['value']
#RefFSRA2 = fRA(fSTims2)*np.pi/180. +  (FCRefP[0] - rCrpi1) * rCdel1['value']
#RefFSDe2 = fDE(fSTims2)*np.pi/180. +  (FCRefP[1] - rCrpi2) * rCdel2['value']

RefRADEC1 = aU.rad2radec(ra = RefFSRA1, dec = RefFSDe1, prec=1, hmsdms=True, delimiter=',', verbose=False)
#RefRADEC2 = aU.rad2radec(ra = RefFSRA2, dec = RefFSDe2, prec=1, hmsdms=True, delimiter=',', verbose=False)

#----- Calculate the Brightness Temperature at the Refercen region in Full-Sun
FStat1 = imstat(FSfile1, region='circle[['+RefRADEC1+']],'+str(rBmaj1['value'])+'arcsec]')
print('circle[['+RefRADEC1+']],'+str(rBmaj1['value'])+'arcsec]')
#FStat2 = imstat(FSfile2, region='circle[['+RefRADEC2+'],'+str(rBmaj1['value'])+'arcsec]')
#print('circle[['+RefRADEC2+'],'+str(rBmaj1['value'])+'arcsec]')
#BTRef = interpolate.interp1d([fSTims1,fSTims2],[FStat1['mean'][0],FStat2['mean'][0]], kind='linear')
BTRef = FStat1['mean'][0]


#---- setup for the summary plots
Fimtime = [datetime.strptime(fSTime1,'%Y/%m/%d/%H:%M:%S'),datetime.strptime(fSTime1,'%Y/%m/%d/%H:%M:%S') ]
FulSBT = [FStat1['mean'][0],FStat1['mean'][0]]
Rimtime = []
RorgTem = []
RestTem = []
CorRatV = []

#======= Manipulate each regional map===============

for inpfile in orgfile:
    rslfile = inpfile[:-6]+'.rescl'
    print('Now processing: [Org Map] '+inpfile+' / [Mod Map]'+rslfile)

#----- Delet old file
    if os.path.exists(rslfile):
        shutil.rmtree(rslfile)

#------ Get the time of the reginal maps
    rSTime = imhead(inpfile, mode='get', hdkey='date-obs')[:-7]
    rsTimD = datetime.strptime(rSTime,'%Y/%m/%d/%H:%M:%S')
    Rimtime.append(rsTimD)

#------ Calculate the Brightness Temperature at the Refercen region in a Regional map
    RStat1 = imstat(inpfile, region='circle[['+str(FCRefP[0])+'pix,'+str(FCRefP[1])+'pix],'+str(rBmaj1['value'])+'arcsec]')
    RorgTem.append(RStat1['mean'][0])
#    BTEst = BTRef(qa.convert(qa.quantity(rSTime), outunit='s')['value'])
    BTEst = BTRef
    RestTem.append(BTEst)

#------ Calculate the correction ratio and apply it to the image
    CorRat = BTEst / RStat1['mean'][0]
    CorRatV.append(CorRat)
    immath(imagename=inpfile, expr = 'IM0*'+str(CorRat), outfile = rslfile)

#------ Create FITS file
    if os.path.exists(inpfile[:-6]+'_rscl.fits'):
        os.remove(inpfile[:-6]+'_rscl.fits')
    exportfits(imagename = rslfile, fitsimage = inpfile[:-6]+'_rscl.fits')

#========plot the averaged Brightness Temperature of the reference region=========

fig = plt.figure(figsize=(8., 9.), tight_layout=True)
ax01 = fig.add_subplot(2, 1, 1,ylabel='Ave. Brightness Temp. @ Ref. Point', title='Solar TP Regional Map/Flux Cal.: '+vis)
ax02 = fig.add_subplot(2, 1, 2,ylabel='Correction Factor')

ax01.plot(Rimtime, RorgTem, marker = 'p', color='g', label='Estimate from Tsys Cal. Only')
ax01.plot(Rimtime, RestTem, marker = 'o', color='b', label='Estimate from Full-Sun maps')
ax01axes = ax01.xaxis
ax01axes.set_major_formatter(DateFormatter('%H:%M'))
ax01.legend()

ax02.plot(Rimtime, CorRatV, marker = '<', color='r')
ax02axes = ax02.xaxis
ax02axes.set_major_formatter(DateFormatter('%H:%M'))

fig.savefig("FluxCalPlot.png")

plt.close()
