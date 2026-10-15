import os
import numpy as np
import math
import analysisUtils as aU
from scipy import interpolate

#---------Input Parameters------------------------
# Name of Calibrated MS
asdm ='uid___A002_X117bb05_Xac00'
#--------- Measurement Set Name
mso = asdm + '.ms'
vis = mso + '.split.cal'
#-----------
#the directory name of ephemris for target region 
#ephem = vis+'/FIELD/EPHEM2_Sun_reg_59872.5.tab'
#head of output file : "member."+MOUS ID+'regSun'
ohead = 'member.uid___A001_X3666_X2b0.regSun'
#Size of Fov [arcsec] (Double-Circle)
sFoV= 5.*60.
#image size = sFoV/cell[6"(B3), 3.3"(B5), 3"(B6), 2.8"(B7)] * 1.5
imsize=[80,80]
#Spectrum Window ID
spw = 1
# Antenna name
antenna = 'PM04'
# Field ID and Field namefor regional mapping
fieldID = '2'
# Scan# of Regional mapping
scanNums = [8, 10, 12]

#Scan Patttern
pattern = 'DC'         # Lis: Lissajous, DC:Double-Circle
#sampling length /Band3 = 25 , Band5 = , Band6 = 15, Band7 = 12,
samplingLength = 25    # arcseconds

#-----------------------------
##--Scanning parameters## Don't change###
azimuth = 0.           # degree (Do not need to change)
elevation = 60.         # degere (Do not need to change)

##--Spec of Antenna
areaAspectMax = 9.0
limitFrequencyMax = 1.6 #Hz
defaultFreqHigh = 1.0 # Hz
limitFrequency = 1.6 # Hz
limitAccelerationMax = math.radians(3.0) # deg/s^2
limitVelocityMax     = math.radians(1.0) # deg/s
timeStep = 0.048

#-----------------------------
def maxSafeScaleFactor(freq, vel, accel, limitFrequencyMax):
        return min([limitFrequencyMax / freq, limitVelocityMax / vel, math.sqrt(limitAccelerationMax / accel)])

def generateCenterdPatternSpec(width, height, samplingLength, limitFrequencyMax, subPattern = "DOUBLE_CIRCLE"):

    width = math.radians(width/3600.)
    height = math.radians(height/3600.)
    samplingLength = math.radians(samplingLength/3600.)

    if (limitFrequencyMax > limitFrequency):
        limitFrequencyMax = limitFrequency

	# These can also be done in two terms
    nc = 2
    # Find resolution needed
    size = max(width, height)

    nbig = int(math.pi*size/(2*samplingLength) + 1)

    # Set some limits
    if(nbig < 2):
         nbig = 2
    elif(nbig > 500):
        print('Unreasonably large map - nbig = %d' % nbig)
        nbig = 501

    if subPattern == "DOUBLE_CIRCLE":
        lbig = 1
        kbig = nbig
    elif subPattern == "STAR":
        lbig = nbig - 1
        kbig = lbig + 1
    elif subPattern == "CLOVER":
        lbig = (nbig-2) / 2 + 1
        kbig = lbig*2 + 1
    elif subPattern == "LUCKY":
        lbig = (nbig-2)/3 + 1
        kbig = lbig*3 + 1
    else:
        raise Exception("Unsupported subPattern: '%s'" % str(subPattern))

    f_high = float(defaultFreqHigh)
    duration = kbig / f_high
    f_low = lbig / duration

    #print("original F_low  = %.6f" % f_low)
    #print("original F_high = %.6f" % f_high)
    print("original duration = %.6f" % duration)

    # Max speed is when the two circular motions add
    v_max = 2.0*math.pi*(f_low + f_high) * size/4.0
    # Similarly max acceleration is when these add.
    a_max = math.pi*math.pi*(f_low*f_low + f_high*f_high) * size
    # Check these against limits
    scale = maxSafeScaleFactor(f_high, v_max, a_max, limitFrequencyMax)
    #print("original scale = %.6f" % scale)
    # Adjust to make sure time for pattern is an integer number of TEs, round up.
    n_points = int((duration / scale) / timeStep) + 1
    scale = duration / (n_points * timeStep)

    print("scale = %.6f" % scale)

    f_high   *= scale
    f_low    *= scale
    v_max    *= scale
    a_max    *= scale*scale
    #duration /= scale

    #print("Max velocity:     %9.6f deg/s"   % math.degrees(v_max))
    #print("Max acceleration: %9.6f deg/s^2" % math.degrees(a_max))
    #print("generateCenteredPatternSpec() result: K = %d, L = %d, f_high = %.6f Hz, f_low = %.6f Hz, duration = %.4f s" % (kbig, lbig, f_high, f_low, duration))

    return duration, 3600.0*math.degrees(0.25*width), 3600.0*math.degrees(0.25*height), f_low, f_high
#----------------------------------
obsPar = generateCenterdPatternSpec(sFoV, sFoV, samplingLength, limitFrequencyMax, subPattern = "DOUBLE_CIRCLE")
patDur = obsPar[0]


#============================================================================

#--- What band are we observing (3, 5, 6, or 7)?

# Observing Band not yet determined.
oBand=-1

# Use the (always present) 1st EB to estimate the LO frequency
msc0 = vis

tb.open(msc0+'/SPECTRAL_WINDOW', nomodify=True)
# Ref. frequencies for all SPWs
obsFreqs = tb.getcol("REF_FREQUENCY")
tb.close()
# Average of ref. frequencies for science SPWs 0,1,2, and 3 = approx. LO freq.
# rescaled to GHz.
freqLO=0.25*(obsFreqs[0]+obsFreqs[1]+obsFreqs[2]+obsFreqs[3])*1.0e-9

# Fixed LO freqs. (in GHz) of Solar ObsMode and corresponding Bands:
# A dictionary
solarLOs={100.0:3, 194.0:5, 239.0:6, 346.6:7}

# Determine the observing band
for solFreq,solBand in solarLOs.items():

    # Calculate relative diff. to nominal solar LO frequency
    relDiff=abs(freqLO-solFreq)/solFreq
    
    # Are we within +/- 10% off the nominal LO freq?
    if (relDiff < 0.1): 
        oBand=solBand
        print("Your observation is in Band "+str(oBand)+".")
        break

# This 'else' is an alternative to 'for'
else:
    sys.exit("Error: Observing frequency does not correspond to Solar ObsMode.")


#--- Set imaging parameters specific for the Band

if (oBand == 3):
        
    # Imaging parameters for Band 3
    cell=['6.0arcsec','6.0arcsec']
#    imsize=[400,400]
    convsupport = 6
    qs_level = 7300.0 # For immath

    # _SADA_: Shift the circle a bit off the center to avoid the bright areas
    # from the statistics, if needed.
    region='circle[[200pix,200pix],40pix]' #imstat

elif (oBand==5):
    
    # Imaging parameters for Band 5

    # TO DO: Set the proper values, so far they are just roughly estimated
    #        by interpolation between Band 3 and 6 (MB).        
    
    cell=['3.3arcsec','3.3arcsec']
#    imsize=[720,720]
    convsupport = 6
    qs_level = 6300.0 # For immath

    # _SADA_: Shift the circle a bit off the center to avoid the bright areas
    # from the statistics, if needed.
    region='circle[[360pix,360pix],72pix]' #imstat

elif (oBand==6):
    
    # Imaging parameters for Band 6
    cell=['3.0arcsec','3.0arcsec']
#    imsize=[800,800]
    convsupport = 6
    qs_level = 5900.0 # For immath

    # _SADA_: Shift the circle a bit off the center to avoid the bright areas
    # from the statistics, if needed.
    region='circle[[400pix,400pix],80pix]' #imstat

else:
    # Imaging parameters for Band 7
    cell=['2.8arcsec','2.8arcsec']
#    imsize=[860,860]
    convsupport = 6
    qs_level = 5500.0 # For immath

    # _SADA_: Shift the circle a bit off the center to avoid the bright areas
    # from the statistics, if needed.
    region='circle[[430pix,430pix],86pix]' #imstat


#--- Common parameters (not Band-specific).

# Frequency for selected spw in MHz as integer
imFreq=int(round(obsFreqs[spw]*1.0e-6))


# Actual Durations are calculated from the data obtained on 2 July 2021.
# "Simulation" are the "original duration" calculated with "generateCenterdPatternSpec".
#------------------
# Actual Duration for scanning of 100" with double-circle: 11 sec
# Simulation: original duration = 11.000000
# ActualDuration for scanning of 200" with double-circle: 21 sec
# Simulation: original duration = 21.000000
# ActualDuration for scanning of 300" with double-circle: 32 sec
# Simulation: original duration = 32.000000
# ActualDuration for scanning of 600" with double-circle: 63 sec
# Simulation: original duration = 63.000000
#---------------------------------

#----------------------------------
for j in range(0, len(scanNums)):
    scanNum = scanNums[j]
    fileHead = ohead+'.scn'+'{:02d}'.format(scanNum)+'_num'

    scanTim = aU.computeDurationOfScan(scanNum, vis=vis, returnSubscanTimes=True)
    imgNum = int(np.fix(scanTim[0]/patDur))

    for i in range(0, imgNum):
        stScanS = qa.quantity(scanTim[2][1][0]+i*patDur, 's')
        stScan = qa.time(qa.quantity(scanTim[2][1][0]+i*patDur, 's'), form='ymd',prec=8)
        edScan = qa.time(qa.quantity(scanTim[2][1][0]+(i+1)*patDur, 's'), form='ymd',prec=8)
        timerange= stScan[0]+' ~ '+edScan[0]
        outfile = fileHead + '{:04d}'.format(i)+'_sci.spw'+str(spw)+'_'+str(imFreq)+'MHz.sd.TP.cont.I.manual.sd'
        print('Now processing: ' + '{:04d}'.format(i) +'/'+'{:04d}'.format(imgNum) + ':' +timerange+':'+ outfile)

        os.system('rm -rf '+outfile+'*')
        tsdimaging(infiles = vis,
                   spw = str(spw),
                   antenna =antenna,
                   field = fieldID,
                   brightnessunit='K',
                   nchan = 1,
                   cell = cell,
                   outfile = outfile,
                   scan = str(scanNum),
                   imsize=imsize,
                   gridfunction = 'SF',
                   convsupport = 6,
                   intent ='*ON_SOURCE*',
                   timerange = timerange,
                   stokes = 'I',
                   specmode='cube',
                   outframe='lsrk')

