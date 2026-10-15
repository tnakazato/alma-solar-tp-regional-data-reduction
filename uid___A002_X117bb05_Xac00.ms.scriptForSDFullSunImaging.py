# ALMA Data Reduction Script for CASA5
#
# Original template by Masumi Shimojo (EA ARC)
# Modified by Miro Barta <barta@asu.cas.cz> (EU ARC, Czech node) on 2018/07/14
# Naming scheme updated to convention for Cy6+, 2019/05/14, barta@asu.cas.cz
# Included Band7, naming scheme updated to convention for Cy7+,
#         on 2020/07/02, masumi.shimojo@nao.ac.jp
# Generalisation for Band 5 (and more), on 2021/03/14, barta@asu.cas.cz
#
#
# 2024-03-04 Masumi Shimojo <masumi.shimojo@nao.ac.jp>, EA ARC - NAOJ
#   "tsdimaging" task is used instead of "sdimaging" task. 
#
# Recommended CASA version: 6.0 and higher
#-------------------------------------------

#============================================================================
# _SADA_: This part to be changed by data analyst

# Project ID
projectID='2023.1.01040.S'

# SB name
sbName='Sun_10_e_03_TP'

# EB/ASDM uid - upper part common for entire SB.
asdm_commonPart_ID='uid___A002_X117bb05_Xac00'

# MOUS uid - find it in AQUA as ObsUnitSet Status UID
mousID='uid___A001_X3666_X2b0'

# 2nd full-sun scan 
SecondScan = 'none'

# Antenna. Select that one with good Tsys for all EBs (check Tsys.plots)
antenna='PM04'

# SPW with the best Tsys values.
spw=2

#============================================================================

#============================================================================
# 'applyonly=False' (QA2 mode) is usefull for just one thing:
# Make the 1st image extra in the QA2 mode for fast checks. QA2 analyst
# should set it (applyonly=False) before executing the script. If he/she forgets
# nothing serious happens, just all the EBs are being processes which may take
# a long time to finaly find that the, e.g., region parameter is wrongly
# selected.

if 'applyonly' not in globals(): applyonly = True

#============================================================================

#--- What band are we observing (3, 5, 6, or 7)?

# Observing Band not yet determined.
oBand=-1

# Use the (always present) 1st EB to estimate the LO frequency
asdm0 = asdm_commonPart_ID
msc0 = asdm0+'.ms.split.cal'

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
    imsize=[400,400]
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
    imsize=[720,720]
    convsupport = 6
    qs_level = 6300.0 # For immath

    # _SADA_: Shift the circle a bit off the center to avoid the bright areas
    # from the statistics, if needed.
    region='circle[[360pix,360pix],72pix]' #imstat

elif (oBand==6):
    
    # Imaging parameters for Band 6
    cell=['3.0arcsec','3.0arcsec']
    imsize=[800,800]
    convsupport = 6
    qs_level = 5900.0 # For immath

    # _SADA_: Shift the circle a bit off the center to avoid the bright areas
    # from the statistics, if needed.
    region='circle[[400pix,400pix],80pix]' #imstat

else:
    # Imaging parameters for Band 7
    cell=['2.8arcsec','2.8arcsec']
    imsize=[860,860]
    convsupport = 6
    qs_level = 5500.0 # For immath

    # _SADA_: Shift the circle a bit off the center to avoid the bright areas
    # from the statistics, if needed.
    region='circle[[430pix,430pix],86pix]' #imstat


#--- Common parameters (not Band-specific).

# Frequency for selected spw in MHz as integer
imFreq=int(round(obsFreqs[spw]*1.0e-6))
   
#--- Loop over calibrated EBs, do imaging for all in the list

if (SecondScan == 'none'):
    rgNum = 1
else: 
    rgNum = 2

for imgNum in range(0, rgNum):

    asdm=asdm_commonPart_ID
    msc=asdm+'.ms.split.cal'

    if (imgNum == 0):
        scanNum = '6'
    else: 
        scanNum= SecondScan

    print('Now processing: '+msc)

    # Construct the image name proper for the ALMA archive
    # The SD naming convention since Cy7 (2019/20) has the following pattern: 
    # <member | group>.<m/gousID>.<sourcename>[-num#]_<intent>.
    #    <freq_info>.<obspatt>.<arrays>.<imagetype>.
    #    <pol>.<modifier>.<datatype>.<format>
    # Example:
    # member.uid___A001_X1465_X1363.Sun-num2_sci.spw1_245000MHz.sd.TP.cont.
    #    I.manual.sd.fits
    # The <format> section, i.e. the FITS-file suffix '.fits', is appended
    # later, at exportfits().

    sdImgName='member.'+mousID+'.fullSun-num'+str(imgNum)+'_sci.spw'+str(spw)
    sdImgName+='_'+str(imFreq)+'MHz.sd.TP.cont.I.manual.sd'

    print("About to create an image "+sdImgName)

    # Do SD imaging
    os.system('rm -rf '+sdImgName+'*')
    tsdimaging(infiles = msc,
               spw = str(spw),
               antenna =antenna,
               phasecenter = 'Sun',
               scan = scanNum,
               field='1',
               nchan = 1,
               cell = cell,
               outfile = sdImgName,
               imsize=imsize,
               gridfunction = 'SF',
               convsupport = convsupport,
               intent ='*ON_SOURCE*',
               stokes = 'I',
               brightnessunit='K',
               specmode='cube',
               outframe='lsrk')

    # Final image to be brightness temparature map.
    # Rescaling to Kelvins (See White et al. 2017)
    os.system('rm -rf '+sdImgName+'.rescl')
    res=imstat(sdImgName+'.image',region=region)

    immath(imagename=sdImgName+'.image',
           expr = 'IM0*'+str(qs_level/res['mean'][0]),
           outfile = sdImgName+'.rescl')

    # Export an image to the FITS file, append <format> section to product name.
    sdFITSname=sdImgName+'.fits'
    exportfits(imagename=sdImgName+'.rescl',fitsimage=sdFITSname)

    # Increment the image counter
    imgNum+=1

    print('ScanNum: '+str(scanNum)+' imaging finished.')
    print('')

print('SB '+sbName+' imaging completed.')
