# ALMA Data Reduction Script for CASA5
#
# Original template by Masumi Shimojo (EA ARC)
# Modified by Miro Barta <barta@asu.cas.cz> (EU ARC, Czech node) on 2018/07/14
# Modified by Masumi Shimojo <masumi.shimojo@nao.ac.jp> (EA ARC) on 2020/06/24
# Generalisation for more bands, barta@asu.cas.cz, on 2021.03.14
# Modified for Casa6 (Python3) by Masumi Shimojo on 2023.02.16
#          for Regional Mapping data
#
# Recommended CASA version: 6.0 and higher
#-------------------------------------------

#Switch for each process
thesteps = []
step_title = {0: 'Import of the ASDM',
              1: 'Create the files of the observing infromation',
              2: 'A priori flagging1',
              3: 'Generation of the Tsys cal table',
              4: 'Do initial flagging',
              5: 'Calibration of the data into Kelvins',
              6: 'Correction of the forward scattering & spillover'}

# The Python variable 'mysteps' will control which steps
# are executed when you start the script using
#   execfile('scriptForSDCalibration.py')
# e.g. setting
#   mysteps = [2,3,4]# before starting the script will make the script execute
# only steps 2, 3, and 4
# Setting mysteps = [] will make it execute all steps.

#==================================================
#--- ASDM Directory Name

# Compose the full ASDM name for this EB.
asdm='uid___A002_X117bb05_Xa2d1'

#SPW for solar observations
spwSun = '0,1,2,3'

#==================================================
#Check CASA ver.

from casatools import ctsys
version = ctsys.version_string()
print("You are using " + version)
if (version < '6.0.0'):
    print("YOUR VERSION OF CASA IS TOO OLD FOR THIS SCRIPT.")
    print("PLEASE UPDATE IT TO 6.0 OR HIGHER BEFORE PROCEEDING.")

#==================================================

# Since Cycle 5 / CASA 5.3 we _DO_ have scriptForPI.py for the SD data, too.
# It cares about setting applyonly=True.

if 'applyonly' not in globals(): applyonly = False

try:
  print('List of steps to be executed ...', mysteps)
  thesteps = mysteps
except:
  print('global variable mysteps not set.')

if (thesteps==[]):
  thesteps = range(0,len(step_title))
  print('Executing all steps: ', thesteps)


#==================================================
#Parameters Setup

# Measurement Set names
mso = asdm + '.ms'
msc = mso + '.split.cal'

#--- Import libraries
import re
import os
import csv
import numpy as np

#=========================================================================
# Import of the ASDM
mystep = 0

if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  if os.path.exists(mso) == False:
      importasdm(asdm = asdm, vis = mso, lazy = True, bdfflags=False,
                 with_pointing_correction=True,
                 asis='Antenna Station Receiver Source CalAtmosphere CalWVR CorrelatorMode SBSummary CalDevice')

#=========================================================================
#   This shall be done _always_, regardless of the 'mysteps' selection!

#--- What band are we observing (3, 5, 6, or 7)?

# Observing Band not yet determined.
oBand=-1

if os.path.exists(mso) == True:

    tb.open(mso+'/SPECTRAL_WINDOW', nomodify=True)
    # Ref. frequencies for all SPWs
    obsFreqs = tb.getcol("REF_FREQUENCY")
    tb.close()
    # Average of ref. frequencies for science SPWs 0,1,2, and 3 = approximate
    # LO frequency, rescaled to GHz.
    freqLO=0.25*(obsFreqs[0]+obsFreqs[1]+obsFreqs[2]+obsFreqs[3])*1.0e-9

    # Fixed LO freqs. (in GHz) of Solar ObsMode and corresponding Bands:
    # A dictionary
    solarLOs={100.0:3, 194.0:5, 239.0:6, 346.6:7}

    #-- Try to identify the observing band
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
        sys.exit("Error: Observing frequency not allowed in Solar ObsMode.")
    
# ASDM has not yet been imported - we need to run 'Step 0' first
else:
  print("Observing band and pointing parameters will be set after ASDM import.")
  print("Run the Step 0 first!")
  sys.exit('Exiting.')



#=========================================================================
# Create the files of the observing infromation
mystep = 1

if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print('Step ', mystep, step_title[mystep])

    # Dump of the observing information
    #
    # NB: So far, as of 2018-08-06, the correct file suffix for SD data is
    # really '.listobs'. The suffix 'listobs.txt' consistent with the INT data
    # processing is NOT accepetd by the packager

    '''
    # To be - perhaps - used later, after packager update
    os.system('rm -rf '+ mso + '.listobs.txt')
    listobs(mso, listfile = mso + '.listobs.txt')
    '''

    os.system('rm -rf '+ mso + '.listobs')
    listobs(mso, listfile = mso + '.listobs')

#=========================================================================
# A priori flagging1
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

#=========================================================================
# Generation of the Tsys cal table

mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  os.system('rm -rf ' + mso + '.tsys')
  gencal(vis = mso, caltable = mso + '.tsys', caltype = 'tsys')

  os.system('rm -rf '+ mso + '.tsys.plots')
  os.system('mkdir '+ mso + '.tsys.plots')
  plotms(vis=mso+'.tsys',xaxis='freq',yaxis='tsys', gridrows=2, gridcols=2,
          iteraxis='antenna', showgui=False,
          plotfile=mso + '.tsys.plots/'+mso+'.tsys.plots.png')


#=========================================================================
# Do initial flagging
mystep = 4

if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

#=========================================================================
# Calibration of the data into Kelvins

mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

#  sdcal(mso, calmode='ps', spw='0,1,2,3', fluxunit='K',
#        outfile=msc, outform='MS2')

  sdcal(infile=mso, calmode = 'ps,tsys,apply', spw=spwSun)

  os.system('rm -rf '+msc+'*')
  os.system('rm -rf '+msc+'*.flagversions')
  split(vis = mso,
        spw = spwSun,
        outputvis = msc,
        datacolumn = 'corrected',
        keepflags = True)

#=========================================================================
# Forward & Spillover efficiency correction
mystep = 6
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print('Step ', mystep, step_title[mystep])

    # --- Correction factors (White et al, Sol Phys, 2017)

    if(oBand == 3):
      ampScale2=0.893  # for Band 3: 0.893
    elif (oBand == 6):
      ampScale2=1.0    # for Band 5: 1.0 - not known yet
    elif (oBand == 6):
      ampScale2=0.862  # for Band 6: 0.862
    else:
      ampScale2=1.0    # for Band 7: 1.0 - not known yet


    # Generate the scaling table.
    gencal(vis=msc,
           caltable=msc+'.eff.tbl',
           caltype='amp', spw='',
           parameter=[np.sqrt(ampScale2)])

    # Apply scaling.
    applycal(vis=msc,
             gaintable=msc+'.eff.tbl')
