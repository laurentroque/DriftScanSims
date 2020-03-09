# REACH Drift Scan Simulations
# Import stuff
import numpy
from scipy.io import loadmat as load
import time

# Time and Frequency Range
obsFreq = numpy.arange(50, 151)*1E+6 # Observation frequency in Hz
LST = numpy.arange(0, 360, 2.5) # LST in degrees (2.5deg = 10min resolution)
ploton = 1 # 1==plot sky (slower), 0==no plot
cf = 1 # 1==calculate chromaticity function, 0==do not calculate
fref = 50E+6 # Reference frequency to calculate chromaticity function
rot = 0 # Antenna rotation in degrees

# Sky maps and antenna beam directories in l-m plane
pdir = 'npy_beam_dipoleandbalun' # Beam directory
mdir = 'ta_REACHmaps' # Map directory

# Load sky maps (LFSM) 40 - 250MHz
maps = load('maps_LFSM_40Mto250M_3deg.mat') # High resolution LFSM maps
maps = maps['maps']
mapf_index = numpy.arange(40, 250)*1E+6 # Map frequency index (Shouldn't have to change)

# Initial values
d2r = numpy.pi / 180
r2d = 1.0 / d2r
res = 512 # Resolution for all maps and beams (512x512), 512==0.35deg
ra = numpy.linspace(360, 0, maps.shape[1]) # Right Ascension in degrees
dec = numpy.linspace(-90, 90, maps.shape[0]) # Declination in degrees
RAi, DECi = numpy.meshgrid(ra, dec)

# Coordinates of REACH 5 site
Lat = -30.836511
Long = 21.370872

Tsky = numpy.zeros((len(obsFreq), len(LST)))

if cf == 1:
  CF = numpy.zeros((len(obsFreq), len(LST)))

# Frequency index
for fidx in range(0, len(obsFreq)):
  # Load antenna beams
  beam_name = pdir + '/reach_hexdip_pat_' + str(int(obsFreq[fidx]/1E+6)) + '.mat'
  Pat = load(beam_name)['Pat']
  Pat[Pat < 0.001] = 0 # set floor to zero (outside of l-m beam)

  if cf == 1: # Load reference beam
    pp = pdir + '/reach_hexdip_pat_' + str(int(fref/1E+6)) + '.mat'
    beamref = load(pp)['Pat']
    beamref[beamref < 0.001] = 0 # Set floor to zero (outside of l-m beam)

  # Time index
  for tidx in range (0, len(LST)):
    tic = time.time()
    LSTi = LST[tidx]
    
    # Sky map (UV) - needs to generate and load first time around
    # Creating strings of paths to maps for loading
    if str(numpy.around(LST[tidx]*1000)/1000).endswith('.0'):
      LStime = str(int(numpy.around(LST[tidx]*1000)/1000))
    else:
      LStime = str(numpy.around(LST[tidx]*1000)/1000)
    map_name = mdir + '/sky_map_uv_f' + str(int(numpy.around(10*obsFreq[fidx]/1E+6/10))) + 'MHz_LST' + LStime + 'deg.mat'
