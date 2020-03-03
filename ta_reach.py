# REACH Drift Scan Simulations
# Import stuff
import numpy
import scipy.io.loadmat as load

# Time and Frequency Range
obsFreq = numpy.arange(50, 151)*1E+6 # Observation frequency in Hz
LST = numpy.arange(0, 360, 2.5) # LST in degrees (2.5deg = 10min resolution)
ploton = 1 # 1==plot sky (slower), 0==no plot
cf = 1 # 1==calculate chromaticity function, 0==do not calculate
fref = 50E+6 # Reference frequency to calculate chromaticity function
rot = 0 # Antenna rotation in degrees

# Sky maps and antenna beam directories in l-m plane
pdir = DriftScanSims.npy_beam_dipoleandbalun # Beam directory
mdir = DriftScanSims.ta_REACHmaps # Map directory

# Load sky maps (LFSM) 40 - 250MHz
maps = load('maps_LFSM_40to250M_3deg.mat') # High resolution LFSM maps
mapf_index = numpy.arange(40, 250)*1E+6 # Map frequency index (Shouldn't have to change)

# Initial values
d2r = numpy.pi / 180
r2d = 1.0 / d2r
res = 512 # Resolution for all maps and beams (512x512), 512==0.35deg
ra = numpy.linspace(360, 0, maps.shape[1]) # Right Ascension in degrees
dec = numpy.linspace(-90, 90, maps.shape[0]) # Declination in degrees
RAi, DECi = numpy.meshgrid(ra, dec) 
^^^^ check this is correct
