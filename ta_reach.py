# REACH Drift Scan Simulations
# Import stuff
import numpy

# Time and Frequency Range
obsFreq = numpy.arange(50, 151)*1E+6 # Observation frequency in Hz
LST = numpy.arange(0, 360, 2.5) # LST in degrees (2.5deg = 10min resolution)
ploton = 1 # 1==plot sky (slower), 0==no plot
cf = 1 # 1==calculate chromaticity function, 0==do not calculate
fref = 50E+6 # Reference frequency to calculate chromaticity function
rot = 0 # Antenna rotation in degrees

# Sky maps and antenna beam directories in l-m plane
pdir = REACHDriftScan.npy_beam_dipoleandbalun # Beam directory
mdir = REACHDriftScan.ta_REACHmaps # Map directory
