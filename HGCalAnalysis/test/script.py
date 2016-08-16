# Kai Chang - Caltech CMS-CERN 2016
#
# Program runs a cmsRun program on a large scale. Outputs a series of ROOT
# files.
#
#
# Needs to have a chosen file selected, and the filesname variable untampered,
# Has to run on CMS environment (cmsenv)
# =============================================================================

import sys, os
import fileinput

## PARAMETERS YOU FILL IN ##
file_dir	= './deep-learning/data/neutrinos.list'
name		= 'pdg12_pt35'

# open and create neutrino list file
nfile = open(file_dir)
lines = nfile.readlines()
nfile.close()

for line in lines:
	f = open('particle_nosmear_calib.py', 'r')
	oldfile = f.read()
	f.close()

	newfile = oldfile.replace("eos/to/reco/root", line)

	if line[-9] != '_': #single digit 
		newfile = newfile.replace("rootname", name + "_" + str(line[-7]) + ".root")
	else:
		newfile = newfile.replace("rechit_unformatted.npy", name + "_" + str(line[-8:-6]) + ".root")

	f = open('particle_nosmear_calib.py', 'w')
	f.write(newfile)
	f.close()

	# runs the root processor
	os.system('cmsRun particle_nosmear_calib.py')

	# return to original format
	f = open('particle_nosmear_calib.py', 'w')
	f.write(oldfile)
	f.close()