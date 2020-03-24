--- WrLINE: A script to extract helical axis and calculate 'writhe' from an AMBER trajectory of circular DNA ---

-- REQUIREMENTS --
1) AMBERTOOLS with the PTRAJ module
2) PYTHON 2.6 or 2.7 with the NUMPY library installed

-- INPUT --
1) An AMBER .mdcrd trajectory of a circular DNA (you could try it on the linear DNA but the safety is not guaranteed!)
2) The corresponding AMBER .prmtop parameter file

-- USAGE --
python WrLINE.py [jobname] [prmtop] [mdcrd] [n basepairs] [n timesteps]

-- EXAMPLE --
Try running this command line to see if the script works
>> python WrLINE.py test test.prmtop test.mdcrd 336 8

-- OUTPUT --
1)** jobname/C1.xyz: xyz trajectory of the helical axis 
2)** jobname/writhe.ser: the calculated writhe time series
and
3) jobname/tw.ser: 2D array of twist at each basepair and time step 
4) jobname/sinreg.ser: 2D array of bending register angle at each basepair and time step 

-- TRANSFORM FILES TO AMBER FORMAT --

On Amber folder there is all scripts needed to adapt WrLINE helical axis trajectory to AMBER format

1) 1st transform trajectory of the helical axis from xyz format to AMBER mdcrd. We use vmd but you can use other programs
>> vmd -e xtctomdcrd.tcl
2) Obtain the related AMBER topology
>> cpptraj -i ptraj.in

Trick: change 1st line of axis trajectory in AMBER format for further analysis in ptraj

