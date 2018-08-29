# synth_spect.py

Python routine that generates a synthetic sub-mm/mm rotational spectrum from selected species using the CDMS/JPL databases and, optionally, compares it with an observed spectrum.

This python script reads online spectroscopic data through the JPL and CDMS databases. The synth_spect.py routine needs 2 to 3 input ASCII files: 1) input.in for the list of input parameters; 2) the file listing the properties of molecules needed to generate the synthetic spectrum; 3) optionally, the file describing the observed spectrum. 

## input.in

The input parameters needed by synth_spect.py are the following:

- choice_obs: 'yes' or 'no' for the use of an observed spectrum. If yes, specify fileobs
- filemod: Name of input model file 
- fileobs: Name of input observed file 
- choice_y: Unit of intensity value ('Tpeak' for K, 'Fpeak' for mJy)
- prefix: Prefix of output files
- beamsize: Beam size [arcsec]
- choice_tau: Include opacity in calculations (yes or no)
- numin: Minimum frequency [GHz]
- numax: Maximum frequency [GHz]
- dnu: Spectral resolution [MHz]
- Vmin: Minimal value for selected transitions (K or mJy depending on choice_y)
- rms: rms noise of synthetic spectrum (K or mJy depending on choice_y)
- choice_plot: Output for the plot: 1) pdf, 2) python window


## modeled spectrum file

The modeled spectrum file contains the list of species to include in the synthetic spectrum together with their properties: 

- total column density Ntot
- excitation temperature Tex
- source size ss
- offset velocity relative to the source velocity v_off
- the FWHM width of transitions width
- the spectroscopic database (jpl or cdms)
- vibrational correction factor
- reference (optional) of the properties 

Please respect the number of characters as specified in the template. Lines starting with a ! or a # will be considered as comments. The file ends with a *.


## observed spectrum file

The observed spectrum file is a 2-column ASCII file giving the frequency (in MHz) in the 1st column and the intensity in the 2nd one, with the unit depending on choice_y specified in input.in.


## How to run the script

Install python3 and the usual packages numpy, pandas, matplotlib, and scipy. Once python3 is install, type: `python3 synth_spect.py`