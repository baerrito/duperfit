# Duperfit (Duplicated Superfit in Python) version 0.1

Spectral classification software, adapting Superfit by Dr. D. Andrew Howell into Python, implementing Python's object-oriented capabilities with aims for a more generally accessible user experience. Utilizes least-squares objective fitting of supernova templates plus galaxy flux and relative reddening via a Cardelli law.

## System Requirements

At present, Duperfit is only tested to work on Linux systems. Mac OS X and Windows testing to come down the road. Duperfit is written in Python 3.8.8, and utilizes the following packages:

 - Astropy 4.3.1
 - extinction 0.4.6
 - NumPy 1.21.5
 - PyDL 0.7.0
 - SciPy 1.7.3

## Installation

Installation is as simple as cloning the repository to your machine. You will also need to set up your OS environment by setting the variable `DFDIR` to your installation directory, e.g. `export DFDIR=/home/mbaer/duperfit` with bash on my machine. One may add such a line to their `.bashrc` file if they so desire.

The `picklejar` directory contains pickle files with all the templates for efficient processing -- i.e., all data is loaded in one object rather than loaded on each iteration. However, these files may or may not work correctly when pulled to your machine. If this occurs, run `tempsetup.py` to set up these files for running the program. The old files will be overwritten, and should work on your machine once this process completes. Note that this script can also be used to add new spectra to your template libraries relatively smoothly.

## Execution

Duperfit can be run a couple of different ways; with the main function, with a GUI, with a Class, or with a provided script. The latter 3 options will run, but are not yet finished with development; expect major updates to come to those in the next year or so.

### Main Function

This is the least user-friendly option, but the most flexible. The function with a full preamble containing what input is expected is provided in the module `duperfit.py`.

### GUI

This is being developed as the most user-friendly option, but may have the least flexibility. This has not been bug-tested since adding the environment variable, so do be warned if it presently does not run as expected.

### Class/script

The Class and script utilize the same input structure, as the script utilizes the Class to run. The provided `params.json` file is structured the same way as a `params` dictionary would be expected to follow in the Class. This structure is as follows:

 - `"IO"`: Handles input/output
   - `"object_dir"`: Directory containing input spectrum
   - `"object_file"`: Input spectrum file
   - `"SN_templates"`: Choice of SN template selections. Options are generically the names of a Pickle file without the extension, barring "All SNe" and "SNe <=10d", which the program can interpret. Default Pickle files are included in the `picklejar` archive.
   - `"user_SN_template"`: Full path to a user-provided pickle file containing SN templates to evaluate (CURRENTLY UNTESTED)
   - `"gal_templates"`: List of galaxy templates to evaluate
   - `"user_gal_template"`: Full path to a user-provided ASCII file containing a galaxy template to evaluate (CURRENTLY UNTESTED)
   - `"output_path"`: Directory for saving the output file
   - `"output_file"`: Name of the output file, defaults to `"object_file"` with the input extension replaced with `".dfo"`
 - `"fit_params"`: Adjusts fitting parameters
   - `"use_exact_z"`: Set `true` to use `"z_min"` as an exact redshift
   - `"z_min"`: Minimum redshift
   - `"z_max"`: Maximum redshift
   - `"delta_z"`: Redshift grid size
   - `"Av_min"`: Minimum A_V
   - `"Av_max"`: Maximum A_V
   - `"Rv"`: Relative extinction
   - `"max_template_scale"`: Upper bound for SN scaling
   - `"max_galaxy_scale"`: Upper bound for galaxy scaling
 - `"fit_weight"`: Handles the weights of the fit
   - `"weight_source"`: Source of weigths for fit  - recommend leaving at either `"incl"` (included) or `"uw"` (unweighted); `"tell"` still should be tested
   - `"estimate_error"`: Set `true` to estimate flux uncertainties with an iterative B-spline fit
 - `"sigma_clipping"`: Handles options for sigma-clipping of the input spectrum
   - `"sigma_clip"`: Set `true` to sigma-clip the signal (recommended)
   - `"sigma_source"`: Source of sigma; can either be the flux errors (`"incl"`) or calculated (`"calc"`)
   - `"n_iterations"`: Number of iterations
   - `"n_grow"`: Grow parameter
   - `"n_sigma"`: Multiple of sigma past which to clip
 - `"wavelength_range"`: Handles the wavelength range and required coverage for the spectrum
   - `"min_wavelength"`: Minimum wavelength for fitting in Angstroms
   - `"max_wavelength"`: Maximum wavelength for fitting in Angstroms
   - `"wavelength_bin"`: Wavelength bin window width in Angstroms
   - `"minimum_wavelength_fraction"`: Minimum fractional wavelength coverage for fitting F, where 0 < F <= 1 (>=0.7 recommended)
 - `"options"`: Additional options
   - `"silence_messages"`: Set `true` to silence warning messages during run (setting `false` is primarily for debugging)
   - `"save_output"`: Set `true` to save the output file
   - `"optimizer`": Choice of optimization algorithm, recommended to use L-BFGS-B, as my tests have shown it to return the same output as TRF in shorter time. This option may get updated or removed in future versions.
