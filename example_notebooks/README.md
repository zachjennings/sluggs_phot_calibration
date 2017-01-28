## Example Work Notebooks:

Currently incomplete, more notebooks will be added going forward.

Here are a collection of example work notebooks. These notebooks correspond to the calibrations necessary to go from three reduced images in different filters
to a full 3-filter catalog. When starting a new galaxy, just copy all these notebooks into a new directory. The user will have to make a copy of each notebook corresponding
to each filter, i.e. usually the user will have three copies of each notebook for a typical 3-filter photometric catalog (the execption being the final catalog notebook).

All current notebooks are taken from the examples on the g-filter photometry on NGC 3115, so there will be situations where certain values, files, etc. used for NGC 3115 will need to be changed for other galaxies and filters. 

### Current Notebooks:
aperture_photometry_calibrations_example.ipynb: Example notebook demonstrating how to select out correct aperture corrections and generate the initial list of photometry. Copy this over to working directory for future photometry work, recommend naming convention of: galaxy_filter_aperture_calibrations.ipynb, e.g. n3115_i_aperture_calibrations.ipynb

### Current Files:
apertures.param: Source extractor parameter file listing out correct output parameters for the apertures catalog. Also can be used for the completeness catalog generation, although with warnings for too many apertures listed.

example_completeness.seconfig: Source extractor config file for completeless corrections. Appropriate values will need to be modified.

example_aper.seconfig: Source extractor config file for aperture corrections. Values will need to be modified for appropriate catalog names, aperture values, etc. 

panstarrs_casjobs_query.sql: Sample sql query for PanSTARRS casjobs server, used for zero-pointing images. Coordinates need to be modified for particular galaxies.



