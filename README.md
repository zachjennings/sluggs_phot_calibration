# SLUGGS Survey Photometric Calibration
This repository is storage for SLUGGS survey photometric calibration work. All work from initial creation of raw catalogs for each image to final all-filter calibrated photometric catalogs will be stored here. 

## Directory Structure
Each pointing has it's own directory. Some galaxies have been imaged with multiple pointings and therefore will have more than one directory. Each directory will contain all photometric work for all filters at that particular pointing. 

To do photometric work on a new pointing, create a new directory and copy over the example "cookbooks" from the main directory. Each filter will get its own version of the calibration notebooks. The idea behind these cookbooks is not a "one-size-fits-all" approach. Rather these give the appropriate steps that a user can iterate on to determine the correct calibrations to make final photometric catalogs for each image. However, each image may have its own quirks to account for.

## General Workflow:
Calibrating photometry on an image has a few general steps. Use the "master notebook" to keep track of all the calibration files. Also keep track of all calibration values in the google doc [located here](https://docs.google.com/spreadsheets/d/1ONeyCJNgF9Db3GR8uYT5LBEJpJ3ylQE8tbSavvmmz_U/edit?usp=sharing) (ask Zach for edit access).

#### 1. Initial Catalog Creation and Aperture Corrections:
The first notebook makes the initial photometric catalog, then calculates aperture corrections. We are performing aperture photometry on our sources, i.e. summing all the light within a fixed aperture around the source. These measurements then need to be calibrated for the additional light from the source that falls outside the specified aperture. We do this by using bright, unsaturated point-sources in the image to determine what fraction of light from a point source falls outside the measurement aperture. These corrections are of course invalid for objects that are not point-like (i.e. galaxies and partially-resolved clusters) as their light distributions will not match that of the bright point-like stars we have measured the correction on.

The two values to be determined here are: 1. The photometric aperture that maximizes the S/N, and 2. The correction from this aperture to the total flux from the point source.

Follow the steps in the aperture photometry cookbook to make initial SExtractor files, then use the remaining steps to select the correct measurement aperture and calibrate the aperture correction. These steps need to be performed for each image, and the corrections are different for each filter.

#### 2. Calibrating Zero-Points:
After we have calibrated our aperture photometry, we need to shift our measurements on to the correct photometric systems. While we previously used SDSS photometry for this purpose, some of our survey targets don't overlap with SDSS and therefore need to be calibrated using other fields throughout the night. Now that Pan-STARRS has finalized its first photometric data release, we will be using it for photometric calibration since it has similar sky coverage to our survey's photometry (all data are taken from Hawaii).

Zero-pointing is done by measuring bright, unsaturated stars in our field and comparing those to published Pan-STARRS photometry. Currently Pan-STARRS photometry needs to be pulled down from the CasJobs site [located here](http://mastweb.stsci.edu/ps1casjobs/). 


#### 3. Calculating Completeness Curves:
Our ability to detect sources decreases as the objects get fainter and fainter, approaching 100% at bright magnitudes and approaching 0% at faint magnitudes. To correctly perform our inference on our photometric populations, we need to quantify this completeness. We do this using fake star tests. We insert thousands of fake sources at random locations in our image, then try to recover them with the same methods that we use to create our initial catalogs. Our rate of successful recovery as a function of magnitude is how we quantify our completeness. 

The notebook.

#### 4. Matching catalogs across filters.
Once zero-points and aperture corrections have been performed, we are free to match catalogs across our filters. Just follow the final steps in the example notebook, and complete the final photometric catalog for the galaxy. At this point we are free to feed our final catalogs to our GC inference software.