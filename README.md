# Solar-Meridional-circulation
Code done to measure meridional flow of the sun from SOLIS synoptic data

## Computational Physics Research Portfolio

This repository contains codes used to calculate differential rotation and meridional circulation of the sun through solar cycle 23 done under the supervision of  of Prof. Bidya Binay Karak at IIT BHU.

## Project Summary
### Status: Ongoing

### Objective: To measure and analyse the meridional flow of the sun over cycle 23 and study the patterns over the rising and falling phases.

### Methodology & Implementations:-
The project employs image cross-correlation to estimate the displacement of magnetic structures between successive magnetograms. This method uses the collective magnetic-field pattern rather than tracking individual features, which may emerge, evolve, merge, or disappear.

For differential rotation, longitudinal shifts are measured across latitude strips after accounting for the reference rotation of the Carrington coordinate system. These shifts are converted into angular rotation rates.

For meridional circulation, the analysis measures latitudinal shifts after correcting for differential rotation. The angular displacement is converted into a physical velocity using the solar radius and the time interval between observations.

Measurements from accepted magnetogram pairs are averaged at each latitude, and their standard deviations are used to estimate the scatter in the results. The analysis is also repeated for the complete dataset and the rising and declining phases of Solar Cycle 23.

### Data
The study uses photospheric line-of-sight magnetic-field maps obtained from the Vector Spectromagnetograph (VSM) instrument at the Synoptic Optical Long-term Investigations of the Sun (SOLIS) observatory. The magnetograms are stored as FITS files and represented as Carrington synoptic maps using a cylindrical equal-area (CEA) projection, with a resolution of (360 times 180) pixels.

Successive maps, separated by approximately one Carrington rotation, are paired to track the displacement of magnetic structures. The analysis focuses on Solar Cycle 23, from 1996 to 2008, and uses latitudes between -50 degrees and +50 degrees, where the measurements are considered more reliable. Maps containing missing data, insufficient magnetic structure, or poor correlations are excluded.

## Next steps
*   Higher Resolution Cross-Examination: Transitioning the pipeline to ingest space-based **SDO/HMI** data to exploit higher spatial resolution and instrument quality.
*   Cadence Scaling: Planning to ingest full-disk magnetogram data for a 1-month test period from the **JSOC** database with a high-frequency 6-hour cadence.

##  Languages & Libraries Used
- Python 3
- astropy,sunpy
- matplotlib (for plotting)
- glob, csv (for handling data)
- numpy, datetime, Scipy (for dataprocessing)
