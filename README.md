LenSNpop
======

LenSNpop uses the SCOTCH catalogues (Lokken et al. 2023) and LensPop (Collet et al. 2015) to simulate lensed transients and their host galaxies, as described in Sainz de Murieta et al. 2024

Before starting
------
* Download the SCOTCH catalogue file (*scotch_Z3_1.1.hdf5*) from [here](https://zenodo.org/records/7563623#.Y8_2my9h2yA) and save it in /data/SourcePop
* Download *redshift_splines.pkl*, *lenspop_splines.pkl* and the folder /2dpdfs/ from [LensPop](https://github.com/tcollett/LensPop) and save it in /data/LensPop
* Download the PLASTICC templates you are interested in from [here](https://zenodo.org/records/2612896) and save them in /data/Templates
* Due to GitHub's limits on file size, the folder */Simulation_catalogues* and the file *lsstobs.csv* are not complete. If you want the full sample, feel free to email me at [ana.sainz-de-murieta@port.ac.uk](ana.sainz-de-murieta@port.ac.uk)
  
How to simulate glSNe + hosts
-----
* host_transient.py : generate small catalogue from SCOTCH
* lenspop.py : generate lens population
* combine_lens_and_source.py : combine output from *host_transient.py* with output from *lenspop.py* to generate N lensing systems per source
* lensSN.py : Compute positions for SNe, create catalogue for simsurvey
* simsurvey_glSNIa.py : Simulate catalogue of SNe using LSST observing strategy (uploaded for glSNe Ia), modifying sncosmo Source allows to simulate other subtypes
* galaxy_galaxy_lensing.py : Lens host galaxy of glSN and compute magnification
