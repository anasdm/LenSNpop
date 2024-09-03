LenSNpop
======

LenSNpop uses the SCOTCH catalogues (Lokken et al. 2023) and LensPop (Collet et al. 2015) to simulate lensed transients and their host galaxies.

Before starting
------
* Download the SCOTCH catalogue file (*scotch_Z3_1.1.hdf5*) from [here](https://zenodo.org/records/7563623#.Y8_2my9h2yA) and save it in /data/SourcePop
* Download *redshift_splines.pkl*, *lenspop_splines.pkl* and the folder /2dpdfs/ from [LensPop](https://github.com/tcollett/LensPop) and save it in /data/LensPop
* 
How to simulate your own glSNe
-----
* host_transient.py : generate small catalogue from SCOTCH
* lenspop.py : generate lens population
* combine_and_lens.py : combine output from *host_transient.py* with output from *lenspop.py* to generate N lensing systems per source
