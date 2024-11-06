# World Ocean Model of Biogeochemistry and Trophic-dynamics (WOMBAT)


     (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.
     / o o \  : :.-.: :: ,. :: `' :: .; :: .; :`-. .-'
    (   "   ) : :: :: :: :: :: .. ::   .':    :  : :  
     \__ __/  : `' `' ;: :; :: :; :: .; :: :: :  : :  
               `.,`.,' `.__.':_;:_;:___.':_;:_;  :_;     

Welcome to the development area of WOMBAT. Here I am developing both WOMBAT-lite (bio_v3.inc) and WOMBAT-mid (bio_v4.inc) versions.

**Note: WOMBAT-lite (bio_v3.inc) is currently the only active version**

The World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT) is based on a NPZD (nutrient–phytoplankton–zooplankton–detritus) model.
The "lite" version of WOMBAT includes one class each of phytoplankton, zooplankton and sinking detritus, as well as nitrate (NO3), bio-available iron (Fe), 
dissolved inorganic carbon (DIC), calcium carbonate (CaCO3), alkalinity (ALK), and oxygen (O2). Fe is carried through the zooplankton and detrital pools as well.
Gas exchange follows OCMIP2 protocols, but in the very near future will be updated in line with MOCSY.

Version "mid" is under beta development and involves 3 phytoplankton, 2 zooplankton and 3 detrital pools.

## Getting up and running (ACCESS-OM2-BGC)

These versions of WOMBAT are currently designed to be run in the ACCESS-OM2-BGC configuration at 1 deg horizontal resolution. 
Follow these instructions to get it running on Gadi (Australia's supercomputer with the National Computing Infrastructure).

First, visit https://github.com/COSIMA/access-om2 and follow the instructions to download the source code for ACCESS-OM2-BGC. This involves:

1. join the "vk83" project on myNCI and load in payu using _module use /g/data/vk83/modules; module load payu_
2. git cloning the src code
3. git checkout of the bgc branch
4. modifying the install.sh script to enable the BGC component
5. run the _install.sh_ command

Second, set up an appropriate run directory to run your model. This involves:

5. git cloning the run directory from https://github.com/COSIMA/access-om2/wiki/Getting-started#quick-start 
6. following the instructions there to checkout the BGC configurations
7. modifying the config.yaml file to point towards your newly compiled executable that you compiled in step 4

Third, check that the default model runs (no changes to the code yet!)

8. go to your run directory, and type _payu sweep; payu setup_
9. check in the newly create _work_ directory that all the necessary files and executables are in there
10. type _payu run -f_ to run the model

Fourth, now we set up the new WOMBAT.

11. copy bio_v3.inc, bio_v4.inc, csiro_bgc.F90 and ocmip2_co2calc.F90 from this repository into your source code under the /access-om2/src/mom/src/mom5/ocean_csiro_bgc/ directory
12. _install.sh_
13. make a new run directory
14. download the input_files by emailing me (see bottom) and ocean_files included here. Put the input_files in your home directory. Put the ocean_files in your rundirectory/ocean/ directory.
15. modify the config.yaml file in the run directory to point to your new executable (step 12) and the WOMBAT-lite inputs that were placed somewhere in your home directory (step 14)
16. _payu sweep; payu setup_
17. check that all files are present in the newly created _work_ directory
18. _payu run -f_

And Voilà. You should be up and running with WOMBAT-lite in the ACCESS-OM2-BGC global 1 degree resolution configuration.

Any questions or concerns or suggestions, please contact Pearse.Buchanan@csiro.au.
