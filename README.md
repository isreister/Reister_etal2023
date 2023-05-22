# Reister_etal2023
This is all the code used for the calculation made in Reister et al 2023 paper 'The dispersal of low salinity over the northeastern Gulf of Alaska shelf with a focus on the Copper River plume' 
Code files in this repo are:
Faster_Globec_and_NGA_LTER_stero.m %%%%%% this code combines data from the historical Seward Line with the more recent NGA LTER seward line and cacluates the freshwater height for the water colulmn for all seasons and all years.

LTERcruisedata_makerV2.m %%%%%%takes in the L2 or L3 (if we can get it) or L1 if that's all we have. Here is a description of L1 and L2:

L1 via SBE data processing:
1) Conversion
2) Filter
3) Align
4) Thermal Mass Correction
5) Derive to check
5) Seaplot to check descent rate, surface soak, and salinity spikes
7) Loop Edit
8) Wild Edit - added to the workflow last year because of bad termination on SKQ that caused a lot of noise
9) Derive
10) Copies of this result to Brita (and SUNA processing?)
11) Bin Averaging to 1 dbar bins
12) Bottle Summary

L2 via MATLAB and Python:
The matlab part of this is Liz's MATLAB version of Seths's code.
1)Something gets sucked into conductivity cell
2)“Jitter” from a bad connection included in the average values
3)Surface values on downcast very different from upcast
4)Waves pull the CTD up and down which throws of the Thermal Mass correction
    i)Usually worse at shallow depths 
    ii) Seen in several casts in row, typically before a gap in sampling when they run for cover.
5)Downcast is completely a mess but upcast is OK, so use that instead. (Make note in .hdr file)

Basically, this code checks to make sure downcasts are the appropirate choice. Its a hands on process. 


