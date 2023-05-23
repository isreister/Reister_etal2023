# Reister_etal2023
This is all the code used for the calculation made in Reister et al 2023 paper 'The dispersal of low salinity over the northeastern Gulf of Alaska shelf with a focus on the Copper River plume' 
Code files in this repo are:

Faster_Globec_and_NGA_LTER_stero.m %%%%%% this code combines data from the historical Seward Line with the more recent NGA LTER seward line and cacluates the freshwater height for the water colulmn for all seasons and all years. Practical Salinity is used in this calculation. Is this OK because it's just a ratio? Note that the Seward Line data here goes from the start of GLOBEC 1997 all the way to NGA LTER 2020 and includes the Gulf watch alaska years (2012-2020) and some other years that go from 2005-2012 (who was that?)

LTERcruisedata_makerV2.m %%%%%% a matlab script where pathnames for L2 or L3. See the "Description_of_ProcessingFiles_for_NGA_LTER_station_data" for more details. L2 or L3 are then run through deep_parseLTERV2 (see below). LTERcruisedata_makerV2 takes the parsed structure and makes a larger structure type dataset that is a compliation of all the NGA LTER cruises. 

deep_parselLTERV2.m %%%%%% is a function that takes in an table (produced by LTERcruise_makerV2.m) and identifies data columns for a desired set of variables that I need by matching names. The columns are rearranged into a standard format. Date time is formalized to be a datenum. Practical Salinity and Absolute Salinity are calculated. There is some small amount of QCing that removes casts that regerister no value anywhere in the upper 300 meters. Casts of this type are removed since they are not real casts. Freshwater height is also cacluated here. Note that at this time freshwater height is calculated with Practical Salinity.


