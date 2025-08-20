# Reister_etal2023
This is all the code used for the calculation made in Reister et al 2023 paper 'The dispersal of low salinity over the northeastern Gulf of Alaska shelf with a focus on the Copper River plume' 
Code files in this repo are:

Faster_Globec_and_NGA_LTER_stero.m %%%%%% this code combines data from the historical Seward Line with the more recent NGA LTER seward line and cacluates the freshwater height for the water colulmn for all seasons and all years. Practical Salinity is used in this calculation. Note that the Seward Line data here goes from the start of GLOBEC 1997 all the way to NGA LTER 2020 and includes the Gulf watch alaska years (2012-2020) as well as cruises of opportunity (2005-2012).


HourlyDataProcessing_2002_2022_windanalysis.m %%%% This script is how wind is files from copernicus are averaged to form our estimates of wind direction and speed for the a region of interest at the mouth of the copper river. historgram roses are also made with this script.


LTERcruisedata_makerV2.m %%%%%% a matlab script where pathnames for L2 or L3. See the "Description_of_ProcessingFiles_for_NGA_LTER_station_data" for more details. L2 or L3 are then run through deep_parseLTERV2 (see below). LTERcruisedata_makerV2 takes the parsed structure and makes a larger structure type dataset that is a compliation of all the NGA LTER cruises. 


WMSsatretrieval.m %%%%% this is a function that is used within the code script 'mydates_satretrieval_v2.m'. Inputs for the function are a land mask (land), a list of dates (t), and an index value of 1,2,3,4,...length(t). The function uses an API to recover the cloud_fraction for aqua and terra sattillite images of the geospatial location defined by latlim and lonlim. Future users need to modify latlim and lonlim to match their region of interest. The function uses a user defined region of interest (ROI) within the latlim and lonlim to decide which satillite image is better between the aqua and terra images for a given date. The ROI is specified by the user on the first run through the code. ROI was useful for the Copper River because I cared most about if the region near the copper river mouth was clear of clouds, even though the satillite image was about 4 times as big as that ROI. The best satillite image (if there is one) is saved to a directory that is automatically created. Beware if you choose lots of dates to consider, this directory can get very large.


deep_parselLTERV2.m %%%%%% is a function that takes in an table (produced by LTERcruise_makerV2.m) and identifies data columns for a desired set of variables that I need by matching names. The columns are rearranged into a standard format. Date time is formalized to be a datenum. Practical Salinity and Absolute Salinity are calculated. There is some small amount of QCing that removes casts that regerister no value anywhere in the upper 300 meters. Casts of this type are removed since they are not real casts. Freshwater height is also cacluated here. Note that at this time freshwater height is calculated with Practical Salinity.


mydates_satretrieval_v2.m %%%% %This script is where you can select what satillite dates you want to retreive from the aqua and terra MODIS satillites. IT then runs the function WMSsatretrieval.m which does the heavy lifting of deciding which if the satillite images are worth keeping given the cloud cover, as well as creating a directory and saving the image in that directory. Search for root in WMSsatretrieval to see the directory path.

