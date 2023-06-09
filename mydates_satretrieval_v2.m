clear all
%all of this is 250 m resolution.
%note there was an error in aqua from october 22 to october 31 of 2005.
%These dates were skipped. Terra does have some good data for some of those
%dates.
cd('C:\Users\cfosadmin\Documents\MATLAB\Research\code\API'); 
load('dates.mat'); %we load up our remaining dates here;
tseg=t; %we take the first segment being index 1:36 (days))
%finding the land mask first


landmaskdate=datetime(2018,06,1,'Format','yyyy-MM-dd'); %arbitrary choice of landmaskdate. This retreives one instance of landmask that we apply to the analysis of each image. Land doesn't change to much here and perspective is identical in all images.
strlandmaskdate=string(landmaskdate);% gets a string version of our date time. So we can use strcat later.
url = 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi?'; 
landmaskinfo = wmsinfo(url);
landmask = landmaskinfo.Layer.refine('osm');
landmask = landmask.refine('OSM_Land_Mask','SearchField','LayerName');
latlim = [59.5, 60.6];
lonlim = [-147,-144];
imagelength = 2048;
[LA,LR] = wmsread(landmask,'ImageHeight',imagelength,'ImageWidth',imagelength,...
    'Latlim',latlim,'Lonlim',lonlim,'ImageFormat','image/png');
land=find(LA~=0);
clear landmaskinfo;
clear LA;

for i=1:length(tseg)
    

   cd('C:\Users\cfosadmin\Documents\MATLAB\Research\code\API'); 
   WMSsatretrieval(land,tseg,i); 
    
end

t(1:36)=[]; %we delete the days we just calculated from our date time values;
cd('C:\Users\cfosadmin\Documents\MATLAB\Research\code\API'); 
save('dates.mat','t');%resave t, replacing the one we loaded at the beginning, in preperation
%for the next loading