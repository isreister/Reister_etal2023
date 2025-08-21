%This script is saved to isreister/Reister_etal2023 github account as
%mydates_satretrieval_v2.m

%This script is where you can select what satillite dates you want to
%retreive from the aqua and terra MODIS satillites. IT then runs the
%function WMSsatretrieval.m which does the heavy lifting of deciding which
%if the satillite images are worth keeping given the cloud cover, as well
%as creating a directory and saving the image in that directory. Search for
%root in WMSsatretrieval to see the directory path.

clear all
%all of this is 250 m resolution.
%note there was an error in aqua from october 22 to october 31 of 2005.
%These dates were skipped. Terra does have some good data for some of those
%dates.
cd('C:\Users\funkb\Documents\MATLAB\Research\code\API'); % Note that this location is not where the satillite images will be saved, you'll see that in the WMSsatretrieval function.

t1 = datetime(2004,05,21,'Format','yyyy-MM-dd');% Make datetime vector for date span you are intersted in.
t2 = datetime(2010,12,31,'Format','yyyy-MM-dd');
tseg = t1:t2;

%finding the land mask first


landmaskdate=datetime(2018,06,1,'Format','yyyy-MM-dd'); %arbitrary choice of landmaskdate. This retreives one instance of landmask that we apply to the analysis of each image. 
strlandmaskdate=string(landmaskdate);% gets a string version of our date time. So we can use strcat later.

url1 = 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi?SERVICE=WMS&REQUEST=GetMap&VERSION=1.3.0&LAYERS=MODIS_Terra_MODIS_Terra_Cloud_Fraction_Day&STYLES=&FORMAT=image%2Fpng&TRANSPARENT=true&HEIGHT=256&WIDTH=256&TIME=';
url2= strcat(strlandmaskdate, '&CRS=EPSG:4326&BBOX=-22.5,0,0,22.5'); 
url= char(strcat(url1,url2));
%url = 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi?'; 
landmaskinfo = wmsinfo(url);
landmask = landmaskinfo.Layer.refine('osm');
landmask = landmask.refine('OSM_Land_Mask','SearchField','LayerName');
latlim = [59.5, 60.6];
lonlim = [-147,-144];
imagelength = 2048;
[LA,LR] = wmsread(landmask,'ImageHeight',imagelength,'ImageWidth',imagelength,...
    'Latlim',latlim,'Lonlim',lonlim,'ImageFormat','image/png');
imshow(LA)
land=find(double(LA(:,:,1))==75);
% testing=zeros(2048,2048);
% testing(land)=1;
% imshow(testing);
clear landmaskinfo;
clear LA;

for i=1:length(tseg)
    

   cd('C:\Users\funkb\Documents\MATLAB\Research\code\API'); 
   WMSsatretrievalv2(land,tseg,i); 
    
end


cd('C:\Users\cfosadmin\Documents\MATLAB\Research\code\API'); 
save('dates.mat','t');%resave t