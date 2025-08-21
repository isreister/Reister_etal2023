%This function is saved to isreister/Reister_etal2023 github account as
%WMSsatretrieval.m

%this is a function that is used within the code script 'mydates_satretrieval_v2.m'. 
% Inputs for the function are a land mask (land), a list of dates (t), and an index 
% value of 1,2,3,4,...length(t). The function uses an API to recover the cloud_fraction for 
% aqua and terra sattillite images of the geospatial location defined by latlim and lonlim. 
% Future users need to modify latlim and lonlim to match their region of interest. The function 
% uses a user defined region of interest (ROI) within the latlim and lonlim to decide which satillite 
% image is better between the aqua and terra images for a given date. The ROI is specified by the user
% on the first run through the code. ROI was useful for the Copper River because I cared most about 
% if the region near the copper river mouth was clear of clouds, even though the satillite image was 
% about 4 times as big as that ROI. The best satillite image (if there is one) is saved to a directory
% that is automatically created. Beware if you choose lots of dates to consider, this directory can get very large.



function WMSsatretrievalv2(land,t,i)

clearvars -global -except land latlim lonlim imagelength url1 url2 t i
currentdate=t(i);
writtenmonth=datetime(t(i),'Format','MMMM d, yyyy');
strcurrentdate=string(currentdate);
disp(strcurrentdate)
charcurrentdate=char(strcurrentdate);
charmonth=char(writtenmonth);
year=charcurrentdate(1:4);
month=charmonth(1:4);
day=charcurrentdate(9:10);




url1 = 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi?SERVICE=WMS&REQUEST=GetMap&VERSION=1.3.0&LAYERS=MODIS_Terra_MODIS_Terra_Cloud_Fraction_Day&STYLES=&FORMAT=image%2Fpng&TRANSPARENT=true&HEIGHT=256&WIDTH=256&TIME=';
url2= strcat(strcurrentdate, '&CRS=EPSG:4326&BBOX=-22.5,0,0,22.5'); 
url= char(strcat(url1,url2));
%url = 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi?'; %in theory this should work..it's the same as the url used to the land mask and that works fine.
      

info = wmsinfo(url);
modis = info.Layer.refine('modis'); %this is all the MODIS satilite layers.
Cloud_Fraction = modis.refine('Cloud_Fraction','SearchField','LayerTitle');%all of the cloud_fraction layers
Cloud_Fraction_aqua = Cloud_Fraction.refine('Aqua','SearchField','LayerTitle'); %all of the aqua cloud fraction layers
Cloud_Fraction_terra = Cloud_Fraction.refine('Terra','SearchField','LayerTitle');
Cloud_Fraction_aqua_day = Cloud_Fraction_aqua.refine('Day','SearchField','LayerTitle');% the day portion of the cloud fraction layers?
Cloud_Fraction_terra_day = Cloud_Fraction_terra.refine('Day','SearchField','LayerTitle'); %
latlim = [59.5, 60.6];
lonlim = [-147,-144];
imagelength = 2048;
[A] = wmsread(Cloud_Fraction_aqua_day,'ImageHeight',imagelength,'ImageWidth',imagelength,... %now read the part that actually matters for us.
    'Latlim',latlim,'Lonlim',lonlim);
[B]= wmsread(Cloud_Fraction_terra_day,'ImageHeight',imagelength,'ImageWidth',imagelength,...
    'Latlim',latlim,'Lonlim',lonlim);



truecolortesting = modis.refine('MODIS_Aqua_CorrectedReflectance_TrueColor','SearchField','LayerTitle');
[T]= wmsread(truecolortesting,'ImageHeight',imagelength,'ImageWidth',imagelength,...
    'Latlim',latlim,'Lonlim',lonlim);



    figure(1) %uncomment to see what you are looking at.
    clf
    imshow(T)%original cloud fraction image

    figure(2) %uncomment to see what the heck your looking at.
    clf
    imshow(A)%original cloud fraction image
%%
clear info;
%% Decide betweeen Terra and Aqua.


tempCell=cell(1,2);

tempCell{1,1}=A;
tempCell{1,2}=B;


%---------

tempCell{1,1} = A;
tempCell{1,2} = B;

CustomCloudPercent = nan(1,2);

if i==1
        C = tempCell{1};

    % Step 1: apply land mask first
    C(repmat(land,1,1,3)) = 0;

    % Step 2: show image and let user draw polygon
    % User should draw the polygon around the region of interest within the
    % satillite image. Don't worry about cloud cover or land. Land will be
    % masked even if it is part of the region of interest. cloud cover is
    % delt examined within  each image, and a cloud fraction for the region
    % of interest is calculated to determin if this satillite image is
    % worth anything to you.
    figure(1); clf;
    imshow(C);
    title('Draw polygon around region of interest, finish by clicking on your starting point again');
    h = drawpolygon; % user draws ROI
    BW = createMask(h); % logical mask of polygon
    save('polygon_of_interest.mat',"BW")
else
    load('polygon_of_interest.mat','"BW"')
end

for ix = 1:length(tempCell)
    C = tempCell{ix};

    % Step 1: apply land mask first
    C(repmat(land,1,1,3)) = 0;

    
    % Step 2: apply polygon mask
    C(~repmat(BW,1,1,3)) = 0;

    % Step 3: count red (cloud) pixels inside polygon
    RED = uint8([230,50,50]);
    inRED = (C(:,:,1) > RED(1)) & (C(:,:,2) < RED(2)) & (C(:,:,3) < RED(3));

    % Step 4: count valid (non-black) pixels inside polygon
    BLACK = uint8([0,0,0]);
    inBLACK = (C(:,:,1)==BLACK(1)) & (C(:,:,2)==BLACK(2)) & (C(:,:,3)==BLACK(3));

    pixelsthatareNOTblack = numel(inBLACK) - nnz(inBLACK);

    % Step 5: calculate cloud fraction
    CustomCloudPercent(ix) = (nnz(inRED) / pixelsthatareNOTblack) * 100;
end

%disp(CustomCloudPercent);

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%


   
aquaorterra=["Aqua" "Terra"]; 
[mymin,myindi]=min(CustomCloudPercent); %Here the we just say "ok given aqua and terra images, which has less clouds".
if mymin<50 %Here we say Ok given that we are looking at the best image between aqua and terra, is it actually good enough (aka, less than 50% clouds).
    
    bestofaquaandterra=char(aquaorterra(myindi));
    
    true_color = modis.refine(bestofaquaandterra,'SearchField','LayerTitle');
    true_color = true_color.refine('TrueColor','SearchField','LayerTitle');
    [A,R] = wmsread(true_color,'ImageHeight',imagelength,'ImageWidth',imagelength,...
        'Latlim',latlim,'Lonlim',lonlim);
    clear modis;
    f = figure('visible','off');
    geoshow(A,R);
    
    
    root='C:\Users\funkb\Documents\MATLAB\Research\Figures\Testing'; %@ this is where the images will be saved. These images will be 50 percent clouds or less for the section of the image you specify while cropping the image. If you just want to get 50% or less clouds for a square image, just remove the cropping.
    
    path1=[root,'\',year]; 
    if ~exist(path1, 'dir')
       mkdir(path1) % beware! this will make folders in your directory specified by root. Which is nice for organization but just be aware of that.
    end
    month=strrep(month,' ',[]); %put this in to take care of the month of may
    %which has a space at it's end.
    path2=[root,'\',year,'\',month];
    if ~exist(path2, 'dir')
       mkdir(path2)
    end
    cd(path2);
    filename=[path2,'\',day,'.png'];
    filename2=[path2,'\',day,'AR.mat'];
    
    
    saveas(f,filename) %this is a .png for the portion of globe you are interested saved according to the projection used here which is  EPSG:4326. you can always change the projection later, but that is of course another can of beans.
    save(filename2,'A','R'); %the data that is save here is an matrix of colors (A) and a raster file (R). With those the world is pretty much your oyster.  
end

end

