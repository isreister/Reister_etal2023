
function WMSsatretrieval(land,t,i)

clearvars -global -except land latlim lonlim imagelength url1 url2 t i
currentdate=t(i);
writtenmonth=datetime(t(i),'Format','MMMM d, yyyy');
strcurrentdate=string(currentdate)
charcurrentdate=char(strcurrentdate);
charmonth=char(writtenmonth);
year=charcurrentdate(1:4);
month=charmonth(1:4);
day=charcurrentdate(9:10);





url = 'https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi?';
info = wmsinfo(url);
modis = info.Layer.refine('modis');
Cloud_Fraction = modis.refine('true color','SearchField','LayerTitle');
Cloud_Fraction_aqua = Cloud_Fraction.refine('Aqua','SearchField','LayerTitle');
Cloud_Fraction_terra = Cloud_Fraction.refine('Terra','SearchField','LayerTitle');
Cloud_Fraction_aqua_day = Cloud_Fraction_aqua.refine('Day','SearchField','LayerTitle');
Cloud_Fraction_terra_day = Cloud_Fraction_terra.refine('Day','SearchField','LayerTitle');
latlim = [59.5, 60.6];
lonlim = [-147,-144];
imagelength = 2048;
[A] = wmsread(Cloud_Fraction_aqua_day,'ImageHeight',imagelength,'ImageWidth',imagelength,...
    'Latlim',latlim,'Lonlim',lonlim);
[B]= wmsread(Cloud_Fraction_terra_day,'ImageHeight',imagelength,'ImageWidth',imagelength,...
    'Latlim',latlim,'Lonlim',lonlim);

clear info;
%% Decide betweeen Terra and Aqua.


tempCell=cell(1,2);

tempCell{1,1}=A;
tempCell{1,2}=B;


%---------
%In this section we rotate the image a few times and make black sections of the cloud faction calculation that
%we do not need (this applies the land mask as well the east of kayak
%island and west of hinchinbrook island. The seaward edge is approximately
%the blanked out at the shelf edge.
    % figure(1)
    % clf
    % imshow(A)%original cloud fraction image
    CustomCloudPercent=nan(1,2);
for ix=1:length(tempCell)

    C=tempCell{1,ix};

    C(land)=0;%land mask applied
%     
%         figure(2)
%         clf
%         imshow(C)
    
    C = imrotate(C,27,'crop');
%         figure(4)
%         imshow(C)
    bottomquarter=imagelength-(imagelength/4);
    C(bottomquarter:imagelength,:,:)=0; %mask of "beyond the shelf edge" applied.
%         figure(5)
%         imshow(C)
    
    C = imrotate(C,7,'crop');
%         figure(6)
%         imshow(C)
    rightninth=imagelength-(imagelength/9);
    C(:,rightninth:imagelength,:)=0; %mask of east of kayak island applied
%         figure(7)
%         imshow(C)
    
    C = imrotate(C,25,'crop');%the crop here is enough to remove "west of hinchinbrook".
%         figure(8)
%         imshow(K)
    %this completes the specification of our bounding box.
    
    %----------------
    Re=C;
    RED=[uint8(230),uint8(50),uint8(50)].';
    inRED=find(Re(:,:,1)>RED(1) & Re(:,:,2)<RED(2) & Re(:,:,3)<RED(3));
    clear Re;
    % R1=R(:,:,1);%designates the different layers.
    % R2=R(:,:,2);
    % R3=R(:,:,3);
    % R1(inRED)=uint8(255);
    % R2(inRED)=uint8(255);
    % R2(inRED)=uint8(255);
    % R(:,:,1)=R1;
    % R(:,:,2)=R2;
    % R(:,:,3)=R3;%makes all the 100% clouds within our "bounding box that is actually created by just removing
    % %everything we don't want" white
    % clear R1;
    % clear R2;
    % clear R3;
    %     figure(9)
    %     clf
    %     imshow(R)
    
    BLACK=[uint8(0),uint8(0),uint8(0)].';
    inBLACK=find(C(:,:,1)==BLACK(1) & C(:,:,2)==BLACK(2) & C(:,:,3)==BLACK(3));
    clear C;
    
    pixelsthatareblack=length(inBLACK); %this is our total number of pixels that we DO NOT want to consider. Note that
    %the last value (and max value) here designates the bottom right corner of our image
    %and is 4194304.Its imporant to know that find is giveing back a vector of
    %indice values that correspond to the row column space ONLY. it does not
    %correspond to the row column level space (3D space of the matrix).
    %Therefor if we want the pixels that are not black, we take the total
    %number of pixels in row column space and subtract black colored pixels,
    %leaving behing other colored pixesl, which is, the "bounding box" we want
    %to consider as 100% for our cloud faction analysis. 
    clear inBLACK;
    pixelsthatareNOTblack=(2048.*2048)-pixelsthatareblack;
    disp(pixelsthatareNOTblack); %this value should not change.
    disp(length(inRED)); %this value should change
    CustomCloudPercent(ix)=(length(inRED)./pixelsthatareNOTblack).*100;
    %This commented section will show that the 100% cloud faction pixels are 
    %being dectected correctly. They are correct if the red pixels in figure 1
    %become black pixels in figure 2.The length of InRED gives us how many 100
    %pixels there there are. Divided by the total number of pixels, gives us
    %the percentage of the image that registers as 100% cloud faction. This
    %percentage is what we use to filter out cloudy days.
%------------
end
aquaorterra=["Aqua" "Terra"]; 
[mymin,myindi]=min(CustomCloudPercent);
if mymin<50

bestofaquaandterra=char(aquaorterra(myindi));

true_color = modis.refine(bestofaquaandterra,'true color','SearchField','LayerTitle');
true_color = true_color.refine(bestofaquaandterra,'SearchField','LayerTitle');
[A,R] = wmsread(true_color,'ImageHeight',imagelength,'ImageWidth',imagelength,...
    'Latlim',latlim,'Lonlim',lonlim);
clear modis;
f = figure('visible','off');
geoshow(A,R);


root='C:\Users\Isaac\Documents\MATLAB\Research\Figures\API';

path1=[root,'\',year];    
if ~exist(path1, 'dir')
       mkdir(path1)
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


saveas(f,filename)
save(filename2,'A','R');
end

end

% some text
% teoms dsf