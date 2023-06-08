%load wind.
%now we need the wind activity for those dates that expressed case 1, and
%we need the wind activity for those dates that expressed case 2. See if
%there is any notable difference. 
%concerning wind variabiles u and v, these are eastward and northward
%components of wind. IOW they are the horizontal speed of air moving towards the east and towards north respectively,
%at a height of ten metres above the surface of the Earth, in metres per second.

%verbatium for v from copernicus website:
% This parameter is the northward component of the 10m wind. 
% It is the horizontal speed of air moving towards the north, 
% at a height of ten metres above the surface of the Earth, in
%  metres per second. Care should be taken when comparing this 
% parameter with observations, because wind observations vary on 
% small space and time scales and are affected by the local terrain, 
% vegetation and buildings that are represented only on average in the
%  ECMWF Integrated Forecasting System (IFS). This parameter can be combined 
% with the U component of 10m wind to give the speed and direction of 
% the horizontal 10m wind.
clear all 
close all
load('sampletimes.mat');
load('C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\50Ktest\v4_less_strict128_EOFallresultsPP.mat');
sT=T;
ncsource1='C:\Users\cfosadmin\Documents\MATLAB\Research\data\WindCopperRiver\hourlyData2002-2022\adaptor.mars.internal-1686094391.3367488-3810-18-653f8a5d-8d78-49a3-a1b5-eb25911fd2fa.nc'; %uv 2000-2007
ncsource2='C:\Users\cfosadmin\Documents\MATLAB\Research\data\WindCopperRiver\hourlyData2002-2022\adaptor.mars.internal-1686099790.2569082-17318-11-064e734d-f946-4bb2-84c7-ffa634b441f2.nc'; %uv 2008-2013
ncsource3='C:\Users\cfosadmin\Documents\MATLAB\Research\data\WindCopperRiver\hourlyData2002-2022\adaptor.mars.internal-1686105502.5966501-751-8-086cb908-aa97-4523-a44b-28781413bf36.nc'; %uv 2014-2020
ncsource4='C:\Users\cfosadmin\Documents\MATLAB\Research\data\WindCopperRiver\hourlyData2002-2022\adaptor.mars.internal-1686109298.838593-14340-1-46022ebd-4ace-434f-afd3-0142e971dc47.nc'; %uv 2021-2022

%first need to load the wind data.
info1=ncinfo(ncsource1);
info2=ncinfo(ncsource2);
info3=ncinfo(ncsource3);
info4=ncinfo(ncsource4);


%next make a compass rose histogram for the wind, stacked with wind speeds
%as well. Do this for all. Then do this for the dates of case1 and case2.
%Also do this the specific 3 day period that you identified earlier. 

ncdisp(ncsource1)
ncdisp(ncsource2)
ncdisp(ncsource3)
ncdisp(ncsource4)


v10_1=ncread(ncsource1,'v10'); %this is all days from 2002 to 2007
u10_1=ncread(ncsource1,'u10');
ncT_1=ncread(ncsource1,'time');

v10_2=ncread(ncsource2,'v10'); %this all from 2008 to 2013
u10_2=ncread(ncsource2,'u10');
ncT_2=ncread(ncsource2,'time');

v10_3=ncread(ncsource3,'v10'); %this all from 2014 to 2020
u10_3=ncread(ncsource3,'u10');
ncT_3=ncread(ncsource3,'time');

v10_4=ncread(ncsource4,'v10'); %this all from 2021 to 2022
u10_4=ncread(ncsource4,'u10');
ncT_4=ncread(ncsource4,'time');



singlepoint_u10=squeeze(u10_3(2,4,:)); %so this gets the uv for all available data.Hourly. 
singlepoint_v10=squeeze(v10_3(2,4,:));
figure
 scatter(singlepoint_u10,singlepoint_v10)

convertedT1=datetime(1900,1,1)+hours(ncT_1);
convertedT2=datetime(1900,1,1)+hours(ncT_2);
convertedT3=datetime(1900,1,1)+hours(ncT_3);
convertedT4=datetime(1900,1,1)+hours(ncT_4);

pre_u10=cat(3,u10_1,u10_2,u10_3,u10_4);
pre_v10=cat(3,v10_1,v10_2,v10_3,v10_4); %combined on the 3rd dimension.
% u10=u10_1;
% v10=v10_1;

pre_wT=[convertedT1;convertedT2;convertedT3;convertedT4]; %here is our complete set of dates and hours.


%To get a 13UTC to 13UTC 24 we generated the a daily bin index, and then deleted
%the first 13 and last 9 indices in the time, and deleted the last 22
%indices in the daily bin index. This effectively makes the first 23
%indices correspond to the our first "day" the next 23 indicies correspond
%to the next "day" and so on. If you doubt me, just display the wT and Y
%and check it for yourself.

%As we are going from 13UTC to 13UTC the next day, two days are involved.
% We will qualify the daily averages by when the avearge ends as that will be what is being matched 
% with satillite data. So for
% example an average that bins from january 1 13 UTC to january 2 13 UTC
% will be called January 2.

%then we discretize.
[Y,E] = discretize(pre_wT,'day');
%clip the ends of data and time.
Y(end-22:end)=[];

E(1)=[];
E(end)=[];

pre_wT(1:13)=[];
pre_wT(end-9:end)=[];

pre_u10(:,:,1:13)=[];
pre_u10(:,:,end-9:end)=[];
pre_v10(:,:,1:13)=[];
pre_v10(:,:,end-9:end)=[];

[m,n,z]=size(pre_v10);
u10=nan(m,n,length(E));
v10=nan(m,n,length(E));
for ix=1:max(Y)
    inditoavgover=find(Y==ix);
u10(:,:,ix)=mean(pre_u10(:,:,inditoavgover),3);
v10(:,:,ix)=mean(pre_v10(:,:,inditoavgover),3);
% if ~isempty(inditoavgover) %this bit of code is handy if there are
% missing indcies, but in this case we know everything is accounted for.
% %Ymat(1:length(discharge(inditoavgover)),ix)=discharge(inditoavgover);
% end

end

wT=E.';
%now we have our complete set of binned wind velocityes for the 24 hours previous to the image.
% corresponding u10 and v10 daily averages.




%finding the wind data times that match our sample dates. Our original data set was averages.Thus these times all
%correspond to the daily averages.
[value,P]=intersect(wT,sT);

sampleu10=u10(:,:,P);
samplev10=v10(:,:,P);
%we now have the the u and v wind data for our imagery samples. 

%now we choose the middle point.<----6/8/2022 why are we doing this? Why
%are we not averaging all of the points together?
singlepoint_sampleu10=squeeze(sampleu10(2,3,:)); %this gets the uv for just our sample. Which in this case is the 500 image dataset, that is good. 
singlepoint_samplev10=squeeze(samplev10(2,3,:));

singlepoint_u10=squeeze(u10(2,3,:)); %so this gets the uv for all available data.Hourly. 
singlepoint_v10=squeeze(v10(2,3,:));
figure
 scatter(singlepoint_u10,singlepoint_v10)

for i=1:length(squeeze(u10(2,3,:)))
averaged_u10(i)=mean(u10(:,:,i),'all');
averaged_v10(i)=mean(v10(:,:,i),'all');
end

for i=1:length(squeeze(sampleu10(2,3,:)))
averaged_sampleu10(i)=mean(sampleu10(:,:,i),'all');
averaged_samplev10(i)=mean(samplev10(:,:,i),'all');
end


%run princ as on Winddata.
singlepoint_sampleuv10=[singlepoint_sampleu10,singlepoint_samplev10]; %this is a :x2 data set. uv pairs for single point

singlepoint_uv10=[singlepoint_u10,singlepoint_v10];

averaged_singlepoint_uv10=[averaged_u10.',averaged_v10.'];

averaged_singlepoint_sampleuv10=[averaged_sampleu10.',averaged_samplev10.'];

% Run data through pca
X = singlepoint_sampleuv10;
XX = singlepoint_uv10;
AX=averaged_singlepoint_sampleuv10;
AXX=averaged_singlepoint_uv10;
% De-mean (MATLAB will de-mean inside of PCA, but I want the de-meaned values later)
% Do I want wind anomolies or actual wind? We are trying to figure out if
% the direction of the plume anomolies can be connected to the wind
% direction. With the original wind data, we can say that the wind was
% blowing in this direction at this time. What we suspect is that the main
% reason behind the plume heading west along the coast is coriolis force.
% Then it seems like the most important consideration is, "is the wind
% blowing against this, or with this, and thus might it expalin the
% anomlies we see in the plume". Using demeaned wind magnitudes, we are no
% longer looking at the true wind direction. If we plotted wind anomlies
% verus time, we are looking at when the wind speed slowed down. 

%lets say we wave a rotating tank, and have some dye and we can see a
%coriolis force moving the dye to the right. And then lets say we have
%multiple fans and we can turn them on and off to simulate diferent wind
%directions. If i'm interesteing in plume variations on the surface.
%plotting those variations (anomlies) against the true wind direction is a
%way of finding out if the wind directions are causing those variations, or
%if something else is causeing them. If I have a large fan that is always
%pointing in one direction, lets say at slight angle to the downstream
%flow, if I remove that I'm doing myself a disservice. Wind that truely
%changes it's direction is going to have an impact on the flow of the
%surface dye and on variations. But if I have wind anomlies, well now I'm
%saying if the wind less than or greater than some mean state, it's going
%to have an impact on the surface water and that's just not neccessarily
%true in a first order sense. That's really investigating if wind shear has
%an impact on the variations, and maybe it does, or maybe it doesn't. But
%in the first order the true wind direction intuitively makes the most
%sense to consider

%then again, if there is a mean wind flow, the plume would have already
%adjusted to it right? No because the wind flipping back and forth more
%than the curren. The current is always going in one direction more or
%less and the surface goes along with it. The wind is...more variable so
%it's mean is not felt by the surface? Ugh. I'm still not sure. 

mean_uvall=mean(XX);
mean_uvsample=mean(X);
mean_avguvall=mean(AXX);
mean_avguvsample=mean(AX);
X = X - mean(X); 
XX= XX - mean(XX); %repeat for all data, not just sample data.
AX=AX-mean(AX); %above is single point. Now this is the averaged point.
AXX=AXX-mean(AXX);

singlepointalldirection_for_originaluv=atan2(XX(:,1),XX(:,2));
%save('C:\Users\cfosadmin\Documents\MATLAB\Research\data\WindCopperRiver\my_results\Winddirection_originaluv_demeaned.mat','singlepointalldirection_for_originaluv')
% Do the PCA
[coeff,score,latent,~,explained] = pca(X);
[coeff2,score,latent,~,explained] = pca(XX);%repeat for all data, not just sample data.
[coeff3,score,latent,~,explained] = pca(AX);
[coeff4,score,latent,~,explained] = pca(AXX);
% Calculate eigenvalues and eigenvectors of the covariance matrix
covarianceMatrix = cov(X);
covarianceMatrix2 = cov(XX);
[V,D] = eig(covarianceMatrix);
[V2,D2] = eig(covarianceMatrix2);
% "coeff" are the principal component vectors.
% These are the eigenvectors of the covariance matrix.
% Compare "coeff" and "V". Notice that they are the same,
% except for column ordering and an unimportant overall sign.


% Multiply the original data by the principal component vectors to get the
% projections of the original data on the principal component vector space.
% % This is also the output "score". Compare ...
sampleuv10_InPCSpace_singlepoint = X*coeff;
uv10_InPCSpace_singlepoint = XX*coeff2;

sampleuv10_InPCSpace = X*coeff3;
uv10_InPCSpace = XX*coeff4;


%rotation angles
alpha = atan2d(coeff(1,2), coeff(1,1));
alpha2 = atan2d(coeff2(1,2), coeff2(1,1));
alpha3 = atan2d(coeff3(1,2), coeff3(1,1));
alpha4 = atan2d(coeff4(1,2), coeff4(1,1));

figure
scatter(AXX(:,1),AXX(:,2));
hold on
scatter(XX(:,1),XX(:,2));
scatter(AX(:,1),AX(:,2));
scatter(X(:,1),X(:,2));

%now we are talking about alongshore and cross shore winds. Great!

%find out what wind data matches.
%%
%divdie wind and images into Oct-April call that winter. call May-September
%Summer
%then divide winter and summer each into strong east, strong west, weak east,
%weak west

%then make composite images of the all the associated images.


%call weakwind 0-5 strong wind <5 m/s

%splitting up summer and winter samples.
indisummertime=sort(find(month(sT)<=9 & month(sT)>=5));
summertime=sT(indisummertime);%these are all of our summer month times for the sample. 


indiwintertime=sort(find(month(sT)<=4 | month(sT)>=10));
wintertime=sT(indiwintertime);

%P is a list of indices for the sample day within the span of days that go
%from 2002 01 01 to 2021 05 22 (I think it's 05 22, donest really matter
%since our last sample date is in feburary).
winterwinds=sampleuv10_InPCSpace(indiwintertime,1); %by choosing column one we are looking at along shore winds.
summerwinds=sampleuv10_InPCSpace(indisummertime,1);

%find indi of all our winds.
%idx = A<-15 | A>15;

indisummereaststrong=find(summerwinds>5);
indisummereastweak=find(summerwinds<=5 & summerwinds>0);
indisummerweststrong=find(summerwinds<-5);
indisummerwestweak=find(summerwinds>=-5 & summerwinds<0);

indiwintereaststrong=find(winterwinds>5);
indiwintereastweak=find(winterwinds<=5 & winterwinds>0);
indiwinterweststrong=find(winterwinds<-5);
indiwinterwestweak=find(winterwinds>=-5 & winterwinds<0);


%find actual wind

summereaststrong=summerwinds(indisummereaststrong);
summereastweak=summerwinds(indisummereastweak);
summerweststrong=summerwinds(indisummerweststrong);
summerwestweak=summerwinds(indisummerwestweak);

wintereaststrong=winterwinds(indiwintereaststrong);
wintereastweak=winterwinds(indiwintereastweak);
winterweststrong=winterwinds(indiwinterweststrong);
winterwestweak=winterwinds(indiwinterwestweak);

%
titles=["Winter West Weak", "Winter West Strong", "Winter East Weak", "Winter East Strong", "Summer West Weak",...
    "Summer West Strong", "Summer East Weak", "Summer East Stong"];
alongshorewind_varnames=["winterwestweak", "winterweststrong", "wintereastweak", "wintereaststrong", "summerwestweak",...
"summerweststrong", "summereastweak",  "summereaststrong"];

figure(1);
for i=1:8
subplot(4,2,i)
charwind=char(alongshorewind_varnames(i));
histbins=[];
eastwestcheck=strfind(alongshorewind_varnames(i),"west");
if ~isempty(eastwestcheck) %figures out if we have an east or a west wind(air moving towards the west!)
    direction_sig=-1;
else
    direction_sig=1;
end

strongweakcheck=strfind(alongshorewind_varnames(i),"strong");
if ~isempty(strongweakcheck) %figure out if we have a strong or weak wind category.
    %weak is greater than zero but less than or equal to 5. Strong is
    %greater than 5 and we assume there are none greater than 10 for this.
    edges=[5 6 7 8 9 10];
else
    edges=[0 1 2 3 4 5];
end


evalstatement=['histogram(direction_sig.*' charwind ',edges);'];
eval(evalstatement);
evalstatement2=['N=length(' charwind ');'];
eval(evalstatement2);
strN=string(N);
chrN=char(strN);
chrtitle=char(titles(i));
currenttitle= [ chrtitle ':N=' chrN];
title(currenttitle);
ylabel('Samples')
xlabel('Wind Speed (m/s)')
end
sgtitle('Histograms of East/West Strong/Weak Winds')
%%
%find images that correspond
set(gcf, 'Position', get(0, 'Screensize'));
saveas(figure(1),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\Differencing\HistEastWestStrongWeakWinds.png');
figure(2)

long=[-147:0.0236:-144];
lat=[59.5:0.0086:60.6];
%images corresponding to winter winds.

avgedimages_for_winterWeWWeS_EaWEaS=nan(128, 128, 4);
for i=1:4
    subplot(2,2,i)
    charwind=char(alongshorewind_varnames(i));
    evalstatement=['avgcurrentwind=H(:,:,indi' charwind ');'];
    eval(evalstatement);
    compositeavg=mean(avgcurrentwind,3);
    avgedimages_for_winterWeWWeS_EaWEaS(:,:,i)=compositeavg;
    m_proj('Lambert','longitudes',[-147.4 -143.8],'latitudes',[59.3 60.8], 'rectbox', 'off'); 
    m_pcolor(long,lat,flipud(compositeavg));
    shading flat;
    hold on
    cb=colorbar;
    set(get(cb,'label'),'string','Green Channel Shade','fontsize',9);
    hold on
    m_gshhs_h('patch',[.5 .5 .5]); %adds coastline
    m_grid('tickdir','out','linewi',2,'fontsize',9);
    hold on
    clevels=[0:15:255];

    m_contour(long,lat,flipud(compositeavg),clevels,'k');
    hold on
    
    evalstatement2=['N=length(' charwind ');'];
    eval(evalstatement2);
    strN=string(N);
    chrN=char(strN);
    chrtitle=char(titles(i));
    currenttitle=[ chrtitle ':N=' chrN];
    title(currenttitle);
   
    hold off
end
sgtitle('Mean pixel shade of Summer Winds (May-September)')

figure(3)
avgedimages_for_summerWeWWeS_EaWEaS=nan(128, 128, 4);
% images corresponding to summer winds
for i=1:4
    subplot(2,2,i)
    charwind=char(alongshorewind_varnames(i+4));
    evalstatement=['avgcurrentwind=H(:,:,indi' charwind ');'];
    eval(evalstatement);
    compositeavg=mean(avgcurrentwind,3);
    avgedimages_for_summerWeWWeS_EaWEaS(:,:,i)=compositeavg;
    m_proj('Lambert','longitudes',[-147.4 -143.8],'latitudes',[59.3 60.8], 'rectbox', 'off'); 
    m_pcolor(long,lat,flipud(compositeavg));
    shading flat;
    hold on
    cb=colorbar;
    set(get(cb,'label'),'string','Green Channel Shade','fontsize',9);
    hold on
    m_gshhs_h('patch',[.5 .5 .5]); %adds coastline
    m_grid('tickdir','out','linewi',2,'fontsize',9);
    hold on
    clevels=[0:15:255];

    m_contour(long, lat, flipud(compositeavg),clevels,'k');
    hold on
    
    evalstatement2=['N=length(' charwind ');'];
    eval(evalstatement2);
    strN=string(N);
    chrN=char(strN);
    chrtitle=char(titles(i+4));
    currenttitle=[chrtitle ':N=' chrN];
    title(currenttitle);
    hold off
    hold off
end
sgtitle('Mean pixel shade of Winter Winds (October-April)')

figure(4)
m_pcolor(long,lat,S2) 

shading flat
hold on

cb=colorbar('Location', 'southoutside');
    set(get(cb,'label'),'string','{\bf {\sigma}} of pixel intensity','fontsize',14);
   

    hold on
    m_gshhs_h('patch',[.5 .5 .5]); %adds coastline
    m_grid('tickdir','out','linewi',2,'fontsize',14);
    hold on
    clevels=[25:2.5:55];

    m_contour(long, lat, S2,clevels,'k');
    hold on
    
    title('Standard Deviation ({\sigma})');
    
hold off

figure(5)
[CS,CH]=m_etopo2('contour',[-250 -100],'edgecolor', 'k'); %adds bathymetry.

figure(4)
hold on

%%
%Contour line simplification
%     The contour matrix C is a two row matrix of contour lines. Each
%     contiguous drawing segment contains the value of the contour,
%     the number of (x, y) drawing pairs, and the pairs themselves.
%     The segments are appended end-to-end as
%  
%         C = [level1 x1 x2 x3 ... level2 x2 x2 x3 ...;
%              pairs1 y1 y2 y3 ... pairs2 y2 y2 y3 ...]

cl=[-250 -100];
for inn=1:length(cl)
currentlevel=cl(inn);
indispecificlevel=find(CS(1,:)==currentlevel);
for in=1:(length(indispecificlevel)-1)
    currentindi=indispecificlevel(in);
if CS(2,currentindi)<50 %removes contours less than 50 positions long
   CS(1:2,indispecificlevel(in)+1:indispecificlevel(in+1)-1)=NaN; 
end
end
end

%takes in coutour matrix and gives back a structure that i can work with
s=contourdata(CS);
% S(k).level contains the contour level height of the k-th line.
% S(k).numel contains the number of points describing the k-th line.
% S(k).isopen is True if the k-th contour is open and False if it is closed.
% S(k).xdata contains the x-axis data for the k-th line as a column vector.
% S(k).ydata contains the y-axis data for the k-th line as a column vector.
%
% For example: PLOT(S(k).xdata,S(k).ydata)) plots just the k-th contour.

%plots contours with jet colors

figure(4);
for k=1:length(s)
if s(k).level==-100
  color=[0.9412    0.6000    0.8275]; %light pink
end

if s(k).level==-250
  color=[1.0000    0.7373    0.3412]; %peach
end


plot(s(k).xdata,s(k).ydata,'color',color,'LineWidth',3);

hold on
end



m_text(-144.4,59.720,'100m','fontsize',14,'fontweight','bold', 'fontangle', 'italic','color', [0.9412    0.6000    0.8275]);
m_text(-144.4,59.57,'250m','fontsize',14,'fontweight','bold', 'fontangle', 'italic','color', [1.0000    0.7373    0.3412]);
set(gca,'fontsize',16);

samplews10=((sampleuv10_InPCSpace(:,1).^2)+(sampleuv10_InPCSpace(:,2).^2)).^(1/2);
direction=atan2(sampleuv10_InPCSpace(:,1),sampleuv10_InPCSpace(:,2));

% obtain wind speed and direction for each location for all wind data
% including non sample days. 
saveas(figure(4),'C:\Users\cfosadmin\Documents\SiteReveiw\satilliteUTC13_special_v3_std_dev_of_image.png')
avgws10=((uv10_InPCSpace(:,1).^2)+(uv10_InPCSpace(:,2).^2)).^(1/2);
alldirection=atan2(uv10_InPCSpace(:,1),uv10_InPCSpace(:,2));
%%

%THIS SECTION IS ONLY REQUIRED IF WE ARE NOT DOING SINGLE POINT.
%orentating our matrix to be like our images.
%EDIT NO roation needed since we are now just doing a single point.
% samplews10=permute(samplews10,[2,1,3]);
% direction=permute(direction,[2,1,3]);
% sampleu10=permute(sampleu10,[2,1,3]);
% samplev10=permute(samplev10,[2,1,3]);
% 
% %all wind including non sample days.
% avgws10=permute(avgws10,[2,1,3]);
% alldirection=permute(alldirection,[2,1,3]);


%EDIT only needed if we want to average the area.
%averages for sample wind.
% for i=1:length(samplews10(1,1,:))
% totalarea_avg_sample_ws10(i)=mean(samplews10(:,:,i),'all');
% totalarea_dir(i)=mean(direction(:,:,i),'all');
% end

%gives us the averages for all the wind including non-sample days.
% for i=1:length(avgws10(1,1,:))
% totalarea_avg_ws10(i)=mean(avgws10(:,:,i),'all');
% totalarea_alldir(i)=mean(alldirection(:,:,i),'all');
% end

%%
%this is all wind anomlies in the time period:

char_sample=char(string(length(averaged_samplev10)));
char_all=char(string(length(averaged_v10)));
wind_rose(rad2deg(alldirection),avgws10);
char_yearstart=char(string(max(year(wT))));
char_yearend=char(string(min(year(wT))));
title(['Demeaned  All Wind anomolies for ',char_yearstart,'-',char_yearend,' (N=',char_all,')'])
saveas(figure(6),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\annualcycle\satilliteUTC13_special_allwind.png')

%this is all sample wind


wind_rose(rad2deg(direction),samplews10);
title(['Demeaned All Wind Speed and Direction anomlies for all samples(N=',char_sample,')'])
saveas(figure(7),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\annualcycle\satilliteUTC13_special_allsamplewindP.png')

%this is z>-5000 

char_samplesummer=char(string(length(samplews10(indisummertime))));
wind_rose(rad2deg(direction(indisummertime)),samplews10(indisummertime));
title(['Demeaned Wind anmolies for Summer mode of E1 (N=',char_samplesummer,')'])
saveas(figure(9),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\annualcycle\satilliteUTC13_special_summerwindP.png')

%this is z<-10000 
char_samplewinter=char(string(length(samplews10(indiwintertime))));
wind_rose(rad2deg(direction(indiwintertime)),samplews10(indiwintertime));
title(['Demeaned Wind Anomlies for Winter mode of E1 (N=',char_samplewinter,')'])

saveas(figure(10),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\annualcycle\satilliteUTC13_special_winterwindP.png')


%adding the mean back in:

avgws10_true=(((uv10_InPCSpace(:,1)+mean_uvall(:,1)).^2)+(uv10_InPCSpace(:,2)+mean_uvall(:,2)).^2).^(1/2);
alldirection_true=atan2(uv10_InPCSpace(:,1)+mean_uvall(:,1),uv10_InPCSpace(:,2)+mean_uvall(:,2));

samplews10_true=(((sampleuv10_InPCSpace(:,1)+mean_uvsample(:,1)).^2)+(sampleuv10_InPCSpace(:,2)+mean_uvsample(:,2)).^2).^(1/2);
direction_true=atan2(sampleuv10_InPCSpace(:,1),sampleuv10_InPCSpace(:,2));



%this is all wind in the time period:
wind_rose(rad2deg(alldirection_true),avgws10_true);
title(['True Wind for ',char_yearstart,'-',char_yearend,' (N=',char_all,')'])
saveas(figure(11),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\annualcycle\satilliteUTC13_special_allwind_true.png')

%this is all sample wind


wind_rose(rad2deg(direction_true),samplews10_true);
title(['True Wind Speed and Direction anomlies for all samples(N=',char_sample,')'])
saveas(figure(12),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\annualcycle\satilliteUTC13_special_allsamplewindP_true.png')

%this is z>-5000 


wind_rose(rad2deg(direction_true(indisummertime)),samplews10_true(indisummertime));
title(['True Wind anmolies for Summer mode of E1 (N=',char_samplesummer,')'])
saveas(figure(13),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\annualcycle\satilliteUTC13_special_summerwindP_true.png')

%this is z<-10000 


wind_rose(rad2deg(direction_true(indiwintertime)),samplews10_true(indiwintertime));
title(['True Wind Anomlies for Winter mode of E1 (N=',char_samplewinter,')'])

saveas(figure(14),'C:\Users\cfosadmin\Documents\MATLAB\Research\Figures\API\sample\results\annualcycle\satilliteUTC13_special_winterwindP_true.png')

samplecoeff_singlepoint=coeff;
allwindcoeff_singlepoint=coeff2;

samplecoeff=coeff3;
allwindcoeff=coeff4;

%single point
rotangle_sample_singlepoint=alpha;
rotangle_all_singlepoint=alpha2;

%averaged 6 grid point
rotangle_sample_6pointavg=alpha3;
rotangle_all_6pointavg=alpha4;


save('C:\Users\cfosadmin\Documents\MATLAB\Research\Data\WindCopperRiver\my_results\satilliteUTC13_special_windanalysis_results.mat','direction','P','wT','singlepoint_sampleuv10','singlepoint_uv10'...
,'mean_uvsample_singlepoint','mean_uvall_singlepoint','samplecoeff_singlepoint','allwindcoeff_singlepoint','rotangle_all_6pointavg','rotangle_sample_6pointavg','rotangle_sample_singlepoint','rotangle_all_singlepoint','sampleuv10_InPCSpace','uv10_InPCSpace','mean_uvall','mean_uvsample','direction_true','samplews10_true','alldirection_true','avgws10_true','allwindcoeff','samplecoeff','sampleuv10_InPCSpace_singlepoint','uv10_InPCSpace_singlepoint');
%load('C:\Users\cfosadmin\Documents\MATLAB\Research\Data\WindCopperRiver\my_results\windanalysis_results.mat')


