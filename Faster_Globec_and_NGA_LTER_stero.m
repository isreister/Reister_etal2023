%parse seward matlab.
clear all
%Don't forget to change your path name here!!!!!!
load('C:\Users\cfosadmin\Documents\MATLAB\Research\Data\NGALTER\SewardLineStations');

%think of S and T as 3 dimensional matrix of S and T.
%there are 13 stations.
%There are 80 cruises.
%There are 301 depths.
%So a salinity measuremnt is given for a particular time cruise, at a
%particular depth, at a particular station. We only consider the top 300 as
%this is the GAK line.

%We want to put this into our other structure, which at it's surface level
%is structured first by season, or cruise. parsed(1) is spring 2018.
%we can load our other structure and start there.
strddd=["25","50","75","100","200","270","300"];
ddd=[25,50,75,100,200,270,300];
for currentdepth=1:length(ddd)
    load('C:\Users\cfosadmin\Documents\MATLAB\Research\Data\NGALTER\SewardLineStations');
 strddd=["25","50","75","100","200","270","300"];
ddd=[25,50,75,100,200,270,300];
figure(1)
clf
figure(2)
clf
depth=strddd(currentdepth);
depthnum=ddd(currentdepth);
load('C:\Users\cfosadmin\Documents\MATLAB\Research\Data\NGALTER\ALL_NGALTERdata_dec2020');
seasons=80;
higheststation=13;
stations=1:13;
Jan=find(TGD(:,2)==1);% compiling month indices.
Feb=find(TGD(:,2)==2);
Mar=find(TGD(:,2)==3);
Apr=find(TGD(:,2)==4);
May=find(TGD(:,2)==5);
Jun=find(TGD(:,2)==6);
Jul=find(TGD(:,2)==7);
Aug=find(TGD(:,2)==8);
Sep=find(TGD(:,2)==9);
Oct=find(TGD(:,2)==10);
Nov=find(TGD(:,2)==11);
Dec=find(TGD(:,2)==12);
YYYY1=1997;%year limits
YYYY2=2020;
years=(YYYY1:YYYY2);
years=string(years);%string of years

FWinarray=NaN(80,higheststation);% this creates a FW height for each station between 1 and 13, in time and space.
for i= 1:13
    for ii=1:seasons
    top50S=S(1:depthnum,i,ii);
    top50T=T(1:depthnum,i,ii);
    avg_Sal50=mean(top50S);
    avg_Tem50=mean(top50T);
    Ref_SP=33.8; % This is Reference salinty in practical salinity. Absolute salinty would be 34.0 and is calucalted 
    % using Gibbs seawater toolkit sw_SA_from_SP. The SP
    % used was the reference of 33.8 from Weingartner 2005 which in turn comes
    % from Dodimead. This is SP!!! OK to use the same ref salinity for all depth??
    Ref_p=depthnum/2; %Must change if depth changes. This is our half way point in our uppper 50 meters
    Ref_rho=gsw_rho(Ref_SP,avg_Tem50,Ref_p); %We don't have a reference CT to use so we are using the CT in the data.
    conc1=Ref_SP*Ref_rho;
    %depth=C(:,7); This line and next is method to get the measured volume
    %(-1)*depth(end); It doesn't handle the NaNs well/might not be needed... 
    V1= depthnum; %Just using 50 for volume as we are assuming a 1x1x(~50) column here
    Current_rho=gsw_rho(avg_Sal50,avg_Tem50,Ref_p); %Using the same pressure (25 db) for both ref and current.
    conc2=avg_Sal50*Current_rho;
    FWinput=((conc1*V1)/conc2)-V1;
    FWinarray(ii,i)=FWinput;
    end
    
end


dt_TGD=datetime(TGD); %makes a datetime array for use later.



%multiple GAK1s sometimes. In the first cruise it looked like we did one at
%the beginning of the cruise and then 2 at the end of the cruise. Choosing
%the 1st of the ones at the end is probably best? But for now just taking
%the average of the three.

%update, this needs to be looked at closer. We must use the combed over
%data set that Seth provided. That should do away with this requirement.

%This adds in GAK14 an GAK 15 to our 1-13 stations
test=NaN;
for i=1:6
    numbersonly=parsed(i).GAK.Numbersonly;
    highstation=max(parsed(i).GAK.Numbersonly);
    
    GAKFW2_high=NaN(highstation,1);
    %making space for lat longs of GAK1
    if i==1
        GAKLAT=NaN(higheststation,1);
        GAKLONG=NaN(higheststation,1);
    end
    integerstation_indi=find(floor(numbersonly(:))==numbersonly(:));
    integerstation=parsed(i).GAK.Station(integerstation_indi);
    integerFW=parsed(i).GAK.LTERFW50(integerstation_indi);
    
    %Grabbing the lat and long for the spring cruise in 2018 (which was all
    %the stations for GAK
    if i==1
    integerLat=parsed(i).GAK.LTERlat(integerstation_indi);
    integerLong=parsed(i).GAK.LTERlong(integerstation_indi);
    end
%     for ii=1:length(parsed(i).GAK.Station)
%     test=parsed(i).GAK.Station(ii);
%     
%     test2=char(test);
%     if test2(end)=='I'%skip intermediate
%         continue
%     end
%     end
for ii=1:highstation-1
    strstat=string(ii+1);
    fullstatstr=strcat("GAK",strstat);
    consolidated_station_indi=eval(strcat('find(integerstation=="',fullstatstr,'")'));
    FWGAK=integerFW(consolidated_station_indi);
    GAKFW2_high(ii+1)=mean(FWGAK);
    %averaging lat and long for the 1st cruise only
    if i==1
    LATGAK=integerLat(consolidated_station_indi);
    LONGGAK=integerLong(consolidated_station_indi);
    GAKLAT(ii+1)=mean(LATGAK);
    GAKLONG(ii+1)=mean(LONGGAK);
    end
end    

FWGAK1=mean(parsed(i).GAK1.LTERFW50);%MAYBE the right thing to do here is take the first value...assuming it to be the prod cast.
    if i==1%again this is just gathering up lats and longs for GAK 1 integer stations.
    LATGAK1=mean(parsed(i).GAK1.LTERlat);
    LONGGAK1=mean(parsed(i).GAK1.LTERlong);
    GAKLAT(1)=LATGAK1;
    GAKLONG(1)=LONGGAK1;
    end
GAKFW2_high(1)=FWGAK1;
allGAK=GAKFW2_high;
%GAKline=vertcat(FWGAK1,FWGAK);
%S = timerange('04/01','05/30');
% '2018-02-01','months'
%month(startPeriod) <= month(rowTime) and month(rowTime) <= month(endPeriod).
% dt_TGD=datetime(TGD);
for iii=1:highstation 
%FWinarray(71+i,iii)=allGAK(iii);%REMOVED AS WE ARE KEEPING MY ORIGINAL
%PARSED DATA SET OUT OF THIS, EXCPET WE GRAB LAT AND LONG WITH IT
end
end

T2=array2table(FWinarray);
T1=array2table(dt_TGD);
Y=TGD(:,1); %year column.used later for filling in the gaps. 
T=horzcat(T1,T2);
FWinarray_plus_y=horzcat(FWinarray,Y);
%T now represents an 80 by 15 array that, when there were multiple
%measurements for a station at a given time, those measurements are
%averaged.

GAKLAT(14:15)=[];%REmoving last two stations in the line which we are not plotting this time.
GAKLONG(14:15)=[];

JanT=T(Jan,:);%not acutaly used yet
FebT=T(Feb,:);
MarT=T(Mar,:);
AprT=T(Apr,:);
MayT=T(May,:);
JunT=T(Jun,:);
JulT=T(Jul,:);
AugT=T(Aug,:);
SepT=T(Sep,:);
OctT=T(Oct,:);
NovT=T(Nov,:);
DecT=T(Dec,:);

%here we make arrays for different months and pad everything so dimensions
%are the same.
monthabrv=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
for i=1:12
carray=eval(strcat(monthabrv(i),"T"));
if ~isempty(carray)
temp=eval(strcat("FWinarray_plus_y(",monthabrv(i),",:);"));
eval(strcat("array",monthabrv(i),"=temp.';"));
end
if isempty(carray)
    eval(strcat("array",monthabrv(i),"=[];"))
end

end

numyears=str2double(years);

%This loop creates a NaN temp 16 by 24, bottomost row populated by years.
%it then looks at each year in the current array, and maps the year with
%it's data to the Temp. This leaves years where no values were recorded as
%NaN columns. The result is arrays of equal dimension for the time span of
%24 years and 15 stations and a bottomost row that is the years. At the
%very end the year row is removed and the rows are flipped so station 15 is
%the topmost row. 
for i=1:length(monthabrv)
    temp=NaN(higheststation+1,24);%this forloop makes arrays with equal dimensions for each month.
    temp(higheststation+1,:)=numyears;
    currentarray=eval(strcat("array",monthabrv(i)));
if ~isempty(currentarray)

for ii=1:length(numyears)
    indx=find(currentarray(higheststation+1,:)==numyears(ii));
    if ~isempty(indx)
    temp(:,ii)=currentarray(:,indx);
    
    end
   
end
end
eval(strcat("array",monthabrv(i),"=flipud(temp(1:higheststation,:));"))
end





yaxis=(1:higheststation).'; %x and y axis of our space time conception

xaxis=(1:24).';

montharrays=NaN(length(yaxis),length(xaxis),12);
%all month arrays are put into this larger array.
for i=1:12
    eval(strcat("montharrays(:,:,i)=array",monthabrv(i),";"))
end

%getting matrix of averages (station by month)
avg_stationmonth=NaN(higheststation,12); %rows are station, columns are months
std_stationmonth=NaN(higheststation,12);
CI_stationmonth=NaN(higheststation,12);%rows are station, columns are months
rng_stationmonth=NaN(higheststation,12);
for i=1:12
    for ii=1:higheststation
        
        avg_stationmonth(ii,i)=mean(montharrays(ii,:,i),'omitnan'); %the average for a particular station and a particular month, of all the years, of freshwater height caculated from a specific integration depth. 
        %since we do not expect these values to be different, that is why
        %we determine the 95 % confidence interval for these measurements
        %that should be the same.
        std_stationmonth(ii,i)=std(montharrays(ii,:,i),'omitnan');
        rng_stationmonth(ii,i)=range(montharrays(ii,:,i),'omitnan');
        %this is the range of the of integrated freshwater heighe
        %measurents that occur for a particular station and a particular
        %year. We expect that these ranges would be the same, but of course
        %they are not. That is why we need to determine the 95% confidence
        %interval. Our measuremnet in this case is the range. just as our
        %measurement in the last one was the average. 
        N=length(montharrays(ii,~isnan(montharrays(ii,:,i)),i)); %gives us our sample size, N. y and err are y axis values and standard deviation of the sample
        ySEM=(std(montharrays(ii,:,i),'omitnan'))/sqrt(N); %this is the standard error.
    
    %then we apply the formula for the t distribution for 95% confidence interval. 
        CI95 = tinv([0.025 0.975], N-1);
        yCI95 = bsxfun(@times, ySEM, CI95(:));
        CI_stationmonth(ii,i)=yCI95(1,1);
    end
end

%Now we start to construct the figures
% colormap(turbo)
% labels=string([]); %for us in legend
% stations_avg_allyears=NaN(higheststation,12);
% for month=1:12
% arraymonth=montharrays(:,:,month);
% 
% 
% 
% number_of_colors = 1200; %Just make sure this is higher than your number of data points so that you get enough colors to represent the spectrum of your data.
% 
% 
% figure(1);
% Station_avg_allyears=NaN(length(yaxis),1);%makes a space for the averages for each station for a particular month.
% for i=1:length(Station_avg_allyears)
%    Station_avg_allyears(i)=mean(arraymonth(i,:),'omitnan');%if the whole year is Nan, meaning no data for that station was ever recorded for that month, e.g. there are no January data values for station 15 ever, that becomes a Nan.
% end
% stations_avg_allyears(:,month)=Station_avg_allyears;
% %assuming we have all the GAK locations for the month we can just use
% %locations from a fully sampled GAK cruise. Spring 2018 is one such cruise.
% 
% clf
% m_proj('Stereographic','longitudes',[-147.5],'latitudes',[59], 'radius', [3.5], 'rectbox', 'on'); %put map on our figure.
% m_gshhs_h('patch',[.5 .5 .5]);
% m_grid('xlabeldir','end','fontsize',10);
% 
% Station_avg_allyears=flipud(Station_avg_allyears); %flip the stations so that GAK1 is at the top, thus matching our latlong
% 
% Max2=max(avg_stationmonth,[],'all');%this is the max of averages of all months
% 
% Min2=min(avg_stationmonth,[],'all'); %this is the min of averages of all months
% %Generating color matrix for averaged values
% idx_in_colorbar2 = floor(1+(Station_avg_allyears - Min2) / (Max2 -Min2) * (number_of_colors-1)); %if a column is NaN, idx_in_colorbar is all NaN. 
% %if part of a column is NaN, that part of the idx_in_colorbar is also NaN.
% cm2 = turbo(number_of_colors);
% 
% matrix_with_rgb=NaN(length(idx_in_colorbar2),3);
% for i=1:length(idx_in_colorbar2)%this deals with the NaN values by replacing them with white in our matrix_with_rbg
%     if isnan(idx_in_colorbar2(i))
%        matrix_with_rgb(i,:)=[1,1,1];
%        continue
%     end
%     
%     if ~isnan(idx_in_colorbar2(i))
%     matrix_with_rgb(i,:)=cm2(idx_in_colorbar2(i),:); %this assigns the correct color for values that are not NaN.   
%     end
% 
% end
% 
% 
% for i=1:length(GAKLAT)
%     m_line(GAKLONG(i),GAKLAT(i),'marker','o','markersize',5,'color',matrix_with_rgb(i,:), 'markerfacecolor',matrix_with_rgb(i,:));
%     hold on
% end
% c=colorbar;
% c.TicksMode='manual';
% c.LimitsMode='manual';
% c.Limits=[0,1];
% c.Ticks=[0,.5,1];
% middle=Min2+(Max2-Min2)/2;
% Max2=round(Max2,1);
% Min2=round(Min2,1);
% middle=round(middle,1);
% c.TickLabels=[Min2,middle,Max2];
% c.Label.String = 'Freshwater Height (m)';
% title(strcat(monthabrv(month)," Map:Seward Line Average(97-20) Freshwater Height(m) in Upper ",depth,"m"))
% fileplace=strcat("C:\Users\Isaac\Documents\MATLAB\Research\Figures\FWcontent\Final\",depth,"\",monthabrv(month),depth,"Stero_avg_1997_2020_GAK_FWin.png");
% saveas(gcf,fileplace)
% hold off
% clf
% plot (stations,Station_avg_allyears)
% ax=gca;
% ax.XTickMode='manual';
% ax.XTickLabelMode='manual';
% ax.YTickMode='manual';
% ax.YTickLabelMode='auto';
% ax.XTick=1:1:higheststation;
% strstat=string(stations);
% ax.XTickLabel = strstat;
% 
% ax.YLimMode='manual';
% ax.YLim=[0,10];
% ax.YTick=1:1:10;
% 
% ax.XLimMode='manual';
% ax.XLim=[0,higheststation+1];
% ax.XTick=1:1:higheststation;
% 
% %strstation=string(1:10);
% %ax.YTickLabel=strstation;
% xlabel("Stations");
% ylabel("Avg (1997-2020) Freshwater height (m)");
% title(strcat(monthabrv(month),":Seward Line Average(97-20)Freshwater Height (m) in Upper ",depth, "m"))
% fileplace=strcat("C:\Users\Isaac\Documents\MATLAB\Research\Figures\FWcontent\Final\",depth,"\",monthabrv(month),depth,"plotted_avg_1997_2020_GAK_FWin.png");
% saveas(gcf,fileplace)
% 
% figure(2)
% 
% if sum(isnan(avg_stationmonth(:,month)))==higheststation %a score of 15 means that all stations are NaN so there is no usable data for the month.Don't plot that month, plot the rest an hold
% continue
% end
% hold on   
% 
% plot (stations,Station_avg_allyears)
% labels(end+1)=monthabrv(month);
% hold off
% 
% 
% end
% figure(2)
% ax=gca;
% ax.XTickMode='manual';
% ax.XTickLabelMode='manual';
% ax.YTickMode='manual';
% ax.YTickLabelMode='auto';
% ax.XTick=1:1:higheststation;
% strstat=string(stations);
% ax.XTickLabel = strstat;
% 
% ax.YLimMode='manual';
% ax.YLim=[0,10];
% ax.YTick=1:1:10;
% 
% ax.XLimMode='manual';
% ax.XLim=[0,higheststation+1];
% ax.XTick=1:1:higheststation;
% legend(labels)
% %strstation=string(1:10);
% %ax.YTickLabel=strstation;
% xlabel("Stations");
% ylabel("Avg (1997-2020) Freshwater Height (m)");
% title(strcat("All Months:Seward Line Average(97-20)Freshwater Height (m) to Upper ",depth, "m"))
% fileplace=strcat("C:\Users\Isaac\Documents\MATLAB\Research\Figures\FWcontent\Final\",depth,"\",depth,"allmonths_plotted_avg_1997_2020_GAK_FWin.png");
% saveas(gcf,fileplace)


% monthss=1:12;
% for i=1:15
%     
%     Station_avg_allyears=NaN(length(yaxis),1);%makes a space for the averages for each station for a particular month.
% for ii=1:length(Station_avg_allyears)
%    Station_avg_allyears(ii)=mean(arraymonth(ii,:),'omitnan');%if the whole year is Nan, meaning no data for that station was ever recorded for that month, e.g. there are no January data values for station 15 ever, that becomes a Nan.
% end
%    plot(monthss,a
%     
% end
fileplace=strcat("C:\Users\cfosadmin\Documents\MATLAB\Research\Data\NGALTER\",depth,"avgmonths.mat");
save(fileplace,'avg_stationmonth','std_stationmonth','CI_stationmonth','montharrays','rng_stationmonth')
 clear all
end