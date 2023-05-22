clear all
close all
%%%%%%%%%

cd 'C:\Users\cfosadmin\Documents\MATLAB\Research\code\SMAP_LTER_comp'
Single_Station_timeline='GAK15';
numfiles=12;
%%%%%%
fpath='C:\Users\cfosadmin\Documents\MATLAB\Research\Data\NGALTER\Raw_LTER\';
%below is prep for loading files
names=["NGA_SKQ201810S_ctd_L3_v2";...
   "NGA_WSD201807_ctd_L3_v2";...
   "NGA_TGX201809_ctd_L3_v2";...
   "NGA_TGX201904_ctd_L3_v1";...
   "NGA_SKQ201915S_ctd_L3_v1";...
   "NGA_TGX201909_ctd_L3_v1";...
   "NGA_SKQ202006S_ctd_L2_v1";...
   "NGA_SKQ202010S_ctd_L2_v1";...
   "NGA_SKQ202012S_ctd_L2_v1";...
   "NGA_SKQ202106S_ctd_L2_v1";...
   "NGA_SKQ202110S_ctd_L2_v1";...
   "NGA_TGX202109_ctd_L2_v1"];...

%below is creation of structure to be filled
for i=1:numfiles
 struct_array(i)=struct('Data',[],'name',[],'SMAPcoord',[],'Station',[],'SMAPtime_indx',[],'SMAPSSSindxed',[],'cruisename',[],...
     'GAK',struct('LTER',[],'LTERSal50',[],'LTERlat',[],'LTERlong',[],'LTERFW50',[],'SMAP',[],'Station',strings,'Numbersonly',[],'CastNumbersonly',[]),...
     'GAK1',struct('LTER',[],'LTERSal50',[],'LTERlat',[],'LTERlong',[],'LTERFW50',[],'SMAP',[],'Station',strings,'Numbersonly',[],'CastNumbersonly',[]),...
     'MID',struct('LTER',[],'LTERSal50',[],'LTERlat',[],'LTERlong',[],'LTERFW50',[],'SMAP',[],'Station',strings,'Numbersonly',[],'CastNumbersonly',[]),...
     'KOD',struct('LTER',[],'LTERSal50',[],'LTERlat',[],'LTERlong',[],'LTERFW50',[],'SMAP',[],'Station',strings,'Numbersonly',[],'CastNumbersonly',[]),...
     'CS',struct('LTER',[],'LTERSal50',[],'LTERlat',[],'LTERlong',[],'LTERFW50',[],'SMAP',[],'Station',strings,'Numbersonly',[],'CastNumbersonly',[]),...
     'RES',struct('LTER',[],'LTERSal50',[],'LTERlat',[],'LTERlong',[],'LTERFW50',[],'SMAP',[],'Station',strings,'Numbersonly',[],'CastNumbersonly',[]),...
     'GEO',struct('LTER',[],'LTERSal50',[],'LTERlat',[],'LTERlong',[],'LTERFW50',[],'SMAP',[],'Station',strings,'Numbersonly',[],'CastNumbersonly',[]),...
     'PWS',struct('LTER',[],'LTERSal50',[],'LTERlat',[],'LTERlong',[],'LTERFW50',[],'SMAP',[],'Station',strings,'Numbersonly',[],'CastNumbersonly',[]));
 
 parsed=struct_array;
end
 %assigning approriate seasons to each filename, while also capturing
 %original file name
 cruises=["Spring 2018","Summer 2018","Fall 2018","Spring 2019","Summer 2019","Fall 2019","Spring 2020","Summer 2020","Fall 2020","Spring 2021","Summer 2021","Fall 2021"];
 for i = 1:numfiles
    str = strcat(fpath,' ', names(i),'.csv');
    
    
    myFilenames{i} = str;
    parsed(i).name=names(i);
    parsed(i).cruisename=cruises(i);
 end
myfilenames2=string(myFilenames);

%below enters the LTER station data and station text, into structure. The
%raw data must be in the form of a matlab table, stored in the same
%directory or the path for it must be specified somehow. parseLTER (which also must be either in the same directory or 
%have a path that goes to it) takes
%that table and extracts avg lat, avg long, calculated FW input, Surface
%salinity, time and the text for station number(eg. GAK5I).

for i=1:numfiles
    
    LTERtable=myfilenames2(i);
    [parsed(i).Data,parsed(i).Station,parsed(i).fullcolumnData,parsed(i).castnumbersall]=deep_parseLTERV2(LTERtable);
end


%below takes the lat long for the station and finds the closest match for
%the SMAP data, using the findSMAPpoint function. findSMAPpoint must be in
%the same directory or otherwise be specified.Also SMAP data for the area must 
%be downloaded. This file assumes that the SMAP data is contained in an SMAP folder, but the path can be changed within
%the file. SMAPcoord contains long and
%lat for the SMAP point, and the associated sea surface salinty through all
%time, as well as the time stamps themselves. All of which goes into the
%structure as another structure.

% for i=1:numfiles
%     datalat=parsed(i).Data(:,2);
%     datalong=parsed(i).Data(:,1);
%     
%    parsed(i).SMAPcoord=findSMAPpoint(datalong,datalat);
%   
%     
% end

%below looks though the SMAP time and finds that timestamp that is closest
%to the LTER time, and enters it into the structure
% %%
% for i=1:numfiles
% t2=datetime(2015,1,1,0,0,parsed(i).SMAPcoord.times);
% temp=NaN(length(parsed(i).Data(:,5)),1);
% for ii=1:length(parsed(i).Data(:,5))
%     
%     t1 = datetime(parsed(i).Data(ii,5),'ConvertFrom','datenum');
%     nearest_index = interp1(t2, 1:length(t2), t1, 'nearest');
%     temp(ii)=nearest_index;
% 
% end
% parsed(i).SMAPtime_indx=temp;
% end
%%
%--------------------------------uncomment if working with SMAP again.
%below then finds the SMAP SSS assoicated to that SMAP time selected above
%and enters it into the structure.
% for i=1:numfiles
%    temp=NaN(length(parsed(i).Data(:,5)),1); 
%    for ii=1:length(parsed(i).Data(:,5))
%     
%     temp(ii)=parsed(i).SMAPcoord.tsSSS(ii,parsed(i).SMAPtime_indx(ii));
%    end
%     parsed(i).SMAPSSSindxed=temp;
% end
%-------------------------------------
%make plots

%     F(inn,1)=avg_Lat;
%     F(inn,2)=avg_Long;
%     F(inn,3)=FWinput;
%     F(inn,4)=SurfaceSP;
%     F(inn,5)=time;
%scatter(parsed(1).Data(:,4),parsed(1).SMAPtime

% figure(1)
% hold on
% for i=1:numfiles
%     if isempty(rmmissing(parsed(i).SMAPSSSindxed)) %don't plot data wherein the entire satallite data is missing
%         continue
%     end
% scatter(parsed(i).Data(:,4),parsed(i).SMAPSSSindxed)
% 
% end
%lsline

% hold off


% subplot(2,2,1)
% xlabel('Station SSS (PSU)')
% ylabel('Matched SMAP SSS (PSU)')
% legend({'Spring 2018','Summer 2018','Fall 2018','Spring 2019','Fall 2019'},'Location','northwest')
% title('Correlation:NGA LTER to SMAP, zoom') 
% xlim([0 35])
% ylim([0 35])
% 
% subplot(2,2,2)
% xlabel('Station SSS (PSU)')
% ylabel('Matched SMAP SSS (PSU)')
% legend({'Spring 2018','Summer 2018','Fall 2018','Spring 2019','Fall 2019'},'Location','northwest')
% title('Correlation:NGA LTER to SMAP, zoom') 
% xlim([20 35])
% ylim([20 35])
%mystery why we don't see parsed(4) which is summer 2019.
%Solved:the summer 2019 mystery. Beginning of July of 2019 is a gap in the
%SMAP data for the region at least. Unclear as to why. 
% 
% temp=NaN(numfiles,2);
% recordinx=NaN(numfiles,2);
% for i=1:numfiles
% 
% inx=strmatch(Single_Station_timeline,parsed(i).Station,'exact');
% if isempty(inx)
%     continue
% end
% if length(inx)>1
%     inx=inx(1);
% end
% recordinx(i,1)=i;
% recordinx(i,2)=inx;
% temp(i,1)=parsed(i).Data(inx,4);
% temp(i,2)=parsed(i).Data(inx,5);
% end

%this is code to fix the scenario where the last run doesn't contain an inx
%value.
% 
% recordinx=num2cell(recordinx);
% m = any(cellfun('isempty', recordinx), 2);  % Any cell in row is empty
% recordinx(m, :) = [];
% 
% recordinx=cell2mat(recordinx);
% recordinx=rmmissing(recordinx);
% inx=recordinx(end,2);
% i=recordinx(end,1);

% figure(2)
% hold on
% %assuming geographic location of Gak 9 cast has not changed much we choose
% %the SMAP coordinate data for location that aligned best with GAK 9 in the
% %most recent cruise and then create the SMAP timeseries for that mark.
% %Also we have a bug where if the last cruise did not contain the searched
% %for parameter then inx is zero. We want the recent record of the dataset.
% %Also the parsed(i) needs to match with the inx.
% 
% TotalSattimes=datetime(2015,1,1,0,0,parsed(5).SMAPcoord.times);
% plot(TotalSattimes,parsed(i).SMAPcoord.tsSSS(inx,:))
% scatter(datetime(temp(:,2),'ConvertFrom','datenum'),temp(:,1), 'r')
% hold off
% xlabel('time(date)')
% ylabel('Salinity (PSU)')
% titlestring=strcat('Comparison of SMAP and NGA LTER SSS of ',Single_Station_timeline);
% title(titlestring)









%This is the start of the per cruise figure of subplots
%important note, for individual stations, the SMAPSSS may have a NAN value,
%in which case it will not show up. this is the case for GAK1 SMAP
%stations.


% figure(3)
% %look for paticular station by changing this:
% pstation="GAK2";
% 
% for i=1:numfiles
%     if isempty(rmmissing(parsed(i).SMAPSSSindxed)) %don't plot data wherein the entire satallite data is missing, this is Summer 2019
%         continue
%     end
% h(i)=subplot(3,2,i);
% tlen=length(parsed(i).Data(:,4));
% hold on
% scatter(parsed(i).Data(:,4),parsed(i).SMAPSSSindxed,'w');
% lsline
% for ii=1:tlen
% scatter(parsed(i).Data(ii,4),parsed(i).SMAPSSSindxed(ii),'w');
% text(parsed(i).Data(ii,4),parsed(i).SMAPSSSindxed(ii),parsed(i).Station(ii))
% if parsed(i).Station(ii)==pstation
% text(parsed(i).Data(ii,4),parsed(i).SMAPSSSindxed(ii),parsed(i).Station(ii),'Color','red')
% end
% end
% hold off
% %lsline
% 
% %titlestring=strcat(shname(i),':Fresh water content (m^3) to upper 50 m');
% title(parsed(i).cruisename);
% xlabel('Station SSS (PSU)')
% ylabel('Matched SMAP SSS (PSU)')
% 
% end

% saveas(figure(3),'C:\Users\Isaac\Documents\MATLAB\Research\Figures\SMAP_LTER compare\sttext_cruisesplit_alltransects.png')
% 
% %This figure divides up by transect (aka line) as well as cruise.




%Below takes what is already in the structure, restructures it, and puts it
%into another part of the structure.  It does not use an outside function. 

for i=1:numfiles

F2=parsed(i).Data;
% SMAP=parsed(i).SMAPSSSindxed;
Stationstrg=parsed(i).Station;
Lines=["Kodiak";"Seward";"Middleton";"Prince William Sound";"Kayak Island"];

for ix= 1: length(F2) %this is the number of casts in the cruise.
    Gcheck= strmatch("GAK1", Stationstrg(ix), 'exact'); %first checks if GAK1 is there.
if contains(Stationstrg(ix),"GAK")==1
    if ~isempty(Gcheck)
%     if isempty(parsed(i).GAK1.LTER(end))
%     parsed(i).GAK1.LTER(end)=F2(ix,4);
%     parsed(i).GAK1.SMAP(end)=SMAP(ix);
%     parsed(i).GAK1.Station(end)=Stationstrg;
%     continue
%     end
        parsed(i).GAK1.LTER(end+1)=F2(ix,4);
        parsed(i).GAK1.LTERSal50(end+1)=F2(ix,6);
        parsed(i).GAK1.LTERlat(end+1)=F2(ix,2);
        parsed(i).GAK1.LTERlong(end+1)=F2(ix,1);
        parsed(i).GAK1.LTERFW50(end+1)=F2(ix,3);
%         parsed(i).GAK1.SMAP(end+1)=SMAP(ix);
        parsed(i).GAK1.Station(end+1)=Stationstrg(ix);
        parsed(i).GAK1.Station=string(rmmissing(cellstr(parsed(i).GAK1.Station)));
        parsed(i).GAK1.CastNumbersonly(end+1)=parsed(i).castnumbersall(ix);
        continue
    end
%     if isempty(parsed(i).GAK.LTER(end))
%     parsed(i).GAK.LTER(end)=F2(ix,4);
%     parsed(i).GAK.SMAP(end)=SMAP(ix);
%     parsed(i).GAK.Station(end)=Stationstrg;
%     continue
%     end
parsed(i).GAK.LTER(end+1)=F2(ix,4);
parsed(i).GAK.LTERSal50(end+1)=F2(ix,6);
parsed(i).GAK.LTERlat(end+1)=F2(ix,2);
parsed(i).GAK.LTERlong(end+1)=F2(ix,1);
parsed(i).GAK.LTERFW50(end+1)=F2(ix,3);
% parsed(i).GAK.SMAP(end+1)=SMAP(ix);
parsed(i).GAK.Station(end+1)=Stationstrg(ix);
parsed(i).GAK.Station=string(rmmissing(cellstr(parsed(i).GAK.Station)));
parsed(i).GAK.CastNumbersonly(end+1)=parsed(i).castnumbersall(ix);

continue
end

if contains(Stationstrg(ix),"MID")==1
%     if isempty(parsed(i).MID.LTER(end))
%     parsed(i).MID.LTER(end)=F2(ix,4);
%     parsed(i).MID.SMAP(end)=SMAP(ix);
%     parsed(i).MID.Station(end)=Stationstrg;
%     continue
%     end
parsed(i).MID.LTER(end+1)=F2(ix,4);
parsed(i).MID.LTERSal50(end+1)=F2(ix,6);
parsed(i).MID.LTERlat(end+1)=F2(ix,2);
parsed(i).MID.LTERlong(end+1)=F2(ix,1);
parsed(i).MID.LTERFW50(end+1)=F2(ix,3);
% parsed(i).MID.SMAP(end+1)=SMAP(ix);
parsed(i).MID.Station(end+1)=Stationstrg(ix);
parsed(i).MID.Station=string(rmmissing(cellstr(parsed(i).MID.Station)));
parsed(i).MID.CastNumbersonly(end+1)=parsed(i).castnumbersall(ix);
continue
end

if contains(Stationstrg(ix),"CS")==1
%     if isempty(parsed(i).CS.LTER(end))
%     parsed(i).CS.LTER(end)=F2(ix,4);
%     parsed(i).CS.SMAP(end)=SMAP(ix);
%     parsed(i).CS.Station(end)=Stationstrg;
%     continue
%     end
parsed(i).CS.LTER(end+1)=F2(ix,4);
parsed(i).CS.LTERSal50(end+1)=F2(ix,6);
parsed(i).CS.LTERlat(end+1)=F2(ix,2);
parsed(i).CS.LTERlong(end+1)=F2(ix,1);
parsed(i).CS.LTERFW50(end+1)=F2(ix,3);
% parsed(i).CS.SMAP(end+1)=SMAP(ix);
parsed(i).CS.Station(end+1)=Stationstrg(ix);
parsed(i).CS.Station=string(rmmissing(cellstr(parsed(i).CS.Station)));
parsed(i).CS.CastNumbersonly(end+1)=parsed(i).castnumbersall(ix);
continue
end

if contains(Stationstrg(ix),"KOD")==1
%     if isempty(parsed(i).KOD.LTER(end))
%     parsed(i).KOD.LTER(end)=F2(ix,4);
%     parsed(i).KOD.SMAP(end)=SMAP(ix);
%     parsed(i).KOD.Station(end)=Stationstrg;
%     continue
%     end
parsed(i).KOD.LTER(end+1)=F2(ix,4);
parsed(i).KOD.LTERSal50(end+1)=F2(ix,6);
parsed(i).KOD.LTERlat(end+1)=F2(ix,2);
parsed(i).KOD.LTERlong(end+1)=F2(ix,1);
parsed(i).KOD.LTERFW50(end+1)=F2(ix,3);
% parsed(i).KOD.SMAP(end+1)=SMAP(ix);
parsed(i).KOD.Station(end+1)=Stationstrg(ix);
parsed(i).KOD.Station=string(rmmissing(cellstr(parsed(i).KOD.Station)));
parsed(i).KOD.CastNumbersonly(end+1)=parsed(i).castnumbersall(ix);
continue
end

if contains(Stationstrg(ix),"RES")==1
%     if isempty(parsed(i).RES.LTER(end))
%     parsed(i).RES.LTER(end)=F2(ix,4);
%     parsed(i).RES.SMAP(end)=SMAP(ix);
%     parsed(i).RES.Station(end)=Stationstrg;
%     continue
%     end
parsed(i).RES.LTER(end+1)=F2(ix,4);
parsed(i).RES.LTERSal50(end+1)=F2(ix,6);
parsed(i).RES.LTERlat(end+1)=F2(ix,2);
parsed(i).RES.LTERlong(end+1)=F2(ix,1);
parsed(i).RES.LTERFW50(end+1)=F2(ix,3);
% parsed(i).RES.SMAP(end+1)=SMAP(ix);
parsed(i).RES.Station(end+1)=Stationstrg(ix);
parsed(i).RES.Station=string(rmmissing(cellstr(parsed(i).RES.Station)));
parsed(i).RES.CastNumbersonly(end+1)=parsed(i).castnumbersall(ix);
continue
end

if contains(Stationstrg(ix),"TEST")==1
continue
end

if contains(Stationstrg(ix),"GEO")==1
%     if isempty(parsed(i).GEO.LTER(end))
%     parsed(i).GEO.LTER(end)=F2(ix,4);
%     parsed(i).GEO.SMAP(end)=SMAP(ix);
%     parsed(i).GEO.Station(end)=Stationstrg;
%     continue
%     end
parsed(i).GEO.LTER(end+1)=F2(ix,4);
parsed(i).GEO.LTERSal50(end+1)=F2(ix,6);
parsed(i).GEO.LTERlat(end+1)=F2(ix,2);
parsed(i).GEO.LTERlong(end+1)=F2(ix,1);
parsed(i).GEO.LTERFW50(end+1)=F2(ix,3);
% parsed(i).GEO.SMAP(end+1)=SMAP(ix);
parsed(i).GEO.Station(end+1)=Stationstrg(ix);
parsed(i).GEO.Station=string(rmmissing(cellstr(parsed(i).GEO.Station)));
parsed(i).GEO.CastNumbersonly(end+1)=parsed(i).castnumbersall(ix);
continue
end

%     if isempty(parsed(i).PWS.LTER(end))
%     parsed(i).PWS.LTER(end)=F2(ix,4);
%     parsed(i).PWS.SMAP(end)=SMAP(ix);
%     parsed(i).PWS.Station(end)=Stationstrg;
%     continue
%     end
parsed(i).PWS.LTER(end+1)=F2(ix,4);
parsed(i).PWS.LTERSal50(end+1)=F2(ix,6);
parsed(i).PWS.LTERlat(end+1)=F2(ix,2);
parsed(i).PWS.LTERlong(end+1)=F2(ix,1);
parsed(i).PWS.LTERFW50(end+1)=F2(ix,3);
% parsed(i).PWS.SMAP(end+1)=SMAP(ix);
parsed(i).PWS.Station(end+1)=Stationstrg(ix);
parsed(i).PWS.Station=string(rmmissing(cellstr(parsed(i).PWS.Station)));
parsed(i).PWS.CastNumbersonly(end+1)=parsed(i).castnumbersall(ix);





%m_line(Sat(ix,2),Sat(ix,1),'marker','o','markersize',5,'color','k');
%This works because everything has stayed in order. Each lat and long as we
%go from 1 to 88 corresponds with the 1 to 88 colors (representing our
%data) in matrix_with_rgb
end


%saveas(figure(2),'Gak15timeseries.png')

end




%Below we take the unorganized GAK line stations and organize them from
%near shore to off shore. We do that that for all values (lat, long, FW
%content, etc) of GAK.
Transects=["GAK1","GAK","MID","KOD","CS","RES","GEO","PWS"];



for i=1:numfiles
for ix=1:length(Transects)
    if isempty(eval(strcat("parsed(i).",Transects(ix),".LTER")))==1 %this makes it so we skip for example parsed(4).CS since it contines no values.
        continue
    end
    
    if Transects(ix)=="PWS" %Skips the prince William sound section since ordering that doesn't really matter. More of just a scattershot there.
        continue
    end
    if Transects(ix)=="GEO"
        continue
    end
    Tran=Transects(ix);
    Unorganized=eval(strcat("parsed(i).",Tran,".Station"));
    OrgLTER=eval(strcat("parsed(i).",Tran,".LTER"));
    OrgLTERSal50=eval(strcat("parsed(i).",Tran,".LTERSal50"));
    OrgLTERFW50=eval(strcat("parsed(i).",Tran,".LTERFW50"));
%     OrgSMAP=eval(strcat("parsed(i).",Tran,".SMAP"));
    OrgLat=eval(strcat("parsed(i).",Tran,".LTERlat"));
    OrgLong=eval(strcat("parsed(i).",Tran,".LTERlong"));
    OrgCastNumbers=eval(strcat("parsed(i).",Tran,".CastNumbersonly"));
    Organized=Unorganized;
    Numbersonly=cellstr(Unorganized);
    index=NaN(1,length(Unorganized)); %making a blank vector to fill.
for ii=1:length(Unorganized)
    Station=Unorganized(ii);
    test=char(Station);
    if test(end)=='I'
        test(end)=[];
        Staion=string(test);
        Station=strcat(Station,'.5');
    end
    s=Station;
    Station=cellfun(@(x) x(regexp(x,'[^a-zA-Z]')),s,'un',0);
    Numbersonly(ii)=Station;
    index(ii)=ii;
end
Numbersonly=cellfun(@str2num, Numbersonly);%turns out cell array into a double array
Num_and_Indi=[Numbersonly;index]; %combines the index and the station number
Num_and_Indi=Num_and_Indi.';%transpose
SortedNumbers=sortrows(Num_and_Indi);   %sorts according to station. the index is sorted alongside so we know where the station number used to be.
index2=[1:1:length(SortedNumbers(:,1))]; %list of indices consectutive from 1 to number of stations by 1.
index2=index2.'; %transposed
Num_Indi_Indi2=horzcat(SortedNumbers,index2);% so now in columns we have station numbers, changed index, and original index.
Organized(Num_Indi_Indi2(:,3))=Unorganized(Num_Indi_Indi2(:,2));
%before the above line, oranganized and unorganized are the same. They are
%the stations unstripped. If we replace the entire set of stations numbers unstripped
%such that index 1 (the left most postion) is equal to indice 20, we have changed the station number from GAK15 to be GAK2
%. In this way indice 2 for (formerly GAK15) become indice 19 which is
%GAK3. Thus our stations get organized, and we obtain a map that can be
%used on our other variables.
%
%Applying the same changing over to the rest of that year's station's
%variable. On the left, is the final organized product. On the right, it's
%labelled organized but SECRETLY it is unorganized. That's just so I don't
%have to have an organized and unorganized for each variable, only the
%ones above where we do the dirty work.

OrgLTER(Num_Indi_Indi2(:,3))=OrgLTER(Num_Indi_Indi2(:,2));
OrgLTERSal50(Num_Indi_Indi2(:,3))=OrgLTERSal50(Num_Indi_Indi2(:,2));
OrgLTERFW50(Num_Indi_Indi2(:,3))=OrgLTERFW50(Num_Indi_Indi2(:,2));
% OrgSMAP(Num_Indi_Indi2(:,3))=OrgSMAP(Num_Indi_Indi2(:,2));
OrgLat(Num_Indi_Indi2(:,3))=OrgLat(Num_Indi_Indi2(:,2));
OrgLong(Num_Indi_Indi2(:,3))=OrgLong(Num_Indi_Indi2(:,2));
OrgCastNumbers(Num_Indi_Indi2(:,3))=OrgCastNumbers(Num_Indi_Indi2(:,2));
    
eval(strcat("parsed(i).",Tran,".LTER=OrgLTER"));
eval(strcat("parsed(i).",Tran,".LTERSal50=OrgLTERSal50"));
eval(strcat("parsed(i).",Tran,".LTERFW50=OrgLTERFW50"));
% eval(strcat("parsed(i).",Tran,".SMAP=OrgSMAP"));
eval(strcat("parsed(i).",Tran,".LTERlat=OrgLat"));
eval(strcat("parsed(i).",Tran,".LTERlong=OrgLong"));
eval(strcat("parsed(i).",Tran,".Station=Organized"));
eval(strcat("parsed(i).",Tran,".Numbersonly=SortedNumbers(:,1).'"));
eval(strcat("parsed(i).",Tran,".CastNumbersonly=OrgCastNumbers"));
end

end


valuesrange=[];
for i=1:numfiles-1
valuesrange = vertcat(valuesrange,parsed(i).Data(:,3)); % combining all fw data for all cruises. 
end
values_min = min(valuesrange); %Alternatively 1 % Finding the min and max of my data. This will
%be the range of the colorbar
values_max = max(valuesrange); %Alternatively 5.


save('C:\Users\cfosadmin\Documents\MATLAB\Research\Data\NGALTER\ALL_NGALTERdata_feb2023.mat','parsed','values_max','values_min');
%---- the R-square statistic, the F statistic and p value
  %  for the full model, and an estimate of the error variance.