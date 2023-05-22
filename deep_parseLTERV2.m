function [F2,statnum5,gridded,castnumbers]= deep_parseLTERV2(LTERtable)
%Record. Revisions seperated by %--------%
%2/7/2023 Isaac revised the code to extract PAR.
%---------------
%important: If surface volume changes,  V1 and ref_p must change.

%S=%Dataset
%Calculate absolute salinity of each casts
%Parse data into seperate casts??
%Average the upper 50 meters of salinity of each cast

%Reference salinty of NGA via Weingartner is 33.8
%How much freshwater had to be added to this 1 by 1 meter 50 meter deep
%column of water to get the current salinity?
%Assume  that the water column is well mixed. 
%working through a problem set: lets say ref salinity of a 1 x 1 x50 water
%column is 33.8. Measured is 32
%How much fresh water is needed to dilute 33.8 to 32

%-----------------------

%We have 34678 g of salt dispersed in 1 cubic meters of fresh water (derived from salinity 33.8) 
%This is in fact salt water. This is our initial concentration
%How much water do we add to get to 32828.8 g of salt dispersed in 1 cubic
%meter of fresh water (aka fresher salt water, also dervied from final salinity of 32)
%? This is our final concentration.

%Known: use the fact that c1V1 = c2V2
%Known:33.8 g/kg * 1026 kg/m^3 = 34678 g/ m^3
%Known: 32 g/ kg * 1025.9 kg/m^3 = 32828.8 g/ m^3
%Known: V1 = 50
%Known: V2 = 50 + x
%Unknown: x (the water added)
%Solving:
%(c1V1/c2) - 50 = x

%x=2.8164 cubic meters of fresh water is added to the water column to
%reduce the salinty from 33.8 to 32.

%---------------

% Our code will follow a similar logic. 


%Task 1 import the dataset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=readtable(LTERtable);
%/Users/isaacr/Documents/MATLAB/LTERdata
%Task 2 Calculate SA ect %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Numcasts=table2array(S(:,8));
DateStrings=table2array(S(:,4));
ttime = datetime(DateStrings,'InputFormat','yyyy-MM-dd''T''HH:mm:SS');
tttime=datenum(ttime);
statnum=table2array(S(:,2));

%finding the variables according to their names. Some names differ by a
%case, hence strcmpi.
salcol = find(strcmpi(S.Properties.VariableNames, 'Salinity__psu_'), 1);
loncol = find(strcmpi(S.Properties.VariableNames, 'Longitude__decimal_degrees_east_'), 1);
latcol = find(strcmpi(S.Properties.VariableNames, 'Latitude__decimal_degrees_north_'), 1);
castcol = find(strcmpi(S.Properties.VariableNames, 'Cast'), 1);
presscol = find(strcmpi(S.Properties.VariableNames, 'Pressure__dbar_'), 1);
tempcol = find(strcmpi(S.Properties.VariableNames, 'Temperature__C_'), 1);
temp2col = find(strcmpi(S.Properties.VariableNames, 'Temperature2__C_'), 1);
condcol = find(strcmpi(S.Properties.VariableNames, 'Conductivity__S_m_'), 1);
florcol = find(strcmpi(S.Properties.VariableNames, 'Fluorescence__mg_m3_'), 1);
sigcol = find(strcmpi(S.Properties.VariableNames, 'Density__kg_m3_sigmat_'), 1);
parcol = find(strcmpi(S.Properties.VariableNames, 'par__umol_photons_m2_sec_'),1);
%Develop a search algorithm.



AA=table2array(S(:,[loncol,latcol,castcol,presscol,tempcol,temp2col,condcol,salcol,florcol,sigcol]));
AA2=table2array(S(:,parcol));
%AA=table2array(S(:,[5,6,8,10,11,12,13,salcol]));
AA=[AA tttime AA2]; %appending time.


%long ,lat ,cast, pressure_dbar, temperature_C, temperature2_C,
%conductivity_Sm, salinity psu, flouresence, sigdensity, time
% 1      2    3         4            5                 6             7,
% 8,9,10,11
lines=length(AA);
B = eye( length(AA) , 13 ); 


SP = AA(:,8); %this is the processed psu data.
SA = gsw_SA_from_SP(SP,AA(:,4),AA(:,1),AA(:,2));%think about removing due to high lats.
z = gsw_z_from_p(AA(:,4),AA(:,1));
CT = gsw_CT_from_t(SA,AA(:,5),AA(:,4));
B(:,1)=AA(:,1); %Long
B(:,2)=AA(:,2); %Lat
B(:,3)=SP;     %SP
B(:,4)=SA;     %SA
B(:,5)=AA(:,4); %Pressure
B(:,6)=AA(:,3); %Cast
B(:,7)=z;      %depth
B(:,8)=AA(:,5);      %normal temp
B(:,9)=AA(:,9); %flouresence.
B(:,10)=AA(:,10);%sigt
B(:,11)=AA(:,11); %TIME
B(:,12)=gsw_rho(SA,CT,AA(:,4))-1000; %better sigt

%B is the full cruise data, yet to be split into transects.
%making the station numbers be enterable into the array.
for ix=1:size(statnum,1)
    test=char(statnum(ix));
    if test(end)=='I'
        test(end)=[]; %takes away the I
        
        test=strcat(test,'.5'); %adds a .5 to represent it.
        statnum(ix)=cellstr(test);
    end
end
Station=cellfun(@(x) x(regexp(x,'[^a-zA-Z]')),statnum,'un',0);
Station=cellfun(@str2num,Station,'UniformOutput',false);


A=cellfun(@isempty,Station); %finds out what cells are empty (they are because that station signifer doesn't contain a number. GEO for example.
Station(A)={NaN}; %assigns all the empty cells to be nans


Station=cell2mat(Station);
B(:,13)=Station; %station numbers are assigned to the thirteenth row.
B(:,14)=AA(:,12);%PAR is tagged onto the end 2/7/2023
%{'t','s','dens','chl', 'particle', 'CDOM' };


Bsize=size(B,2);

%  define the grid
dp = 1; % grid spacing 0.5 meters
maxp = 300; %ceil( nanmax( p ));
pgrid = 0:dp:maxp;

% make pressure matrices
%gridded.p = pgrid';

%B parse
%Long ,   lat,    SP,       SA,        Pressure,   Cast,    depth, CT,time,
%sgt
%Task 3 Determine cast subsets. Also find what is in the upper 50 meters
F=NaN(Numcasts(end),6);
statnum4=string(NaN(Numcasts(end),1));%making a place for station number.
gridded=NaN(length(pgrid'),Bsize,Numcasts(end)); %this is our twelve variables, down to a 65 db, and each stack is a different cast.
count=0;
    for inn = 1: Numcasts(end) %for each cast of this cruise...
        Castin= find(B(:,6)==inn & B(:,5)<=300);%finds all the rows that correspond to the first cast AND are less than 65. We choose 65 because that is the max for acrobat profiles.
        
        if isempty(Castin) %skip cast if there is no >50 record for a cast, aka the cast doesn't exist.
            continue
        end
        count=count+1;
        datin=B(Castin,:); %it's already bin averaged to 1 db before it even gets to matlab....
        %pin=B(pin,5);
        %binout = binaverage( pin, datin, pgrid' );
        castnumbers(count)=inn;
        
        maxdepth=max(abs(B(Castin(:),7)));%calculate the max depth.
       
        gridslab=NaN(length(pgrid'),Bsize);
        Map=~isnan(datin);
        for i=1:Bsize
     
        currentcolumn=Map(:,i);
        gridslab(currentcolumn,i)=datin(currentcolumn,i);
        
        
        end
        if isempty(gridslab) %sometimes the entire cast is just empty, in which case we skip.
            continue
        end
        gridded(:,:,inn)=gridslab;
        
        statnum2=string(NaN(length(Castin),1));
     
        for innn = 1:length(Castin)
       
           C(innn,:)=B(Castin(innn),:); %this is putting the data for the cast for all values into a matix c
           statnum2(innn,:)=string(statnum(Castin(innn),:)); %put the full data
        end
        %Task 4 Calculate needed averages
     
        avg_SP=mean(C(:,3));
        avg_Lat=mean(C(:,2));
        avg_Long=mean(C(:,1));
        avg_CT=mean(C(:,8));
        SurfaceSP=C(1,3);
        time=C(1,9);
        statnum3=statnum2(1,1);
        %Task 5 Calcuating additional FW input to this region (in cubic m)
        
        Ref_SP=33.8; %Reference salinty of NGA via Weingartner is 33.8% This is SP!!!
        Ref_p=maxdepth./2; %max depth is determined by either being ~50 (technically the closest pressure to 50,usualy comming in at 49.somthing.Ref p is then just the mid point of max depth. in shallow areas max depth adjusts.prevening overestimates of say the copper river.
        Ref_rho=gsw_rho(Ref_SP,avg_CT,Ref_p); %We don't have a reference CT to use so we are using the CT in the data.
        conc1=Ref_SP*Ref_rho;
        %depth=C(:,7); This line and next is method to get the measured volume
        %(-1)*depth(end); It doesn't handle the NaNs well/might not be needed... 
        V1= maxdepth; %see comments for Ref_p
        Current_rho=gsw_rho(avg_SP,avg_CT,Ref_p); %Using the same pressure (25 db) for both ref and current.
        conc2=avg_SP*Current_rho;
        FWinput=((conc1*V1)/conc2)-V1;
        
        %Task 6 Putting FWinput and Lat and Long into a matrix
        F(inn,1)=avg_Long;
        F(inn,2)=avg_Lat;
        F(inn,3)=FWinput;
        F(inn,4)=SurfaceSP;
        F(inn,5)=time;
        F(inn,6)=avg_SP;
        statnum4(inn,1)=statnum3;
       
        
          
    end


F2=rmmissing(F);
statnum5=rmmissing(statnum4);
castnumbers=rmmissing(castnumbers);
%B is the full cast.
%F2 is single value for a full cast. 
%statnum5 is the string of stations. 



end