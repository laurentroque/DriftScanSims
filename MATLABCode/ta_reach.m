Drift Scan Sim
% REACH drift scan smiulations

%% Time and Frequency range
obsFreq = (50:150).*1e6; %% Observation frequency in Hz
LST = 0:2.5:357.5; %% LST in degrees (2.5deg = 10min resolution)
ploton = 1; %% 1==plot sky (slower with plotting), 0==no plot
cf = 1; %% 1==calculate CF, 0==do not calculate
fref = 50e6; %% reference frequency to calculate chromaticity function (CF) 
%rot = 0; %% antenna rotation in degrees

%% sky maps and antenna beams directories in l-m plane (pattern converted using pat_load3)
%pdir = 'C:\Users\nraza\Documents\MATLAB\xarray17_hera\Dipole_REACH\'; % beams directory
%pdir = 'C:\Users\nraza\Dropbox\REACH_sims\npy_beam_dipoleandbalun\'; % beam directory
pdir = '/media/ian/ssd2020/REACH/REACH_sims/npy_beam_dipoleandbalun/'; %Linux beam directory
%mdir = 'C:\Users\nraza\Dropbox\v7data\'; % maps directory 
mdir = '/media/ian/ssd2020/REACH/REACH_sims/ta_REACHmaps/'; %Linux map directory

%% load sky maps (LFSM) - 40:1:250MHz
load('maps_LFSM_40Mto250M_3deg.mat'); % High reolsution LFSM maps
mapf_indx = (40:250).*1e6; % map frequency index (do not change)

%% Initial Values
d2r = pi/180;
r2d = 1./d2r;
res = 512; %resolution for all maps and beams (512x512), note 512==0.35deg
ra = linspace(360,0,size(maps,2)); % Right Ascension (RA) in degrees
dec = linspace(-90,90,size(maps,1)); % Declination (DEC) in degrees
[RAi, DECi] = meshgrid(ra,dec);

% REACH 5 site
Lat = -30.836511; %-30.7612; %-30.721389;
Long = 21.370872; %21.4082; %21.428333;

% obsTime = linspace(0,23.59,142); % Observation time in decimal hours (e.g. ~10min resolution)
% for 0:360 LST range use these e.g. values but otherwise just set LST
% Day = 30;
% Month = 8;
% Year = 2018;
% Dec = Dec_a + Dec_b./60;        % Declination (DEC) in degrees
% Time = Time_h + Time_m./60;     % Time in decimal hours
% var=1;
% if Lat_a < 0
%     sss = -1;
% else
%     sss = 1;
% end
% Lat = Lat_a + sss.*Lat_b./60;       % Latitude in degrees
% Long = val*(Long_a + Long_b./60);   % Longitude in degress (West == -1*)

% %% Tables
% Tab_A = [ % Month Year(normal) Year(Leap)
%     1   0   0; 2   31  31; 3   59  60; 4   90  91; 5   120 121; 6   151 152;
%     7   181 182; 8   212 213; 9   243 244; 10  273 274; 11  304 305; 12  334 335];
%
% Tab_B = [ % Year Days
%     1998    -731.5; 1999    -366.5; 2000    -1.5; 2001    364.5; 2002    729.5;
%     2003    1094.5; 2004    1459.5; 2005    1825.5; 2006    2190.5; 2007    2555.5;
%     2008    2920.5; 2009    3286.5; 2010    3651.5; 2011    4016.5; 2012    4381.5;
%     2013    4747.5; 2014    5112.5; 2015    5477.5; 2016    5842.5; 2017    6208.5;
%     2018    6573.5; 2019    6938.5; 2020    7303.5; 2021    7669.5];
% %% Days from J2000
%     if ismember(Year,[2000 2004 2008 2012 2016 2020]) == 1
%         Days = obsTime(tidx)./24 + Tab_A(Month,3) + Day + Tab_B(Year-1997,2);
%     else
%         Days = obsTime(tidx)./24 + Tab_A(Month,2) + Day + Tab_B(Year-1997,2);
%     end
% %% Local Siderial Time (LST) and Hour Angle (HA)
%     LSTi = 100.46 + 0.985647*Days + Long + 15*obsTime(tidx); % Approximate formula
%     LSTi = mod(LSTi,360);% LST in degrees
%     % setappdata(0,'LST',LSTi);
%%

Tsky = zeros(length(obsFreq),length(LST));

if cf==1
    CF = zeros(length(obsFreq),length(LST));
end

if exist('rot','var')==0
    rot = 0;
end

if exist('fref','var')==0
    fref = 50e6; 
end

for fidx = 1:length(obsFreq)
    %% load antenna beams
    %beam_name = [pdir 'reach_dipole_pat_' num2str(obsFreq(fidx)/1e6) '.mat'];
    beam_name = [pdir 'reach_hexdip_pat_' num2str(obsFreq(fidx)/1e6) '.mat'];
    load(beam_name);
    if rot~=0
        imT = imTransform(Pat,rot*pi/180,1,'cubic');
        if ismember(res,size(imT,1))==1
            Pat = imT;
        else
            cen = floor(size(imT,1)/2);
            lim = round(res/2); clear Pat;
            Pat = imT(cen-lim+1:cen+lim,cen-lim+1:cen+lim);
        end
    end
    Pat(Pat<0.001) = 0; % set floor to zero (outside of l-m beam)
    
    if cf==1 % load ref beam
        %pp = load([pdir 'reach_dipole_pat_' num2str(fref/1e6) '.mat']);
        pp = load([pdir 'reach_hexdip_pat_' num2str(fref/1e6) '.mat']);
        beamref = pp.Pat;
        if rot~=0
            imT = imTransform(beamref,rot*pi/180,1,'cubic');
            if ismember(res,size(imT,1))==1
                beamref = imT;
            else
                cen = floor(size(imT,1)/2);
                lim = round(res/2); clear beamref;
                beamref = imT(cen-lim+1:cen+lim,cen-lim+1:cen+lim);
            end
        end
        beamref(beamref<0.001) = 0; % set floor to zero (outside of l-m beam)
    end
        
    
    for tidx = 1:length(LST)
        tic
        LSTi=LST(tidx);
        
        %% Sky map (UV) - need to generate and load first time around
        map_name = [mdir 'sky_map_uv_f' num2str(round(10*obsFreq(fidx)/1e6)/10) 'MHz_LST'  num2str(round(LST(tidx)*1000)/1000) 'deg.mat'];
        if exist(map_name,'file')==2 % check if map exists, then load
            sky = matread(map_name);
        else % otherwise generate map and save it
            %% Conversion to theta/phi
            HA = LSTi - RAi;
            HA = mod(HA,360);    % HA in degrees
            
            sin_ALT = sin(DECi.*d2r).*sin(Lat*d2r) + cos(DECi.*d2r).*cos(HA.*d2r).*cos(Lat*d2r);
            ALT = asin(sin_ALT).*r2d;
            ALT = real(ALT);
            
            cos_A = (sin(DECi.*d2r) - sin(ALT.*d2r).*sin(Lat*d2r))./(cos(ALT.*d2r).*cos(Lat*d2r));
            A = acos(cos_A).*r2d;
            A = real(A);
            clear sin_ALT cos_A
            
            temp = sin(HA.*d2r);
            AZ = zeros(size(temp));
            AZ(temp<0) = A(temp<0);
            AZ(temp>=0) = 360 - A(temp>=0);
            
            TH = 90-ALT;
            PH = AZ;
            clear AZ ALT A temp
            
            Y = fliplr(squeeze(maps(:,:,mapf_indx==obsFreq(fidx)))); % use LFSM map
            
            %choose top hemisphere (one way to asign imag to other plane)
            PH(TH>90) = 1j;
            Y(TH>90) = 1j;
            TH(TH>90) = 1j;
            
            u = sin(TH*d2r).*cos(PH*d2r);
            v = sin(TH*d2r).*sin(PH*d2r);
            
            sky = gridder(real(Y),real(u),real(v),res,res);
            matwrite(sky,map_name);

        end
        
        if ploton==1 % plot maps to check (when ploton==1)
            pcolor(linspace(-1,1,res),linspace(-1,1,res),log(sky)); shading interp;
            title(sprintf('LST = %0.3fh, freq = %0.1fMHz',LST(tidx)/15,obsFreq(fidx)/1e6),'fontsize',12); 
            hcb=colorbar; title(hcb,'log(K)'); drawnow;
        end
        
        %% main convolution code
        % sky = sky.*(sky_temp(f/1e9,0)/sky_temp(0.408,0));
        step = 2/length(Pat);
        corr = wcorr(length(Pat)); %cos(th) term;
        T2 = sum(sum(corr.*Pat)).*step^2;
        T1 = sum(sum(corr.*Pat.*sky)).*step^2;
        %     T2 = trapz(abs(trapz(Pat,2).*step)).*step;
        %     T1 = trapz(abs(trapz(Pat.*Yg,2).*step)).*step;
        Tsky(fidx,tidx) = T1/T2;
        
        %% Chromaticity function 
        if cf==1
            % load ref map (for a given LST)
            map_name2 = [mdir 'sky_map_uv_f' num2str(round(10*fref/1e6)/10) 'MHz_LST'  num2str(round(LST(tidx)*1000)/1000) 'deg.mat'];
            if exist(map_name2,'file')==2 % check if map exists, then load
                skyref = matread(map_name2);
            else
                disp('turn cf function off (cf=0), then run to generate the map at the reference frequency')
            end
            
            CF1 = sum(sum(corr.*Pat.*skyref)).*step^2;
            CF2 = sum(sum(corr.*beamref.*skyref)).*step^2;
            CF(fidx,tidx) = CF1/CF2;
        end
        tloop = toc;
        fprintf('Freq = %0.1fMHz, LST = %0.3fh, Loop Time = %0.2f\n',obsFreq(fidx)/1e6,LST(tidx)/15,tloop)
    end
end

%% plot Tsky as a function of LST and frequency
figure;
plot(LST./15,Tsky); grid on; shg
xlabel('LST /hours','fontsize',12)
ylabel('Antenna Temperature /K','fontsize',12)

if cf==1
    figure;
    plot(LST./15,CF); grid on; shg
    xlabel('LST /hours','fontsize',12)
    ylabel('Chromaticity Function','fontsize',12)
end
