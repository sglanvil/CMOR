% June 16, 2023

clear; clc; close all;

% ---------------------- USER SPECIFY ----------------------
inVarName='OMEGA';
outVarName='Wzm';
inDir="/glade/scratch/strandwg/QBOI/waccm-SC.QBOi.EXP2.LA.001/atm/proc/tseries/month_1/";
outDir="/glade/scratch/sglanvil/QBOi/data/";
fileName=inDir+"waccm-SC.QBOi.EXP2.LA.001.cam.h0."+inVarName+".197901-208001.nc";
outFile="waccm-SC.QBOi.EXP2.LA.001.cam.h0."+outVarName+".197901-208001.nc";
PSfile="waccm-SC.QBOi.EXP2.LA.001.cam.h0.PS.197901-208001.nc";
PSLfile="waccm-SC.QBOi.EXP2.LA.001.cam.h0.PSL.197901-208001.nc";

lon=ncread(fileName,'lon');
lat=ncread(fileName,'lat');
lev=ncread(fileName,'lev');
ilev=ncread(fileName,'ilev');
time=ncread(fileName,'time');
date=ncread(fileName,'date');
hyam=permute(repmat(ncread(fileName,'hyam'),1,length(lon),length(lat)),[2 3 1]);
hybm=permute(repmat(ncread(fileName,'hybm'),1,length(lon),length(lat)),[2 3 1]);
PSin=ncread(inDir+PSfile,'PS');
PSin=permute(repmat(PSin,1,1,1,length(lev)),[1 2 4 3]);
P0=1000;

varIn=ncread(fileName,inVarName);
varZM=NaN(length(lat),length(ilev),length(time)); % allocate space
for itime=1:length(time)
    itime
    var0=squeeze(varIn(:,:,:,itime));
    PS0=squeeze(PSin(:,:,:,itime));    
    hybp=(hyam.*P0+hybm.*PS0)./100;
    varOut=NaN(length(lon),length(lat),length(ilev));
    for ilat=1:length(lat)
        for ilon=1:length(lon)
            varOut(ilon,ilat,:)=interp1(squeeze(hybp(ilon,ilat,:)),...
                squeeze(var0(ilon,ilat,:)),ilev);
        end
    end
    varZM(:,:,itime)=squeeze(mean(varOut,1,'omitnan'));
end

if strcmp(inVarName,'OMEGA')==1
    ilevRep=100*permute(repmat(ilev,1,length(lat),length(time)),[2 1 3]);
    varZM=-varZM.*7000./ilevRep;
end

%%

ncName=sprintf(outDir+outFile);
cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
ncid = netcdf.create(ncName,cmode);
%Define the dimensions
dimidlon = netcdf.defDim(ncid,'lon',length(lon));
dimidlat = netcdf.defDim(ncid,'lat',length(lat));
dimidlev = netcdf.defDim(ncid,'lev',length(lev));
dimidilev = netcdf.defDim(ncid,'ilev',length(ilev));
dimidtime = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
%Define IDs for the dimension variables (pressure,time,varsitude,...)
lon_ID=netcdf.defVar(ncid,'lon','double',[dimidlon]);
lat_ID=netcdf.defVar(ncid,'lat','double',[dimidlat]);
lev_ID=netcdf.defVar(ncid,'lev','double',[dimidlev]);
ilev_ID=netcdf.defVar(ncid,'ilev','double',[dimidilev]);
time_ID=netcdf.defVar(ncid,'time','double',[dimidtime]);
var_ID=netcdf.defVar(ncid,outVarName,'float',[dimidlat dimidilev dimidtime]);
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid,lon_ID,lon);
netcdf.putVar(ncid,lat_ID,lat);
netcdf.putVar(ncid,lev_ID,lev);
netcdf.putVar(ncid,ilev_ID,ilev);
netcdf.putVar(ncid,time_ID,0,length(time),time);
netcdf.putVar(ncid,var_ID,varZM);
netcdf.close(ncid)

% date_ID=netcdf.defVar(ncid,'date','int',[dimidtime]);
% netcdf.putVar(ncid,date_ID,date);


