load('hexdip.mat');
fmhz = 50:200;
tidx = find(squeeze(dat(:,1,1)) <= 90); % find range inside hemisphere (ignoring back beam)
th = squeeze(dat(tidx,1,1));
ph = squeeze(dat(tidx,2,1));
TH = reshape(th,91,360);
PH = reshape(ph,91,360);
theta = linspace(0.0001,89.9999,91);
phi = linspace(0.00001,359.9999,360);
u = sin(theta*pi/180)'*cos(phi*pi/180);
v = sin(theta*pi/180)'*sin(phi*pi/180);
pnts = 512; %grid points
[xnew,ynew] = meshgrid(linspace(-1,1,pnts),linspace(-1,1,pnts));

for fidx = 1:101
dd = squeeze(dat(tidx,3,fidx));
D = reshape(dd,91,360);
Pat = griddata(u,v,D,xnew,ynew,'cubic');
Pat = Pat./max(Pat(:));
Pat(isnan(Pat)) = min(Pat(:))+0.000001;

pcolor(10*log10(Pat)); shading interp; colorbar;
title(sprintf('Frequency %dMHz',fmhz(fidx))); shg
pause(0.1);

save(['reach_hexdip_pat_' num2str(fmhz(fidx)) '.mat'],'Pat');
    
[fmhz(fidx) stats(Pat)]
end