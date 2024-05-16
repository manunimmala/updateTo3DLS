% script to test-run the LS stable model as I make changes


%% 5:47 PM on 5/15/24
% modified LS_makeGrid to take and output n, e, and z grid info to account
% for 3D. BUT relabeled all the inputs, so LS_stableModel should run as
% before. 
% Yes it works. 

clear all 
close all

% make up some input params
ustar = 1;
wstar = 0;
L = .1;
z_i = 100;
z0 = 0.03;
xmin = -10;
xmax = 1000;
zmin = 1;
zmax = 100;
np = 10000;
vs = .03;
x0 = 0;
h0 = 1;
nxgrid = 500;
nzgrid = 100;
C0 = 3;

% run it
[xgrid, zgrid, cgrid, depgrid] = LS_stableModel(ustar, wstar, L, ...
    z_i, z0, xmin, xmax, zmin, zmax, np, vs, x0, h0, nxgrid, nzgrid, C0);

% plot it 
figure
levels = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.1,1.2,1.4,1.6,2]/10000;
%contourf(squeeze(cgrid))
contourf(xgrid,zgrid,squeeze(cgrid)',levels)
ylabel('height (m)')
xlabel('downwind distance (m)')
colorbar
grid(gca,'on')

%% 
