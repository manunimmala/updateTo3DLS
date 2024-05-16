function [egrid, ngrid, zgrid, egridConstant, ngridConstant,... 
    zgridConstant, pgrid, depgrid, egridCellSize, ngridCellSize,zgridCellSize]... 
    = LS_makeGrid(emin, emax, nmin, nmax, zmin, zmax, negrid, nngrid, nzgrid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LS_MAKEGRID Creates grids for particle residence times and deposition.
% [ngrid, egrid, zgrid, ngridConstant, egridConstant, zgridConstant, pgrid,
% depgrid, ngridCellSize, egridCellSize, zgridCellSize] =
% LS_makeGrid(nmin, nmax, emin, emax, zmin, zmax, nngrid, negrid, nzgrid)
%
% Inputs:
% nmin - Minimum north-coordinate of the domain (m)
% nmax - Maximum north-coordinate of the domain (m)
% emin - Minimum east-coordinate of the domain (m)
% emax - Maximum east-coordinate of the domain (m)
% zmin - Minimum z-coordinate of the domain (m)
% zmax - Maximum z-coordinate of the domain (m)
% nngrid - Number of grid cells in the north-direction (must be greater
% than 1)
% negrid - Number of grid cells in the east-direction (must be greater
% than 1)
% nzgrid - Number of grid cells in the z-direction (must be greater
% than 1)
%
% Outputs:
% ngrid - North-coordinates of the grid centers (m)
% egrid - East-coordinates of the grid centers (m)
% zgrid - Z-coordinates of the grid centers (m)
% ngridConstant - Constant used for computing grid indices in north-direction
% egridConstant - Constant used for computing grid indices in east-direction
% zgridConstant - Constant used for computing grid indices in z-direction
% pgrid - Grid for particle residence times (s)
% depgrid - Grid for particle deposition (particles/m)
% ngridCellSize - Size of each grid cell in the north-direction (m)
% egridCellSize - Size of each grid cell in the east-direction (m)
% zgridCellSize - Size of each grid cell in the z-direction (m)
%
% This function creates grids for tracking particle residence times and
% deposition. It computes the grid cell sizes, grid coordinates, and
% constants used for efficiently computing grid indices. The function
% also initializes the particle residence time grid (pgrid) and the
% deposition grid (depgrid).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ngridCellSize = (nmax-nmin)/(nngrid-1);
egridCellSize = (emax-emin)/(negrid-1);
zgridCellSize = (zmax-zmin)/(nzgrid-1);

ngrid = linspace(nmin, nmax-ngridCellSize, nngrid-1);
egrid = linspace(emin, emax-egridCellSize, negrid-1);
zgrid = linspace(zmin, zmax-zgridCellSize, nzgrid-1);

ngridConstant = ngrid(1)/ngridCellSize-1;
egridConstant = egrid(1)/egridCellSize-1;
zgridConstant = zgrid(1)/zgridCellSize-1;

pgrid = zeros(length(egrid),length(ngrid),length(zgrid));
depgrid = pgrid(:,1);

end