function [up0, wp0] = LS_stablev0(ustar, z_i, z, np)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LS_STABLEV0 Computes initial velocity fluctuations for particles in a stable 
% boundary layer.
%   [up0, wp0] = LS_unstablev0(ustar, wstar, L, z_i, z, np)
%
%   Inputs:
%       ustar   - Friction velocity (m/s)
%       z_i     - Boundary layer height (m)
%       z       - Initial height of particles (m)
%       np      - Number of particles
%
%   Outputs:
%       up0     - Initial horizontal velocity fluctuations (m/s)
%       wp0     - Initial vertical velocity fluctuations (m/s)
%
%   This function computes the initial velocity fluctuations for particles
%   in an a stable boundary layer using the vertical and horizontal velocity 
%   variance from Kantha and Clayson (2000).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% From Kantha and Clayson (2000)
sigu = 2*ustar*(1 - z/z_i).^(3/4);
sigw = 1.732*ustar*(1 - z/z_i).^(3/4); 

% Select from gaussian distribution with mean 0 and variance sigu/sigw
up0 = sigu*randn(np,1);
wp0 = sigw*randn(np,1);
end
