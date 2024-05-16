function [x, z, t, dt, up, wp] = LS_stableStep(ustar, L, z_i, z0, C0, vs, x, ...
    z, t, up, wp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LS_STABLESTEP Computes particle velocities and positions in a stable 
% boundary layer using the Langevin equation.
%   [x, z, t, dt, up, wp] = LS_stableStep(ustar, L, z_i, z0, C0, vs, x, z, 
% t, up, wp)
%
%   Inputs:
%       ustar   - Friction velocity (m/s)
%       L       - Obukhov length (m)
%       z_i     - Boundary layer height (m)
%       z0      - Roughness length (m)
%       C0      - Lagrangian structure function constant
%       vs      - Particle settling velocity (m/s)
%       x       - Current x-coordinate of the particle (m)
%       z       - Current z-coordinate of the particle (m)
%       t       - Current time (s)
%       up      - Current horizontal velocity of the particle (m/s)
%       wp      - Current vertical velocity of the particle (m/s)
%
%   Outputs:
%       x       - Updated x-coordinate of the particle (m)
%       z       - Updated z-coordinate of the particle (m)
%       t       - Updated time (s)
%       dt      - Time step (s)
%       up      - Updated horizontal velocity of the particle (m/s)
%       wp      - Updated vertical velocity of the particle (m/s)
%
%   This function computes the particle velocities and positions in a
%   stable boundary layer using the Langevin equation. It calculates the
%   necessary wind statistics and Langevin coefficients based on the input
%   parameters.
%
%   The wind statistics are computed as follows:
%   - Mean horizontal velocity (ubar) is calculated using the Brost &
%     Hotslag stable profile (eq. 28 & 32).
%   - Horizontal velocity variance (sigu2) and its derivative (dsigu2) are
%     computed using Kantha (2000) eq. 3.5.6.
%   - Vertical velocity variance (sigw2) and its derivative (dsigw2) are
%     computed using a modified version of Kantha (2000) eq. 3.5.6.
%   - Covariance between horizontal and vertical velocities (uw) and its
%     derivative (duw) are computed using eq. 12.31 in Rodean (1996) and
%     eq. 36 in Nieuwstadt (1984).
%   - Dissipation rate (eps) is calculated using eq. 12.18 in Rodean (1996).
%
%   The Langevin coefficients (au, aw, bu, bw) are computed using the
%   formulation given by Aylor (2001). The time step (dt) is determined as
%   a fraction of the Lagrangian time scale (tau), which is computed using
%   the vertical velocity variance and dissipation rate.
%
%   Finally, the function updates the particle velocities and positions
%   using the computed increments based on the Langevin equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stable wind statistics
% Brost & Hotslag Stable ubar eq. 28 & 32, % a = 1; b = 2/3; c = 5; d = 0.35;
psi = 1*z/L + 2/3*(z/L - 5/.35).*exp(-.35*z/L) + 2/3*5/.35;                 
ubar = ustar/0.4*(log(z/z0) + psi);                         
ubar = max(ubar,0);

% Kantha (2000) eq. 3.5.6
sigu2 = 4*ustar^2*(1 - z/z_i).^(3/2);                       
dsigu2 = -6*ustar^2*sqrt(1 - z/z_i)/z_i;                    
sigw2 = 3*ustar^2*(1 - z/z_i).^(3/2);                      
dsigw2 = -9*ustar^2*sqrt(1 - z/z_i)/(2*z_i);                    

% Eq. 12.31 in Rodean 1996 and eq. 36 in Nieuwstadt 1984
q = 0;  % for stable conditions
uw = -ustar^2*(1 - z/z_i)^(3/2 - q);                           
duw = -(ustar^2*(1 - z/z_i)^(1/2 - q)*(q - 3/2))/z_i;

% Rodean 1996 eq. 12.18
eps = ustar^3./(0.4*z).*(1 + 3.5*z/z_i*z_i/L).*(1 - 0.85*z/z_i).^(3/2); 

% Luhar 1996
tau = 2*sigw2/(C0*eps);
dt = 0.02*tau;


%% Langevin coefficients
% Aylor (2001)
A = 2*(sigu2*sigw2 - uw^2);
bu = (C0*eps)^.5;
bw = bu;
au = 1/A*bu^2*(uw*wp - sigw2*up) + 0.5*duw + 1/A*(sigw2*dsigu2*wp*up - ...
    uw*dsigu2*wp^2 - uw*duw*up*wp + sigu2*duw*wp^2);
aw = 1/A*bw^2*(uw*up - sigu2*wp) + 0.5*dsigw2 + 1/A*(sigw2*duw*up*wp - ...
    uw*duw*wp^2 - uw*dsigw2*up*wp + sigu2*dsigw2*wp^2);

%%  find u,w velocity increments (langevin equation)

% Find u,w velocity increments (langevin equation)
dup = au*dt + bu*randn*sqrt(dt);
dwp = aw*dt + bw*randn*sqrt(dt);

% Increment u,w velocities
up = up + dup;
wp = wp + dwp;

% Find x,z position increment
dx = (up + ubar)*dt;
dz = (wp - vs)*dt;

% Increment t, x, and z
t = t + dt;
x = x + dx;
z = z + dz;




end





