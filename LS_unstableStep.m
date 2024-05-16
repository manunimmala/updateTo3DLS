function  [x, z, t, dt, up, wp] = LS_unstableStep(ustar, wstar, L, z_i, z0, ...
    C0, vs, x, z, t, up, wp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LS_UNSTABLESTEP Computes particle velocities and positions in an unstable 
% boundary layer using the Langevin equation.
%   [x, z, t, dt, up, wp] = LS_unstableStep(ustar, wstar, L, z_i, z0, C0, 
% vs, x, z, t, up, wp)
%
%   Inputs:
%       ustar   - Friction velocity (m/s)
%       wstar   - Convective velocity scale (m/s)
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
%   This function computes the particle velocities and positions in an
%   unstable boundary layer using the Langevin equation. It implements the
%   closure scheme described in Luhar (1996) and extends it to include
%   particle settling velocity as in Boehm et al. (2005).
%
%   The function first computes the vertical velocity variance (sigw2) and
%   its derivative (dsigw2) using a merged parameterization from Boehm &
%   Aylor (2008). The parameterization combines the convective boundary
%   layer (CBL) and neutral boundary layer contributions to the variance.
%
%   Next, the horizontal wind statistics are computed. The mean horizontal
%   velocity (ubar) is calculated using the logarithmic wind profile with a
%   stability correction term (phi). The horizontal velocity variance
%   (sigu2) is obtained from Luhar (2002).
%
%   The function then computes the necessary parameters for the bi-Gaussian
%   vertical velocity distribution and the Langevin coefficients. It
%   calculates the skewness (S), the closure scheme parameters (m, r, A, B),
%   and the Gaussian distribution parameters (sigA, sigB, wA, wB) for
%   updrafts and downdrafts.
%
%   Finally, the function updates the particle velocities and positions
%   using the computed increments based on the Langevin equation. The time
%   step (dt) is determined as a fraction of the Lagrangian time scale (tau).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute wind statistics
% Vertical velocity variance 
% Merged CBL/SL parameterization from Boehm & Aylor 2008, Eq. 3
sigw2_CBL = 1.7*wstar^2*(z/z_i)^(2/3).*(1-0.9*(z/z_i)).^(4/3);
dsigw2_CBL = (17*wstar^2*(1 - (9*z)/(10*z_i))^(4/3))/(15*z_i*(z/z_i)^(1/3))...
    - (51*wstar^2*(1 - (9*z)/(10*z_i))^(1/3)*(z/z_i)^(2/3))/(25*z_i);

sigw2_neutral = ustar^2*(1.7-z/z_i);
dsigw2_neutral = -ustar^2/z_i;

sigw2 = ((1 - exp(z/L))*wstar^3.*sigw2_CBL + 5*exp(z/L)*ustar^3.*sigw2_neutral)./...
    ((1 - exp(z/L))*wstar^3 + 5*exp(z/L)*ustar^3);
dsigw2 = - (wstar^3*(exp(z/L) - 1)*dsigw2_CBL - 5*ustar^3*exp((z/L))*dsigw2_neutral...
    - (5*ustar^3*exp((z/L))*sigw2_neutral)/L + (wstar^3*exp((z/L))*sigw2_CBL)/L)/(5*ustar^3*exp((z/L))...
    - wstar^3*(exp(z/L) - 1)) - ((5*ustar^3*exp((z/L))*sigw2_neutral...
    - wstar^3*sigw2_CBL*(exp(z/L) - 1))*((5*ustar^3*exp((z/L)))/L - ...
    (wstar^3*exp((z/L)))/L))/(5*ustar^3*exp((z/L)) - wstar^3*(exp(z/L) - 1))^2;

% Mean horizontal wind velocity and horizontal velocity variance
ex = (1-15*z/L)^0.25;
phi = -2*log((1+ex)/2) - log((1+ex^2)/2) + 2*atan(ex) - pi/2;
ubar = max(ustar/.4*(log(z/z0) + phi),0);  
sigu2 = (0.6*wstar)^2;   %Luhar 2002, eq. 15a


% Vertical velocity skewness
% Parameterization from Luhar 2002 eq. 15 d
w3 = 1.2*(z./z_i).*(1 - (z./z_i)).^(3/2)*wstar^3;
dw3 = (6*wstar^3*(1 - z/z_i)^(3/2))/(5*z_i) - (9*wstar^3*z*sqrt(1 - z/z_i))/(5*z_i^2);

% Merged SL/CBL turbulence dissipation rate from Boehm & Aylor 2008 eq. 2
eps = 0.4*wstar^3/z_i + ustar^3 * (1 - z/z_i) * (1 - 15*z/L)^(-1/4)/(.4*z); % eq. 2

% Luhar 1996
tau = 2*sigw2/(C0*eps);
dt = 0.02*tau; 

%% Computes parameters for the bi-gaussian vertical velocity pdf
% Using the closure scheme described in Luhar 1996.
%
% Then compute parameters necessary to calculate the langevin coefficients. The
% first are the gaussian distributions for updrafts and downdrafts and the
% combined bi-gaussian pdf. See Luhar 1989, 1996 and Rodean 1996.
% Boehm et al. 2005 extends Luhar 1996 to include a particle settling
% velocity vs in the phi term.


S = w3./sigw2.^(3/2); % skewness definition, stated immediately before sec 2.2
dS = dw3/sigw2^(3/2) - (3*w3*dsigw2)/(2*sigw2^(5/2));


% new closure scheme
m = 2/3*S^(1/3); % eq. 5, Luhar 1996
dm = (2*dS)/(9*S^(2/3));

r = (1 + m^2)^3*S^2/((3 + m^2)^2*m^2);  % eq. 6e, Luhar 1996
dr = (2*S*(m^2 + 1)^3*dS)/(m^2*(m^2 + 3)^2) +...
    (6*S^2*(m^2 + 1)^2*dm)/(m*(m^2 + 3)^2)...
    - (4*S^2*(m^2 + 1)^3*dm)/(m*(m^2 + 3)^3)...
    - (2*S^2*(m^2 + 1)^3*dm)/(m^3*(m^2 + 3)^2);

A = .5*(1-(r/(4+r))^.5);    % eq. 6c, Luhar 1996
dA = -(dr/(r + 4) - (r*dr)/(r + 4)^2)/(4*sqrt(r/(r + 4)));

B = 1 - A;   % eq. 6d, Luhar 1996
dB = -dA;

sigA = (sigw2*B/(A*(1+m^2)))^.5; % eq. 6a, Luhar 1996
dsigA = (B*dsigw2/(A*(m^2 + 1)) + sigw2*dB/(A*(m^2 + 1))...
    - B*sigw2*dA/(A^2*(m^2 + 1))...
    - 2*B*m*sigw2*dm/(A*(m^2 + 1)^2))/(2*sqrt((B*sigw2)/(A*(m^2 + 1))));

sigB = (sigw2*A/(B*(1+m^2)))^.5; % eq. 6b, Luhar 1996
dsigB = ((A*dsigw2)/(B*(m^2 + 1)) + (sigw2*dA)/(B*(m^2 + 1))...
    - (A*sigw2*dB)/(B^2*(m^2 + 1)) - ...
    (2*A*m*sigw2*dm)/(B*(m^2 + 1)^2))/(2*sqrt((A*sigw2)/(B*(m^2 + 1))));


wB = m*sigB;   % eq. 4b, Luhar 1996
dwB = dsigB*m + dm*sigB;

wA = m*sigA;    % eq. 4a, Luhar 1996
dwA = dsigA*m + dm*sigA;


P_A = (2*pi*sigA^2)^(-1/2) * exp(-(wp - wA)^2/(2*sigA^2));
P_B = (2*pi*sigB^2)^(-1/2) * exp(-(wp + wB)^2/(2*sigB^2));
P_E = A*P_A + B*P_B;

Q = A*(wp - wA)/sigA^2 * P_A + B*(wp + wB)/sigB^2 * P_B;

% Boehm 2005, phi for particles with a settling velocity
phi = -.5 * (A * dwA + (wA - vs)*dA) * erf((wp - wA)/(sqrt(2) * sigA)) + ...
    sigA * (A * dsigA * (wp*(wp - vs)/sigA^2 + 1)...
    + A * (wp - vs)/sigA^2 * (sigA * dwA - wA * dsigA) + sigA * dA) * P_A ...
    +.5 * (B * dwB + (wB + vs)*dB) * erf((wp + wB)/(sqrt(2) * sigB)) + ...
    sigB * (B * dsigB * (wp*(wp - vs)/sigB^2 + 1)...
    + B * (wp - vs)/sigB^2 * (-sigB * dwB + wB * dsigB) + sigB * dB) * P_B;

%% Langevin coefficients
% Equation 11, Luhar 1996, section 5
aw =  (-C0*eps*Q/2 + phi)/P_E;
bw = (C0*eps)^.5;

% Extend to 2D using luhar 2002
au = -up*C0*eps/(2*sigu2);   %Luhar 2002, eq. 13a
bu = bw;

%% Langevin equation

% Find u,w velocity increments (langevin equation)
dup = au*dt + bu*randn*sqrt(dt);
dwp = aw*dt + bw*randn*sqrt(dt);

% Increment u,w velocities
up = up + dup;
wp = wp + dwp;

% Find x,z position increment
dx = (up + ubar)*dt;
dz = (wp-vs)*dt;

% Increment t, x, and z
x = x + dx;
z = z + dz;
t = t + dt;


end

