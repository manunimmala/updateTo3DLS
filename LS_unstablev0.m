function [up0, wp0] = LS_unstablev0(ustar, wstar, L, z_i, z, np)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LS_UNSTABLEV0 Computes initial velocity fluctuations for particles in an 
% unstable boundary layer.
%   [up0, wp0] = LS_unstablev0(ustar, wstar, L, z_i, z, np)
%
%   Inputs:
%       ustar   - Friction velocity (m/s)
%       wstar   - Convective velocity scale (m/s)
%       L       - Obukhov length (m)
%       z_i     - Boundary layer height (m)
%       z       - Initial height of particles (m)
%       np      - Number of particles
%
%   Outputs:
%       up0     - Initial horizontal velocity fluctuations (m/s)
%       wp0     - Initial vertical velocity fluctuations (m/s)
%
%   This function computes the initial velocity fluctuations for particles
%   in an unstable boundary layer using the vertical velocity variance and
%   skewness derived from Boehm et al. 2008 and the bi-Gaussian velocity
%   distribution from Luhar 1996.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% From Boehm et al. 2008 
sigw2_CBL = 1.7*wstar^2*(z/z_i)^(2/3).*(1-0.9*(z/z_i)).^(4/3);
sigw2_neutral = ustar^2*(1.7-z/z_i);
sigw2 = ((1 - exp(z/L))*wstar^3.*sigw2_CBL + 5*exp(z/L)*ustar^3.*sigw2_neutral)./...
    ((1 - exp(z/L))*wstar^3 + 5*exp(z/L)*ustar^3);
w3 = 1.2*(z./z_i).*(1 - (z./z_i)).^(3/2)*wstar^3;

% From Luhar 1996 closure 
S = w3./sigw2.^(3/2); % skewness definition, stated immediately before sec 2.2
m = 2/3*S^(1/3); % eq. 5, Luhar 1996
r = (1 + m^2)^3*S^2/((3 + m^2)^2*m^2);  % eq. 6e, Luhar 1996
A = .5*(1-(r/(4+r))^.5);    % eq. 6c, Luhar 1996
B = 1 - A;   % eq. 6d, Luhar 1996
sigA = (sigw2*B/(A*(1+m^2)))^.5; % eq. 6a, Luhar 1996
sigB = (sigw2*A/(B*(1+m^2)))^.5; % eq. 6b, Luhar 1996
wB = m*sigB;   % eq. 4b, Luhar 1996
wA = m*sigA;    % eq. 4a, Luhar 1996


mu = [wA; -wB];

% Concatenate two(=k) 1x1 covariance matrices via the third dimension 
% Sigma(:,:,k) contains a covariance matrix for each gaussian component k
sigma = cat(3, sigA^2, sigB^2);    
p = [A,B];

% Horizontal initial velocity
up0 = randn(np,1)*.6*wstar; % luhar 2002, eq 15a * randn for sigu

% Vertical initial velocity
% Make a bi-gaussian gmdistribution using the above statistics
% Select np random numbers from it to initialize each particle
wp0 = random(gmdistribution(mu,sigma,p),np);


end

