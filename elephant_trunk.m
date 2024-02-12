% Compute the acoustic field in a bent waveguide with varying cross-sectional area
%
% This example is based on the article:
%
% Propagation in waveguides with varying cross section and
% curvature: a new light on the role of supplementary modes
% in multi-modal methods
%
% Agnès Maurel, Jean-François Mercier and Simon Félix3
%
% Proceedings of the Royal Society A: Mathematical, 
% Physical and Engineering Sciences, 2014
%
% http://dx.doi.org/10.1098/rspa.2014.0008
% 
% However, contrarilly to the original article, a 3D instead of a 2D 
% waveguide is considered here.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimension of the rectangle (m)
a = 0.05;   % width 
b = 0.03;   % height 

% number of propagation modes 
nModes = 50;

% admittance at the boundary
Y = 0.;

% sound speed
c = 340;

% wavenumber at which the field is computed
k = 3/2/b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the geometrical transformation matrices from analytical expression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[C, D, E, kn2, mn] = buildMatricesCDE_analyticalRectangle (a/b, 1, nModes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve multimodal problem with the elephant's trunk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r0 = 1.25*2*b;  % radius of the central part of the bent 
thetaL = 2.62;  % angle of the bent
L = thetaL*r0;  % length of curvilinear abscissa of the center of the guide
nX = 800;       % number of discretization points
X = linspace(-1*b, L+b, nX);    % generate curvilinear abscissa points
% compute scaling factor l, derivative of the scaling factor dl, and 
% curvatur kappa
[l, dl, kappa] = trunkElephant(X, L, b);    

%% compute the admittance matrices from the end to the beginning using an 
%  order 4 Magnus-Moebius scheme

Ad = zeros(nX, nModes, nModes);         % preallocate admittance matrices
K = diag(1i*sqrt((k*l(end))^2 - kn2));  % compute K matrix at the end
% initialize admittance matrix at the end
Ad(end,:,:) = K*l(end)^2;   

for ii = 1:nX-1
  xn = X(end - (ii-1));     % starting curvilinear abcissa
  xnp1 = X(end - (ii));     % ending curvilinear abcissa
  dx = xnp1 - xn;           % distance between abscissa
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute the first abscissa of the Magnus-Moebius scheme
  xH1 = xn + (0.5 - sqrt(3)/6)*dx;
  % compute the corresponding scaling l, scaling derivative dl, and
  % curvature kappa
  [l, dl, kappa] = trunkElephant(xH1, L, b);
  % build matrix K
  K = diag(1i*sqrt((k*l)^2 - kn2));
  % build matrix H
  H1 = [dl*E/l , (eye(nModes) - kappa*l*C)/(l^2);
 (K.^2 + kappa*l*(C*(k*l)^2 - D)) , -dl*transpose(E)/l];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute the second abscissa of the Magnus-Moebius scheme
  xH2 = xn + (0.5 + sqrt(3)/6)*dx;
  % compute the corresponding scaling l, scaling derivative dl, and
  % curvature kappa
  [l, dl, kappa] = trunkElephant(xH2, L, b);
  % build matrix K
  K = diag(1i*sqrt((k*l)^2 - kn2));
  % build matrix H
  H2 = [dl*E/l , (eye(nModes) - kappa*l*C)/(l^2);
 (K.^2 + kappa*l*(C*(k*l)^2 - D)) , -dl*transpose(E)/l];  
  
  % compute the exponential propagator
  On = expm(0.5*dx*(H1 + H2) + sqrt(3)*(dx^2)*(H2*H1 - H1*H2 )/12);
  
  % extract the submatrices of the propagator
  E1 = On(1:nModes, 1:nModes);
  E2 = On(1:nModes, nModes+1:end);
  E3 = On(nModes+1:end, 1:nModes);
  E4 = On(nModes+1:end, nModes+1:end);
  
  % Compute the admittance matrix at the next point
  Ad(nX-ii, :, :) = (E3 + E4*squeeze(Ad(nX-(ii-1),:,:)))/...
  (E1 + E2*squeeze(Ad(nX-(ii-1),:,:)));

end

%% compute the pressure amplitude from the beginning to the end using 
% an order 4 Magnus-Moebius scheme

p = zeros(nModes, nX);  % preallocate the pressure amplitude vectors
% intialize the first vector to generate a plane wave at the input
p(1,1) = 1;             

for ii = 1:nX-1
    
  xn = X(ii);       % starting curvilinear abcissa
  xnp1 = X(ii+1);   % ending curvilinear abcissa
  dx = xn - xnp1;   % distance between abscissa
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute the first abscissa of the Magnus-Moebius scheme
  xH1 = xn + (0.5 - sqrt(3)/6)*dx;
  % compute the corresponding scaling l, scaling derivative dl, and
  % curvature kappa
  [l, dl, kappa] = trunkElephant(xH1, L, b);
  % build matrix K
  K = diag(1i*sqrt((k*l)^2 - kn2));
  % build matrix H
  H1 = [dl*E/l , (eye(nModes) - kappa*l*C)/(l^2);
  (K.^2 + kappa*l*(C*(k*l)^2 - D)) , -dl*transpose(E)/l];  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute the second abscissa of the Magnus-Moebius scheme
  xH2 = xn + (0.5 + sqrt(3)/6)*dx;
  % compute the corresponding scaling l, scaling derivative dl, and
  % curvature kappa
  [l, dl, kappa] = trunkElephant(xH2, L, b);
  % build matrix K
  K = diag(1i*sqrt((k*l)^2 - kn2));
  % build matrix H
  H2 = [dl*E/l , (eye(nModes) - kappa*l*C)/(l^2);
  (K.^2 + kappa*l*(C*(k*l)^2 - D)) , -dl*transpose(E)/l];  
  
  % compute the exponential propagator
  On = expm(0.5*dx*(H1 + H2) + sqrt(3)*(dx^2)*(H2*H1 - H1*H2 )/12);
  
  % extract the submatrices of the propagator
  E1 = On(1:nModes, 1:nModes);
  E2 = On(1:nModes, nModes+1:end);
  
  % compuute the pressure amplitude vector at the next abscissa
  p(:,ii+1) = (E1 + E2*squeeze(Ad(ii+1,:,:)))\p(:,ii);
end

%% plot acoustic field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute analytical eigenmodes indexes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx = 800;   % nummber along x axis
% generate the x coordinate of the points
x = linspace(-b/4, 6*b, nx);
ny = 800;   % nummber along y axis
% generate the y coordinate of the points
y = linspace(-b, 3.5*b, ny);
% preallocate the corresponding transformed coordinates X and Y
Xp = zeros(nx,ny);
Yp = zeros(nx,ny);
% set the transverse transformed coordinate to 0 toextract the field on 
% the central (X, Y) plane
Zp = 0.;
% preallocate the acoustic pressure
Pp = zeros(nx,ny);
% preallocate Magnus-Moebius point index
idx = 0;


for ii = 1:nx
  for jj = 1:ny
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the (X, Y) transformed coordinates corresponding 
    % to the cartesian (x, y) coordinates 
    
    % if the points is in the small straight tube section
    if and(y(jj)<0, abs(x(ii)) <= 0.125)
      Xp(ii,jj) = y(jj);
      Yp(ii,jj) = 2*x(ii)/b+0.5;
    % if the points is in the bent
    else
      % compute the distance to the center of curvature (r0, 0)
      d = sqrt((x(ii) - r0)^2 + y(jj)^2);
      % compute the corresponding curvilinear abscissa in the bent
      Xp(ii,jj) = acos((r0 - x(ii))/d)*r0;
      % compute the corresponding scaling l, scaling derivative dl, and
      % curvature kappa
      [l, dl, kappa] = trunkElephant(Xp(ii,jj), L, b);
      % compute the corresponding transformed coordinate Yp in the bent
      Yp(ii,jj) = (r0-d)/l+0.5;
      
      % if the point is outside of the bent compute the transformed 
      % coordinates in the large straight tube section
      if Xp(ii,jj) > L
        if y(jj)>0
          theta =acos((r0-x(ii))/d) - L/r0;
        else
          theta =  pi - L/r0 + acos((x(ii) - r0)/d);
        end
        Xp(ii,jj) = L + d*sin(theta);
        Yp(ii,jj) = (r0 - d*cos(theta))/2/b+0.5;
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the acoustic pressure in the waveguide
    
    % check if the point is inside the waveguide, and if so,
    % identify the index of the closest Magnus-Moebius scheme point
    pointInside = false;
    if and(and(Xp(ii,jj) <= L+b, Xp(ii,jj) > 0), and(Yp(ii,jj)<1, Yp(ii,jj)>0))
      [toto, idx] = min(abs(X - Xp(ii,jj)));
      pointInside = true;
    elseif and(and(Xp(ii,jj) <0, Xp(ii,jj) > -b), and(Yp(ii,jj)<1, Yp(ii,jj)>0))
      idx = find(X>Xp(ii,jj), 1, "first");
      pointInside = true;
    end
    
    % if the point is inside compute the acoustic pressure
    if pointInside
      for mm = 1:nModes
        if mn(mm,1) == 0, Mm = 1; else, Mm = sqrt(2); end
        if mn(mm,2) == 0, Nn = 1; else, Nn = sqrt(2); end
        Pp(ii,jj) = Pp(ii,jj) + p(mm,idx)*(Mm*cos(mn(mm,1)*pi*Zp) + ...
        Nn*cos(mn(mm,2)*pi*Yp(ii,jj)));
      end
    end
  end
end
Pp(Pp == 0) = nan;

% compute the contour of the waveguide to display them
theta = linspace(0, thetaL, 100);
R1 = 3*b*((theta./thetaL).^3)/2 - 9*b*((theta./thetaL).^2)/4 + r0 -b/4;
translat = [sqrt(r0^2 + b^2)*cos(pi-thetaL-atan(b/r0)) - r0*cos(pi-thetaL),...
sqrt(r0^2 + b^2)*sin(pi-thetaL-atan(b/r0)) - r0*sin(pi-thetaL)];
R2 = -3*b*((theta./thetaL).^3)/2 + 9*b*((theta./thetaL).^2)/4 + r0 +b/4;
xLowerCont = [b/4, b/4, r0-R1.*cos(theta), r0-R1(end).*cos(theta(end)) + translat(1)];
yLowerCont = [-b, 0, R1.*sin(theta), R1(end).*sin(theta(end)) + translat(2)];
xUpperCont = [-b/4, -b/4, r0-R2.*cos(theta), r0-R2(end).*cos(theta(end)) + translat(1)];
yUpperCont = [-b, 0, R2.*sin(theta), R2(end).*sin(theta(end)) + translat(2)];

% plot the acoustic pressure
h = figure;
% imagesc(x,y,real(Pp)')
imagesc(x,y,20*log10(abs(Pp)'))
colorbar
hold on
plot(xUpperCont, yUpperCont, "k")
plot(xLowerCont, yLowerCont, "k")
axis xy
axis equal
print(h, "-dpng", "colormap.png");

% plot the contour of the acoustic pressure
h = figure; contour(x,y,real(Pp.'),100), axis equal
hold on
plot(xUpperCont, yUpperCont, "k")
plot(xLowerCont, yLowerCont, "k")
xlim([-b/2 6.5*b])
print(h, "-dpng", "contour.png");


%% Function to compute the scaling factor l, the derivative of the 
% scaling fator and the curvature at a specific survilinear abscissa X

function [l, dl, kappa] = trunkElephant(X, L, b)
  nX = length(X);       % number of requested abscissa
  
  % preallocate the output variables
  l = zeros(1,nX);      
  dl = zeros(1,nX);
  kappa = zeros(1,nX);
  
  % if the point is in the small tube
  l(X<0) = b/2;
  dl(X<0) = 0;
  kappa(X<0) = 0;
  
  % if the point is inside the bent
  l(and(X>0, X<L)) = 0.5*b*(1 + 9*(X(and(X>0, X<L))/L).^2 - ...
  6*(X(and(X>0, X<L))/L).^3);
  dl(and(X>0, X<L)) = 9*b*X(and(X>0, X<L)).*(1-X(and(X>0, X<L))/L)/(L^2);
  kappa(and(X>0, X<L)) = 1/2/1.25/b;
  
  % if the point is in the large tube
  l(X>L) = 2*b;
  dl(X>L) = 0;
  kappa(X>L) = 0;
end


%% A function to compute the propagation matrices C, D, and E, the squared 
% wavenumbers of the transverse modes kn2 and the modes indexes mn
%
% The expressions are given in the appendix of (blandin et al, 2022)

function [C, D, E, kn2, mn] = buildMatricesCDE_analyticalRectangle (a, b, nModes)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Compute analytical eigenmodes indexes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % sort modes indexes
  k2 = zeros(10000,1);
  mn = zeros(10000,2);
  idx = 0;
  for m=1:100
    for n=1:100
      idx = idx + 1;
      k2(idx,1) = ((m-1)/a)^2+((n-1)/b)^2;
      mn(idx,:) = [m,n];
    end
  end
  [k2, idx] = sort(k2);
  mn = mn(idx,:);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Build the geometrical transformation matrices from analytical expression
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% remove the +1 shift of the index (necessary for the fonction computing the analytical modes)
  mn = mn - 1;

  % *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  %% assembly of matrix C
  % *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 

  % preallocation 
  C = zeros(nModes, nModes);

  for ii=1:nModes
    for jj = 1:nModes
      m = mn(ii,1);
      n = mn(ii,2);
      o = mn(jj,1);
      p = mn(jj,2);
      
      if m == o
        if n == p
          C(ii,jj) = b/2;
        elseif and(n ~= 0, p ~= 0)
          C(ii,jj) = b*((cos((n+p)*pi)-1)/((n+p)^2) + ...
          (cos((n-p)*pi)-1)/((n-p)^2))/(pi^2);
        elseif and(n == 0, p ~= 0)
          C(ii,jj) = sqrt(2)*b*(cos(p*pi)-1)/((p*pi)^2);
        else
          C(ii,jj) = sqrt(2)*b*(cos(n*pi)-1)/((n*pi)^2);
        end
      end
      
    end
  end

  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  %% assembly of matrix D
  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 

  % preallocation 
  D = zeros(nModes, nModes);

  for ii = 1:nModes
    for jj = 1:nModes
      m = mn(ii,1);
      n = mn(ii,2);
      o = mn(jj,1);
      p = mn(jj,2);
      
      % Compute the first part of the integral
      I1 = 0;
      if and( m == o, and( m ~= 0, o ~= 0))
        if n == p
          I1 = 0.5*b*((m*pi/a)^2);
        elseif and(n ~= p, and(n ~= 0, p ~= 0))
          I1 = b*((cos((n+p)*pi)-1)/((n+p)^2) + ...
          (cos((n-p)*pi)-1)/((n-p)^2))*((m/a)^2);
        elseif and(n ~= p, and(n == 0, p ~= 0))
          I1 = sqrt(2)*b*(cos(p*pi)-1)*((m/a/p)^2);
        elseif and(n ~= p, and(n ~= 0, p == 0))
          I1 = sqrt(2)*b*(cos(n*pi)-1)*((m/a/n)^2);
        end
      end
      
      % Compute the second part of the integral
      I2 = 0;
      if m == o
        if and(n == p, and(n ~= 0, p ~= 0))
          I2 = ((n*pi)^2)/2/b;
        elseif and(n ~= p, and(n ~= 0, p ~= 0))
          I2 = n*p*((cos((n-p)*pi) - 1)/((n-p)^2) - ...
          (cos((n+p)*pi) - 1)/((n+p)^2))/b;
        end
      end
      
      D(ii,jj) = I1 + I2;
      
    end
  end

  % *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  %% assembly of matrix E
  % *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  
  % preallocation
  E = zeros(nModes, nModes);

  for ii = 1:nModes
    for jj = 1:nModes
      m = mn(ii,1);
      n = mn(ii,2);
      o = mn(jj,1);
      p = mn(jj,2);
      
      % Compute the first part of the integral
      I1 = 0;
      if and(n == p, o ~= 0)
        if m == 0
          I1 = sqrt(2)*cos(o*pi);
        elseif and(m == o, m ~= 0)
          I1 = 0.5;
        elseif and(m ~= o, m ~= 0)
          I1 = o*(cos((m+o)*pi)/(m+o) - cos((m-o)*pi)/(m-o));
        end
      end
      
      % Compute the second part of the integral
      I2 = 0;
      if and(m == o, p ~= 0)
        if n == 0
          I2 = sqrt(2)*cos(p*pi);
        elseif and(n == p, n ~= 0)
          I2 = 0.5;
        elseif and(n ~= p, n ~= 0)
          I2 = p*(cos((n+p)*pi)/(n+p) - cos((n-p)*pi)/(n-p));
        end
      end
      
      E(ii,jj) = I1 + I2;

    end
  end
  
  % *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  %% Wavnumbers
  % *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  
  kn2 = (mn(1:nModes,1)*pi/a).^2 + (mn(1:nModes,2)*pi/b).^2;

end

