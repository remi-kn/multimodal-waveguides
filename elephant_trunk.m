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
a = 0.05;
b = 0.03;

% number of propagation modes 
nModes = 20;

% admittance at the boundary
Y = 0.;

% sound speed
c = 340;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Build the geometrical transformation matrices from analytical expression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ca, Da, Ea, kn2] = buildMatricesCDE_analyticalRectangle (a/b, 1, nModes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve multimodal problem with the elephant's trunk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r0 = 1.25*2*b;
thetaL = 2.62;
L = thetaL*r0;%3.27*sqrt(b*a);
nX = 800;
X = linspace(-1*b, L+b, nX);
[l, dl, kappa] = trunkElephant(X, L, b);

Ad = zeros(nX, nModes, nModes);
k = 3/1/b;
K = diag(1i*sqrt((k*l(end))^2 - kn2));
Ad(end,:,:) = K*l(end)^2;

%% loop over the section to compute the admittance
for ii = 1:nX-1
  xn = X(end - (ii-1));
  xnp1 = X(end - (ii));
  dx = xnp1 - xn;
  
  xA1 = xn + (0.5 - sqrt(3)/6)*dx;
  [l, dl, kappa] = trunkElephant(xA1, L, b);
  K = diag(1i*sqrt((k*l)^2 - kn2));
 A1 = [dl*Ea/l , (eye(nModes) - kappa*l*Ca)/(l^2);
 (K.^2 + kappa*l*(Ca*(k*l)^2 - Da)) , -dl*transpose(Ea)/l];
  
  xA2 = xn + (0.5 + sqrt(3)/6)*dx;
  [l, dl, kappa] = trunkElephant(xA2, L, b);
  K = diag(1i*sqrt((k*l)^2 - kn2));
 A2 = [dl*Ea/l , (eye(nModes) - kappa*l*Ca)/(l^2);
 (K.^2 + kappa*l*(Ca*(k*l)^2 - Da)) , -dl*transpose(Ea)/l];  
  
  On = expm(0.5*dx*(A1 + A2) + sqrt(3)*(dx^2)*(A2*A1 - A1*A2 )/12);
  
  E1 = On(1:nModes, 1:nModes);
  E2 = On(1:nModes, nModes+1:end);
  E3 = On(nModes+1:end, 1:nModes);
  E4 = On(nModes+1:end, nModes+1:end);
  
  Ad(nX-ii, :, :) = (E3 + E4*squeeze(Ad(nX-(ii-1),:,:)))/...
  (E1 + E2*squeeze(Ad(nX-(ii-1),:,:)));

end

p = zeros(nModes, nX);
p(1,1) = 1;

%% loop over the section to compute the pressure amplitude
for ii = 1:nX-1
  xn = X(ii);
  xnp1 = X(ii+1);
  dx = xnp1 - xn;
  
    xA1 = xn + (0.5 - sqrt(3)/6)*dx;
  [l, dl, kappa] = trunkElephant(xA1, L, b);
  K = diag(1i*sqrt((k*l)^2 - kn2));
 A1 = [dl*Ea/l , (eye(nModes) - kappa*l*Ca)/(l^2);
 (K.^2 + kappa*l*(Ca*(k*l)^2 - Da)) , -dl*transpose(Ea)/l];  
  
  xA2 = xn + (0.5 + sqrt(3)/6)*dx;
  [l, dl, kappa] = trunkElephant(xA2, L, b);
  K = diag(1i*sqrt((k*l)^2 - kn2));
 A2 = [dl*Ea/l , (eye(nModes) - kappa*l*Ca)/(l^2);
 (K.^2 + kappa*l*(Ca*(k*l)^2 - Da)) , -dl*transpose(Ea)/l];  
  dx = -dx;
  
  On = expm(0.5*dx*(A1 + A2) + sqrt(3)*(dx^2)*(A2*A1 - A1*A2 )/12);
  
  E1 = On(1:nModes, 1:nModes);
  E2 = On(1:nModes, nModes+1:end);
  
  p(:,ii+1) = (E1 + E2*squeeze(Ad(ii+1,:,:)))\p(:,ii);
end

%% plot acoustic field

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

mn = mn - 1;
nx = 800;
x = linspace(-b/4, 6*b, nx);
ny = 800;
y = linspace(-b, 3.5*b, ny);
Xp = zeros(nx,ny);
Yp = zeros(nx,ny);
Pp = zeros(nx,ny);


for ii = 1:nx
  for jj = 1:ny
    if and(y(jj)<0, abs(x(ii)) <= 0.125)
      Xp(ii,jj) = y(jj);
      Yp(ii,jj) = 2*x(ii)/b+0.5;
    else
      d = sqrt((x(ii) - r0)^2 + y(jj)^2);
      Xp(ii,jj) = acos((r0 - x(ii))/d)*r0;
      [l, dl, kappa] = trunkElephant(Xp(ii,jj), L, b);
      Yp(ii,jj) = (r0-d)/l+0.5;
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
    
    if and(and(Xp(ii,jj) <= L+b, Xp(ii,jj) > 0), and(Yp(ii,jj)<1, Yp(ii,jj)>0))
      [toto, idx] = min(abs(X - Xp(ii,jj)));
      for mm = 1:nModes
        if mn(mm,1) == 0, Mm = 1; else, Mm = sqrt(2); end
        if mn(mm,2) == 0, Nn = 1; else, Nn = sqrt(2); end
        Pp(ii,jj) = Pp(ii,jj) + p(mm,idx)*(Mm*cos(mn(mm,1)*pi*0.1) + ...
        Nn*cos(mn(mm,2)*pi*Yp(ii,jj)));
      end
    elseif and(and(Xp(ii,jj) <0, Xp(ii,jj) > -b), and(Yp(ii,jj)<1, Yp(ii,jj)>0))
      [idx] = find(X>Xp(ii,jj),1,"first");
      for mm = 1:nModes
        if mn(mm,1) == 0, Mm = 1; else, Mm = sqrt(2); end
        if mn(mm,2) == 0, Nn = 1; else, Nn = sqrt(2); end
        Pp(ii,jj) = Pp(ii,jj) + p(mm,idx)*(Mm*cos(mn(mm,1)*pi*0.1) + ...
        Nn*cos(mn(mm,2)*pi*Yp(ii,jj)));
      end
    end
  end
end
Pp(Pp == 0) = nan;

h = figure;
imagesc(x,y,fliplr(real(Pp)'))
colormap("jet")
axis xy
axis equal

print(h, "-dpng", "colormap-elephant-analytical.png");

theta = linspace(0, thetaL, 100);
R1 = 3*b*((theta./thetaL).^3)/2 - 9*b*((theta./thetaL).^2)/4 + r0 -b/4;
translat = [sqrt(r0^2 + b^2)*cos(pi-thetaL-atan(b/r0)) - r0*cos(pi-thetaL),...
sqrt(r0^2 + b^2)*sin(pi-thetaL-atan(b/r0)) - r0*sin(pi-thetaL)];
R2 = -3*b*((theta./thetaL).^3)/2 + 9*b*((theta./thetaL).^2)/4 + r0 +b/4;

h = figure; contour(x,y,real(Pp.'),100), axis equal
colormap("jet")
hold on
plot([b/4, b/4, r0-R1.*cos(theta), r0-R1(end).*cos(theta(end)) + translat(1)],...
 [-b, 0, R1.*sin(theta), R1(end).*sin(theta(end)) + translat(2)], "k")
plot([-b/4, -b/4, r0-R2.*cos(theta), r0-R2(end).*cos(theta(end)) + translat(1)],...
 [-b, 0, R2.*sin(theta), R2(end).*sin(theta(end)) + translat(2)], "k")
xlim([-b/2 6.5*b])
print(h, "-dpng", "contour-elephant-analytical.png");



function [l, dl, kappa] = trunkElephant(X, L, b)
  nX = length(X);
  
  l = zeros(1,nX);
  dl = zeros(1,nX);
  kappa = zeros(1,nX);
  
  l(X<0) = b/2;
  dl(X<0) = 0;
  kappa(X<0) = 0;
  
  l(and(X>0, X<L)) = 0.5*b*(1 + 9*(X(and(X>0, X<L))/L).^2 - ...
  6*(X(and(X>0, X<L))/L).^3);
  dl(and(X>0, X<L)) = 9*b*X(and(X>0, X<L)).*(1-X(and(X>0, X<L))/L)/(L^2);
  kappa(and(X>0, X<L)) = 1/2/1.25/b;
  
  l(X>L) = 2*b;
  dl(X>L) = 0;
  kappa(X>L) = 0;
end



function [Ca, Da, Ea, kn2] = buildMatricesCDE_analyticalRectangle (a, b, nModes)
  
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

  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  %% assembly of matrix C
  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 

  %% preallocation 
  Ca = zeros(nModes, nModes);

  for ii=1:nModes
    for jj = 1:nModes
      m = mn(ii,1);
      n = mn(ii,2);
      o = mn(jj,1);
      p = mn(jj,2);
      
      if m == o
        if n == p
          Ca(ii,jj) = b/2;
        elseif and(n ~= 0, p ~= 0)
          Ca(ii,jj) = b*((cos((n+p)*pi)-1)/((n+p)^2) + ...
          (cos((n-p)*pi)-1)/((n-p)^2))/(pi^2);
        elseif and(n == 0, p ~= 0)
          Ca(ii,jj) = sqrt(2)*b*(cos(p*pi)-1)/((p*pi)^2);
        else
          Ca(ii,jj) = sqrt(2)*b*(cos(n*pi)-1)/((n*pi)^2);
        end
      end
      
    end
  end



  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  %% assembly of matrix D
  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 

  %% preallocation 
  Da = zeros(nModes, nModes);

  for ii = 1:nModes
    for jj = 1:nModes
      m = mn(ii,1);
      n = mn(ii,2);
      o = mn(jj,1);
      p = mn(jj,2);
      
  %%    if and(ii == 6, jj == 5), keyboard, end
      
      %% Compute the first part of the integral
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
      
      %% Compute the second part of the integral
      I2 = 0;
      if m == o
        if and(n == p, and(n ~= 0, p ~= 0))
          I2 = ((n*pi)^2)/2/b;
        elseif and(n ~= p, and(n ~= 0, p ~= 0))
          I2 = n*p*((cos((n-p)*pi) - 1)/((n-p)^2) - ...
          (cos((n+p)*pi) - 1)/((n+p)^2))/b;
        end
      end
      
      Da(ii,jj) = I1 + I2;
      
    end
  end

  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  %% assembly of matrix E
  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 

  Ea = zeros(nModes, nModes);

  for ii = 1:nModes
    for jj = 1:nModes
      m = mn(ii,1);
      n = mn(ii,2);
      o = mn(jj,1);
      p = mn(jj,2);
      
      %% Compute the first part of the integral
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
      
      %% Compute the second part of the integral
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
      
      Ea(ii,jj) = I1 + I2;

    end
  end
  
  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  %% Wavnumbers
  %% *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   * 
  
  kn2 = (mn(1:nModes,1)*pi/a).^2 + (mn(1:nModes,2)*pi/b).^2;

end

