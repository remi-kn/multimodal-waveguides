% implement admittance computation with Mobius/Magnus scheme from 
%
% Multimodal admittance method in waveguides and singularity
% behavior at high frequencies
%
% Pagneux 2010

clc

% parameters to reproduce Fig. 4a
nModes = 60;
k = 40.3*pi;

% generate space discretisation (with decreasing x)
nSections = 334;
dx = 2.4/nSections;
x = (nSections:-1:1)*dx -1.2;

% set radiation boundary condition on the right of the guide Yc = jkn
Yc = 1i*diag(sqrt(k^2 - (((0:nModes-1).')*pi).^2));

% initialise admittance matrix
Y = zeros(nSections, nModes, nModes);
Y(1,:,:) = Yc;
E1 = zeros(nSections, nModes, nModes);
E2 = zeros(nSections, nModes, nModes);

% loop over the sections
for ii = 2:nSections
  
  % x_n
  xn = x((ii-1));
  % x_n+1
  xnp1 = x(ii);
  % compute spacing
  dx = xnp1 - xn;
  
  H1 = computeH(xn + dx*(3-sqrt(3))/6, k, nModes);
  H2 = computeH(xn + dx*(3+sqrt(3))/6, k, nModes);

  % compute exp( omega_n )
  On = expm(0.5*dx*(H1 + H2) + sqrt(3)*(dx^2)*(H2*H1 - H1*H2 )/12);

  
  E1(ii,:,:) = On(1:nModes, 1:nModes);
  E2(ii,:,:) = On(1:nModes, nModes+1:end);
  E3 = On(nModes+1:end, 1:nModes);
  E4 = On(nModes+1:end, nModes+1:end);

  % compute the admittance matrix of the next step (Eq. 15)
  Y(ii,:,:) = (E3 + E4*squeeze(Y((ii-1),:,:))) ...
  /(squeeze(E1(ii,:,:)) + squeeze(E2(ii,:,:))*squeeze(Y((ii-1),:,:)));
end

%% compute mode amplitude

% initialise the modes amplitude
a = zeros(nModes, nSections);
% set the mode amplitude vector at the beginning of the geometry 
% (end of the vector)
a(10,end) = 1;

% loop over the sections, but in reversed order 
for ii = 1:nSections -1
  % x_n
  xn = x(end - (ii-1));
  % x_n+1
  xnp1 = x(end - ii);
  % compute spacing
  dx = xnp1 - xn;
  
  H2 = computeH(xn + dx*(3-sqrt(3))/6, k, nModes);
  H1 = computeH(xn + dx*(3+sqrt(3))/6, k, nModes);
  dx = -dx;
  % compute exp( omega_n )
  On = expm(0.5*dx*(H1 + H2) + sqrt(3)*(dx^2)*(H2*H1 - H1*H2 )/12);

  
  E1 = On(1:nModes, 1:nModes);
  E2 = On(1:nModes, nModes+1:end);
  
  % compute the amplitude of the previous step (Eq. 16)
  a(:,end - ii) = (E1 + E2*squeeze(Y(end - ii,:,:))) \ a(:,end - (ii-1));
end

%% compute acoustic field

% discretisation over y axis
ny = 200;
y = linspace(0, 1.3, ny);

% initialise field
field = zeros(ny, nSections);

% loop over sections
for ii = 1:nSections
    % compute h
    if and(x(ii) <1, x(ii)>-1)
     % compute height 
     h = 1 + 0.15*(1 + cos(pi*x(ii)));
    else
     % compute height 
     h = 1;
    end
    
    % loop over modes
    for jj = 1:nModes
        field(:,ii) = field(:,ii) + a(jj,ii)*sqrt(2/h)*sin((jj-1)*pi*(y')/h);
    end
    
    % replace the filed outside of the guide by nans
    field(y>h,ii) = nan;
end

% plot field
figure,
imagesc(x,y,log10(abs(field)))
caxis([-1, 0.5])
axis xy
