function[Ey,neff,alpha]=TM_solve_f2(y,eps,lambda,nmodes,neff_min,neff_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% One problem here is the operator (d/z)*(1/eps(z))*(d/z) inside the Wave equation
% Without any modification, the results are not accurate and show discontinuity 
% at the material interface due to the difference of epsilon
% One very nice approach is to use the mid-point mass:

% Paul Harrisson
% Quantum Wells, Wires and Dots.
% 4th edition (2016),
% chap 3 : "Numerical Solutions"
% 3.13: "Extention to variable effective mass"
% page 102

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c    = 2.99792458e8;      %% m/s
eps0 = 8.85418782e-12;    %% F/m
mu0  = 4*pi*1e-7;         %% H/m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 2*pi/lambda;
w  = c*k0;
%eps(1) = eps(end) = 1;
dy=y(2)-y(1);
Ny=length(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Building of the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsp = [         (eps(1:end-1) + eps(2:end)) / 2   eps(end) ];
epsm = [ eps(1) (eps(1:end-1) + eps(2:end)) / 2              ];

b = (1./epsp + 1./epsm) .* ones(1,Ny) ;
a = 1./epsm(2:end) .* ones(1,Ny-1) ;
c = 1./epsm(2:end) .* ones(1,Ny-1) ;

DY2 = (-1)*diag(b)  +  (1)*diag(a,-1)  +  (1)*diag(c,+1) ;

DY2 = DY2 / dy^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the equation is correct (which I am not sure), 
% the solver is for sure vey well accurate

H = diag(eps(:)) * DY2 + diag(eps(:)) * k0^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building and solving of the Hamiltonien %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%H=sparse(H);
%[Hx,Beta] = eigs(H,nmodes,'LR');
[Hx,Beta] = eig(H);
Beta = diag(Beta);

neff=sqrt(Beta)/k0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Filtering and reshaping the Wavefunction %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx1=real(neff)>neff_min;
idx2=real(neff)<neff_max;
idx=logical( idx1.*idx2);

neff=neff(idx);
Hx=Hx(:,idx);

alpha=2*k0*imag(neff);


% Ey = -Beta/(w*eps0*eps) * Hx

for i=1:length(neff)
    Ey(:,i)=Beta(i)*Hx(:,i)./(w*eps0*eps(:));       %transformation of H in electric field
    Ey(:,i)=Ey(:,i)/sqrt(sum(abs(Ey(:,i)).^2)*dy);  % normalisation of the wave function Ey
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here is a small patch due to differences between Octave and Matlab
% Matlab order the eigen values while Octave reverse it

if length(neff)>1
if neff(2)>neff(1)
  Ey=Ey(:,end:-1:1);
  neff=neff(end:-1:1);
end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
