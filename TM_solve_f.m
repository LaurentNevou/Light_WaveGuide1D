function[Ey,neff,alpha]=TM_solve_f(y,eps,lambda,nmodes,neff_min,neff_max)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Building of the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AA=ones(1,length(y));
BB=ones(1,length(y)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DY1 = (-0.5)*diag(BB,-1) + (0.5)*diag(BB,1); %symetric derivative

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DY2 = (-2)*diag(AA) + (1)*diag(BB,1) + (1)*diag(BB,-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building and solving of the Hamiltonien %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H0 = DY2/dy^2 + diag(eps(:)) * k0^2;

% Those 2 writtings give the same result!!!
H = H0 + diag(eps)    *  diag(  DY1*(1./(eps(:))) ) * DY1  /dy^2 ;
%H = H0 - diag(1./eps) *  diag(  DY1*eps(:) )        * DY1  /dy^2 ;

% Like in 2D, it gives different results with the folowing writting
%H = H0 -   diag(  DY1*log(eps(:)) )        * DY1  /dy^2 ;


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
