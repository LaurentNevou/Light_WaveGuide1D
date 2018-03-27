function[Ex,neff,alpha]=TE_solve(y,eps,lambda,nmodes,neff_min,neff_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k0 = 2*pi/lambda;
dy=y(2)-y(1);
%eps(1) = eps(end) = 1;

%%%%%%%%%%%%%%%%Building of the operators%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
AA=ones(1,length(y));
BB=ones(1,length(y)-1);

DY2=(-2)*diag(AA) + diag(BB,1) + diag(BB,-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building and solving of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = DY2/dy^2 + diag(eps) * k0^2;

H=sparse(H);
%[Ex,Beta] = eigs(H,nmodes,'LR');
[Ex,Beta] = eig(H);

neff=sqrt(diag(Beta))/k0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Filtering and reshaping the Wavefunction %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx1=real(neff)>neff_min;
idx2=real(neff)<neff_max;
idx=logical( idx1.*idx2);

neff=neff(idx);
Ex=Ex(:,idx);

alpha=2*k0*imag(neff);


for i=1:length(neff)
    Ex(:,i)=Ex(:,i)/sqrt(sum(abs(Ex(:,i)).^2)*dy);  % normalisation of the wave function Ex
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%