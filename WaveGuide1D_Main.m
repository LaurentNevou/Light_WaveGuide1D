%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 27March2018, lne %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Solving Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy=5e-8;              %% resolution in y
lambda=6e-6;          %% Wavelength (meter)

neff_min=2;           %% filter the solutions where the effective index is superior than
neff_max=4;           %% filter the solutions where the effective index is inferior than
nmodes=4;             %% number of solutions asked (not activated)

TE=1;                 %% activate the calcul for the TE mode
TM=1;                 %% activate the calcul for the TM mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Optical index definition %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% real index / imag index / thickness µm / doping cm-3 / eff mass / Drude time ps

M=[

2    0  3    0  0.067  0.1
3.6  0  1.5  0  0.067  0.1
2    0  3    0  0.067  0.1

];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Computation of the imaginary part from the doping %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=2.99792458E8;                 %% speed of light (m/s)
e=1.602176487E-19;              %% electron charge (C)
eps0=8.854187817620E-12;        %% vaccum dielectric constant (F/m)
mu0=1/(eps0*c^2);               %% vaccum permeabiliy
me=9.10938188E-31;              %% electron mass (kg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=2*pi*c/lambda;
M(:,4)=M(:,4)*1E6;          %% conversion of the doping from cm-3 to m-3
M(:,5)=M(:,5)*me;
M(:,6)=M(:,6)*1E-12;        %% conversion of tau from ps to second

for i=1:length(M(:,1))
    
    Plasma_2=M(i,4)*e^2 / ( eps0 * M(i,5) );  %%% Plasma frequency
    Ki= - Plasma_2/(w^2+ 1i*w/M(i,6));        %%% suceptibility
    
    if isnan(Ki)==1
        Ki=0;
    end
    
    epst(i)=M(i,1)^2 + Ki ;
    
    nn=real( sqrt(epst(i)) );
    kk=imag( sqrt(epst(i)) );
    
    M(i,1)=nn;
    if kk~=0
        M(i,2)=kk;
    end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Discretisation of Epsilon %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yt=M(:,3)*1E-6;  % conversion of the length from um to meter

y=[];eps=[];
y(1)=0;eps(1)=epst(1);

for i=1:length(yt)
    t=yt(i);
    yv= (y(end)+dy): dy : (y(end)+dy)+t;
    y=[y yv];
    epsv=ones(size(yv))*epst(i);
    eps=[eps epsv];
end    

n=sqrt(eps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TE==1
    [Ex,neffTE,alphaTE]=TE_solve(y,eps,lambda,nmodes,neff_min,neff_max);
    n1=length(neffTE);
    neffTE
    alphaTE=alphaTE*0.01   %% losses in cm-1
    ScF=0.5; 
    for i=1:n1
         Ex(:,i)=abs(Ex(:,i)).^2/max(abs(Ex(:,i)).^2)*ScF + real(neffTE(i)); % normalisation for the plotting
    end
end

if TM==1
    [Ey,neffTM,alphaTM]=TM_solve(y,eps,lambda,nmodes,neff_min,neff_max);
    n2=length(neffTM);
    neffTM
    alphaTM=alphaTM*0.01   %% losses in cm-1
    ScF=0.5; 
    for i=1:n2
         Ey(:,i)=abs(Ey(:,i)).^2/max(abs(Ey(:,i)).^2)*ScF + real(neffTM(i)); % normalisation for the plotting
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0fig=50; Y0fig=100;
Wfig=1200;Hfig=1000;

figure('Name','Electrical Field','position',[X0fig Y0fig Wfig Hfig])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(111,'fontsize',20)
hold on;grid on;

%plot(y*1e6,abs(n),'b-','linewidth',2)      % plot the refractive index
plot(y*1e6,real(n),'b-','linewidth',2)


if TE==1
    for i=1:n1
         plot(y*1e6,Ex(:,i),'g','linewidth',2)  % plot Ex for TE-modes
    end
end

if TM==1
    for i=1:n2
         plot(y*1e6,Ey(:,i),'r','linewidth',2)  % plot Ey for TM-modes
    end
end


s={
'\fontsize{20}\color{blue}Optical index'
'\fontsize{20}\color{green}TE-modes: Ex(y)'
'\fontsize{20}\color{red}TM-modes: Ey(y)'
};

xlabel('y (um)');
ylabel('effective index');

x0=0.8*get(gca,'xlim')(end);
y0=0.9*(get(gca,'ylim')(end) - get(gca,'ylim')(1)) + get(gca,'ylim')(1);

text(x0,y0,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%