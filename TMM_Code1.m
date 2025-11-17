%%%%%%%%%%%%%%%%%%%%%%% In The Name of Allah %%%%%%%%%%%%%%%%%%%%%%%%%%
%{
                      Transfer Matrix Method (TMM)
            
                     %% Written by Yaser Khorrami %%
                         % Y.khorrami@gmail.com %
                        %%  version 1.0  %%
						 %%  7 Feb 2022 %%
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Written by Yaser Khorrami %%%%%%%%%%%%%%%%%%%
tic
close all;
clc;
clear all;

%% Step 1: Define the Units %%%
deg = pi/180;
meters = 1;
micrometers = 1e-6 * meters;
centimeters = 1e-2 * meters;
c = 2.99792458e8;                                   %light speed
%%% Parameters %%%
pol = 'E';

i = 1;
for f = [0.1:0.1:2.9]*1e9    %input frequency 
  fr(i) = f;
lam0 = c/f ;
theta = 0 * deg;
phi  = 0 * deg;
ur1 = 1.0 ;           %permeability of reflection layer
er1 = 1.0 ;           %permittivity of reflection layer
ur2 = 1.0 ;           %permeability of transmission layer
er2 = 1.0 ;           %permittivity of transmission layer 
ER = [1.1 , 2.35];    %permittivity of each layer
UR = [1.0 , 1.0];     %permeability of each layer
L = [1.999 , 6.5187]*centimeters;
ninc = sqrt(ur1 * er1) ;
ni = [sqrt(ER(1) * UR(1)) sqrt(ER(2) * UR(2))] ;
lambda = lam0/ninc ;
k0 = 2*pi/lam0 ;
switch pol
    case 'E' 
        pte  = 1;                                        %1:TE and 0:TM
        ptm  = 0;                                        %0:TE and 1:TM
    case 'H'
        pte  = 0;                                        %1:TE and 0:TM
        ptm  = 1;                                        %0:TE and 1:TM
end
%% Step 2: Wavevectors
kinc = 1 * ninc * [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)] ; % kx,ky and kz are normalized by k0
kx = kinc(1);
ky = kinc(2);
urh=1;
erh= 1 + (kx^2) + (ky^2);
Qh = (1/urh)*[kx*ky (urh*erh)-(kx^2);(ky^2)-(urh*erh) -kx*ky];
Vh=-1j*Qh;
Wh=eye(2);
%% step 3 : Main loop through layers
SG.S11=zeros(2);
SG.S12=eye(2);
SG.S21=eye(2);
SG.S22=zeros(2);

nlay = length(L);  % Number of layers
SG.S11=zeros(2);
SG.S12=eye(2);
SG.S21=eye(2);
SG.S22=zeros(2);

%% step 4 : Main Loop Iterates Through Layers
for nn=1:nlay
P=zeros(2);
Q=zeros(2);
V(:,:,nn)=zeros(2);
W(:,:,nn)=zeros(2);
LAM2(:,:,nn)=zeros(2);
OMEGA2(:,:,nn)=zeros(2);
A=zeros(2,2,nn);
B=zeros(2,2,nn);
D=zeros(2,2);
F=zeros(2,2);
end

for  ii = 1 : 2
    Q(:,:,ii) = (1/UR(ii))*[kx*ky (ER(ii)*UR(ii))-(kx^2);(ky^2)-(ER(ii)*UR(ii)) -kx*ky];
    P(:,:,ii) = (1/ER(ii))*[kx*ky (ER(ii)*UR(ii))-(kx^2);(ky^2)-(ER(ii)*UR(ii)) -kx*ky];
    kz(ii)    = sqrt((ER(ii)*UR(ii)) - kx^2 - ky^2);
    OMEGA(:,:,ii) = sqrt(P(:,:,ii)*Q(:,:,ii));
    [W(:,:,ii),LAM2(:,:,ii)] = eig(OMEGA(:,:,ii));
    V(:,:,ii) = (Q(:,:,ii)*W(:,:,ii))*inv(LAM2(:,:,ii));
    X(:,:,ii) = expm(LAM2(:,:,ii) * k0 * L(ii));
    A(:,:,ii) = inv(W(:,:,ii))*Wh + inv(V(:,:,ii))*Vh;
    B(:,:,ii) = inv(W(:,:,ii))*Wh - inv(V(:,:,ii))*Vh;
    
    S.S11(:,:,ii) = inv(A(:,:,ii) - X(:,:,ii)*B(:,:,ii)*inv(A(:,:,ii))*X(:,:,ii)*B(:,:,ii))*(X(:,:,ii)*B(:,:,ii)*inv(A(:,:,ii))*X(:,:,ii)*A(:,:,ii) - B(:,:,ii));
    S.S12(:,:,ii) = inv(A(:,:,ii) - X(:,:,ii)*B(:,:,ii)*inv(A(:,:,ii))*X(:,:,ii)*B(:,:,ii))*X(:,:,ii)*(A(:,:,ii) - B(:,:,ii)*inv(A(:,:,ii))*B(:,:,ii));
    S.S21(:,:,ii) = S.S12(:,:,ii);
    S.S22(:,:,ii) = S.S11(:,:,ii);
    %d. Update device scattering matrix
    D = SG.S12*inv(eye(2)-(S.S11(:,:,ii)*SG.S22));
    F = S.S21(:,:,ii)*inv(eye(2)-SG.S22*S.S11(:,:,ii));

    SG.S11 = SG.S11 + D*S.S11(:,:,ii)*SG.S21;
    SG.S12 = D*S.S12(:,:,ii);
    SG.S21 = F*SG.S21;
    SG.S22 = S.S22(:,:,ii) + F*SG.S22*S.S12(:,:,ii);
end
%% step 5 : Transmission and Reflection side
 kzref =-conj(sqrt((ur1*er1) - kx.^2 - ky.^2)); 
 kztrn = conj(sqrt((ur2*er2) - kx.^2 - ky.^2));
% Compute Reflection Side Connection S-Matrix
  %only one layer 
 Qref=(1/ur1)*[kx*ky (ur1*er1-kx.^2);(ky.^2-ur1*er1) -ky*kx];
 Wref = eye(2);
 lamref=[-1j*kzref 0;0 -1j*kzref]; 
 Vref = Qref * inv(lamref);
  
 Aref = inv(Wh)*Wref + inv(Vh)*Vref;
 Bref = inv(Wh)*Wref - inv(Vh)*Vref;

 Sref.S11 = -inv(Aref)*Bref;
 Sref.S12 = 2*inv(Aref);
 Sref.S21 = 0.5*(Aref - Bref*inv(Aref)*Bref);
 Sref.S22 = Bref*inv(Aref); 
% Compute Transmission Side Connection S-Matrix  
  %only one layer 
 Qtrn = (1/ur2)*[kx*ky (er2*ur2-kx.^2);(ky.^2-er2*ur2) -ky*kx];
 Wtrn =eye(2);
 lamtrn = [1j*kztrn 0;0 1j*kztrn];
 Vtrn =Qtrn*inv(lamtrn);
  
 Atrn = inv(Wh)*Wtrn + inv(Vh)*Vtrn;
 Btrn = inv(Wh)*Wtrn - inv(Vh)*Vtrn;

 Strn.S11 = Btrn*inv(Atrn);
 Strn.S12 = 0.5 * (Atrn - Btrn*inv(Atrn)*Btrn);
 Strn.S21 = 2*inv(Atrn);
 Strn.S22 = -inv(Atrn)*Btrn; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 6: Compute Global Scattering Matrix
  
 for nn = 1 : 2
    if nn==1
        SA.S11=Sref.S11;    SB.S11=SG.S11;
        SA.S12=Sref.S12;    SB.S12=SG.S12;
        SA.S21=Sref.S21;    SB.S21=SG.S21;
        SA.S22=Sref.S22;    SB.S22=SG.S22;
    elseif nn==2
        SA.S11=SG.S11;      SB.S11=Strn.S11;
        SA.S12=SG.S12;      SB.S12=Strn.S12;
        SA.S21=SG.S21;      SB.S21=Strn.S21;
        SA.S22=SG.S22;      SB.S22=Strn.S22;
    end
        
 D = SA.S12 * inv(eye(2) - SB.S11*SA.S22);
 F = SB.S21 * inv(eye(2) - SA.S22*SB.S11);

 SG.S11 = SA.S11 + D*SB.S11*SA.S21;
 SG.S12 = D*SB.S12;
 SG.S21 = F*SA.S21;
 SG.S22 = SB.S22 + F*SA.S22*SB.S12;
 end
%% Step 7: Compute Reflected and Transmitted Fields
 n = [0; 0; 1];    %Polarization unit vector (only along z)
 if theta == 0
   ate = [0 1 0]; % unit vector of TE Polarization 
 else
   ate = -1 * cross(n,kinc) ./ norm(cross(n,kinc));
 end
 atm = cross(ate,kinc) / norm(cross(ate,kinc));  % unit vector of TM Polarization
 P = -((pte*ate) + (ptm*atm));
 %Compute Source Field
  Esrcx = P(1);
  Esrcy = P(2);
 %Compute Delta vector
  esrc = [Esrcx ; Esrcy ];  %  Source Field
 %   REFLECTED FIELDS
  eref = SG.S11 * esrc;
  Exref = eref(1,1);
  Eyref = eref(2,1);
 %   TRANSMITTED FIELDS
  etrn = SG.S21 * esrc;
  Extrn = etrn(1,1);
  Eytrn = etrn(2,1);
%%%Calculate Longitudinal Field Components
 Ezref = +(kx*Exref + ky*Eyref)/kzref;
 Eztrn = -(kx*Extrn + ky*Eytrn)/kztrn;
%%%%%
Eref2 = abs(Exref)^2 + abs(Eyref)^2 + abs(Ezref)^2 ;
Etrn2 = abs(Extrn)^2 + abs(Eytrn)^2 + abs(Eztrn)^2 ;

R(i) = Eref2;
T(i) = -(Etrn2) * real((ur1 * kztrn)/(ur2 * kzref));

i = i+1;

end
% CON = R + T; 
% 
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% if (CON > 1)
%     disp('materials have gain')
% elseif (CON ~= 1)
%     disp('materials have no loss and no gain')
% elseif (CON < 1)
%     disp('materials have loss')
% end
%% Step 8: Plot
figure(1)
plot(fr/1e9,R*100,fr/1e9,T*100,'linewidth',1.6)
xlabel('Frequency (GHz)'); ylabel('%');
title('Transfer Matrix Method');
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')

toc





