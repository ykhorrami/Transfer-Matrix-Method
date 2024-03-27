
close all;
clc;
clear all;

tic
%
%% Step 1: Define the Units %%%
deg = pi/180;
meters = 1;
nanometers = 1e-9 * meters;
micrometers = 1e-6 * meters;
centimeters = 1e-2 * meters;
c = 2.99792458e8;                                   %light speed
%% Parameters %%%
pol = 'E';
M =  load("Weaver.txt");
i = 1;
j = 1;
freq = M(:,4).*1e12; 
eps_re = M(:,5);
eps_im = M(:,6);
h1 = (eps_re-1i.*eps_im);
for p = 1:numel(freq)
 for f = [50:1:200]*10^(12)  %input frequency 
     fr(i) = f;
    lam0 = c/f ;  
    theta = 0 * deg;
    phi  = 0 * deg;
    ur1 = 1.0 ;           %permeability of reflection layer
    er1 = 1.0 ;           %permittivity of reflection layer
    ur2 = 1.0 ;           %permeability of transmission layer
    er2 = 10.8+1e-6*1i ;   %permittivity of transmission layer
    perm_W = h1(p);
    ER = [1.0 , 17.0 , perm_W , 17.0 ,perm_W , 17.0 , perm_W , 10.8+1e-6*1i];    %permittivity of each layer
    UR = [1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0];     %permeability of each layer
    L = [1000 , 310 , 20 , 310 , 20 , 310 , 20, 5*1e+07]*nanometers;
    ninc = sqrt(ur1 * er1) ;
    ni = [sqrt(ER(1) * UR(1)) sqrt(ER(2) * UR(2)) sqrt(ER(3) * UR(3)) sqrt(ER(4) * UR(4)) sqrt(ER(5) * UR(5)) sqrt(ER(6) * UR(6)) sqrt(ER(7) * UR(7)) sqrt(ER(8) * UR(8))] ;
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
    

%%Step 2: Wavevectors
    kinc = 1 * ninc * [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)] ; % kx,ky and kz are normalized by k0

    kx = kinc(1);
    ky = kinc(2);
    urh=1;
    erh= 1 + (kx^2) + (ky^2);
    Qh = (1/urh)*[kx*ky (urh*erh)-(kx^2);(ky^2)-(urh*erh) -kx*ky];
    Vh=-1j*Qh;
    Wh=eye(2);

%%step 3 : Main loop through layers
    SG.S11=zeros(2);
    SG.S12=eye(2);
    SG.S21=eye(2);
    SG.S22=zeros(2);
    
    nlay = numel(L);  % Number of layers
    SG.S11=zeros(2);
    SG.S12=eye(2);
    SG.S21=eye(2);
    SG.S22=zeros(2);

%%step 4 : Main Loop Iterates Through Layers
    for nn=1:nlay
        P=zeros(2);
        Q=zeros(2);
        V(:,:,nn)=zeros(2);
        W(:,:,nn)=zeros(2);
        LAM2(:,:,nn)=zeros(2);
        OMEGA2(:,:,nn)=zeros(2);
        A=zeros(2,2,nn);
        B=zeros(2,2,nn);
        D=zeros(2,2,nn);
        F=zeros(2,2,nn);
    end
   
    
    for  ii = 1 : nlay
        Q(:,:,ii) = (1/UR(ii))*[kx*ky (ER(ii)*UR(ii))-(kx^2);(ky^2)-(ER(ii)*UR(ii)) -kx*ky];
        P(:,:,ii) = (1/ER(ii))*[kx*ky (ER(ii).*UR(ii))-(kx^2);(ky^2)-(ER(ii).*UR(ii)) -kx*ky];
        kz(ii)    = sqrt((UR(ii).*ER(ii)) - kx^2 - ky^2);
        OMEGA(:,:,ii) = sqrt(P(:,:,ii)*Q(:,:,ii));
        [W(:,:,ii),LAM2(:,:,ii)] = eig(OMEGA(:,:,ii));
        V(:,:,ii) = (Q(:,:,ii)*W(:,:,ii))./(LAM2(:,:,ii));
        X(:,:,ii) = expm(LAM2(:,:,ii) * k0 * L(ii));
        A(:,:,ii) = Wh./(W(:,:,ii)) + Vh./(V(:,:,ii));
        B(:,:,ii) = Wh./(W(:,:,ii)) - Vh./(V(:,:,ii));
        
        S.S11(:,:,ii) = B(:,:,ii)./(A(:,:,ii) - X(:,:,ii)./(A(:,:,ii))*X(:,:,ii)*B(:,:,ii))*(X(:,:,ii)*B(:,:,ii)./(A(:,:,ii))*X(:,:,ii)*A(:,:,ii) - B(:,:,ii));
        S.S12(:,:,ii) = 1./(A(:,:,ii) - X(:,:,ii)*B(:,:,ii)./(A(:,:,ii))*X(:,:,ii)*B(:,:,ii))*X(:,:,ii)*(A(:,:,ii) - B(:,:,ii)./(A(:,:,ii))*B(:,:,ii));
        S.S21(:,:,ii) = S.S12(:,:,ii);
        S.S22(:,:,ii) = S.S11(:,:,ii);
        %d. Update device scattering matrix
        D = SG.S12./(eye(2)-(S.S11(:,:,ii).*SG.S22));
        F = S.S21(:,:,ii)./(eye(2)-SG.S22.*S.S11(:,:,ii));
        
        SG.S11 = SG.S11 + D*S.S11(:,:,ii)*SG.S21;
        SG.S12 = D*S.S12(:,:,ii);
        SG.S21 = F*SG.S21;
        SG.S22 = S.S22(:,:,ii) + F*SG.S22*S.S12(:,:,ii);
    end

%%step 5 : Transmission and Reflection side
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

%%Step 6: Compute Global Scattering Matrix
    for nn = 1 : 2
        if nn==2*nn -1
            SA.S11=Sref.S11;    SB.S11=SG.S11;
            SA.S12=Sref.S12;    SB.S12=SG.S12;
            SA.S21=Sref.S21;    SB.S21=SG.S21;
            SA.S22=Sref.S22;    SB.S22=SG.S22;
        elseif nn==2*nn
            SA.S11=SG.S11;      SB.S11=Strn.S11;
            SA.S12=SG.S12;      SB.S12=Strn.S12;
            SA.S21=SG.S21;      SB.S21=Strn.S21;
            SA.S22=SG.S22;      SB.S22=Strn.S22;
        end
        
        D = SA.S12 ./(eye(2) - SB.S11*SA.S22);
        F = SB.S21 ./(eye(2) - SA.S22*SB.S11);
        
        SG.S11 = SA.S11 + D*SB.S11*SA.S21;
        SG.S12 = D*SB.S12;
        SG.S21 = F*SA.S21;
        SG.S22 = SB.S22 + F*SA.S22*SB.S12;
    end

%%Step 7: Compute Reflected and Transmitted Fields
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
    
    R(i,j) = Eref2;
    T(i,j) = -(Etrn2) * real((ur1 * kztrn)/(ur2 * kzref));
    
    i = i+1;    
end
j = j+1;

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
end
%%Step 8: Plot
figure(1)
plot(fr/1e12,T*100,'linewidth',1.6)
xlabel('Frequency (THz)'); ylabel('%');
title('Transfer Matrix Method');
set(gca,'linewidth', 1.5,'fontsize',12,'fontname','Times New Roman')

% xlim([70 210])
% ylim([10 70])

toc


