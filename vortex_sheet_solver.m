%
%
%           Simple Vortex sheet solver
%
%
%==========================================================~
clear all; %close all;

%% Paramters ===============================
%St= 0.1:0.001:0.2;   % strouhal numbers to be computed
%St= 0.5:-0.001:0.4;   % strouhal numbers to be computed
St=[0.4, 0.4];

params.Mj    = 1.2;  % jet Mach number
params.gam   = 1.4;  % ratio of specific heats
params.S     = 1;    % density ratio (1 => isothermal jet) =0; % azimuthal mode number
params.m     = 0;    % azimuthal mode number
params.n     = 1;    % radial mode number
params.Rj    = 0.5;  % jet radius

%% initial guess solutions =================
omega = 2*pi*St*params.Mj;
k_th  = omega/(params.Mj) - 2.1j;


%% Computing the dispersion relation roots =
options = optimset('TolX',1e-20,'TolFun',1e-20);
for i1=2:length(St)
    k_th(i1)  = fsolve(@(kk) vs_resid(omega(i1),kk,params), k_th(i1-1), options);
    % if (abs(imag(k_th(i1))) > 1e-12)
    %     %the mode no longer sits on the real axis
    %     break;
    % end
end;
St_th = St(2:i1)
k_th  = k_th(2:i1)


%% Plot the dispersion relations ===========
%sonic waves :
ka = -linspace(0,10,100);
Sta = -ka/(2*pi*params.Mj);

%vonvected waves :
k_conv = linspace(0,10,100);
St_conv = k_conv/(0.9*2*pi*params.Mj);

% figure(101); clf;
% hold on;
% plot(ka,Sta,'k-');
% plot(k_conv,St_conv,'k--');
% plot(real(k_th),St_th,'ks-');
% axis([-20 20 0 2]);
% xlabel('k');
% ylabel('St');
% drawnow;
