%% Do not edit any part of this file except the line SNumber=0;
% Model script for CW2 of ECM3739 autumn 2018
% variables; angle: theta (u1), angular velocity: dot{theta} (u2)
% free parameter; forcing amplitude: a 
% f=OdeRhs(t,u,a) where u(1)=theta, u(2)=dot{theta}
% You will need to first create a map M that integrates OdeRhs forward 
% exactly T=2*pi/omega using MyIVP referred to as the scroboscopic map in
% lecture notes
clear 
format compact

%% Set your personal parameter value for theta and beta
% insert your student number xxxyyyzzz (usually starts with a 6)
% rename this file as forced_pendulum_xxxyyyzzz.m and call this
% in your scripts for each questions. There is no need to modify any 
% part of this cript except the next line
SNumber=660031764;
if SNumber==0
   error(['Enter your 9 digit student number in your local copy of forced_oscillator.m   ',...
            'See instruction sheet and contact James Rankin at j.a.rankin@exeter.ac.uk if unsure.']); 
end
rng(SNumber)
phase = 2*pi*rand;
dm=0.12+0.035*rand;
format longg
disp('Your personal value of parameter dm:');
disp(dm)
disp('Your personal value of parameter phase:');
disp(phase)
format short

%% parameters and forcing terms
plen=0.2;
height=0.01;
g=9.81;
magnetforce=0.05;
damping=0.3;
omega=16;
dist=@(d,theta)(d-plen*sin(theta)).^2+(plen*cos(theta)-(plen+height)).^2;
fproj=@(d,theta)d.*cos(theta)-sin(theta)*(plen+height);
Fm=@(d,theta)magnetforce*fproj(d,theta)./dist(d,theta).^1.5;
%% sum of forces
F_fpm=@(t,theta,thetad,a)...
    -damping*thetad...                      % damping
    +a.*omega^2*sin(omega*t+phase).*cos(theta)... % forcing
    -sin(theta)*g...                        % gravity
    +Fm(dm,theta)+Fm(-dm,theta);            % two magnets (symmetrically placed)
OdeRhs=@(t,u,a)[...
    u(2,:);
    F_fpm(t,u(1,:),u(2,:),a)/plen];
