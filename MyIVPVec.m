function [xend]=MyIVPVec(f,x0,tspan,N,varargin)
h=(tspan(2)-tspan(1))/N;                    %we find the increment to increase t by each time
t=tspan(1):h:tspan(2);                      %creating a matrix with t ranging between tspan(1) and tspan(2) in increment h
x=x0;
for i=1:N
    k1=h*f(t(i),x(:,:));                   %we follow the Dormand-Prince method to explicitly solve the ODE for every i=1:N
    k2=h*f(t(i)+(1/5)*h,x(:,:)+(1/5)*k1);
    k3=h*f(t(i)+(3/10)*h,x(:,:)+(3/40)*k1+(9/40)*k2);
    k4=h*f(t(i)+(4/5)*h,x(:,:)+(44/45)*k1-(56/15)*k2+(32/9)*k3);
    k5=h*f(t(i)+(8/9)*h,x(:,:)+(19372/6561)*k1-(25360/2187)*k2+(64448/6561)*k3-(212/729)*k4);
    k6=h*f(t(i)+h,x(:,:)+(9017/3168)*k1-(355/33)*k2+(46732/5247)*k3+(49/176)*k4-(5103/18636)*k5);
    x(:,:)=x(:,:)+(35/384)*k1+(500/1113)*k3+(125/192)*k4-(2187/6784)*k5+(11/84)*k6;     
end
xend=x(:,:);                             %we find the last value of x and output as xend
end