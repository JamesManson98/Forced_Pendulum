close all
forced_pendulum_660031764
%% finding the fixed point for a=0
u0=[-0.5,0,0.5;0,0,0];  %defining initial guesses for the equilibria
T=2*pi/omega;           %defining the time period T
a=0;                    %this is the case where amplitude=0
N=round(T/0.005,0);
M=@(u,a)MyIVPVec(@(t,u)OdeRhs(t,u,a),u,[0,T],N);   %defining the function M as required
df=@(u)MyJacobian(@(u)(M(u,a)-u),u,1e-6);   %finding df, the Jacobian in terms of u
x=zeros(2,length(u0));  %preallocating x for speed
for i=1:length(u0)
    [x(:,i),converged,J]=MySolve(@(u)(M(u,a)-u),u0(:,i),df,1e-4,10);    %we find all the equilibria for all of the guesses in u0
end
%% finding which points go to which fixed points for a=0
res=201;        %we choose a value of resolution to create a meshgrid
xincrement=(1-(-1))/(res-1);  %we find the increment to increase x1 by in the meshgrid
yincrement=(6-(-6))/(res-1);  %we find the increment to increase y1 by in the meshgrid
x1=-1:xincrement:1;   %we create a vector of x points to put into meshgrid
y1=-6:yincrement:6;   %we create a vector of y points to put into meshgrid
[X,Y]=meshgrid(x1,y1);  %we create a grid of values from x1=-1 to 1 and y1=-6 to 6 of resolution res*res
B=reshape(X,[1,res*res]);   %we reshape all the x values into a 1*(res*res)
C=reshape(Y,[1,res*res]);   %we reshape all the y values into a 1*(res*res)
A=[B;C];    %we pair all the x and y values that correspond to each other in B and C
D=zeros(1,res*res); %we preallocate a row of zeros with length res*res
for i=1:200
    A=M(A,a);   %we iterate the Map M a maximum of 80 times
    for j=1:res*res
        if norm(A(:,j)-x(:,1))<(2*1e-1) %we find which equilibria each initial point tends to as T->infinity and categorise them into 4 groups for the 3 equilibria and for those that do not converge successfully
            D(j)=0;
        elseif norm(A(:,j)-x(:,2))<(2*1e-1)
            D(j)=1;
        elseif norm(A(:,j)-x(:,3))<(2*1e-1)
            D(j)=2;
        else
            D(j)=3;
        end
    end
    if ~ismember(3,D)   %if D(j) is not 3, every point in the matrix A has converged to a fixed point
        break
    end
end
E_1=reshape(D,[res,res]); %we reshape this into a resxres grid that can be plotted by imagesc
xlim=[-1,1];    %we define the x limit of the image we produce
ylim=[-6,6];    %we define the y limit of the image we produce
imagesc(xlim,ylim,E_1)    %we plot using imagesc
set(gca,'YDir','normal')    %we flip the y direction of imagesc
colormap([1 0 0; 1 1 0; 1 1 1; 0 0 1])  %we define the colours used in imagesc
hold on
plot(x(1,1),x(2,1),'x','color',[0 1 0],'markersize',14) %we plot the 3 equilibrium points
plot(x(1,2),x(2,2),'+','color',[0 1 0],'markersize',14)
plot(x(1,3),x(2,3),'*','color',[0 1 0],'markersize',14)
title('a=0')    %we define the titles
xlabel('u_1 (=\theta)')     %we define the x-axis label
ylabel('u_2 (=d\theta/dt)') %we define the y-axis label
legend('u_-','u_c','u_+');  %we define the legend for the 3 equilibrium points
%% finding the fixed point for a=0.007
u0=x;    %defining initial guesses for the equilibria
a=0.007;                    %this is the case where amplitude=0.0005
x=zeros(2,length(u0));  %preallocating x for speed
for i=1:length(u0)
    [x(:,i),converged,J]=MySolve(@(u)(M(u,a)-u),u0(:,i),df,1e-4,10);    %we find all the equilibria for all of the guesses in u0
end
%% finding which points go to which fixed points for a=0.007
A=[B;C];    %we pair all the x and y values that correspond to each other in B and C
D=zeros(1,res*res); %we preallocate a row of zeros with length res*res
for i=1:200
    A=M(A,a);   %we iterate the Map M 25 times
    for j=1:res*res
        if norm(A(:,j)-x(:,1))<(2*1e-2) %we find which equilibria each initial point tends to as T->infinity and categorise them into 4 groups for the 3 equilibria and for those that do not converge successfully
            D(j)=0;
        elseif norm(A(:,j)-x(:,2))<(2*1e-2)
            D(j)=1;
        elseif norm(A(:,j)-x(:,3))<(2*1e-2)
            D(j)=2;
        else
            D(j)=3;
        end
    end
    if ~ismember(3,D)    %if D(j) is not 3, every point in the matrix A has converged to a fixed point
        break
    end
end
E_2=reshape(D,[res,res]); %we reshape this into a resxres grid that can be plotted by imagesc
figure(2)
imagesc(xlim,ylim,E_2)    %we plot using imagesc
set(gca,'YDir','normal')    %we flip the y direction of imagesc
colormap([1 0 0; 1 1 0; 1 1 1; 0 0 1])  %we define the colours used in imagesc
hold on
plot(x(1,1),x(2,1),'x','color',[0 1 0],'markersize',14) %we plot the 3 equilibria
plot(x(1,2),x(2,2),'+','color',[0 1 0],'markersize',14)
plot(x(1,3),x(2,3),'*','color',[0 1 0],'markersize',14)
title('a=0.007')        %we define the title, labels and axis
xlabel('u_1 (=\theta)')
ylabel('u_2 (=d\theta/dt)')
legend('u_-','u_c','u_+');
%% finding the fixed point for a=0.008
a=0.008;                    %this is the case where amplitude=0.007
x=zeros(2,length(u0));  %preallocating x for speed
for i=1:length(u0)
    [x(:,i),converged,J]=MySolve(@(u)(M(u,a)-u),u0(:,i),df,1e-4,10);    %we find all the equilibria for all of the guesses in u0
end
%% finding which points go to which fixed points for a=0.008
A=[B;C];    %we pair all the x and y values that correspond to each other in B and C
D=zeros(1,res*res); %we preallocate a row of zeros with length res*res
for i=1:200
    A=M(A,a);   %we iterate the Map M 25 times
    for j=1:res*res
        if norm(A(:,j)-x(:,1))<(2*1e-2) %we find which equilibria each initial point tends to as T->infinity and categorise them into 4 groups for the 3 equilibria and for those that do not converge successfully
            D(j)=0;
        elseif norm(A(:,j)-x(:,2))<(2*1e-2)
            D(j)=1;
        elseif norm(A(:,j)-x(:,3))<(2*1e-2)
            D(j)=2;
        else
            D(j)=3;
        end
    end
    if ~ismember(3,D)   %if D(j) is not 3, every point in the matrix A has converged to a fixed point
        break
    end
end
E_3=reshape(D,[res,res]); %we reshape this into a resxres grid that can be plotted by imagesc
figure(3)
imagesc(xlim,ylim,E_3)    %we plot using imagesc
set(gca,'YDir','normal')    %we flip the y direction of imagesc
colormap([1 0 0; 1 1 0; 1 1 1; 0 0 1])  %we define the colours used in imagesc
hold on
plot(x(1,1),x(2,1),'x','color',[0 1 0],'markersize',14) %we plot all of the equilibrium points
plot(x(1,2),x(2,2),'+','color',[0 1 0],'markersize',14)
plot(x(1,3),x(2,3),'*','color',[0 1 0],'markersize',14)
title('a=0.008')        %we define the title, labels and axis
xlabel('u_1 (=\theta)')
ylabel('u_2 (=d\theta/dt)')
legend('u_-','u_c','u_+');