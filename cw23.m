close all
forced_pendulum_660031764
%% finding the fixed point uc for a=0.008
u0=[0;0];    %defining initial guesses for the equilibria
T=2*pi/omega;           %defining the time period T
N=round(T/0.005,0);
a=0.008;                    %this is the case where amplitude=0.0005
M=@(u)MyIVPVec(@(t,u)OdeRhs(t,u,a),u,[0,T],N);   %defining the function M as required
df=@(u)MyJacobian(@(u)(M(u)-u),u,1e-6);   %finding df, the Jacobian in terms of u
[uc,~,~]=MySolve(@(u)(M(u)-u),u0,df,1e-4,10)        %we find all the equilibria for all of the guesses in u0
MInv=@(u0)MyIVPVec(@(t,u0)OdeRhs(t,u0,a),u0,[0,-T],N); %we define the inverse function MInv in terms of MyIVPVec which is to be integrated backwards in time
%% computing the M'(uc) and its eigenvectors/eigenvalues
A=@(u)MyJacobian(@(u)(M(u)),u,1e-6);
[V,~]=eig(A(uc));          %we compute the eigenvectors of partial M by partial u
[~,ind]=min(abs(eig(A(uc))));   %we find the minimum of the absolute value of the eigenvalues of A, which corresponds to the stable vector
[~,ind2]=max(abs(eig(A(uc))));  %we find the maximum of the absolute value of the eigenvalues of A, which corresponds to the unstable vector
Vstable=V(:,ind);   %we find the stable eigenvector
Vunstable=V(:,ind2);    %we find the unstable eigenvector
h=1e-2;
ynew=[uc-h*Vunstable,uc+h*Vunstable];   %we find a short line segment Wu
ynew1=[uc-h*Vstable,uc+h*Vstable];      %we find a short line segment Ws
snew=[0,1];     %we find s to put into MapLine
snew1=[0,1];
while length(ynew)<4000
    [ynew,xnew,snew]=MapLine(M,ynew,snew,0.01,5);   %we run MapLine and iterate it until ynew becomes a long curve segment
end
while length(ynew1)<5000
    [ynew1,xnew1,snew1]=MapLine(MInv,ynew1,snew1,0.01,5);   %we run MapLine and iterate it until ynew1 becomes a long curve segment
end
load('cw21vars.mat','E_1','E_2','E_3','xlim','ylim')
imagesc(xlim,ylim,E_3)    %we plot using imagesc
set(gca,'YDir','normal')    %we flip the y direction of imagesc
colormap([1 0 0; 1 1 0; 1 1 1; 0 0 1])  %we define the colours used in imagesc
hold on
plot(uc(1),uc(2),'+','color',[0 1 0],'markersize',14)   %plotting
plot(ynew(1,:),ynew(2,:),'color','c')
plot(ynew1(1,:),ynew1(2,:),'color', 'y')
title('Plot of stable manifolds, unstable manifolds and basins for a=0.008')
legend('Unstable manifold','Stable manifold')
legend('Saddle point','Unstable manifold','Stable manifold')
xlabel('u_1 (=\theta)')
ylabel('u_2 (=d\theta/dt)')
%% finding the fixed point uc for a=0.007
u0=[0;0];    %defining initial guesses for the equilibria
T=2*pi/omega;           %defining the time period T
N=round(T/0.005,0);
a=0.007;                    %this is the case where amplitude=0.0005
M=@(u)MyIVPVec(@(t,u)OdeRhs(t,u,a),u,[0,T],N);   %defining the function M as required
df=@(u)MyJacobian(@(u)(M(u)-u),u,1e-6);   %finding df, the Jacobian in terms of u
[uc,~,~]=MySolve(@(u)(M(u)-u),u0,df,1e-4,10)        %we find all the equilibria for all of the guesses in u0
MInv=@(u0)MyIVPVec(@(t,u0)OdeRhs(t,u0,a),u0,[0,-T],N); %we define the inverse function MInv in terms of MyIVPVec which is to be integrated backwards in time
%% computing the M'(uc) and its eigenvectors/eigenvalues
A=@(u)MyJacobian(@(u)(M(u)),u,1e-6);
[V,~]=eig(A(uc));          %we find Vm the eigenvector of M partially differentiated with respect to u
[~,ind]=min(abs(eig(A(uc))));   %we find the index of the minimum of maximum eigenvalues
[~,ind2]=max(abs(eig(A(uc))));
Vstable=V(:,ind);  %we define Vstable and Vunstable as the corresponding minimum and maximum eigenvalues
Vunstable=V(:,ind2);
h=1e-2; %this is the displacement either side of the saddle point to create the line segment for
ynew2=[uc-h*Vunstable,uc+h*Vunstable];  %we define the line segments ynew2 and ynew3 to put into MapLine
ynew3=[uc-h*Vstable,uc+h*Vstable];
snew2=[0,1];    %snew2 and snew3 are [0,1]
snew3=[0,1];
while length(ynew2)<4000
    [ynew2,xnew2,snew2]=MapLine(M,ynew2,snew2,0.01,5);  %we run MapLine and iterate it until ynew2 becomes a long curve segment
end
while length(ynew3)<5000
    [ynew3,xnew3,snew3]=MapLine(MInv,ynew3,snew3,0.01,5);   %we run MapLine and iterate it until ynew3 becomes a long curve segment
end
figure(2)
imagesc(xlim,ylim,E_2)    %we plot using imagesc
set(gca,'YDir','normal')    %we flip the y direction of imagesc
colormap([1 0 0; 1 1 0; 1 1 1; 0 0 1])  %we define the colours used in imagesc
hold on
plot(uc(1),uc(2),'+','color',[0 1 0],'markersize',14)   %we plot the equilibrium points
plot(ynew2(1,:),ynew2(2,:),'color','c')
plot(ynew3(1,:),ynew3(2,:),'color', 'y')
title('Plot of stable manifolds, unstable manifolds and basins for a=0.007')    %plotting figure
legend('Unstable manifold','Stable manifold')
legend('Saddle point','Unstable manifold','Stable manifold')
xlabel('u_1 (=\theta)')
ylabel('u_2 (=d\theta/dt)')
