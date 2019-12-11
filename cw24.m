close all
forced_pendulum_660031764
%% step 1 period 1
T=2*pi/omega;           %defining the time period T
N=round(T/0.005,0);     %we round N to the nearest T/0.005 to put into the Map M
M=@(u,a)MyIVPVec(@(t,u)OdeRhs(t,u,a),u,[0,T],N);   %defining the function M as required
h=1e-4;     
u0=[0.5;0];
dif=@(u)MyJacobian(@(u)(M(u,0)-u),u0,1e-5);   %finding df, the Jacobian in terms of u
[uinitial,~,~]=MySolve(@(u)(M(u,0)-u),u0,dif,1e-4,10);    %we find all the equilibria for all of the guesses in u0
ytan=[0;0;1];   %we set MyTrackCurve going with a increasing
inistep=1e-3;   %we take a large initial step for j=1
j=1             %we output the step
f{j}=@(y)(M(y(1:2),y(3))-y(1:2)); %we define the equation f for the j=1 term
ytrack=@(y)[f{:}(y)];   %we track the curve at f
df=@(y)MyJacobian(ytrack,y,h);  %we find the jacobian of ytrack, to put into MyTrackCurve
ylist=MyTrackCurve(ytrack,df,[uinitial;0],ytan,j,M,'nmax',50,'smax',0.005,'stepsize',inistep); %we track the curve at f
ystore=cell(1,8);
ystore{1}=ylist;    %we store the full orbits in ystore
plot(ylist(3,:),ylist(1,:),'r');    %we plot the amplitude against theta
hold on
%% converging to the period doubling point pd1
differential=@(u,a)MyJacobian(@(u)M(u,a),u,1e-5);       %defining the input variable y2foldin
[V,D]=eig(differential(ylist((1:2),end),ylist(3,end))); %we find the eigenvectors and eigenvalues     
[~,ind]=min(real(eig(differential(ylist((1:2),end),ylist(3,end))))+1);  %we find where the eigenvector where eigenvalues cross at -1
yguess=[ylist(:,end);V(:,ind)]; %we set the guess to input into MySolve to find the period doubling point
g=@(y)[f{:}(y);differential(y(1:2),y(3))*y(4:5)+y(4:5);(y(4:5))'*(y(4:5))-1];   %we define the system of equations to find the period doubling point 
df=@(y)MyJacobian(g,y,1e-5);    %we find the jacobian of g to put into MySolve
[x,converged,~]=MySolve(g,yguess,df,1e-8,10);   %we solve the system to find the period doubling point
apd=zeros(1,8); %we preallocate an array to store apd, the period doubling points, from 1 to 8.
apd(1)=x(3);    %the first period doubling point is inserted into apd
format long     %we set format to long so we can output the period doubling points to 5 significant figures
period_doubling_point=round(apd(1),5,'significant') %we output the first period doubling point to 5 significant figures
plot(x(3),x(1),'kx')    %we plot the period doubling point
title('bifurcation diagram of forcing amplitude vs. \theta')
xlabel('forcing amplitude') 
ylabel('\theta')
legend('branches of periodic orbits','period doubling points');
%% branching to periodic
for j=2:8
j   %we output what step of the orbit it's on the orbit is 2^(j-1)
r=1e-2/2^(j-1); %we define the initial step to take away from the period doubling point and scale it, such that as j increases, the initial stepsize decreases
uguess1=x(1:2^(j-1))+r*x(2^(j-1)+2:end);    %we define the first guess for all initial u, so we displace it by r from the period doubling point along the vector
uguess2=x(1:2^(j-1))-r*x(2^(j-1)+2:end);    %we define the second guess for all initial u, so we displace it by r from the period doubling point in the opposite direction to uguess1
y1=[uguess1;uguess2;x(2^(j-1)+1)];  %we use the value of a from the previous period doubling point
ytan=[x(2^(j-1)+2:end);-x(2^(j-1)+2:end);0];    %we set MyTrackCurve with an initial tangent [v1;-v1;0]
f=cell(2^(j-1),1);  %we produce an array to store the equations in
for z=1:2^(j-1)-1   
    f{z}=@(y)(M(y(2*z-1:2*z),y(2^j+1))-y(2*z+1:2*z+2)); %we store all of the equations in a cell
end
f{2^(j-1)}=@(y)(M(y(2^j-1:2^j),y(2^j+1))-y(1:2));   %we store the last equation in the cell
evalArray = @(y)cellfun(@(f)f(y),{f{:}},'Uni',false);   
ytrack=@(y)reshape(cell2mat(evalArray(y)),[2^j,1]); %we define ytrack by manipulating the cell arrays, so that it can be input into MyTrackCurve, by evaluating it and converting it to a matrix, then reshaping it
df=@(y)[Gu(y,-1,j,M),Ga(y,j,M)];    %we use the definition of the approximation of the jacobian of ytrack to put into MyTrackCurve
ylist=MyTrackCurve(ytrack,df,y1,ytan,j,M,'nmax',10000,'smax',0.01,'stepsize',inistep/2^(j-1),'smin',1e-1*inistep/2^(j-1)); %we run MyTrackCurve for the system and again scale the stepsize and smin, so that it's inversely proportional to the orbit 2^(j-1)
ystore{j}=ylist;    %we store the full orbits in ystore
sz=size(ylist); %we find the size of ylist
%% plotting
for J=1:2^(j-1)
plot(ylist(2^j+1,:),ylist(2*J-1,:),'r') %we plot all of the amplitudes against the thetas for every value in ylist and every branch
end
%% converging to the period doubling point pd1
hjac=1e-5;  %we define the h to be the deviation in the splitting difference method
Jacobian={@(y)Gu(y(1:2^j+1),-1,j,M),@(y)Ga(y(1:2^j+1),j,M),@(y)zeros(2^j);@(y)(1/(2*hjac))*(Gu(y(1:2^j+1)+hjac*[y(2^j+2:end);0],1,j,M)-Gu(y(1:2^j+1)-hjac*[y(2^j+2:end);0],1,j,M)),@(y)(1/(2*hjac))*(Gu(y(1:2^j+1)+[zeros(2^j,1);hjac],1,j,M)-Gu(y(1:2^j+1)-[zeros(2^j,1);hjac],1,j,M))*y(2^j+2:end),@(y)Gu(y(1:2^j+1),1,j,M);@(y)zeros(1,2^j),@(y)0,@(y)2*y(2^j+2:end)'};  %we define the Jacobian as required
[V,D]=eig(Gu(ylist(:,end),1,j,M));  %we find the eigenvectors 
[~,ind]=min(abs(real(eig(Gu(ylist(:,end),1,j,M)))));    %we find the corresponding eigenvector to put into our initial guess
yguess=[ylist(:,end);real(V(:,ind))];   %we find the yguess to put into MySolve to find the period doubling point
g=@(y)[ytrack(y);Gu(y(1:2^j+1),1,j,M)*y(2^j+2:end);(y(2^j+2:end))'*(y(2^j+2:end))-1];   %%we define the system of equations to find the period doubling point 
evalArray = @(y)reshape(cellfun(@(f)f(y),{Jacobian{:,:}},'Uni',false),[3,3]);   %we evaluate the system of functions and reshape it
C=@(y)(cell2mat(evalArray(y))); %we then convert it into an array
[x,converged,J]=MySolve(g,yguess,C,1e-7,10);    %we then solve the system of equations to find the period doubling point
apd(j)=x(2^j+1);    %then we save the period doubling point in an array apd
period_doubling_point=round(apd(j),5,'significant') %we then round and output the period doubling point to 5 significant figures
if converged==0 %if the system doesn't converge then we break, as something will have gone wrong
    break
end
for J=1:2^(j-1)
plot(x(2^j+1),x(2*J-1),'kx')    %we then plot all the period doubling points, against the theta values of the branches at the period doubling point
end
end
ratio=zeros(1,6);
for i=1:6
    ratio(i)=(apd(i+1)-apd(i))/(apd(i+2)-apd(i+1)); %we find all the ratios of successive distances between period doubling points in the parameter r_i
end    
r_1=ratio(1)    %we output all of the ratios of successive distances between period doublings r_1,...,r_6
r_2=ratio(2)
r_3=ratio(3)
r_4=ratio(4)
r_5=ratio(5)
r_6=ratio(6)