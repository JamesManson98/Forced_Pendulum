close all
forced_pendulum_660031764
T=2*pi/omega;           %defining the time period T
N=round(T/0.005,0);
a=0.019                    %we pick a value of a that causes chaos in the system
M=@(u,a)MyIVPVec(@(t,u)OdeRhs(t,u,a),u,[0,T],N);   %defining the function M as required
u=[-0.5096;0];    %defining initial guesses for the equilibria
[lambda,Rdiag,A]=LyapunovQR(@(u)M(u,a),u,5000); %we call LyapunovQR for 5000 steps
figure(1)
plot((1:5000),A(2,:))               %we plot and label, title and provide legends for the graphs required
hold on
plot((1:5000),A(1,:))
title('Trajectory vs. Iterate Number')
xlabel('Iterate Number') 
ylabel('Trajectory')
legend('u_2 (=d\theta/dt) trajectory','u_1 (=\theta) trajectory');
figure(2)
plot(A(1,(101:5000)),A(2,(101:5000)),'g.')
title('Trajectory in the (\theta,d\theta/dt) plane')
xlabel('u_1 (=\theta)') 
ylabel('u_2 (=d\theta/dt)')
figure(3)
plot((1:5000),cumsum(log(Rdiag),2)./repmat((1:5000),2,1))
title('Lyapunov Exponent vs. Iterate Number')
xlabel('Iterate Number') 
ylabel('Lyapunov Exponent')
legend('Lyapunov Exponent 1', 'Lyapunov Exponent 2')
lambda  %we print the output lambda after the 5000 iterations