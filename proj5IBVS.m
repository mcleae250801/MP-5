%
% miniproject 5 part 4
%
% generate estimate of T_0B and T_TC
%
hand_eye_cal; % from part 3

% load target image
load targetS_part4
uv_des=uv;

% set up IBVS call
N=200;
irb1200.MaxIter=N; % # of iteration for kinematic control
irb1200.StepSize=.03;   % controller gain

% target camera frame with respect to reference frame 0
T_0C_des=[-ey ez -ex [1.8;0;0.55];0 0 0 1];
% PBVS to reduce error
%[q,uv,T_0C]=IBVS(irb1200,cam,pS,uv_des,T_0B,T_TC,T_0B_est,T_TC_est);
[q,uv,T_0C,JI,J,erruv]=IBVS_unlimited_range(irb1200,cam,pS,uv_des,T_0B,T_TC,T_0B_est,T_TC_est);

disp(T_0C{N+1})
disp(T_0C_des)

% show image space progression
for i=1:N 
    H=figure(50);plot(uv{i}(1,:),uv{i}(2,:),'x','linewidth',2);
    xlabel('image plane horizontal axis');
    ylabel('image plane vertical axis'); 
    axis([0 cam.uv0(1)*2 0 cam.uv0(2)*2]);
    set ( gca, 'xdir', 'reverse' );%pause
    %drawnow;
    drawnow;M_Image(i)=getframe(H);
end

for i=1:irb1200.MaxIter;minJI(i)=min(svd(JI{i}));minJ(i)=min(svd(J{i}));end
figure(61);plot((1:irb1200.MaxIter),minJI,'x','linewidth',2);
title('Minimum singular value of the image Jacobian in each step');
figure(62);plot((1:irb1200.MaxIter),minJ,'x','linewidth',2);
title('Minimum singular value of the robot Jacobian in each step');
figure(63);plot((1:irb1200.MaxIter),vecnorm(erruv),'x','LineWidth',2);
title('Image pixel error norm in each step');
H=figure(60);movie(H,M_Image,1,20);
