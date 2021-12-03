clear all;close all;

%% define a pinhole camera

% intrinsic camera parameters
W = 1280; % width in pixels (1280)  uo=W/2 (center of screen)
H = 1024; % height in pixels (1024) vo=H/2 (center of screen)
rhow = 1e-5; % width per pixel (10um)
rhoh = 1e-5; % height per pixel (10um)
rho=[rhow;rhoh];
f = .015; % focal length (0.015m)

phih=atan2(H*rhoh/2,f); % vertical viewing angle range
phiw=atan2(W*rhow/2,f); % horizontal viewing angle range

u0=W/2; %center of image plane in pixel coordinate
v0=H/2; %center of image plane in pixel coordinate
uv0=[u0;v0];

cam=pinholecam(f,rho,uv0);


%% define target points in the inertial frame

% planar target points in the reference frame (total Nx*Ny points)
z0=0; % z-coordinate of the points (planar, so 0)
%nx1=-.2;kx=.1;nx2=.2; % arrays in x
%ny1=-.2;ky=.1;ny2=.2; % arrays in y
nx1=-.2;kx=.2;nx2=.2; % arrays in x
ny1=-.2;ky=.2;ny2=.2; % arrays in y
px=(nx1:kx:nx2);Nx=length(px);
py=(ny1:ky:ny2);Ny=length(py);
Pxy=kron([zeros(1,Nx);ones(1,Nx)],py)+...
  kron(px,[ones(1,Ny);zeros(1,Ny)]);
N = Nx*Ny;
P0=[Pxy;z0*ones(1,Nx*Ny)]; % pack into a 2x(Nx*Ny) vector

%% Specify # of images
M=20;

%% calibration fit matrix
L=zeros(2*M,6);

%% build up the calibration fit matrix camcalib
for i= 1:M
    
    % set image generation flag to 1
    image_gen_flag = 1;
    
    while image_gen_flag > 0
        
        % Randomly generate R_{oc} (p_{0c} set so targets are in front of camera)
        theta=rand*pi/10;k=randn(3,1);k=k/norm(k);
        p=[.8*(rand-.5);.8*(rand-.5);-(rand+1)*3];
        Toc{i} = [rot(k,theta) p; 0 0 0 1];
        
        % projection onto camera image plane
        [uv_true{i},uvw_true{i},Pc_true{i}]=cam_image(cam,Toc{i},P0);
        %[uv_true{i},uvw_true{i},P1]=cam_image_round(cam,Toc{i},P0);
        
        if size(uv_true{i},2)<N;
            disp('some patterns are out of view');
        else
            image_gen_flag = 0;
        end
    
    end
    
    % ns pixel error
    ns=.5;
    uv_meas{i}=uv_true{i}+randn(size(uv_true{i}))*ns;
    
    % find the Homography matrix
    
    Hmat{i}=homographymat(uv_meas{i},Pxy);
    hx=Hmat{i}(:,1);
    hy=Hmat{i}(:,2);
    L(2*i-1,:)=[hx(1)^2-hy(1)^2 hx(2)^2-hy(2)^2 hx(3)^2-hy(3)^2 ...
        2*(hx(1)*hx(2)-hy(1)*hy(2)) 2*(hx(1)*hx(3)-hy(1)*hy(3)) ...
        2*(hx(2)*hx(3)-hy(2)*hy(3))];
    L(2*i,:)=[hx(1)*hy(1) hx(2)*hy(2) hx(3)*hy(3) ...
        hx(1)*hy(2)+hx(2)*hy(1) hx(1)*hy(3)+hx(3)*hy(1) ...
        hx(2)*hy(3)+hx(3)*hy(2)];    
end

[U,S,V]=svd(L);
svalues=sprintf(',  %0.4g',diag(S));
disp(['**** Singular Values ****: [',svalues(2:end),'  ]'])
Bvec=V(:,6);
Best=[Bvec(1) Bvec(4) Bvec(5);Bvec(4) Bvec(2) Bvec(6);Bvec(5) Bvec(6) Bvec(3)];
K_est=inv(chol(Best*Best(3,3)));K_est=K_est/K_est(3,3);
% Define the identified camera
cam_est=cam;
cam_est.K=K_est;
% compare the intrinsic camera matrix
disp('K');
disp(cam.K);
disp('Estimated K');
disp(cam_est.K);
disp('percentage error in K');
cc=sprintf(' %0.2g, %0.2g, %0.2g\n',abs((cam.K-cam_est.K)./cam.K)*100);
cc=strrep(cc,'NaN','-');
cc=strrep(cc,'Inf','-');
disp(cc);
%  reproject plot
% 
% plot color array
plot_color=jet(M);zmax=0;
% 
for i=1:M    
  % estimate T_co
  
  Ti=inv(cam_est.K)*Hmat{i};
  ai=mean([norm(Ti(:,1)) norm(Ti(:,2))]);
  Rpi=Ti/ai;
  Rpi=sign(Rpi(3,3))*Rpi;
  Tco_est=[Rpi(:,1) Rpi(:,2) cross(Rpi(:,1),Rpi(:,2)) Rpi(:,3);0 0 0 1];
  Toc_est{i}=inv(Tco_est);
  [uv_est{i},uvw_est{i},Pc_est{i}]=cam_image(cam_est,Toc_est{i},P0);
 
  % error 
  erruv{i}=uv_true{i}(:,1:size(uv_est{i},2))-uv_est{i};
  erruvnorm(i)=norm(erruv{i});
  erruvmax(i)=max(max(abs(erruv{i})));
  errPc{i}=Pc_true{i}-Pc_est{i};
  errPcnorm(i)=norm(errPc{i});
  errPcmax(i)=max(max(abs(errPc{i})));
end

%  plot in 3D 
figure(1000);
for i=1:M
  plot3(Pc_true{i}(1,:),Pc_true{i}(2,:),Pc_true{i}(3,:),'x','linewidth',1,...
      'color',plot_color(i,:));
  hold on;
  plot3(Pc_est{i}(1,:),Pc_est{i}(2,:),Pc_est{i}(3,:),'o','linewidth',1,...
      'color',plot_color(i,:));
end
xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');
title('3D reprojection based on camera calibration');
grid on;
zz=zlim(gca);
zmax=max(zz(2),zmax);
view([-96,-43])
axis([xlim(gca) ylim(gca) 0 zmax])
hold off

% plot in image plane
figure(2000);
for i=1:M
  plot(uv_true{i}(1,:),uv_true{i}(2,:),'x',...
      uv_est{i}(1,:),uv_est{i}(2,:),'o','linewidth',3);
  hold on
end
title('image plane reprojection based on camera calibration');
grid on;
view(180,-90);
axis([0 W 0 H]);
hold off
% max reprojection error
mxerr=sprintf(',  %4.3g',max(erruvmax));
%mxerr=sprintf(',  %4.3g',max(errPcmax));
disp(['**** Maximum Reprojection Error:  ',mxerr(2:end)])
hold off

%% camcalib_depreciated

% calibration matrix
L=zeros(2*M,6);

% build up the calibration matrix
for i= 1:M
    % Randomly generate R_{oc} (p_{0c} set so targets are in front of camera)
    theta=rand*pi/10;k=randn(3,1);k=k/norm(k);
    p=[.8*(rand-.5);.8*(rand-.5);-(rand+1)*3];
    Toc{i} = [rot(k,theta) p; 0 0 0 1];
    
    % projection onto camera image plane
    %[uv{i},uvw{i},P1]=cam_image(cam,Toc{i},P0);
    [uv{i},uvw{i},P1]=cam_image_round(cam,Toc{i},P0);
    
    if size(uv{i},2)<N;
        disp('some patterns are out of view');
        disp(uv{i});i=i-1;
        return;
    end
    
    % ns pixel error
    ns=.5;
    uv{i}=uv{i}+randn(size(uv{i}))*ns;
    
    % find the Homography matrix
    
    Hmat{i}=homographymat(uv{i},Pxy);
    hx=Hmat{i}(:,1);
    hy=Hmat{i}(:,2);
    L(2*i-1,:)=[hx(1)^2-hy(1)^2 hx(2)^2-hy(2)^2 hx(3)^2-hy(3)^2 ...
        2*(hx(1)*hx(2)-hy(1)*hy(2)) 2*(hx(1)*hx(3)-hy(1)*hy(3)) ...
        2*(hx(2)*hx(3)-hy(2)*hy(3))];
    L(2*i,:)=[hx(1)*hy(1) hx(2)*hy(2) hx(3)*hy(3) ...
        hx(1)*hy(2)+hx(2)*hy(1) hx(1)*hy(3)+hx(3)*hy(1) ...
        hx(2)*hy(3)+hx(3)*hy(2)];    
end

[U,S,V]=svd(L);
Bvec=V(:,6);
Best=[Bvec(1) Bvec(4) Bvec(5);Bvec(4) Bvec(2) Bvec(6);Bvec(5) Bvec(6) Bvec(3)];
KKtr=inv(Best);
KKtr=KKtr/KKtr(3,3);
%
tt=eye(3,3);
tt=tt(3:-1:1,:);
X=chol(tt*KKtr*tt);
Kest=tt*X'*tt;
camest=cam;
camest.K=Kest;
disp('K vs. K_est');
disp(cam.K);
disp(Kest);
for i=1:M
  [uv1{i},uvw1{i},P1]=cam_image(camest,Toc{i},P0);
  figure(M);
  plot(uv{i}(1,:),uv{i}(2,:),'x',uv1{i}(1,:),uv1{i}(2,:),'o','linewidth',3);
  title('reprojection based on camera calibration');
  grid on;
  view(180,-90);
  axis([0 W 0 H]);
  erruv{i}=uv{i}-uv1{i};
  erruvnorm(i)=norm(erruv{i});
  erruvmax(i)=max(max(abs(erruv{i})));
end
hold off


%% functions
% pinholecam.m
%
% generates a pinhole camera object
%
% usage: cam=pinholecam(f,rhox,rhoy,u0,v0)
% 
% input: 
% f: focal length
% rhox: width/pixel
% rhoy: height/pixel
% u0,v0: pixel domain location of optical center
% 
% output:
% cam: camera object with attributes f, rhox, rhoy, u0, v0, and
% K=intrinsic camera matrix
% 

function cam=pinholecam(f,rho,uv0)

cam.f=f;
cam.rho=rho;
cam.uv0=uv0;
cam.K=[f/rho(1) 0 uv0(1);0 f/rho(2) uv0(2);0 0 1];

end


function R = rot(k,theta)
    k=k/norm(k);
    R=eye(3,3)+sin(theta)*hat(k)+(1-cos(theta))*hat(k)*hat(k);
end 


%
% hat.m (converting a vector into a skew-symmetric cross-product matrix
%
% khat = hat(k)
%

function khat = hat(k)
  
  khat=[0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];

end

%
% cam_image_round: 
%
% generating camera images of a set of 3D points on the image plane for a
% specified pinhole camera.   the image coordinates around rounded to the
% nearest integer.
%
% input:
% cam: pinhole camera object
% Toc: camera location in reference frame
% P0: a set of 3D points in reference frame 3xN matrix
%
% output: 
% uv: image coordinate of P0
% uvw: homogenous coordinate of P0
% P1: target points in the camera frame
%

function [uv,uvw,P1]=cam_image_round(cam,Toc,P0)

K=cam.K; % intrinsic camera matrix
Tco=inv(Toc);
C=K*Tco(1:3,:); %full camera matrix

% points in the camera frame
P1=Tco(1:3,1:3)*P0;

% (u',v',w') 
uvw=C*[P0;ones(1,size(P0,2))];

% image plane coordinate in pixels
u=round(uvw(1,:)./uvw(3,:));
v=round(uvw(2,:)./uvw(3,:));
uv=[u;v];

% only keep points within the image plane
uv=uv(:,uv(1,:)<2*cam.uv0(1));
uv=uv(:,uv(2,:)<2*cam.uv0(2));
uv=uv(:,uv(1,:)>0);
uv=uv(:,uv(2,:)>0);

end

%
% cam_image: 
%
% generating camera images of a set of 3D points on the image plane for a
% specified pinhole camera
%
% input:
% cam: pinhole camera object
% Toc: camera location in reference frame
% P0: a set of 3D points in reference frame 3xN matrix
%
% output: 
% uv: image coordinate of P0
% uvw: homogeneous coordinate of P0
% P1: target points in the camera frame
%
function [uv,uvw,P1]=cam_image(cam,Toc,P0)

K=cam.K; % intrinsic camera matrix
Tco=inv(Toc);
C=K*Tco(1:3,:); %full camera matrix

% points in the camera frame
P1=Tco(1:3,1:3)*P0;

% (u',v',w') 
uvw=C*[P0;ones(1,size(P0,2))];

% image plane coordinate in pixels
u=(uvw(1,:)./uvw(3,:));
v=(uvw(2,:)./uvw(3,:));
uv=[u;v];

% only keep points within the image plane
uv=uv(:,uv(1,:)<2*cam.uv0(1));
uv=uv(:,uv(2,:)<2*cam.uv0(2));
uv=uv(:,uv(1,:)>0);
uv=uv(:,uv(2,:)>0);

end

%
% homographymat.m
%
% Find the homography matrix 
%
% usage: H=homographymat(uv,Pxy)
%
% input:
%   uv: image pixel coordinate of N points (2xN matrix)
%   Pxy: planar set of points in the inertial frame
% 
% output:
%   H: Homography matrix
%

function H=homographymat(uv,Pxy)

N=size(Pxy,2);
A=zeros(2*N,9);
z3=zeros(1,3);
alpha=[uv;ones(1,N)];
for i=1:N
    beta=[Pxy(:,i);1];
    A(2*(i-1)+1:2*i,:)=...
        [z3 -alpha(3,i)*beta' alpha(2,i)*beta';...
        alpha(3,i)*beta' z3 -alpha(1,i)*beta'];
end
[U,S,V]=svd(A);h=V(:,9);
%h=-null(A);
H=[h(1:3)';h(4:6)';h(7:9)'];
H=H/norm(H(:,1));

end

%
% pnp_general.m: solution of the general PnP problem
%
% **** Note: this does not work for planar targets (null space is dimension
% 4, and an additional fit is necessary) ****  
% 
% input: 
%       uv = image coordinates in pixels of target points
%       P0 = locations of target points in reference frame
%       K = camera intrinsic calibration matrix
%
% output: 
%       Test = T_C0 estimate
%       Zest = depth estimate of each target point
% 

function [Test,Zest]=pnp_general(uv,P0,K)

% # of target points
N=size(P0,2);

% multiply image by K inverse
y=inv(K)*[uv;ones(1,N)];
% 
%y=reshape([uv;ones(1,N)],3*N,1);
%Pstack=reshape(P0,3*N,1);

% form the regression matrix
C=zeros(3*N,N+3+9);
for i=1:N
   C((i-1)*3+1:3*i,1:3*3)=kron(eye(3,3),P0(:,i)');
   C((i-1)*3+1:3*i,3*3+1:3*(3+1))=eye(3,3);
   C((i-1)*3+1:3*i,3*(3+1)+i)=y(:,i);
end

% use svd to estimate the near null space 
% (check the last few smallest singular values!)
[U,S,V]=svd(C);
% estimate null space
% null space of C is 4 dimensional if target points are in a plane
nC=V(:,12+N);
% extract orientation
R1=reshape(nC(1:9),3,3)';
% use svd to extract orthonormal components
[u,s,v]=svd(R1);
% estimate orientation matrix R_{co}
Rest=u*v'*sign(det(R1));
% find the scaling constant 
a=pinv(reshape(Rest',9,1))*nC(1:9);
% recover p_{co} and depth by scaling appropriately
pest=nC(10:12)/a;
Zest=nC(13:13+N-1)/a;
% This is the estimate T_{co}
Test=[Rest pest;0 0 0 1];
end

