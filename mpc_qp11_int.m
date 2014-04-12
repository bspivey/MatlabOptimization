% SOFC MPC QP Controller

function J = mpc_qp11_int(uaug,x0,upast,A,B,C,yrvec,ph,ch,Qweights,Rweights,O,L,Ct,P,R,Qt,V,nx,ny,nu,UP1,UP2,UP3)

%% System Matrices
%% Initialize variables
% x0 = zeros(nx,1);
yrvec = yrvec(1:ph*ny,1);

%% Input handling
u = uaug(1:nu*ch,1);
slk = uaug(nu*ch+1:length(uaug),1);

%% Solve for the optimal u profile using QP
du = (UP1 - UP2)*u - UP3*upast;

% Regulatory approach
% xdiff = L*(du - dusvec) + O*(x0 - xs)
% Tracking approach
xdiff = L*du + O*x0 - Ct*yrvec;

J = 0.5*xdiff'*Qt*xdiff + 0.5*du'*R*du + 0.5*slk'*V*slk;





