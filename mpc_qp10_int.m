% SOFC MPC QP Controller

function [f,g] = mpc_qp10_int(uaug,x0,upast,A,B,C,yrvec,ph,ch,Qweights,Rweights,O,L,Ct,P,R,Qt,V,W,nx,ny,nu,UP1,UP2,UP3,t,sim_t_final)

%% System Matrices
%% Initialize variables
% x0 = zeros(nx,1);
yrvec = yrvec(1:ph*ny,1);

%% Input handling
u = uaug(1:nu*ch,1);
slk = uaug(nu*ch+1:nu*ch+2*ny*(ph-1),1);
% slk2 = uaug(nu*ch+2*ny*(ph-1)+1:length(uaug)); DEBUG

%% Solve for the optimal u profile using QP
du = (UP1 - UP2)*u - UP3*upast;

% Regulatory approach
% xdiff = L*(du - dusvec) + O*(x0 - xs)
% Tracking approach
xdiff = L*du + O*x0 - Ct*yrvec;

% f = 0.5*xdiff'*Qt*xdiff + 0.5*du'*R*du + 0.5*slk'*V*slk +
% 0.5*slk2'*W*slk2; DEBUG
f = 0.5*xdiff'*Qt*xdiff + 0.5*du'*R*du + 0.5*slk'*V*slk;

if nargout > 1
    g = [[xdiff'*Qt*L*(UP1 - UP2) + du'*R*(UP1 - UP2)], slk'*V];
%     g = [[xdiff'*Qt*L*(UP1 - UP2) + du'*R*(UP1 - UP2)], slk'*V, slk2'*
%     W]; DEBUG
end

% if flag == 1
%     save('TEST_u0.mat','xdiff','du','slk','f');
% elseif flag == 2
%     save('TEST_uaug.mat','xdiff','du','slk','f');  
% end

if t > sim_t_final
% GAMS CODE
size(xdiff)
size(du)
size(slk)
size(x0)
size(yrvec)
size(uaug)

save('mpcdata.mat','Qt','R','V','L','O','Ct','UP1','UP2','UP3','x0','upast','yrvec');
save('debug.mat','xdiff','Qt','yrvec','Ct','J');
end




