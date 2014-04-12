% xseqlse10.m
% Extended Sequential Least Squares Estimator
%
% Ben Spivey, 11/26/10
disp('made it to xseq')
format long e

load('mimo_21ss_10zoh.mat');
[A,B,C,D] = ssdata(mimo_21ss_10zoh);
nA = length(A);
z1A = zeros(nA,ny); % assumes number of disturbances equals ny
z2A = zeros(ny,nA);
IA = eye(ny,ny);

% Input data from mpc algorithm
nx = nA + ny;

Aaug = [A z1A; C*A IA];
Baug = [B; C*B];
IC = eye(ny,ny);
zC = zeros(ny,nA);
Caug = [zC IC];

A = Aaug;
B = Baug;
C = Caug;

% Read observation data
time        = obsdata(:,1);
sofc_obs    = obsdata(:,2);
weight      = obsdata(:,3);
obsnumber   = obsdata(:,4);
tfinal      = max(time);
obsfinal    = length(time);

% Initialization
[X_0,xbar,sigmas,relerr,abserr] = setup10(nx);
Pbar0 = diag(sigmas.^2);
dXbar0 = xbar';
X_0 = X_0';
if exist('x0hist.mat')
    load('x0hist.mat');
    [row,col] = size(x0_mat);
    X_0 = x0_mat(:,col);
end

phi0 = eye(nx);
t_cur = 0;

% Convert phi0 to vector form for integration
for i=1:nx
    first = (i-1)*nx + 1;
    last = first + (nx-1);
    phi0_vec(first:last,1) = phi0(:,i);
end

X(:,1) = X_0;
state0_vec=[phi0_vec; X(:,1)];

P(:,:,1)=Pbar0;
for j=1:obsfinal
    t_obs = time(j);
    
    % Read current observation
    station_obs = obsnumber(j);
    weight_obs = weight(j);
    Y(j,1) = sofc_obs(j);
    % Integrate states and state transition matrix if necessary
    if j==1
        index = 1;
        % StateTemp = [36 Phi values; 6 State values]
        state0_vec=[phi0_vec; X(:,index)];
        StateTemp = phi_and_x_eqns(state0_vec,dupast,A,B,nu);

        % Convert StateTemp phi subvector to phi matrix
        for k=1:nx
            first = (k-1)*nx + 1;
            last = first + (nx-1);
            phi(:,k,j) = StateTemp(first:last,1);
        end

        first = nx^2+1;
        % Integrate the reference trajectory, X
        [row,col] = size(A); % assumes a square matrix
        size(X)
        first
        col
        size(StateTemp)
        X(:,j) = StateTemp(first:first+col-1);

    else
        phi(:,:,j) = eye(nx);
        X(:,j) = X(:,j-1);
    end

    % Compute a priori estimates: propagate state and covariance
    Pbar(:,:,j) = phi(:,:,j)*P(:,:,index)*phi(:,:,j)';
    
    % Compute residual, Htilda, Kalman gain
    G = C*X(:,j);
    Y(j,1);
    dy(j,1) = Y(j,1) - G(station_obs,1);
    Htilda(j,:) = C(station_obs,:);
    R = 1/weight(j,1);
    K(:,j) = Pbar(:,:,j)*Htilda(j,:)'*inv(Htilda(j,:)*Pbar(:,:,j)*Htilda(j,:)' + R);
    yi = dy(j,1);

    % Compute a posterior estimates
    dXhat(:,j) = K(:,j)*yi;
    P(:,:,j) = (eye(nx) - K(:,j)*Htilda(j,:))*Pbar(:,:,j);
    Xprior(:,j) = X(:,j);

    % Update the state vector
    X(:,j) = X(:,j) + dXhat(:,j);
    
    % Calculate postfit residual
    Gpf = C*X(:,j);
    dypf(j,1) = Y(j,1) - Gpf(station_obs,1);

end

% Calculate RMS
pre_rms = 0;
post_rms = 0;
disp('residuals')
pre_rms = sqrt(dy(:,1)'*dy(:,1)/obsfinal);
post_rms = sqrt(dypf(:,1)'*dypf(:,1)/obsfinal);
X(:,obsfinal)

% Calculate covariance and correlation matrices
sigmas_final = sqrt(diag(P(:,:,obsfinal)));
sigmas_final_mat = diag(sigmas_final);
inv_sigmas_mat = inv(sigmas_final_mat);
corr_mat = inv_sigmas_mat*P(:,:,obsfinal)*inv_sigmas_mat - eye(nx) + sigmas_final_mat;

% Output results to screen
% diary('ben.out.extended_seq.03.txt')
% disp('       End of file encountered')
% fprintf('\n\n\n')
% fprintf('No. of data points= %d   ', obsfinal)
% fprintf('RMS before filter= %d   ', pre_rms)
% fprintf('RMS after filter= %d\n\n\n', post_rms)
% fprintf('       XHAT_0\n')
% disp(dXhat(:,obsfinal))
% fprintf('       Time of estimated state')
% disp(tfinal)
% fprintf('       Estimated State\n')
% disp(X(:,obsfinal))
% fprintf('       Covariance Matrix\n')
% disp(P(:,:,obsfinal))
% fprintf('       Correlation Matrix\n')
% disp(corr_mat)
% diary off

% Save results to a matrix file
save('xseqsle0.mat','X_0','xbar','sigmas');

% figure
% plot(dy(:,1))
% figure
% plot(dypf(:,1))

    
        
