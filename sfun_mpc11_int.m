function [sys,x0,str,ts] = sfun_mpc11_int(t,x,u,flag)

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    sizes = simsizes;
    
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 1; % a discrete state is chosen for default purposes but not used
    sizes.NumOutputs     = 2;
    sizes.NumInputs      = 6;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1; % at least one sample time is needed
    
    sys = simsizes(sizes);
    str = [];
    x0  = [1];
    ts  = [10 0];
    
  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,   
    disp('sample time')
    disp(t)
    tic
    disp('mpc case 3 is called')
    
    ny = 1;
    nu = 1;
    nx = 1;
    
    % Assigns inputs
    ymeas   = u(1:ny,1)
    yref    = u(ny+1:2*ny,1);
    umin    = u(2*ny+1:2*ny+nu,1);
    umax    = u(2*ny+nu+1:2*ny+2*nu,1);
    ymin    = u(2*ny+2*nu+1:2*ny+2*nu+ny,1);
    ymax    = u(3*ny+2*nu+1:3*ny+2*nu+ny,1);
    
    % Create obsdata matrix
    obsdata(:,1) = ones(ny,1); % time
    obsdata(:,2) = ymeas; % observations
    obsdata(:,3) = 100*ones(ny,1);  % weight
    obsdata(:,4) = [1:ny]'; % obsnumber
    
    % Estimate the state vector
    if exist('dupast2.mat', 'file')
        dustruc = load('dupast2.mat');
        dupast = dustruc.du;
    else
        dupast = zeros(nu,1);
    end
    xseqlse11_int;
    disp('x0')
    x0 = X(:,obsfinal)
    
    % Solve for the next optimal 'move'
    if exist('upast2.mat', 'file')
        ustruc = load('upast2.mat');
        upast = ustruc.upast;
    else
        upast = zeros(nu,1);
    end
    run_mpc11_int;
    
    ypred = yp;
    disp('ypred size')
    size(ypred)
    ypred
    delete('ypred2.mat');
    save('ypred2.mat','ypred');
    
    yp_temp = yp;
    if exist('yp2.mat')
        load('yp2.mat');
        size(yp)
        size(yp_temp)
        yp = [yp; yp_temp]
        delete('yp2.mat');
        save('yp2.mat','yp');
    else
        save('yp2.mat','yp');
    end
    
    % Print output to screen
    disp('mpc inputs')
    fprintf('ymeas = %f\n\n', ymeas);
    fprintf('yref  = %f\n\n', yref);
    fprintf('umin = %f\n\n', umin);
    fprintf('umax = %f\n\n', umax);
    disp('mpc outputs')
    disp(mv)
    
    sys = [mv(1:nu);
           obj]
    
    % Store calculated MV for future retrieval
    upast = mv(1:nu);
    save('upast2.mat','upast');
    save('dupast2.mat','du');
    yref0=yref;
    save('yref2.mat','yref0');
    try
        y = evalin('base','y');
        u1 = evalin('base','u1');
        u2 = evalin('base','u2');
        u3 = evalin('base','u3');
        u4 = evalin('base','u4');
    catch
        y=[];
        u1=[];
        u2=[];
        u3=[];
        u4=[];
    end
%     save('y.mat','y');
%     save('u1.mat','u1');
%     save('u2.mat','u2');
%     save('u3.mat','u3');
%     save('u4.mat','u4');
    
    disp('MPC time:')
    toc
    
  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case {1 2 4 9}            % 1: continuous
                            % 2: discrete
    sys=[];                 % 4: calcTimeHit
                            % 9: termination
  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

% end sfuntmpl
