function [sys,x0,str,ts] = sfun_mpc10_int(t,x,u,flag)

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    sizes = simsizes;
    
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 1;
    sizes.NumOutputs     = 5;
    sizes.NumInputs      = 28;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;  % at least one sample time is needed
    
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
    
    ny = 5;
    nu = 4;
    nx = 21;
    
    % Assigns inputs
    ymeas   = u(1:ny,1)
    yref    = u(ny+1:2*ny,1);
    umin    = u(2*ny+1:2*ny+nu,1);
    umax    = u(2*ny+nu+1:2*ny+2*nu,1);
    ymin    = u(2*ny+2*nu+1:2*ny+2*nu+ny,1);
    ymax    = u(3*ny+2*nu+1:3*ny+2*nu+ny,1);
    
    if (t > 1000) && (t <= 1500)
        yref(1) = -0.35;
    elseif (t > 1500) && (t <= 2000)
        yref(1) = 0.65;
    elseif (t > 2000) && (t <= 2500)
        yref(1) = -0.1;
    else
        yref(1) = 1.2;
    end
    
    % Create obsdata matrix
    obsdata(:,1) = ones(ny,1); % time
    obsdata(:,2) = ymeas; % observations
    obsdata(:,3) = 100*ones(ny,1);  % weight
    obsdata(:,4) = [1:ny]'; % obsnumber
    
    % Estimate the state vector
    if exist('dupast.mat', 'file')
        dustruc = load('dupast.mat');
        dupast = dustruc.du;
    else
        dupast = zeros(nu,1);
    end
    xseqlse10_int;
    disp('x0')
    x0 = X(:,obsfinal)
    
    % Solve for the next optimal 'move'
    if exist('upast.mat', 'file')
        ustruc = load('upast.mat');
        upast = ustruc.upast;
    else
        upast = zeros(nu,1);
    end
    run_mpc10_int;
    
    ypred = yp;
    disp('ypred size')
    size(ypred)
    ypred
    delete('ypred.mat');
    save('ypred.mat','ypred');
    
    yp_temp = yp;
    if exist('yp.mat')
        load('yp.mat');
        size(yp)
        size(yp_temp)
        yp = [yp; yp_temp]
        delete('yp.mat');
        save('yp.mat','yp');
    else
        save('yp.mat','yp');
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
    save('upast.mat','upast');
    save('dupast.mat','du');
    yref0=yref;
    save('yref.mat','yref0');
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
    save('y.mat','y');
    save('u1.mat','u1');
    save('u2.mat','u2');
    save('u3.mat','u3');
    save('u4.mat','u4');
    
    disp('MPC time:')
    toc
    
    yref
    
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
