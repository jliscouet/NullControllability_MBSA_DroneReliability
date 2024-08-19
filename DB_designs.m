% designs.m
% Define the different designs in this file

designs = struct;

% Design 1: HARV Hexarotor PNPNPN -----------------------------------------
n_d = 1;                                                                    % Design ID 
designs(n_d).name = 'HARV Hexarotor PNPNPN';                                % Design name
designs(n_d).m = 3.877;                                                     % Total mass (in kg)
designs(n_d).d = 0.4;                                                       % Rotor arrangement radius (in m)
designs(n_d).k = 0.03;                                                      % Reactive torque to thrust ratio (in m)
f_min = 0;                                                                  % Minimum rotor thrusts (in N);
f_max = 14.5;                                                               % Maximum rotor thrusts (in N);
I_xx = 7.998e-2;                                                            % Moment of inertia around the x-axis (in kg·m²)
I_yy = 7.884e-2;                                                            % Moment of inertia around the y-axis (in kg·m²)
I_zz = 0.148;                                                               % Moment of inertia around the z-axis (in kg·m²)
% Normalized control effectiveness matrix
th =           [1   1    1   1   1   1 ]; % Thrust
b1 = sqrt(3)/2*[0  -1   -1   0   1   1 ]; % Roll
b2 =           [1  0.5 -0.5 -1 -0.5 0.5]; % Pitch
b3 =           [1  -1    1  -1   1  -1 ]; % Yaw

J = diag([-designs(n_d).m, I_xx, I_yy, I_zz]);                              % Inertia matrix
designs(n_d).A = [zeros(4, 4), eye(4); zeros(4, 4), zeros(4, 4)];           % State matrix
designs(n_d).B = [zeros(4, 4); inv(J)];                                     % Control matrix
designs(n_d).B_f = diag([1 designs(n_d).d designs(n_d).d designs(n_d).k])*[th;b1;b2;b3]; % Control effectiveness matrix
designs(n_d).u_min = f_min*ones(size(designs(n_d).B_f, 2), 1);              % Minimum rotor thrusts (in N);
designs(n_d).u_max = f_max*ones(size(designs(n_d).B_f, 2), 1);              % Maximum rotor thrusts (in N);

% Design 2: HARV Hexarotor PPNNPN -----------------------------------------
n_d = 2;                                                                    % Design ID 
designs(n_d).name = 'HARV Hexarotor PPNNPN';                                % Design name
designs(n_d).m = 3.877;                                                     % Total mass (in kg)
designs(n_d).d = 0.4;                                                       % Rotor arrangement radius (in m)
designs(n_d).k = 0.03;                                                      % Reactive torque to thrust ratio (in m)
f_min = 0;                                                                  % Minimum rotor thrusts (in N);
f_max = 14.5;                                                               % Maximum rotor thrusts (in N);
I_xx = 7.998e-2;                                                            % Moment of inertia around the x-axis (in kg·m²)
I_yy = 7.884e-2;                                                            % Moment of inertia around the y-axis (in kg·m²)
I_zz = 0.148;                                                               % Moment of inertia around the z-axis (in kg·m²)
% Normalized control effectiveness matrix
th =           [1   1    1   1   1   1 ]; % Thrust
b1 = sqrt(3)/2*[0  -1   -1   0   1   1 ]; % Roll
b2 =           [1  0.5 -0.5 -1 -0.5 0.5]; % Pitch
b3 =           [1  1    -1  -1   1  -1 ]; % Yaw

J = diag([-designs(n_d).m, I_xx, I_yy, I_zz]);                              % Inertia matrix
designs(n_d).A = [zeros(4, 4), eye(4); zeros(4, 4), zeros(4, 4)];           % State matrix
designs(n_d).B = [zeros(4, 4); inv(J)];                                     % Control matrix
designs(n_d).B_f = diag([1 designs(n_d).d designs(n_d).d designs(n_d).k])*[th;b1;b2;b3]; % Control effectiveness matrix
designs(n_d).u_min = f_min*ones(size(designs(n_d).B_f, 2), 1);              % Minimum rotor thrusts (in N);
designs(n_d).u_max = f_max*ones(size(designs(n_d).B_f, 2), 1);              % Maximum rotor thrusts (in N);

% Design 3: HARV Quadcopter PNPN ------------------------------------------
n_d = 3;                                                                    % Design ID 
designs(n_d).name = 'HARV Quadcopter PNPN';                                 % Design name
designs(n_d).m = 3.877;                                                     % Total mass (in kg)
designs(n_d).d = 0.4;                                                       % Rotor arrangement radius (in m)
designs(n_d).k = 0.03;                                                      % Reactive torque to thrust ratio (in m)
f_min = 0;                                                                  % Minimum rotor thrusts (in N);
f_max = 14.5;                                                               % Maximum rotor thrusts (in N);
I_xx = 7.998e-2;                                                            % Moment of inertia around the x-axis (in kg·m²)
I_yy = 7.884e-2;                                                            % Moment of inertia around the y-axis (in kg·m²)
I_zz = 0.148;                                                               % Moment of inertia around the z-axis (in kg·m²)
% Normalized control effectiveness matrix
th =           [1  1  1  1 ];   % thrust
b1 = sqrt(2)/2*[1 -1 -1  1 ];   % roll
b2 = sqrt(2)/2*[1  1 -1 -1 ];   % pitch
b3 =           [1 -1  1 -1 ];   % yaw

J = diag([-designs(n_d).m, I_xx, I_yy, I_zz]);                              % Inertia matrix
designs(n_d).A = [zeros(4, 4), eye(4); zeros(4, 4), zeros(4, 4)];           % State matrix
designs(n_d).B = [zeros(4, 4); inv(J)];                                     % Control matrix
designs(n_d).B_f = diag([1 designs(n_d).d designs(n_d).d designs(n_d).k])*[th;b1;b2;b3]; % Control effectiveness matrix
designs(n_d).u_min = f_min*ones(size(designs(n_d).B_f, 2), 1);              % Minimum rotor thrusts (in N);
designs(n_d).u_max = f_max*ones(size(designs(n_d).B_f, 2), 1);              % Maximum rotor thrusts (in N);

% Design 4: HARV hexarotor with 2 coaxial rotors PPNN[PP][NN] -------------
n_d = 4;                                                                    % Design ID 
designs(n_d).name = 'HARV hexarotor+2 PPNN[PP][NN]';                        % Design name
designs(n_d).m = 3.877;                                                     % Total mass (in kg)
designs(n_d).d = 0.4;                                                       % Rotor arrangement radius (in m)
designs(n_d).k = 0.03;                                                      % Reactive torque to thrust ratio (in m)
f_min = 0;                                                                  % Minimum rotor thrusts (in N);
f_max = 14.5;                                                               % Maximum rotor thrusts (in N);
I_xx = 7.998e-2;                                                            % Moment of inertia around the x-axis (in kg·m²)
I_yy = 7.884e-2;                                                            % Moment of inertia around the y-axis (in kg·m²)
I_zz = 0.148;                                                               % Moment of inertia around the z-axis (in kg·m²)
% Normalized control effectiveness matrix
th =           [1   1    1   1   1   1   1   1 ]; % thrust
b1 = sqrt(3)/2*[0  -1   -1   0   1   1   1   1 ]; % roll
b2 =           [1  0.5 -0.5 -1 -0.5 0.5 -0.5 0.5]; % pitch
b3 =           [1  1    -1  -1   1  -1   1  -1 ]; % yaw

J = diag([-designs(n_d).m, I_xx, I_yy, I_zz]);                              % Inertia matrix
designs(n_d).A = [zeros(4, 4), eye(4); zeros(4, 4), zeros(4, 4)];           % State matrix
designs(n_d).B = [zeros(4, 4); inv(J)];                                     % Control matrix
designs(n_d).B_f = diag([1 designs(n_d).d designs(n_d).d designs(n_d).k])*[th;b1;b2;b3]; % Control effectiveness matrix
designs(n_d).u_min = f_min*ones(size(designs(n_d).B_f, 2), 1);              % Minimum rotor thrusts (in N);
designs(n_d).u_max = f_max*ones(size(designs(n_d).B_f, 2), 1);              % Maximum rotor thrusts (in N);

% Design 5: HARV coaxial hexarotor (PNNPPNNPPNNP) -------------------------
n_d = 5;                                                                    % Design ID 
designs(n_d).name = 'HARV coaxial hexarotor PNNPPNNPPNNP';                % Design name
designs(n_d).m = 3.877;                                                     % Total mass (in kg)
designs(n_d).d = 0.4;                                                       % Rotor arrangement radius (in m)
designs(n_d).k = 0.03;                                                      % Reactive torque to thrust ratio (in m)
f_min = 0;                                                                  % Minimum rotor thrusts (in N);
f_max = 14.5;                                                               % Maximum rotor thrusts (in N);
I_xx = 7.998e-2;                                                            % Moment of inertia around the x-axis (in kg·m²)
I_yy = 7.884e-2;                                                            % Moment of inertia around the y-axis (in kg·m²)
I_zz = 0.148;                                                               % Moment of inertia around the z-axis (in kg·m²)
% Normalized control effectiveness matrix
th =           [1  1   1   1    1    1   1  1   1    1    1   1 ]; % thrust
b1 = sqrt(3)/2*[0  0  -1  -1   -1   -1   0  0   1    1    1   1 ]; % roll
b2 =           [1  1  0.5 0.5 -0.5 -0.5 -1 -1 -0.5 -0.5  0.5 0.5]; % pitch
b3 =           [1 -1  -1   1    1   -1  -1  1   1   -1   -1   1 ]; % yaw

J = diag([-designs(n_d).m, I_xx, I_yy, I_zz]);                              % Inertia matrix
designs(n_d).A = [zeros(4, 4), eye(4); zeros(4, 4), zeros(4, 4)];           % State matrix
designs(n_d).B = [zeros(4, 4); inv(J)];                                     % Control matrix
designs(n_d).B_f = diag([1 designs(n_d).d designs(n_d).d designs(n_d).k])*[th;b1;b2;b3]; % Control effectiveness matrix
designs(n_d).u_min = f_min*ones(size(designs(n_d).B_f, 2), 1);              % Minimum rotor thrusts (in N);
designs(n_d).u_max = f_max*ones(size(designs(n_d).B_f, 2), 1);              % Maximum rotor thrusts (in N);

% Design 6: HARV coaxial quadcopter (PNPNNPPN) ----------------------------
n_d = 6;                                                                    % Design ID 
designs(n_d).name = 'HARV coaxial quadcopter PNPNNPPN';                   % Design name
designs(n_d).m = 3.877;                                                     % Total mass (in kg)
designs(n_d).d = 0.4;                                                       % Rotor arrangement radius (in m)
designs(n_d).k = 0.03;                                                      % Reactive torque to thrust ratio (in m)
f_min = 0;                                                                  % Minimum rotor thrusts (in N);
f_max = 14.5;                                                               % Maximum rotor thrusts (in N);
I_xx = 7.998e-2;                                                            % Moment of inertia around the x-axis (in kg·m²)
I_yy = 7.884e-2;                                                            % Moment of inertia around the y-axis (in kg·m²)
I_zz = 0.148;                                                               % Moment of inertia around the z-axis (in kg·m²)
% Normalized control effectiveness matrix
th =           [1  1  1  1  1  1  1   1];   % thrust
b1 = sqrt(2)/2*[1  1 -1 -1 -1 -1  1   1];   % roll
b2 = sqrt(2)/2*[1  1  1  1 -1 -1 -1  -1];   % pitch
b3 =           [1 -1  1 -1 -1  1  1  -1];   % yaw

J = diag([-designs(n_d).m, I_xx, I_yy, I_zz]);                              % Inertia matrix
designs(n_d).A = [zeros(4, 4), eye(4); zeros(4, 4), zeros(4, 4)];           % State matrix
designs(n_d).B = [zeros(4, 4); inv(J)];                                     % Control matrix
designs(n_d).B_f = diag([1 designs(n_d).d designs(n_d).d designs(n_d).k])*[th;b1;b2;b3]; % Control effectiveness matrix
designs(n_d).u_min = f_min*ones(size(designs(n_d).B_f, 2), 1);              % Minimum rotor thrusts (in N);
designs(n_d).u_max = f_max*ones(size(designs(n_d).B_f, 2), 1);              % Maximum rotor thrusts (in N);

% Design 7: Dummy ----------------------------
n_d = 6;                                                                    % Design ID 
designs(n_d).name = 'Dummy';                   % Design name
designs(n_d).m = 3.877;                                                     % Total mass (in kg)
designs(n_d).d = 0.4;                                                       % Rotor arrangement radius (in m)
designs(n_d).k = 0.03;                                                      % Reactive torque to thrust ratio (in m)
f_min = 0;                                                                  % Minimum rotor thrusts (in N);
f_max = 14.5;                                                               % Maximum rotor thrusts (in N);
I_xx = 7.998e-2;                                                            % Moment of inertia around the x-axis (in kg·m²)
I_yy = 7.884e-2;                                                            % Moment of inertia around the y-axis (in kg·m²)
I_zz = 0.148;                                                               % Moment of inertia around the z-axis (in kg·m²)
% Normalized control effectiveness matrix
th =           [1  1  1  1  1  1  1   1];   % thrust
b1 = sqrt(2)/2*[1  1 -1 -1 -1 -1  1   1];   % roll
b2 = sqrt(2)/2*[1  1  1  1 -1 -1 -1  -1];   % pitch
b3 =           [1 -1  1 -1 -1  1  1  -1];   % yaw

J = diag([-designs(n_d).m, I_xx, I_yy, I_zz]);                              % Inertia matrix
designs(n_d).A = [zeros(4, 4), eye(4); zeros(4, 4), zeros(4, 4)];           % State matrix
designs(n_d).B = [zeros(4, 4); inv(J)];                                     % Control matrix
designs(n_d).B_f = diag([1 designs(n_d).d designs(n_d).d designs(n_d).k])*[th;b1;b2;b3]; % Control effectiveness matrix
designs(n_d).u_min = f_min*ones(size(designs(n_d).B_f, 2), 1);              % Minimum rotor thrusts (in N);
designs(n_d).u_max = f_max*ones(size(designs(n_d).B_f, 2), 1);              % Maximum rotor thrusts (in N);