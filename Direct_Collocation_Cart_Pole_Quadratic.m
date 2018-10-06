function Direct_Collocation_Cart_Pole_Quadratic()
% This function is used to use the direct collocation method to address the
% projectile problem. Instead of conduct the numerical integration, here we
% use the direct collocation method to address the numerical integration
% issue,

% The system to be considered is the Cart-Pole system, where the whole
% system has two degress of freedom.

p.l = 1;
p.m1 = 1;
p.m2 = 1;
p.dmax = 2;
p.umax = 20;
% The variables to be optimized are the time, state and control
p.d = 1;
p.T = 2;
N = 40;
p.N = N;

q1_init = 0;
q1dot_init = 0;
q2_init = 0;
q2dot_init = 0;

p.q1_N = p.d;
p.q1dot_N = 0;
p.q2_N = pi;
p.q2dot_N = 0;

q1_array = linspace(q1_init,p.q1_N, p.N);
q1dot_array = linspace(q1dot_init, p.q1dot_N, p.N);
q2_array = linspace(q2_init, p.q2_N, p.N);
q2dot_array = linspace(q2dot_init, p.q2dot_N, p.N);

p.q1_init = q1_init;
p.q1dot_init = q1dot_init;
p.q2_init = q2_init;
p.q2dot_init = q2dot_init;


u_array = zeros(p.N, 1);

x0 = [q1_array, q2_array, q1dot_array, q2dot_array, u_array']';

% This x0 vector is written in a way that the first N is q1 and the last N
% is the u

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @Cart_Pole_Constraint;

x = fmincon(@Cart_Pole_Obj,x0,A,b,Aeq,beq,lb,ub,nonlcon, options, p);


m1 = 1.0; config.dyn.m1 = m1;   %cart mass
m2 = 1.0; config.dyn.m2 = m2;   %pendulum mass
config.dyn.g = 9.81;
config.dyn.l = 1;


tspan = linspace(0, p.T, N);

% Interpolater(tspan, x, 5, p)

% u_array = x(end-N+1:end,:);
% x = x(1:end-N,:);
% x = reshape(x, N, 4);
% q1_array = x(:,1);
% q2_array = x(:,2);
% q1dot_array = x(:,3);
% q2dot_array = x(:,4);

[x_all, xdot_all, u_all, xdot_all_ref, tspan_ref] = Interpolater(tspan, x, 5, p);
q1_array = x_all(1,:);
q2_array = x_all(2,:);
q1dot_array = x_all(3,:);
q2dot_array = x_all(4,:);

for i = 1:length(u_all)
    
    t = tspan_ref(i);
    
    z = [q1_array(i), q2_array(i), q1dot_array(i), q2dot_array(i)]';
    
    drawCartPole(t,z,config.dyn)
    
    pause(0.01)
    
end

%% Dynamics constraint validation
% This part is used to plot the dynamics discrepancy
figure
plot(tspan_ref, q1_array,'LineWidth',1.5);
hold on
plot(tspan_ref, q2_array,'LineWidth',1.5);
hold on
plot(tspan_ref, q1dot_array,'LineWidth',1.5);
hold on
plot(tspan_ref, q2dot_array,'LineWidth',1.5);

legend('q1', 'q2','q1dot','q2dot')


figure
plot(tspan_ref, xdot_all(1,:) - xdot_all_ref(1,:),'LineWidth',1.5);
hold on
plot(tspan_ref, xdot_all(2,:) - xdot_all_ref(2,:),'LineWidth',1.5);
hold on
plot(tspan_ref, xdot_all(3,:) - xdot_all_ref(3,:),'LineWidth',1.5);
hold on
plot(tspan_ref, xdot_all(4,:) - xdot_all_ref(4,:),'LineWidth',1.5);
legend('q1dot', 'q2dot','q1ddot','q2ddot')




end

function [x_all, xdot_all, u_all, xdot_all_ref, tspan_ref] = Interpolater(tspan, x, N, p)

% Here the inputs are
% Time span
% optimized variables
% Number of Points to be inserted between following grids

Grids_No = length(tspan);
Delta_t = tspan(2) - tspan(1);
delta_t = Delta_t/(N+1);

u_array = x(end-Grids_No+1:end,:);
x = x(1:end-Grids_No,:);
x = reshape(x, Grids_No, 4);
q1_array = x(:,1);
q2_array = x(:,2);
q1dot_array = x(:,3);
q2dot_array = x(:,4);

x_all = [];
xdot_all = [];
u_all = [];

for i = 1:Grids_No-1
    % Here is to retrieve the spline function during this interval which
    % remains to be a known function
    
    x_k = [q1_array(i), q2_array(i), q1dot_array(i), q2dot_array(i)]';
    x_kp1 = [q1_array(i+1), q2_array(i+1), q1dot_array(i+1), q2dot_array(i+1)]';
    u_k = u_array(i);
    u_kp1 = u_array(i+1);
    
    qddot_k = Cart_Pole_rhs(x_k, u_k,  p);
    f_k = qddot_k;
    
    qddot_kp1 = Cart_Pole_rhs(x_kp1, u_kp1,  p);
    f_kp1 = qddot_kp1;
    
    if i == Grids_No-1
        Inner_Grids_No = N+2;
    else
        Inner_Grids_No = N+1;
    end
    for j = 1:Inner_Grids_No
        t = (j-1) * delta_t;
        [x_t, xdot_t]= Quadratic_Spline_Eval(x_k, f_k, f_kp1, Delta_t, t);
        x_all = [x_all x_t];
        xdot_all = [xdot_all xdot_t];
        u_t = Linear_Spline_Eval(u_k, u_kp1, Delta_t, t);
        u_all = [u_all u_t];
    end
end

%% The last one is the real dynamics f(t)

Total_State_No = length(u_all);
xdot_all_ref = [];


for i = 1:Total_State_No
    
    x_k = x_all(:,i);
    u_k = u_all(:,i);
    
    x_kp1 = Cart_Pole_rhs(x_k, u_k,  p);
    xdot_all_ref = [xdot_all_ref x_kp1];
end


tspan_ref = linspace(tspan(1), tspan(end), length(u_all));

end

function u_t = Linear_Spline_Eval(u_k, u_kp1, Delta_t, t)

u_t = u_k + t/Delta_t * (u_kp1 - u_k);

end

function [x_t, xdot_t]= Quadratic_Spline_Eval(x_k, f_k, f_kp1, Delta_t, t)

% This function is used to evaluate the value given the Quadratic Spline
% for state

x_t = x_k + t * f_k - t * t/(2 * Delta_t) * (f_k - f_kp1);

xdot_t = f_k - t/Delta_t * (f_k - f_kp1);

end

function  Obj = Cart_Pole_Obj(x, p)

u_array = x(end-p.N+1:end,:);

Obj = 0;

for i = 1:p.N-1
    
    Obj = Obj + u_array(i) * u_array(i) + u_array(i+1) * u_array(i+1);
end

% Obj = dot(u_array,u_array);

end

function [c,ceq] = Cart_Pole_Constraint(x, p)

% The default way to write down the system dynamics is that c(x)<=0

% The first thing to do is to unzip the decision variables into different
% categories
c = [];
ceq = [];

N = p.N;
h_k = p.T/(N-1);

u_array = x(end-N+1:end,:);
x = x(1:end-N,:);
x = reshape(x, N, 4);
q1_array = x(:,1);
q2_array = x(:,2);
q1dot_array = x(:,3);
q2dot_array = x(:,4);

% Now it is time to impose these constraints

%% 1. Dynamics integration through the grid points
for i = 1:N-1
    % Integration from k to k + 1
    
    q_k = [q1_array(i), q2_array(i), q1dot_array(i), q2dot_array(i)]';
    q_kp1 = [q1_array(i+1), q2_array(i+1), q1dot_array(i+1), q2dot_array(i+1)]';
    u_k = u_array(i);
    u_kp1 = u_array(i+1);
    
    f_k = Cart_Pole_rhs(q_k, u_k,  p);
    f_kp1 = Cart_Pole_rhs(q_kp1, u_kp1,  p);
    
    % System state update
    Dynamics_Constraint = q_kp1 - q_k - 1/2 * h_k * (f_k + f_kp1);
    
    ceq = [ceq; Dynamics_Constraint];
    
end

%% 2. Path constraint (Inequality constraint)
for i = 1:N
    % Integration from k to k + 1
    
    q1_k = q1_array(i);
    u_k = u_array(i);
    
    c = [c; q1_k - p.dmax];
    c = [c; -q1_k - p.dmax];
    c = [c; u_k - p.umax];
    c = [c; -u_k - p.umax];
end

%% 3. Boundary condition
% Initial state
ceq = [ceq; q1_array(1) - p.q1_init];
ceq = [ceq; q1dot_array(1) - p.q1dot_init];
ceq = [ceq; q2_array(1) - p.q2_init];
ceq = [ceq; q2dot_array(1) - p.q2dot_init];

% Final state
ceq = [ceq; q1_array(N) - p.q1_N];
ceq = [ceq; q1dot_array(N) - p.q1dot_N];
ceq = [ceq; q2_array(N) - p.q2_N];
ceq = [ceq; q2dot_array(N) - p.q2dot_N];


end

function zdot = Cart_Pole_rhs(z, u,  p)
l = p.l;
m1 = p.m1;
m2 = p.m2;
g = 9.81;
q1 = z(1);              q2 = z(2);
q1dot = z(3);           q2dot = z(4);

q1ddot_num = l * m2 * sin(q2) * q2dot^2 + u + m2 * g * cos(q2) * sin(q2);
q1ddot_den = m1 + m2 * (1 - cos(q2)^2);
q1ddot = q1ddot_num/q1ddot_den;

q2ddot_num = l * m2 * cos(q2) * sin(q2) * q2dot^2 + u * cos(q2) + (m1 + m2 ) * g * sin(q2);
q2ddot_den = l * m1 + l * m2 * (1 - cos(q2)^2);
q2ddot = -q2ddot_num/q2ddot_den;

zdot = [q1dot, q2dot, q1ddot, q2ddot]';

end


function drawCartPole(t,z,dyn)
%
% This function draws the cart pole's current state


cartWidth = 0.3*dyn.l;
cartHeight = 0.2*dyn.l;

extents = [-3.5,1,-cartHeight - 1.2*dyn.l, cartHeight + 1.2*dyn.l];

cartColor = [33,87,177]/256;
poleColor = [183,27,73]/256;
groundColor = [84,64,64]/256;

x = z(1,:);
q = z(2,:);
dx = z(3,:);
dq = z(4,:);

p = autoGen_cartPoleKinematics(x,q,dx,dq,dyn.l);

p1 = [x;0];   %Center of cart
p2 = p;   % Tip of pendulum

clf; hold on;
axis equal;
% axis(extents);

% Draw the cart:
h = rectangle(...
    'Position',[(p1'-0.5*[cartWidth,cartHeight]),[cartWidth,cartHeight]],...
    'Curvature',0.2*[1,1],...
    'EdgeColor',cartColor,...
    'LineWidth',1,...
    'FaceColor',cartColor);

plot(extents(1:2),-0.55*cartHeight*[1,1],'LineWidth',4,'Color', groundColor);

pos = [p1,p2];
plot(pos(1,:),pos(2,:),'Color',poleColor,'LineWidth',3);
plot(p2(1),p2(2),'.','Color',poleColor','MarkerSize',80);

title(['Cart-Pole,  t = ', num2str(t,2)]);

end

function [p,dp] = autoGen_cartPoleKinematics(x,q,dx,dq,l)
%AUTOGEN_CARTPOLEKINEMATICS
%    [P,DP] = AUTOGEN_CARTPOLEKINEMATICS(X,Q,DX,DQ,L)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    14-Jun-2015 09:12:09

t2 = cos(q);
t3 = sin(q);
p = [x+l.*t3;-l.*t2];
if nargout > 1
    dp = [dx+dq.*l.*t2;dq.*l.*t3];
end
end


