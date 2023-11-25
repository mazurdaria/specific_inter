clear all
global N_p;
global N_m;
global phi_0; 

global eps_p;
global eps_m; 
global sigma_p; 
global sigma_m;
global lD;
global hat_chi;
global V;
global mu_m_bulk;
global mu_p_bulk;
global c;
charge = 1.66*10^-19;
kB = 1.38*10^-23;
T = 300;
eps = 40*8.85*10^-12;
Na = 6.022*10^23;
c_bar = 10^3 * Na;

lD = (eps*kB*T/2/charge^2/c_bar)^0.5;
c=c_bar*lD^3;
%c=10^3*Na*lD^3;
N_p = 1;
N_m = 1;
phi_0 = 0.075;
mu_p_bulk = log(N_p * phi_0) + 1 - N_p * (log(1 - phi_0 * (N_p + N_m)) + 1);
mu_m_bulk = log(N_m * phi_0) + 1 - N_m * (log(1 - phi_0 * (N_p + N_m)) + 1);
eps_p = 0;
eps_m = 0;
sigma_p = 0.5;
sigma_m = 0.5;
global H;
X = [1, 0.1; 0.1, 0.5];
hat_chi = inv(X);
H = 5*10^-9;

z = linspace(0, H/lD, 70);

V=0.1 * charge / kB / T;

solinit = bvpinit(z, @initial_guess);

sol = bvp4c(@bvpfcn, @bcfun, solinit);


c_p_sol = sol.y(1, :);
c_m_sol = sol.y(3, :);
u_sol = sol.y(5, :);

e_z = charge * (c_p_sol - c_m_sol);


Q = trapz(sol.x, e_z);

% Compute capacitance
C = (Q / V)*eps/lD

figure;

subplot(3,1,1);
plot(sol.x, c_p_sol, '-');
title('c_p vs. z');
xlabel('z');
ylabel('c_p');

subplot(3,1,2);
plot(sol.x, c_m_sol, '-');
title('c_m vs. z');
xlabel('z');
ylabel('c_m');

subplot(3,1,3);
plot(sol.x, u_sol, '-');
title('u (or V) vs. z');
xlabel('z');
ylabel('u (or V)');



function dydx = bvpfcn(z,y)
disp(size(y))
global H;
global N_p N_m phi_0 mu_m_bulk mu_p_bulk;
global hat_chi;
global eps_p;
global eps_m; 
global sigma_p; 
global sigma_m;
%global c;
global lD;
global U;

U = y(5);
c_p = y(1); 
c_m = y(3); 

mu_p = log(N_p * phi_0 * c_p) + 1 - N_p * (log(1 - phi_0 * (N_p * c_p + N_m * c_m)) + 1);
mu_m = log(N_m * phi_0 * c_m) + 1 - N_m * (log(1 - phi_0 * (N_p * c_p + N_m * c_m)) + 1);
v_plus = U + mu_p - mu_p_bulk - eps_p * (exp(-(z^2)/(2*sigma_p^2)) + exp(-((H/lD - z)^2)/(2*sigma_p^2)));
v_minus = -U + mu_m - mu_m_bulk - eps_m * (exp(-(z^2)./(2*sigma_m^2)) + exp(-((H/lD-z)^2)./(2*sigma_m^2)));
disp(size(y));

%dydx = zeros(6,1);

dydx(1) = y(2);
dydx(3) = y(4);
dydx(5) = y(6);

chi_plus_plus = hat_chi(1,1);
chi_plus_minus = hat_chi(1,2);
chi_minus_plus = hat_chi(2,1);
chi_minus_minus = hat_chi(2,2);
dydx(2) = (chi_plus_plus .* v_plus + chi_plus_minus .* v_minus);
dydx(4) = (chi_minus_plus .* v_plus + chi_minus_minus .* v_minus);
dydx(6) = 0.5 .* (y(3) - y(1));

end
function res = bcfun(ya,yb)

global V;
res = zeros(6,1);
res(1) = ya(1);
%- 1e-5; 
res(2) = yb(1); 
res(3) = ya(3); 
res(4) = yb(3); 
res(5) = ya(5) - V; 
res(6) = yb(5) - V; 
end

function Vec_guess = initial_guess(z)
global V;
global H;
global c;
global lD;
g_c_p =c;
%((exp(-z) + exp(-(H/lD-z))) + 1)*c*2;


g_dot_c_p = 0;
%((-exp(-z) + exp(-(H/lD-z))))*c*2;



g_c_m =c;
%((exp(-z) + exp(-(H/lD-z))) + 1)*c*2;

g_dot_c_m = 0;
%((-exp(-z) +exp(-(H/lD-z))))*c*2;



g_u = V * exp(-z) + V * exp(-(H/lD-z));
g_dot_u = -V*exp(-z)+V * exp(-(H/lD-z));

Vec_guess = [g_c_p; g_dot_c_p; g_c_m; g_dot_c_m; g_u; g_dot_u];
end
