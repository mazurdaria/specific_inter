clear all
global N_p;
global N_m;
global phi_0; 

global eps_p;
global eps_m; 
global sigma_p; 
global sigma_m;

global hat_chi;
global V;
global mu_m_bulk;
global mu_p_bulk;
global d_plus;
global d_minus;
global c;
global lD;
charge = 1.66*10^-19;
kB = 1.38*10^-23;
T = 300;
eps = 78*8.85*10^-12;
Na = 6.022*10^23;
c_bar = 10^3 * Na;

lD = (eps*kB*T/2/charge^2/c_bar)^0.5;
d_plus = 3.29*2*10^-10/lD;
d_minus = 3.29*2*10^-10/lD;
c=1*c_bar*lD^3;
n0 = c + c;
n1 = 0.5 * d_plus * c + 0.5 * d_minus * c;
n2 = pi * d_plus.^2 * c + pi * d_minus.^2 * c;
n3 = pi / 6 * (d_plus.^3 * c + d_minus.^3 * c);
P = n0 * (1 + n3 + n3.^3) / (1 - n3).^3 + (1/12 * pi * n2.^3 + n1 .* n2 * (1 - n3) - 3 * n0 * n3) / (1 - n3).^3;
mu_p_bulk = log(c/c) + n2.*d_plus./2./(1-n3) + pi.*d_plus.^2.*(n1./(1-n3) + n2.^2./(8*pi.*(1-n3).^2)) + pi.*d_plus.^3.*P./6-log(abs(1-n3));
mu_m_bulk = log(c/c) + n2.*d_minus./2./(1-n3) + pi.*d_minus.^2.*(n1./(1-n3) + n2.^2./(8*pi.*(1-n3).^2)) + pi.*d_minus.^3.*P./6-log(abs(1-n3));
eps_p = 0;
eps_m = 0;
sigma_p = 0.5;
sigma_m = 0.5;
global H;
%X=[30,15;15,30];
X=[10,1;1,10];
%X=[60,30;30,60];
hat_chi = inv(X);
%hat_chi = [0.1, -0.09; -0.09, 0.001]; для первого графика

H=10*10^-9;
H_t=H/lD;
z = linspace(0, H_t, 100);

%De=??
%a=??
%re?=?
epsilon=0.001;
V = 0.1 * charge / kB / T;

solinit = bvpinit(z, @initial_guess);

sol = bvp4c(@bvpfcn, @bcfun, solinit);

sigma_initial = -deval(sol, z_eval)(6);

delta_V = 1e-3; 
V = V_initial + delta_V;
sol_perturbed = bvp4c(@bvpfcn, @bcfun, solinit);


sigma_perturbed = -deval(sol_perturbed, z_eval)(6);


C = (sigma_perturbed - sigma_initial) / delta_V;


disp(['The differential capacitance C at z = 0 is: ', num2str(C)]);


function dydx = bvpfcn(z,y)
global d_plus;
global d_minus;
global H;
global mu_m_bulk mu_p_bulk U;
global hat_chi;
global eps_p;
global eps_m; 
global sigma_p; 
global sigma_m;
global c;
global lD;
U = y(5);
c_p = y(1); 
c_m = y(3)
n0 = c_p + c_m;
n1 = 0.5 * d_plus * c_p + 0.5 * d_minus * c_m;
n2 = pi * d_plus.^2 * c_p + pi * d_minus.^2 * c_m;
n3 = pi / 6 * (d_plus.^3 * c_p + d_minus.^3 * c_m);

P = n0 * (1 + n3 + n3.^3) / (1 - n3).^3 + (1/12 * pi * n2.^3 + n1 .* n2 * (1 - n3) - 3 * n0 * n3) / (1 - n3).^3;
mu_p = log(c_p/c)*0 + n2.*d_plus./2./(1-n3) + pi.*d_plus.^2.*(n1./(1-n3) + n2.^2./(8*pi.*(1-n3).^2)) + pi.*d_plus.^3.*P./6-log(1-n3);
mu_m = log(c_m/c)*0 + n2.*d_minus./2./(1-n3) + pi.*d_minus.^2.*(n1./(1-n3) + n2.^2./(8*pi.*(1-n3).^2)) + pi.*d_minus.^3.*P./6-log(1-n3);

v_plus = U + mu_p - mu_p_bulk - eps_p * (exp(-(z^2)/(2*sigma_p^2)) + exp(-((H/lD - z)^2)/(2*sigma_p^2)));
v_minus = -U + mu_m - mu_m_bulk - eps_m * (exp(-(z^2)/(2*sigma_m^2)) + exp(-((H/lD-z)^2)/(2*sigma_m^2)));
%V_morse = De * (1 - exp(-a*(z-re))).^2 - De;
%if any(c_m <=0)
    %disp('C_m<=0')

%end

dydx = zeros(6,1);

dydx(1) = y(2);
dydx(3) = y(4);
dydx(5) = y(6);
chi_plus_plus = hat_chi(1,1);
chi_plus_minus = hat_chi(1,2);
chi_minus_plus = hat_chi(2,1);
chi_minus_minus = hat_chi(2,2);
dydx(2) = (chi_plus_plus * v_plus + chi_plus_minus * v_minus)*c;
dydx(4) = (chi_minus_plus * v_plus + chi_minus_minus * v_minus)*c;

dydx(6) = 0.5 * (y(3) - y(1))/c;

end
function res = bcfun(ya,yb)
global U;
global V;
%res = zeros(6,1);
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


g_dot_c_p =0;
%((-exp(-z) + exp(-(H/lD-z))))*c*2;



g_c_m = c;
%((exp(-z) + exp(-(H/lD-z))) + 1)*c*2;

g_dot_c_m =0;
%((-exp(-z) +exp(-(H/lD-z))))*c*2;



g_u =V * exp(-z) + V * exp(-(H/lD-z));
g_dot_u = -V*exp(-z)+V * exp(-(H/lD-z));

Vec_guess = [g_c_p; g_dot_c_p; g_c_m; g_dot_c_m; g_u; g_dot_u];
end
