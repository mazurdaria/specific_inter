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

% Инициализация физических констант

charge = 1.66*10^-19; % Заряд электрона
kB = 1.38*10^-23; % Постоянная Больцмана
T = 300; % Температура в кельвинах
eps = 40*8.85*10^-12; % Диэлектрическая проницаемость среды
Na = 6.022*10^23; 

c_bar = 10^3 * Na; % концентрация в объемной фазе 

% Расчет дебаевской длины
lD = (eps*kB*T/2/charge^2/c_bar)^0.5;

c = c_bar*lD^3; % безразмерная концентрация в объемной фазе 

% Инициализация параметров системы
% Размеры ионов
N_p = 1;
N_m = 1; 
phi_0 = 0.075; %объемная доля
% Химический потенциал в объемной фазе
mu_p_bulk = log(N_p * phi_0) + 1 - N_p * (log(1 - phi_0 * (N_p + N_m)) + 1);
mu_m_bulk = log(N_m * phi_0) + 1 - N_m * (log(1 - phi_0 * (N_p + N_m)) + 1);

% Параметры взаимодействия со стенкой
eps_p = 0; 
eps_m = 0;
sigma_p = 0.5; 
sigma_m = 0.5; 

global H;
X = [1, 0.1; 0.1, 0.5]; % Матрица структурных взаимодействий
hat_chi = inv(X); 


H = 5*10^-9; % Ширина поры

% Расчет и визуализация решений
z = linspace(0, H/lD, 70); % Интервал координаты z
V = 0.1 * charge / kB / T; % Потенциал

% Инициализация и решение граничной задачи
solinit = bvpinit(z, @initial_guess);
sol = bvp4c(@bvpfcn, @bcfun, solinit);

% Извлечение решений
c_p_sol = sol.y(1, :); % Концентрация катионов
c_m_sol = sol.y(3, :); % Концентрация анионов
u_sol = sol.y(5, :); % Электрический потенциал


% Визуализация результатов
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

% Функции для решения граничной задачи (bvpfcn, bcfun, initial_guess)
%...


% Функция для решения граничной задачи
function dydx=bvpfcn(z,y)
    global hat_chi;
    global m_bulk;
    global p_bulk;
    global d_p;
    global d_m;
    global c_b;
    global lD;
    global H;
    global eps_p;
    global eps_m;
    global sigma_p;
    global sigma_m;

    U=y(5);
    c_p=y(1);
    c_m=y(3);

    % Расчеты химических потенциалов
    
    mu_p = log(N_p * phi_0 * c_p) + 1 - N_p * (log(1 - phi_0 * (N_p * c_p + N_m * c_m)) + 1);
    mu_m = log(N_m * phi_0 * c_m) + 1 - N_m * (log(1 - phi_0 * (N_p * c_p + N_m * c_m)) + 1);

  
    v_p=U+mu_p-p_bulk-eps_p*(exp(-(z^2)/(2*sigma_p^2))+exp(-((H/lD-z)^2)/(2*sigma_p^2)));
    v_m=-U+mu_m-m_bulk-eps_m*(exp(-(z^2)/(2*sigma_m^2))+exp(-((H/lD-z)^2)/(2*sigma_m^2)));

    % Коэффициенты матрицы структурных взаимодействий
    chi_pp=hat_chi(1,1);
    chi_pm=hat_chi(1,2);
    chi_mp=hat_chi(2,1);
    chi_mm=hat_chi(2,2);

    % Формирование системы уравнений
    dydx=zeros(6,1);

    dydx(1)=y(2);
    dydx(3)=y(4);
    dydx(5)=y(6);

    dydx(2)=(chi_pp*v_p+chi_pm*v_m)*c_b;
    dydx(4)=(chi_mp*v_p+chi_mm*v_m)*c_b;
    dydx(6)=0.5*(y(3)-y(1))/c_b;
end

% Функция для установления граничных условий
function res = bcfun(ya, yb, varargin)
    global V;

    % Установка граничных условий для концентраций и потенциала
    res(1) = ya(1);
    res(2) = yb(1);
    res(3) = ya(3);
    res(4) = yb(3);
    res(5) = ya(5) - V; 
    res(6) = yb(5) - V;
end

% Начальное приближение для решения граничной задачи
function v_guess=initial_guess(z)
    global V;
    global c_b;

    % Начальные значения для концентраций и потенциала
    g_c_p=c_b;
    g_dot_c_p=0;

    g_c_m=c_b;
    g_dot_c_m=0;

    g_u=V;
    g_dot_u=0;

    % Возвращение вектора начального приближения
    v_guess=[g_c_p;g_dot_c_p;g_c_m;g_dot_c_m;g_u;g_dot_u];
end
