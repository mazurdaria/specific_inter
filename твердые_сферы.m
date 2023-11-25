clear all;
global hat_chi;
global V;
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

% Инициализация физических констант
charge=1.66*10^-19; % Заряд электрона
kB=1.38*10^-23; % Постоянная Больцмана
T=300; % Температура в кельвинах
eps=40*8.85*10^-12; % Диэлектрическая проницаемость среды
Na=6.022*10^23; % Число Авогадро 
c_bar=10^3*Na; % Концентрация в объемной фазе
lD=(eps*kB*T/2/charge^2/c_bar)^0.5; % Дебаевская длина

% Инициализация параметров частиц
d_p=3.29*2*10^-10/lD; % Диаметр катионов
d_m=3.32*2*10^-10/lD; % Диаметр анионов
c_b=c_bar*lD^3; % Обезразмеренная концентрация в объемной фазе
% Параметры взаимодействия со стенкой
eps_p=0;
eps_m=0; 
sigma_p=0.5;
sigma_m=0.5; 

% Расчет химического потенциала в объемной фазе
n0=2*c_b;
n1=0.5*(d_p+d_m)*c_b;
n2=pi*c_b*(d_p^2+d_m^2);
n3=pi/6*c_b*(d_p^3+d_m^3);
P=n0*(1+n3+n3^3)/(1-n3)^3+(1/12*pi*n2^3+n1*n2*(1-n3)-3*n0*n3)/(1-n3)^3;


p_bulk=n2*d_p/2/(1-n3)+pi*d_p^2*(n1/(1-n3)+n2^2/(8*pi*(1-n3)^2))+pi*d_p^3*P/6-log(1-n3)+log(c_b);
m_bulk=n2*d_m/2/(1-n3)+pi*d_m^2*(n1/(1-n3)+n2^2/(8*pi*(1-n3)^2))+pi*d_m^3*P/6-log(1-n3)+log(c_b);

% Матрица структурных взаимодействий
X=[10,1;1,10]; 
hat_chi=inv(X); % Обратная матрица взаимодействий

H=10*10^-9; % Ширина поры
H_t=H/lD; % Обезразмеренная ширина
z=linspace(0,H_t,100); % Интервал значений, на котором решаются уравнения

V=-0.1*charge/kB/T; % Потенциал

% Инициализация и решение граничной задачи
solinit=bvpinit(z,@initial_guess);
sol=bvp4c(@bvpfcn,@bcfun,solinit);

% Извлечение решений
c_p_sol = sol.y(1, :); % Концентрация катионов
c_m_sol = sol.y(3, :); % Концентрация анионов
u_sol = sol.y(5, :); % Электрический потенциал
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


    n0=c_p+c_m;
    n1=0.5*(c_p*d_p+c_m*d_m);
    n2=pi*(c_p*d_p^2+c_m*d_m^2);
    n3=pi/6*(c_p*d_p^3+c_m*d_m^3);

    P=n0*(1+n3+n3^3)/(1-n3)^3+(1/12*pi*n2^3+n1*n2*(1-n3)-3*n0*n3)/(1-n3)^3;

    % Расчеты химических потенциалов
    mu_p=n2*d_p/2/(1-n3)+pi*d_p^2*(n1/(1-n3)+n2^2/(8*pi*(1-n3)^2))+pi*d_p^3*P/6-log(1-n3)+log(c_p);
    mu_m=n2*d_m/2/(1-n3)+pi*d_m^2*(n1/(1-n3)+n2^2/(8*pi*(1-n3)^2))+pi*d_m^3*P/6-log(1-n3)+log(c_m);

  
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
