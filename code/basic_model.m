%% ODE of MCM 2021 Problem A
%% Start

%% Natural growth rate
r1 = 2.4;
r2 = 4.8;
r3 = 3.6;

%% Initial quantity of fungi
N10 = 1;
N20 = 1;
N30 = 2;
N0 = [N10; N20; N30];

%% Impact of interact between two species
sigma12 = 1.2;
sigma13 = 0.7;
sigma21 = 0.8;
sigma23 = 1.1;
sigma31 = 1.2;
sigma32 = 0.8;

%% Interval of time
tspan = 1:0.1:100;

%% The derivative of function is not continuous at t=0
tspan2 = tspan(2:length(tspan));

%% Relative humidty

%%%%% Modify the hu to change humidity
hu = 75;
%%%%% Random number
ranNum = 0.1*randn(100);
waveOfHumidity = interp1(1:tspan(length(tspan)), ranNum, tspan, 'spline');
waveOfHumidity = waveOfHumidity(1:length(tspan));
% plot it
plot(tspan, hu * (1 + waveOfHumidity), 'k', 'LineWidth', 2);
yline(60,'--k')
axis([1 tspan(length(tspan)) 0 100])
xlabel('days'); ylabel('RH %');
v = @(t) hu / 50 * waveOfHumidity(length(t));
%% The influence of moisture
%  Described like 1+B*cos(t) here

Wv1 = @(t) 1 + 0.20 * v(t);
Wv2 = @(t) 1 + 0.10 * v(t);
Wv3 = @(t) 1 + 0.40 * v(t);

%% Max quantity that environment can sustain
N1max = @(t) 100;
N2max = @(t) 100;
N3max = @(t) 100;

%% Differential equations set
f = @(t, y)[
        r1 * y(1) * (1 - y(1) / N1max(t) - sigma21 * y(2) / N2max(t) - sigma31 * y(3) / N3max(t)) * Wv1(t)
        r2 * y(2) * (1 - y(2) / N2max(t) - sigma12 * y(1) / N1max(t) - sigma32 * y(3) / N3max(t)) * Wv2(t)
        r3 * y(3) * (1 - y(3) / N3max(t) - sigma13 * y(1) / N1max(t) - sigma23 * y(2) / N2max(t)) * Wv3(t)
        ];

%% Slove the above ODE set
[t, y] = ode45(f, tspan, N0);

%% Show the conclusion with images
figure(2)
% Figure one upper
subplot(2, 1, 1)
plot(t, y(:, 1), t, y(:, 2), t, y(:, 3), 'LineWidth', 2)
title('Number of individuals in the population');
xlabel('days'); ylabel('population'); legend(' population 1', ' population 2', ' population 3')

% Figure one lower
subplot(2, 1, 2)
hold on
% Growth rate
growthRate1 = diff(y(:, 1))/0.1;
growthRate2 = diff(y(:, 2))/0.1;
growthRate3 = diff(y(:, 3))/0.1;
plot(tspan2, growthRate1, tspan2, growthRate2, tspan2, growthRate3, 'LineWidth', 2)
title('Population growth rate')
xlabel('days'); ylabel('growth rate'); legend(' population 1', ' population 2', ' population 3')

% Figure two
figure(3)
% Natural decaying consitent
nDecCon = 30;
% Speed of log-decaying
speedOfDecay = growthRate1 + growthRate2 + growthRate3 + nDecCon;
plot(tspan2, speedOfDecay,'Color', '#A2142F',  'LineWidth', 2)
axis([0 tspan2(length(tspan2)) 0 150])
title('Decaying speed')
xlabel('days'); ylabel('speed')

% Figure 4
figure(5)
hold on
rateOfDecay = cumtrapz(tspan2, speedOfDecay);
plot(tspan2, rateOfDecay, 'Color', '#77AC30', 'LineWidth', 2);
yline(2500, '--k');
title('Decomposition rate')
xlabel('days'); ylabel('Decomposition rate')
disp(rateOfDecay(length(tspan2)))

%% End
