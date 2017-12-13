%% Water Hammer Simulation - Ryan Haseman
clear; close all

% Simulates the amplitude of the transient pressure wave in a fluid system
% due to sudden valve closure.

%% physical vars
% pipe properties
    d = .0508; % m pipe outer diameter
    wt = .00165; % m wall thickness of pipe
    L = 20; % m pipe length
    E = 195e9; % pa Young's modulus of pipe material (304 ss)
% fluid properties
    rho = 1000; % kg/m^3 density of the fluid (water)
    bm = 2.15e9; % pa bulk modulus (water)
    mu = 1.307e-3; %N*s/m^2 dynamic viscosity @ 10c
% accumulator properties
    b = 1000; % N*sec/m
    k = 3000; % N/m  Spring constant

%% calculated vars
% dimensional properties
    ID = d-2*wt; % m ID of pipe
    Vpipe = (ID/2)^2*pi * L; % m^3 volume of pipe run
    Vac = .001647; % m^3?? based on 4" round by 8" long accumulator
    a = pi*(ID/2)^2; % m^2 cross-sectional area of flow volume
    aAccum = .008107; % m^2  area of accumulator "piston"
% wave speed 
    sos = sqrt(1/(rho*((1/bm)+((ID)/(E*wt))))); % wave speed
% lumped parameter calculated
    I = rho*L/a; % kg/m^4  inertance
    c1 = Vpipe/bm; % capacitance of pressurized pipe
    c2 = Vac/bm; % capacitance of pressurized accumulator (shearer murphy eq)
    c = c1 + c2; % capacitance of damped system
    m = 1;%rho*sos*a; % mass of water (??) not sure if right
    R = 20e5;%128*mu*L/(pi*ID^4); % pa*s/m^5  resistance from poiseuille law
    
%% Time and frequency arrays
ts = .001; % time domain sampling period
t = 0:ts:100; % time array
N = length(t); % number of "samples"
fsample = 1/ts; % frequency sampling rate
fnyq = fsample/2; % nyquist freq
f =-fnyq:fsample/(N-1):fnyq; % frequency array


%% State equation formulation for accumulator system
Q_int = 13e-3; % m^3/s  initial flow (about 200gpm)
u = 400e3*ones(length(t),1); % pa pressure input

x0 = [...
    0;...
    Q_int;...
    0;...
    0;...
    ];

A1 = [...
    0, 1/c, aAccum/c, 0; ...
    -1/I, -R/I, 0, 0; ...
    -aAccum/m, 0, -b/m, -1/m; ...
    0, 0, k, 0 ...
    ];

B1 = [...
    0;...
    1/I;...
    0;...
    0 ...
    ];

C1 = [...
    1,0,0,0;...
    -aAccum,0,-b,-1;...
    0,0,0,1 ...
    ];

D1= [0;0;0];

sys1 = ss(A1, B1, C1, D1);
sys1.InputName = 'Ps';
sys1.OutputName = {'Pc';'fm';'fk'};
P1 = lsim(sys1,u,t,x0);

%% State equation formulation for system w/o accumulator
x0 = [...
    0;...
    Q_int;...
    ];

A2 = [...
    0, 1/c; ...
    -1/I, -R/I; ...
    ];

B2 = [...
    0;...
    1/I;...
    ];

C2 = [1, 0];

D2= [0];

sys2 = ss(A2, B2, C2, D2)
sys2.InputName = 'Ps';
sys2.OutputName = 'Pc1';
P2 = lsim(sys2,u,t,x0);

%% Frequency specturm
z1 = fftshift(fft(P1(:,1))); % fourier transfor of pressure with accumulator
z2 = fftshift(fft(P2)); % fourier transform for pressure wave w/o accumulator

%% Plot simulations 
figure;
    subplot(2,1,1)
        plot(t,P2/1000, 'b') % simulation without accumulator
        hold on
        plot(t, P1(:,1)/1000, 'r') % simulation with accumulator
        grid on
        title('Trainsient pressure surge')
        xlabel('time (s)')
        ylabel('pressure (kPa)')
        legend('without accumulator','with accumulator')
    subplot(2,1,2)
        plot(f,abs(z2));
        hold on
        plot(f,abs(z2));
        grid on
        title('Two-Sided Frequency Spectrum');
        xlabel('Frequency (Hz)');
        ylabel('Magnatude');

% figure;
%     bodeplot(sys1, sys2);
%     grid on
