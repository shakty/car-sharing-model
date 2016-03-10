clc
close all
clear

% Open debugger if there is an error.
% dbstop if error

%% Simulation stuff.

nRuns = 10;



%% Variables for the game.
%%%%%%%%%%%%%%%%%%%%%%%%%

% Cars.
CAR_NUMBER = 15; % CAR_SHARE = 0.75;

% Payoffs.
PAYOFF_BUS = 50;

%% Learning Variables.
%%%%%%%%%%%%%%%%%%%%%%

% Car vs Bus variables

% Propensity increase for choosing bus.
INCREASE_BUS = 5;

% Propensity increase for the BUS, when the player does not find a car.
INCREASE_CAR_MISSED = 10;

% Propensity incraese for the CAR, when the player finds a car.
INCREASE_CAR_GOT = 140;

% Time variables.

% Propensity increase of neighboring times when getting a car.
INCREASE_TIME = 20;

% Propensity decrease of neighboring times when missing a car.
DECREASE_TIME = 10;

% How many neighboring times are updated by INCREASE_TIME.
TIME_INTERVAL_INCREASE = 15;

% How many neighboring times are updated by DECREASE_TIME.
TIME_INTERVAL_DECREASE = 10;

% How the interval-increase decays in further times.
INCREASE_DECAY = 0.2;

% How the interval-decrease decays in previous times.
DECREASE_DECAY = 0.1;

% If got car, propensities + INCREASE_SHOCK will be updated.
INCREASE_SHOCK = 5;

% If did not get car, propensities + DECREASE_SHOCK will be updated.
DECREASE_SHOCK = 10;

%% Seed.

seed_random = 0;
seed_machinetime = 1;
seed_fixed = 2;

seed = 0;

seedtype = seed_random;         

% Random seed must be initialized for each batch (level of sigma)
if (seedtype ~= seed_fixed)
    s = RandStream('mcg16807','Seed', seed);
    RandStream.setGlobalStream(s);
else
    batchSeed = seed;
end



%% Save it!
%%%%%%%%%%%

simName = 'b50_c75_rl';

save(['mat/' simName]);
