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
CAR_NUMBER = [5 10 15 ]; % CAR_SHARE = 0.75;

% Payoffs.
PAYOFF_BUS = [ 50 70 ];

%% Learning Variables.
%%%%%%%%%%%%%%%%%%%%%%

%% Roth Erev model

INCREASE_SHOCK = 25;

DECREASE_SHOCK = 40;

TIME_INTERVAL_DECREASE = 10;

% Experimentation / Error.
epsilon = [0.2];

% Forgetting (or recency).
phi = [0.001];

% Strength of initial propensities.
S1 = 1;

% Reference point at time 0 (baseline BUS_PAYOFF).
rho1 = 0;

% Weights assigned to positive and negative reinforcement.
% How much new experience is weighted against old. (1 = only new).
wPlus = 0.6;
wMinus = 0.6;

% Positive Constraint.
upsilon = 0.0001;

% Use data from experiment to set initial propensities and probabilities.
INIT_T1 = 1;

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

simName = 'timeright';

save(['conf/' simName]);
