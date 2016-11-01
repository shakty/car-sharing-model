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

% Experimentation / Error.
epsilon = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

% Forgetting (or recency).
phi = [0.001];

% Strength of initial propensities.
S1 = 1;

% Reference point at time 0.
rho1 = [ 0 ];

% Make rho1 a parameter relative to BUS_PAYOFF (added to it).
rho1_relative_to_bus = 0;

% Weights assigned to positive and negative reinforcement.
% How much new experience is weighted against old. (1 = only new).
wPlus = 1; % When my reward exceeds expectation.
wMinus = 1; % When my reward is below expectation.

% Positive Constraint.
upsilon = 0.0001;

%% Adapting Roth Erev to Car-Sharing

% The increase in target departure time if a car is gotten.
INCREASE_SHOCK = [20];

% The decrease in target departure time if the car is missed.
DECREASE_SHOCK = [20];

% The reward for having got a car if chosen car.
REWARD_GOT_CAR = [0 40];

% The propensities of departure times within
% this interval (both + and -) are updated.
TIME_INTERVAL = 10;

%% Init

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

simName = 'new-deal';

save(['conf/' simName]);

fprintf('File saved: %s\n', ['conf/' simName '.mat']);