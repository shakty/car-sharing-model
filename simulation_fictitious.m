function [ output_args ] = simulation( args )
%SIMULATION Summary of this function goes here
%   Detailed explanation goes here


%% Rules of the game.
%%%%%%%%%%%%%%%%%%%%%

% Nr Players.
N = 20;

% Timesteps.
T = 30;

% Strategies indexes.
BUS = 1;
CAR = 2;

% Available strategies.
avail_strategies_carbus = [ BUS CAR ];

% Nr. available strategies each player each round.
nr_strategies = length(avail_strategies_carbus);

% Avail strategies for car (leaving from 10 to 11).
avail_strategies_time = 0:60;

% Nr available time strategies.
nr_strategies_time = length(avail_strategies_time);
% seed = args.seed;

% Payoffs.
PAYOFF_CAR = 30;
SLOPE_CAR = (100 - PAYOFF_CAR) / 60;
SLOPE_CAR_MISS = PAYOFF_CAR / 60;

%% Input Parameters.
%%%%%%%%%%%%%%%%%%%%

simName = args.simName;
simCount = args.simCount;

DUMP = args.DUMP;
dumpDir = args.dumpDir;

DEBUG = args.DEBUG;
INIT_T1 = args.INIT_T1;

% if (DUMP)
%     OUT_DIR = [ args.dumpDir int2str(simCount) '/' ];
%     mkdir(OUT_DIR);
% end

% Game variables.

PAYOFF_BUS = args.PAYOFF_BUS;
CAR_NUMBER = args.CAR_NUMBER;

if (CAR_NUMBER == 15)
    car_level = 75;
elseif (CAR_NUMBER == 10)
    car_level = 50;
else
    car_level = 25;
end


% Erev Roth.

S1 = args.S1;
epsilon = args.epsilon;
phi = args.phi;
wPlus = args.wPlus;
wMinus = args.wMinus;
upsilon = args.upsilon;

% Rho1, bus relative.
if (args.rho1_relative_to_bus)    
    rho1 = PAYOFF_BUS + args.rho1;
else
    rho1 = args.rho1;
end

% Rho1, absolute.
%rho1 = args.rho1;

TIME_INTERVAL = args.TIME_INTERVAL;

REWARD_GOT_CAR = args.REWARD_GOT_CAR;
INCREASE_SHOCK = args.INCREASE_SHOCK;
DECREASE_SHOCK = args.DECREASE_SHOCK;

REPETITIONS = args.nRuns;

%% Fit params.
FIT = args.FIT;
if (FIT)
    FIT_MEAN_COL = 5;
    FIT_SD_COL = 6;
    suffix = [ num2str(PAYOFF_BUS), '_', num2str(car_level), '.csv' ];
    
    fitBus = csvread(['fit/summary_bus_round_', suffix ]);    
    fitDepTime = csvread(['fit/summary_deptime_round_', suffix ]);
    % need to skip first row because it contains non-numbers (NA).
    fitSwitch = csvread(['fit/summary_switch_round_', suffix ], 1, 0);
    fitSwitch = [ 50 25 1 0 0 0 0 0 0 ; fitSwitch ];
end

%% Data structures for all repetitions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rep_nCars = zeros(REPETITIONS, 1);
% rep_avgCarTime = zeros(REPETITIONS, 1);
% rep_depCarTime = zeros(REPETITIONS, N);
% rep_propensities_carbus = zeros(N, 2, REPETITIONS);
% rep_propensities_time = zeros(N, nr_strategies_time, REPETITIONS);

for r = 1 : REPETITIONS


% Clear simulation data.
    
% Init equal for all strategies.

% T+1 so that we do not worry about the last iteration.
timeTargets = zeros((T+1), N);
expected_payoffs_car = zeros((T+1), N);

timeTargets(1,:) = randi(nr_strategies_time, N, 1)';
expected_payoffs_car(1,:) = PAYOFF_CAR + (SLOPE_CAR * timeTargets(1,:)');

% Set initial propensities and probabilities based on experimental data.
if (INIT_T1) 

    if (PAYOFF_BUS == 50)
        % Probability of taking the bus at round 1.
        PBUS = 0.2078721;
        
        % Avg. departure time and standard deviation at round 1.
        % TCAR = 38.77484;
        % TCAR_SD = 20.68923;
        
        % Distribution of departure times from experiment (counts).
        props = [ ...
            85   1   3   0   0   2   1   0   3   2   6   1   0   1   0   8   2   1   2   0 ...
            10   2   0   2   2   6   2   2   7   7  84   1   1   5   2   8   3   0   1 ...
            8  32   4   4   2   7  56   3   2   1   1  28   0   5   1   3  12   2   4 ...
            7   4 195 ...
        ];
    
        % time_distr = [  1 1 1 1 31 31 31 32 41 46 46 51 56 59 61 61 61 61 61 61 ];
        
    else
        % Probability of taking the bus at round 1.
        PBUS = 0.3522868;
        
        % Avg. departure time and standard deviation at round 1.
        % TCAR = 37.57252;
        % TCAR_SD = 21.82269;
        
        % Distribution of departure times from experiment (counts).
        props = [ ...
             95  1   0   0   2   0   1   0   1   0   3   0   0   1   0   6   0   1   0   2 ...
             5   0   1   2   1   2   1   1   1   2  67   9   6   1   1   2   4   1   3 ...
             2  17   1   2   3   4  51   4   4   3   5  21   3   5   0   1   9   1   4 ...
             7   3 151 ...
        ];
    
        
        % time_distr = [ 1 1 1 21 31 31 31 41 46 46 51 56 61 61 61 61 61 61 61 ];
        % time_distr_last3 = [ 16 36 40 ];
        % time_distr = sort([ time_distr time_distr_last3(randi(3)) ]);
    
    end
    
    increase = referencePoints(1) * (1 - epsilon);
    if (increase == 0)
        increase = 10;
    end
    tmp = (rand(N,1) > PBUS) .* increase;
    propensities_carbus(:, CAR) = propensities_carbus(:, CAR).*tmp;
    [tmp, ~] = find(propensities_carbus(:, CAR) == 0);
    propensities_carbus(tmp,CAR) = S1;
    propensities_carbus(tmp,BUS) = propensities_carbus(tmp,BUS).*increase;
    
    % This snippet initializes ALL probabilities of EACH individual based
    % on the GLOBAL distribution of all the population.
        
    % tmp = normpdf(1:nr_strategies_time, TCAR, TCAR_SD).*increase.*100;
    
    propensities_time = repmat(props, N, 1);    
    for i=1:N
        jitter = 1.1 + (0.9-1.1).*rand(nr_strategies_time, 1);
        
        propensities_time(i,:) = props .* jitter';
        % propensities_time(i,:) = propensities_time(i,:) .* tmp;
        
        sumPropensities = sum(propensities_carbus(i, :));
        probabilities(i, BUS) = propensities_carbus(i, BUS) / ...
            sumPropensities;
        probabilities(i, CAR) = propensities_carbus(i, CAR) / ...
            sumPropensities;
        
        sumPropensities = sum(propensities_time(i, :));
        
        for it = 1 : nr_strategies_time
            probabilities_time(i, it) = propensities_time(i, it) / ...
                sumPropensities;
        end        
    end    
    
    
         
end


% Store strategies, payoffs and choices over all rounds.
payoffs = zeros(N, T);
gotCars = zeros(N, T);
strategies_carbus = ones(N, T);
strategies_time = ones(N, T);

if (DEBUG)
    shareBus = zeros(T, 1);
    meanDep = zeros(T, 1);
end
    
if (FIT)
    mseBus = zeros(T, 1);
    mseTime = zeros(T, 1);
    mseSwitch = zeros(T, 1);
end

%% One simulation.
%%%%%%%%%%%%%%%%%%

for t = 1 : T

    if (FIT && t > 1)
        strategy_switches = zeros(N,1);
    end
        
    for player = 1 : N
        
        % Expected payoff of car based on previous history of play.
        expCar = expected_payoffs_car(t, player);
        
        err = (randn)*(T+1-t)
        expCarErr = expCar + err;
        
        if (PAYOFF_BUS > expCar) 
            curStrategy_carbus = BUS;
        elseif (PAYOFF_BUS < expCar)
            curStrategy_carbus = CAR;
        else
            curStrategy_carbus = randsample(avail_strategies_carbus, 1);
        end
                
        % If CAR was chosen, choose the time.
        if (curStrategy_carbus == CAR)
            
            curStrategy_time = timeTargets(t, player);
            
        else
            % Time for Bus (used for sorting).
            curStrategy_time = 0;
        end
        
        
        % Crazy move.
        if (rand < 0.1)
            if (curStrategy_carbus == BUS) 
                curStrategy_carbus = CAR;
                curStrategy_time = randi(nr_strategies_time);
            else
                curStrategy_carbus = BUS;
            end
        end
        
        if (FIT && t > 1)
            if (curStrategy_carbus ~= strategies_carbus(player,(t-1)))
                strategy_switches(player) = 1;
            end
        end
        
        % Store chosen strategies.
        strategies_carbus(player, t) = curStrategy_carbus;
        strategies_time(player, t) = curStrategy_time;
        
    end
    
    % Sort players by departure time.
    [~, sorted_players] = sort(strategies_time(:, t));
    
    % Compute payoffs, propensities and probabilities.
    leftCars = CAR_NUMBER;
    for player = 1 : N
        
        choseCarGotCar = 0;
        choseBus = 0;
        
        idx = sorted_players(player);
        
        % BUS.
        if (strategies_carbus(idx, t) == BUS)
            payoff = PAYOFF_BUS;
            choseBus = 1;
            
            % CAR.
        else
            time = strategies_time(idx, t);
            if (leftCars > 0)
                choseCarGotCar = 1;
                gotCars(idx, t) = 1;
                leftCars = leftCars - 1;
                payoff = PAYOFF_CAR + (SLOPE_CAR * time);
            else
                payoff = PAYOFF_CAR - (SLOPE_CAR_MISS * time);
            end
        end        
        
        payoffs(idx, t) = payoff;        
        
        % Updating Beliefs.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (choseBus == 1)            
            % No update.
            timeTargets(t+1, idx) = timeTargets(t, idx);
            expected_payoffs_car(t+1, idx) = expected_payoffs_car(t, idx);
            
        % Car.
        else
            
            % Expected payoff of car for next round.
            
            if (choseCarGotCar)
                timeTarget = min(time + INCREASE_SHOCK, ...
                    nr_strategies_time);
                %expCarPayoff = PAYOFF_CAR + (SLOPE_CAR * timeTarget);
            else
                timeTarget = max(time - DECREASE_SHOCK, 1);
                %expCarPayoff = PAYOFF_CAR - (SLOPE_CAR_MISS * timeTarget);
            end
            
            expCarPayoff = PAYOFF_CAR + (SLOPE_CAR * timeTarget);
            expected_payoffs_car(t+1, player) = expCarPayoff;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end

    if (DEBUG)
        carPlayers = find(strategies_carbus(:,t) == CAR);
        shareBus(t) = 20 - length(carPlayers);
        meanDep(t) = mean(strategies_time(carPlayers,t));
    end
    
    if (FIT)         
        meanBus = length(find(strategies_carbus(:,t) == BUS)) / N;
        meanSquaredErrBus = 100 * (fitBus(t, FIT_MEAN_COL) - meanBus)^2;
        
        meanTime = mean(strategies_time(find(strategies_carbus(:,t) == CAR)));
        meanSquaredErrTime = 100 * (fitDepTime(t, FIT_MEAN_COL) - meanTime)^2;
        
        if (t ~= 1)
            meanSwitch = mean(strategy_switches);
            meanSquaredErrSwitch = 100 * (fitSwitch(t, FIT_MEAN_COL) - meanSwitch)^2;
        else
            meanSquaredErrSwitch = 0;
        end
        
        mseBus(t) = meanSquaredErrBus;
        mseTime(t) = meanSquaredErrTime;
        mseSwitch(t) = meanSquaredErrSwitch;
        
    end
    
end

if (DEBUG)
    CAR_NUMBER
    PAYOFF_BUS
    carPlayers = find(strategies_carbus(:,t) == CAR);
    nCars = length(carPlayers)
    avgDepTimeCar = mean(strategies_time(carPlayers,t))
    % pause(0.4);
    
    plot(1:30, shareBus);
    plot(1:30, meanDep);
    pause(0.5);
    
end

if (DUMP)
      
    if (FIT)
        
        fileNameMse = [dumpDir 'mse_' int2str(simCount) '.csv'];
        
        rounds = (1:T)';
        car_levels = repmat(car_level, T, 1);
        bus_payoffs = repmat(PAYOFF_BUS, T, 1);
        sessions = repmat(simCount, T, 1);
        repetitions = repmat(r, T, 1);
        
        
        inits = repmat(INIT_T1, T, 1);
        rewardGotCars = repmat(REWARD_GOT_CAR, T, 1);
        heteros = zeros(T, 1); % TODO: update if heterogeneity is used.
        S1s = repmat(S1, T, 1);
        epsilons = repmat(epsilon, T, 1);
        phis = repmat(phi, T, 1);
        rho1s = repmat(rho1, T, 1);
        wPluss = repmat(wPlus, T, 1);
        wMinuss = repmat(wMinus, T, 1);
        upsilons = repmat(upsilon, T, 1);
        increaseShocks = repmat(INCREASE_SHOCK, T, 1);
        decreaseShocks = repmat(DECREASE_SHOCK, T, 1);
        intervals = repmat(TIME_INTERVAL, T, 1);
        
        csvMatrix = [ ...
            inits, ...
            heteros, ...
            rewardGotCars, ...
            S1s, ...
            epsilons, ...
            phis, ...
            rho1s, ...
            wPluss, ...
            wMinuss, ...
            upsilons, ...
            increaseShocks, ...
            decreaseShocks, ...
            intervals, ...
            sessions, ...
            repetitions, ...
            bus_payoffs, ...
            car_levels, ...
            rounds, ...
            mseBus, ...
            mseTime, ...
            mseSwitch
        ];
        
        if (r == 1)
             header = { ...
            'init', ...
            'hetero', ...
            'reward.car', ...
            'S1', ...
            'epsilon', ...
            'phi', ...
            'rho1', ...
            'wPlus', ...
            'wMinus', ...
            'upsilon', ...
            'increase.shock', ...
            'decrease.shock', ...
            'interval', ...
            'session', ...
            'repetition', ...
            'payoff.bus', ...
            'car.level', ...
            'round', ...
            'msd.bus', ...
            'msd.time', ...
            'msd.switch' ...
            };
            
            csvwrite_with_headers(fileNameMse, csvMatrix, header);
        else            
            dlmwrite(fileNameMse, csvMatrix,'-append','delimiter',',');
        end
    
    else
          fileName = [dumpDir 'data_' int2str(simCount) '.csv'];
          
          rounds = repmat((1:T)',N,1);
          players = repmat(1:N, T, 1);
          reshaped_carbus = reshape(strategies_carbus', N*T, 1);
          reshaped_times = reshape(strategies_time', N*T, 1);
          reshaped_gotCars = reshape(gotCars', N*T, 1);
          reshaped_payoffs = reshape(payoffs', N*T, 1);
          car_levels = repmat(car_level, N*T, 1);
          bus_payoffs = repmat(PAYOFF_BUS, N*T, 1);
          sessions = repmat(simCount, N*T, 1);
          repetitions = repmat(r, N*T, 1);
          
          
          inits = repmat(INIT_T1, N*T, 1);
          rewardGotCars = repmat(REWARD_GOT_CAR, N*T, 1);
          heteros = zeros(N*T, 1); % TODO: update if heterogeneity is used.
          S1s = repmat(S1, N*T, 1);
          epsilons = repmat(epsilon, N*T, 1);
          phis = repmat(phi, N*T, 1);
          rho1s = repmat(rho1, N*T, 1);
          wPluss = repmat(wPlus, N*T, 1);
          wMinuss = repmat(wMinus, N*T, 1);
          upsilons = repmat(upsilon, N*T, 1);
          increaseShocks = repmat(INCREASE_SHOCK, N*T, 1);
          decreaseShocks = repmat(DECREASE_SHOCK, N*T, 1);
          intervals = repmat(TIME_INTERVAL, N*T, 1);
          
          csvMatrix = [ ...
              inits, ...
              heteros, ...
              rewardGotCars, ...
              S1s, ...
              epsilons, ...
              phis, ...
              rho1s, ...
              wPluss, ...
              wMinuss, ...
              upsilons, ...
              increaseShocks, ...
              decreaseShocks, ...
              intervals, ...
              sessions, ...
              repetitions, ...
              bus_payoffs, ...
              car_levels, ...
              players(:), ...
              rounds ...
              reshaped_carbus ...
              reshaped_times ...
              reshaped_gotCars ...
              reshaped_payoffs ...
              ];
          
          
          % Headers for repetition 1.
          if (r == 1)
              
              header = { ...
                  'init', ...
                  'hetero', ...
                  'reward.car', ...
                  'S1', ...
                  'epsilon', ...
                  'phi', ...
                  'rho1', ...
                  'wPlus', ...
                  'wMinus', ...
                  'upsilon', ...
                  'increase.shock', ...
                  'decrease.shock', ...
                  'interval', ...
                  'session', ...
                  'repetition', ...
                  'payoff.bus', ...
                  'car.level', ...
                  'player', ...
                  'round', ...
                  'decision', ...
                  'departure.time', ...
                  'got.car', ...
                  'payoff' ...
                  };
              
              csvwrite_with_headers(fileName, csvMatrix, header);
              
              % Do not append headers.
          else
              dlmwrite(fileName, csvMatrix,'-append','delimiter',',');
          end
    end

% carPlayers = find(strategies_carbus(:,t) == CAR);
% avgDepTimeCar = mean(strategies_time(carPlayers,t));
% 
% rep_nCars(r) = length(carPlayers);
% rep_avgCarTime(r) = avgDepTimeCar;
% rep_depCarTime(r,:) = strategies_time(:,t);
% rep_propensities_carbus(:,:,r) = propensities_carbus;
% rep_propensities_time(:,:,r) = propensities_time;


end        
     

end

