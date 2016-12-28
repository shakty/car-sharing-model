function [ output_args ] = simulation_belief( args )
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

% N beliefs.

% Shares of how many people I expect will take the BUS.
bus_shares = [ 0.9 0.7 0.5 0.3 0.1 ];
nr_beliefs_bus = length(bus_shares);

% Upper limits of time interval for taking the CAR (associated to
% likelihood of getting a car in that time interval).
time_intervals_lowbound = [ 0, 12, 24, 36, 48 ];
time_intervals_highbound = [ 12, 24, 36, 48, 60 ];

time_intervals = [ time_intervals_lowbound ; time_intervals_highbound ];
nr_beliefs_time = length(time_intervals);

beliefs_times_distr           = [ 0.1 0.1 0.3 0.2 0.2 ];
beliefs_times_distr_cumBefore = [ 0   0.1 0.2 0.5 0.8 ];
beliefs_times_distr_cumAfter  = [ 0.9 0.8 0.4 0.2 0   ];

%% Input Parameters.
%%%%%%%%%%%%%%%%%%%%

simCount = args.simCount;

DUMP = args.DUMP;
dumpDir = args.dumpDir;

DEBUG = args.DEBUG;
INIT_T1 = args.INIT_T1;

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
    
    fitPayoffAdj = csvread(['fit/summary_payoff-adj_round_', suffix ]); 
    fitPayoffAdjCar = csvread(['fit/summary_payoff-adj-car_round_', suffix ]); 
end

%% Data structures for all repetitions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rep_nCars = zeros(REPETITIONS, 1);
% rep_avgCarTime = zeros(REPETITIONS, 1);
% rep_depCarTime = zeros(REPETITIONS, N);
% rep_belief_carbus = zeros(N, 2, REPETITIONS);
% rep_belief_time = zeros(N, nr_strategies_time, REPETITIONS);

for r = 1 : REPETITIONS


% Clear simulation data.


% Init RANDOM for all strategies.
% T+1 so that we do not worry about the last iteration.

% I believe that: the majority | half | the minority of people choose BUS.
beliefs_bus = zeros(nr_beliefs_bus, (T+1), N);

% Random beliefs at the beginning.
y = rand(nr_beliefs_bus, N);
y = y./repmat(sum(y,1),size(y,1),1);
beliefs_bus(:,1,:) = y;


% Important! Now we are NOT updating beliefs about others' departure times.
% I believe that: if I choose CAR, and pick time: 0-12, 13-24, 25-36,
% 37-48, 49-60 I will get the CAR.
% beliefs_times = zeros(nr_beliefs_time, (T+1), N);
y  = rand(nr_beliefs_time, N);
y = y./repmat(sum(y,1),size(y,1),1);
beliefs_times_distr = y;

beliefs_times_distr_cumBefore = cumsum(beliefs_times_distr) - beliefs_times_distr;
flipped_beliefs_time = flipud(beliefs_times_distr);
beliefs_times_distr_cumAfter = cumsum(flipped_beliefs_time) - flipped_beliefs_time;
beliefs_times_distr_cumAfter = flipud(beliefs_times_distr_cumAfter);

% Set initial belief and probabilities based on experimental data.
if (INIT_T1) 

    if (PAYOFF_BUS == 50)
        % Probability of taking the bus at round 1.
        PBUS = 0.2078721;
        
        % Avg. departure time and standard deviation at round 1.
        TCAR = 38.77484;
        TCAR_SD = 20.68923;
        

        
    else
        % Probability of taking the bus at round 1.
        PBUS = 0.3522868;
        
        % Avg. departure time and standard deviation at round 1.
        TCAR = 37.57252;
        TCAR_SD = 21.82269;
      
    
    end
    
    
    
    % Sort the probabilities of each player so to make a clear belief
    % that most people will take CAR or BUS.
    for i=1:N
        if (rand < PBUS)
            % Probabilities that other people take bus in these shares: 
            % 0.9, 0.7, 0.5, 0.3, 0.1
            beliefs_bus(:,1,i) = sort(beliefs_bus(:,1,i), 1, 'ascend');
            if (rand > 0.5)
               tmp = beliefs_bus(end,1,i);
               beliefs_bus(end,1,i) = beliefs_bus((end-1),1,i);
               beliefs_bus((end-1),1,i) = tmp;
            end
            
        else
            beliefs_bus(:,1,i) = sort(beliefs_bus(:,1,i), 1, 'descend');          
            if (rand > 0.5)
               tmp = beliefs_bus(1,1,i);
               beliefs_bus(1,1,i) = beliefs_bus(2,1,i);
               beliefs_bus(2,1,i) = tmp;
            end
            
        end
        
        
    end    
   
    % Check.
    % sum(beliefs_bus(1,1,:) > beliefs_bus(3,1,:))
end


% Store strategies, payoffs and choices over all rounds.
payoffs = zeros(N, T);
gotCars = zeros(N, T);
strategies_carbus = ones(N, T);
strategies_time = ones(N, T);
strategies_time_interval = ones(N, T);

if (DEBUG)
    shareBus = zeros(T, 1);
    meanDep = zeros(T, 1);
end
    
if (FIT)
    mseBus = zeros(T, 1);
    mseTime = zeros(T, 1);
    mseSwitch = zeros(T, 1);
    msePayoffAdj = zeros(T, 1);
    msePayoffAdjCar = zeros(T, 1);
    
    meanBusShare = zeros(T, 1);
    meanDepTime = zeros(T, 1);
    meanSwitches = zeros(T, 1);
    meanPayoffAdj = zeros(T, 1);
    meanPayoffAdjCar = zeros(T, 1);
end

%% One simulation.
%%%%%%%%%%%%%%%%%%

for t = 1 : T

    if (FIT && t > 1)
        strategy_switches = zeros(N,1);
    end
    
    probabilities_getcar = zeros(nr_beliefs_bus, nr_beliefs_time, N);
    
    for player = 1 : N
        
        %% Make choice of BUS/CAR and TIME.
        
        beliefs_player = beliefs_bus(:, t, player);
        mostProbableBusBeliefIdx = find(beliefs_player == max(beliefs_player));
        
        % Changed to CAR if a profitable time interval for CAR is found.
        curStrategy_carbus = BUS;
        curStrategy_time = 0;
        curTimeInterval = 0;
        expCarPayoffFound = -1;
        
        for b = 1 : nr_beliefs_bus
        
            bus_share = bus_shares(b); 
           
            for j = 1 : nr_beliefs_time
                
                % Probability of getting a car =
                % people who choose bus
                % + 1/2 people who choose same time interval
                % + people who departed after.
                probGetCar = bus_share +  ...
                    (1-bus_share) * (beliefs_times_distr(j, player) / 2 ) + ...
                    (1-bus_share) * beliefs_times_distr_cumAfter(j, player);
                
                probabilities_getcar(b, j, player) = probGetCar;
                
                % We compute payoffs only for the most probable BUS prob.
                if (b ~= mostProbableBusBeliefIdx)
                    continue;
                end
                
                
                % Time target randomly selected in interval.
                timeTarget = randi(time_intervals(:,j));
                
                expCarPayoff = PAYOFF_CAR + (SLOPE_CAR * timeTarget);
                expCarPayoff = expCarPayoff * probGetCar;
                
                expCarPayoff = expCarPayoff + (randi([-30 30]) / t);
                
                if (0 && DEBUG)
                    probGetCar
                    b
                    timeTarget
                    expCarPayoff
                end
                
                % Choose CAR if exp payoff is larger, or if equal with p=0.5.
                % If current payoff is larger than previous found, or if
                % equal with p=0.5, update.
                if ((expCarPayoff > expCarPayoffFound || ...
                        expCarPayoff == expCarPayoffFound && rand > 0.5) && ...
                        (expCarPayoff > PAYOFF_BUS || ...
                        expCarPayoff == PAYOFF_BUS && rand > 0.5))
                    
                    expCarPayoffFound = expCarPayoff;
                    curStrategy_time = timeTarget;
                    curStrategy_carbus = CAR;
                    curTimeInterval = j;
                    
                end
            end
        end
        
        % Crazy move.
        if (rand < (epsilon/t))
            if (curStrategy_carbus == BUS) 
                curStrategy_carbus = CAR;
                curTimeInterval = randi(nr_beliefs_time);
                curStrategy_time = randi(time_intervals(:,curTimeInterval));
                
                % curStrategy_time = randi(nr_strategies_time);
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
        strategies_time_interval(player, t) = curTimeInterval;
        
    end
    
    % Sort players by departure time.
    [~, sorted_players] = sort(strategies_time(:, t));
    
    % Compute payoffs, belief and probabilities.
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
        player_beliefs = beliefs_bus(:, t, idx);
        
        if (choseBus == 1)            
            % No update.
            beliefs_bus(:, t+1, idx) = player_beliefs;
            
        % Car.
        else
            
            % Given the time interval chosen and if we found a car or not
            % we update the beliefs all car/bus levels.
            chosenInterval = strategies_time_interval(idx, t);
            sumDataProb = 0;
            
            for b = 1 : nr_beliefs_bus
                
                probGetCar = probabilities_getcar(b, chosenInterval, idx);
                
                % Bayes Update.
                if (choseCarGotCar)
                    probGivenEvidence = player_beliefs(b) * probGetCar;
                    beliefs_bus(b, t+1, idx) = probGivenEvidence;
                else
                    probGivenEvidence = player_beliefs(b) * (1 - probGetCar);
                    beliefs_bus(b, t+1, idx) = probGivenEvidence;
                end
               
                sumDataProb = sumDataProb + probGivenEvidence;
            end
            
            beliefs_bus(:, t+1, idx) = beliefs_bus(:, t+1, idx) ./ sumDataProb;
            
            if (sum(beliefs_bus(:, t+1, idx)) < 0.99 || ...
                sum(beliefs_bus(:, t+1, idx)) > 1.01)
                
                beliefs_bus(:, t+1, idx)
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end

    if (DEBUG)
        carPlayers = find(strategies_carbus(:,t) == CAR);
        shareBus(t) = 20 - length(carPlayers);
        meanDep(t) = mean(strategies_time(carPlayers,t));
    end
    
    if (FIT)         
        nBusTakers = length(find(strategies_carbus(:,t) == BUS));
        meanBus = nBusTakers / N;
        meanSquaredErrBus = 100 * (fitBus(t, FIT_MEAN_COL) - meanBus)^2;
        
        carTakers = find(strategies_carbus(:,t) == CAR);
        meanTime = mean(strategies_time(carTakers));
        meanSquaredErrTime = 100 * (fitDepTime(t, FIT_MEAN_COL) - meanTime)^2;
        
        carPayoffs = payoffs(carTakers, t);
        curMeanPayoffAdj = (50 * nBusTakers + sum(carPayoffs)) / N ;
        curMeanPayoffAdjCar = mean(carPayoffs);
        
        meanSquaredErrPayoffAdj = 100 * (fitPayoffAdj(t, FIT_MEAN_COL) - curMeanPayoffAdj)^2;
        meanSquaredErrPayoffAdjCar = 100 * (fitPayoffAdjCar(t, FIT_MEAN_COL) - curMeanPayoffAdjCar)^2;
        
        if (t ~= 1)
            meanSwitch = mean(strategy_switches);
            meanSquaredErrSwitch = 100 * (fitSwitch(t, FIT_MEAN_COL) - meanSwitch)^2;
        else
            meanSwitch = 0;
            meanSquaredErrSwitch = 0;
        end
                
        meanBusShare(t) = meanBus;
        meanDepTime(t) = meanTime;
        meanSwitches(t) = meanSwitch;
        meanPayoffAdj(t) = curMeanPayoffAdj;
        meanPayoffAdjCar(t) = curMeanPayoffAdjCar;
        
        mseBus(t) = meanSquaredErrBus;
        mseTime(t) = meanSquaredErrTime;
        mseSwitch(t) = meanSquaredErrSwitch;
        msePayoffAdj(t) = meanSquaredErrPayoffAdj;
        msePayoffAdjCar(t) = meanSquaredErrPayoffAdjCar;
        
    end
    
end

if (DEBUG)
    CAR_NUMBER
    PAYOFF_BUS
    if (CAR_NUMBER == 10 && PAYOFF_BUS == 50)
        pause(0.1);
    end
    
    carPlayers = find(strategies_carbus(:,t) == CAR);
    nCars = length(carPlayers)
    avgDepTimeCar = mean(strategies_time(carPlayers,t))
    % pause(0.4);
    
    plot(1:T, shareBus);
    plot(1:T, meanDep);
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
            meanBusShare, ...
            meanDepTime, ...
            meanSwitches, ...
            meanPayoffAdj, ...
            meanPayoffAdjCar, ...
            mseBus, ...
            mseTime, ...
            mseSwitch, ...
            msePayoffAdj, ...
            msePayoffAdjCar ...
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
            'bus', ...
            'departure.time', ...
            'decision.switch', ...
            'payoff.adjusted', ...
            'payoff.adjusted.car', ...
            'msd.bus', ...
            'msd.time', ...
            'msd.switch', ...
            'msd.payoff.adjusted', ...
            'msd.payoff.adjusted.car'
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

end        
     

end

