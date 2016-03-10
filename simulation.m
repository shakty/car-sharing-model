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

% if (DUMP)
%     OUT_DIR = [ args.dumpDir int2str(simCount) '/' ];
%     mkdir(OUT_DIR);
% end

CAR_NUMBER = args.CAR_NUMBER;

if (CAR_NUMBER == 15)
    car_level = 0.75;
elseif (CAR_NUMBER == 10)
    car_level = 0.5;
else
    car_level = 0.25;
end


PAYOFF_BUS = args.PAYOFF_BUS;
INCREASE_BUS = args.INCREASE_BUS;
INCREASE_CAR_GOT = args.INCREASE_CAR_GOT;                
INCREASE_CAR_MISSED = args.INCREASE_CAR_MISSED;
INCREASE_TIME = args.INCREASE_TIME;
DECREASE_TIME = args.DECREASE_TIME;
TIME_INTERVAL_INCREASE = args.TIME_INTERVAL_INCREASE;
TIME_INTERVAL_DECREASE = args.TIME_INTERVAL_DECREASE;
INCREASE_DECAY = args.INCREASE_DECAY;
DECREASE_DECAY = args.DECREASE_DECAY;
INCREASE_SHOCK = args.INCREASE_SHOCK;
DECREASE_SHOCK = args.DECREASE_SHOCK;


REPETITIONS = args.nRuns;

increase_decay = INCREASE_DECAY * INCREASE_TIME;
decrease_decay = DECREASE_DECAY * DECREASE_TIME;

%% Data structures for all repetitions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rep_nCars = zeros(REPETITIONS, 1);
% rep_avgCarTime = zeros(REPETITIONS, 1);
% rep_depCarTime = zeros(REPETITIONS, N);
% rep_propensities_carbus = zeros(N, 2, REPETITIONS);
% rep_propensities_time = zeros(N, nr_strategies_time, REPETITIONS);

for r = 1 : REPETITIONS


% Clear simulation data.
    
% Initial propensities for car/bus for each strategy of each player.
propensities_carbus = ones(N, nr_strategies);

% Initial propensities for time for each strategy of each player.
propensities_time = ones(N, nr_strategies_time);

% Probabilities to choose car vs bus. (all equally probable at t=0).
probabilities = ones(N, nr_strategies)*(1/nr_strategies);

% Probabilities to choose a time. (4all equally probable at t=0).
probabilities_time = ones(N, nr_strategies_time)*(1/nr_strategies_time);

% Store strategies, payoffs and choices over all rounds.
payoffs = zeros(N, T);
gotCars = zeros(N, T);
strategies_carbus = ones(N, T);
strategies_time = ones(N, T);


%% One simulation.
%%%%%%%%%%%%%%%%%%

for t = 1 : T
    
    for player = 1 : N
        
        
        % Randomly draw a strategy according to the prob distribution.
        curStrategy_carbus = randsample(avail_strategies_carbus, ...
            1, true, ...
            probabilities(player,:));
        
        
        % If CAR was chosen, choose the time.
        if (curStrategy_carbus == CAR)
            
            % Randomly draw a time strategy according to prob distribution.
            curStrategy_time = randsample(avail_strategies_time, ...
                1, true, ...
                probabilities_time(player,:));
        else
            % Time for Bus (used for sorting).
            curStrategy_time = 0;
        end
        
        % Store chosen strategies.
        strategies_carbus(player, t) = curStrategy_carbus;
        strategies_time(player, t) = curStrategy_time;
        
    end
    
    [~, sorted_players] = sort(strategies_time(:, t));
    
    % Compute payoffs, propensities and probabilities.
    leftCars = CAR_NUMBER;
    for player = 1 : N
        
        choseCarGotCar = 0;
        choseCarMissed = 0;
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
                choseCarMissed = 1;
                payoff = PAYOFF_CAR - (SLOPE_CAR * time);
            end
        end
        
        payoffs(idx, t) = payoff;
        
        if (choseCarGotCar == 1)
            
            % Increase Car propensity.
            propensities_carbus(idx, CAR) = ...
                propensities_carbus(idx, CAR) + INCREASE_CAR_GOT;
            % TODO: could do an in increase proportional payoff.
            
            
            increase = INCREASE_TIME;
            upLimit = min(time + TIME_INTERVAL_INCREASE, ...
                nr_strategies_time);
            downLimit = time + INCREASE_SHOCK;
            for i = downLimit : upLimit
                propensities_time(idx, i) = ...
                    propensities_time(idx, i) + increase;
                increase = increase - increase_decay;
                if (increase <= 0)
                    break;
                end
            end
            
            
        elseif (choseCarMissed == 1)
            
            % Increase Bus propensity.
            propensities_carbus(idx, BUS) = ...
                propensities_carbus(idx, BUS) + INCREASE_CAR_MISSED;
            % TODO: could do an in increase proportional to payoff.
            
            increase = DECREASE_TIME;
            downLimit = max(time - TIME_INTERVAL_DECREASE, 1);
            upLimit = (time - DECREASE_SHOCK);
            for i = downLimit : upLimit
                propensities_time(idx, i) = ...
                    propensities_time(idx, i) + increase;
                increase = increase - decrease_decay;
                if (increase <= 0)
                    break;
                end
            end
            
            
            % Bus.
        else
            
            % Increase Bus propensity.
            propensities_carbus(idx, BUS) = ...
                propensities_carbus(idx, BUS) + INCREASE_BUS;            
            
        end
              
        
        % Update probabilities.
        
        sumPropensities = sum(propensities_carbus(idx, :));
        
        probabilities(idx, BUS) = propensities_carbus(idx, BUS) / ...
            sumPropensities;
        probabilities(idx, CAR) = propensities_carbus(idx, CAR) / ...
            sumPropensities;
        
        % Update time probabilities only if car was chosen.
        if (choseBus == 0)
            sumPropensities = sum(propensities_time(idx, :));
            
            for i = 1 : nr_strategies_time
                probabilities_time(idx, i) = propensities_time(idx, i) / ...
                    sumPropensities;
            end
        end
        
    end
    
    
end

if (DUMP)
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
    
    csvMatrix = [ ...
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
    
     header = { ...
        'session', ...
        'repetition', ...
        'payoff.bus', ...
        'car.level.num', ...
        'player', ...
        'round', ...
        'decision', ...
        'departure.time', ...
        'got.car', ...
        'payoff' ...
    };
    
if (r == 1)    
    csvwrite_with_headers(fileName, csvMatrix, header);
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


% % Stats.
% 
% avgCar = mean(rep_nCars)
% % How close to equilibrium?
% EQ = CAR_NUMBER - avgCar
% 
% avgDepTimeCar = mean(rep_avgCarTime)

% fid = fopen(dataFileName, 'a');
% 
% %% Param
% 
% headers = {
%     'simname', ...
%     'simcount', ...
%     'run', ...
%     'car.number', ...
%     'bus.payoff', ...
%     'round', ...
%     'avg.car.n', ...
%     'avg.car.time' ...
% };
% 
% write_csv_headers(fileName, headers);
                


        
end

