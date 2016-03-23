%% Stefano Balietti

%% Initialization of variables 
%Clear workspace
close all
clear
clc

dbstop if error

%% Add other directories to path
path(path,'util/'); % Help functions

%% Simulation stuff.
%%%%%%%%%%%%%%%%%%%%

dumpDir = 'dump/';
confDir = 'conf/';

%% Loading Conf
load([confDir 'simple'])

%% Modifying params locally.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nRuns = 50;
DUMP = 1;
DEBUG = 0;

% wPlus = 1;
% wMinus = 1;
% 
% rho1 = 0;
% % 
% % 
% % INIT_T1 = 0 ;
% % 
% % Cars.
% CAR_NUMBER = [15]; % CAR_SHARE = 0.75;
% % Payoffs.
% PAYOFF_BUS = [50];
% % 
% % S1 = [1];
% % 
% % INCREASE_SHOCK = 40;
% % 
% % DECREASE_SHOCK = 40;
% % 
% % % Experimentation / Error.
% 
% epsilon = 0;
% % 
% % % Forgetting (or recency).
% phi = 0;
% 
% 
% % Weights assigned to positive and negative reinforcement.
% % How much new experience is weighted against old. (1 = only new).
% wPlus = 0.5;
% wMinus = 0.7;
% 
% TIME_INTERVAL_DECREASE = 10;
% 
% rho1 = - 20;


% csvFile = '/home/stefano/Documents/mypapers/kay_car/data/ALL/summary_exp.csv';
% exp = csvread(csvFile);

%% Start Vectorization of Parameters Sets
fprintf('\nStarting...\n');

uniqueSimName = createSimName(simName, DUMP, dumpDir, 0);

folderName = [ dumpDir uniqueSimName '/' ];


simCount = 1;

% Nest several loops to simulate parameter sets.
for i1=1:length(CAR_NUMBER)
    my_CAR_NUMBER = CAR_NUMBER(i1);
    
    for i2=1:length(PAYOFF_BUS)
        my_PAYOFF_BUS = PAYOFF_BUS(i2);        
                    
    for i9=1:length(TIME_INTERVAL_DECREASE)
        my_TIME_INTERVAL_DECREASE = TIME_INTERVAL_DECREASE(i9);           
        
    for i12=1:length(INCREASE_SHOCK)
        my_INCREASE_SHOCK = INCREASE_SHOCK(i12);
         
    for i13=1:length(DECREASE_SHOCK)
        my_DECREASE_SHOCK = DECREASE_SHOCK(i13); 
        
    for i14=1:length(S1)
        my_S1 = S1(i14);
        
    for i15=1:length(epsilon)
        my_epsilon = epsilon(i15);
        
    for i16=1:length(phi)
        my_phi = phi(i16);
        
    for i17=1:length(rho1)
        my_rho1 = rho1(i17);
        
    for i18=1:length(wPlus)
        my_wPlus = wPlus(i18);
        
    for i19=1:length(wMinus)
        my_wMinus = wMinus(i19);
        
    for i20=1:length(upsilon)
        my_upsilon = upsilon(i20);
         
          
            % Defining seed
            if (seedtype == seed_machinetime)
                % Set seed with milliseconds
                seed1 = sscanf(datestr(now, 'FFF'),'%d') * 1000;
                s = RandStream('mcg16807','Seed', seed1);
                RandStream.setGlobalStream(s);
                rng shuffle
                seed2 = randi(1000);
                seed = seed1 + seed2;
            elseif (seedtype == seed_random)
                % Random seed
                seed = randi(1000000);
            elseif (seedtype == seed_fixed)
                % Seed already fixed. seed = seed;                
            end
            
            paramsObj = struct( ...
                'simName', simName, ...
                'simCount', simCount, ...
                'dumpDir', folderName, ...
                'DUMP', DUMP, ...
                'DEBUG', DEBUG, ...
                'INIT_T1', INIT_T1, ...
                'nRuns', nRuns, ...
                'seed', seed, ...
                'CAR_NUMBER', my_CAR_NUMBER, ...
                'PAYOFF_BUS', my_PAYOFF_BUS, ...          
                'TIME_INTERVAL_DECREASE', my_TIME_INTERVAL_DECREASE, ...
                'INCREASE_SHOCK', my_INCREASE_SHOCK, ...
                'DECREASE_SHOCK', my_DECREASE_SHOCK, ...
                'S1', my_S1, ...
                'epsilon', my_epsilon, ...
                'phi', my_phi, ...                
                'rho1', my_rho1, ...
                'wPlus', my_wPlus, ...
                'wMinus', my_wMinus, ...
                'upsilon', my_upsilon ...
            );
            
            fprintf('\nSim count: %i\n',simCount);            
            fprintf('------------------------------------\n');            
            display(paramsObj);        
            fprintf('------------------------------------\n');
        
            simulation(paramsObj);

            fprintf('\n\n');
            
            % Updating the simulations count.
            simCount=simCount+1;
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
    end
end


fprintf('Execution completed correctly\n');
% Exit Matlab when invoked from command line with -r option
%exit