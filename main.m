%% Stefano Balietti

%% Initialization of variables 
%Clear workspace
close all
clear
clc

%% Add other directories to path
path(path,'util/'); % Help functions

confDir = 'conf/mat/';

%% Loading Conf
load([confDir 'b50_c75_rl'])

%% Simulation stuff.
%%%%%%%%%%%%%%%%%%%%

dumpDir = 'dump/';
confDir = 'conf/mat/';

compLOCAL = 0;
compPARALLEL = 1;
compLSF = 2;

%% Modifying params locally.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Start Vectorization of Parameters Sets
fprintf('\nStarting...\n');

simCount = 1;

% Nest several loops to simulate parameter sets.
for i1=1:length(CAR_NUMBER)
    my_CAR_NUMBER = CAR_NUMBER(i1);
    
    for i2=1:length(PAYOFF_BUS)
        my_PAYOFF_BUS = PAYOFF_BUS(i2);   
    
    for i3=1:length(INCREASE_BUS)
        my_INCREASE_BUS = INCREASE_BUS(i3);
        
    for i4=1:length(INCREASE_CAR_MISSED)
        my_INCREASE_CAR_MISSED = INCREASE_CAR_MISSED(i4);
                
    for i5=1:length(INCREASE_CAR_GOT)
        my_INCREASE_CAR_GOT = INCREASE_CAR_GOT(i5);
                        
    for i6=1:length(INCREASE_TIME)
        my_INCREASE_TIME = INCREASE_TIME(i6);
         
    for i7=1:length(DECREASE_TIME)
        my_DECREASE_TIME = DECREASE_TIME(i7);
                
    for i8=1:length(TIME_INTERVAL_INCREASE)
        my_TIME_INTERVAL_INCREASE = TIME_INTERVAL_INCREASE(i8);
                
    for i9=1:length(TIME_INTERVAL_DECREASE)
        my_TIME_INTERVAL_DECREASE = TIME_INTERVAL_DECREASE(i9);
        
    for i10=1:length(INCREASE_DECAY)
        my_INCREASE_DECAY = INCREASE_DECAY(i10);
         
    for i11=1:length(DECREASE_DECAY)
        my_DECREASE_DECAY = DECREASE_DECAY(i11);
        
    for i12=1:length(INCREASE_SHOCK)
        my_INCREASE_SHOCK = INCREASE_SHOCK(i12);
         
    for i13=1:length(DECREASE_SHOCK)
        my_DECREASE_SHOCK = DECREASE_SHOCK(i13);
 
         
          
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
                'simCount', simCount, ...
                'nRuns', nRuns, ...
                'seed', seed, ...
                'CAR_NUMBER', my_CAR_NUMBER, ...
                'PAYOFF_BUS', my_PAYOFF_BUS, ...
                'INCREASE_BUS', my_INCREASE_BUS, ...
                'INCREASE_CAR_GOT', my_INCREASE_CAR_GOT, ...                
                'INCREASE_CAR_MISSED', my_INCREASE_CAR_MISSED, ...
                'INCREASE_TIME', my_INCREASE_TIME, ...
                'DECREASE_TIME', my_DECREASE_TIME, ...
                'TIME_INTERVAL_INCREASE', my_TIME_INTERVAL_INCREASE, ...
                'TIME_INTERVAL_DECREASE', my_TIME_INTERVAL_DECREASE, ...
                'INCREASE_DECAY', my_INCREASE_DECAY, ...
                'DECREASE_DECAY', my_DECREASE_DECAY, ...
                'INCREASE_SHOCK', my_INCREASE_SHOCK, ...
                'DECREASE_SHOCK', my_DECREASE_SHOCK ...              
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
end


fprintf('Execution completed correctly\n');
% Exit Matlab when invoked from command line with -r option
%exit