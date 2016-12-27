%% Stefano Balietti

%% Initialization of variables 
%Clear workspace
close all
clear
clc

dbstop if error
% dbstop in main at 215

%% Add other directories to path
path(path,'util/'); % Help functions

%% Simulation stuff.
%%%%%%%%%%%%%%%%%%%%

dumpDir = 'dump/';
confDir = 'conf/';

%% Loading Conf
load([confDir 'belief'])

%% Computation type.
%%%%%%%%%%%%%%%%%%%%

compLocal = 1;
compLSF = 2;
compParallel = 3;

clusterDir = '/cluster/home/gess/balistef/matlab/car-sharing-model/';
if (exist(clusterDir, 'dir') == 7)
    comp = compLSF;
else
    comp = compLocal;
end

%% Modifying params locally.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nRuns = 1;
DUMP = 1;
DEBUG = 0;
FIT = 1;

dryLSF = 0;

% nRuns = 20;
% epsilon = 0.1;
% phi = 0.1;


% PAYOFF_BUS = [ 70 50 ];
% CAR_NUMBER = [ 15 10 5 ];

% INCREASE_SHOCK = [10];
% DECREASE_SHOCK = [10];
% 
% TIME_INTERVAL_DECREASE = 3;
% rho1 = 0;

% INIT_T1 = 1;

%% Start Vectorization of Parameters Sets
fprintf('\nStarting...\n');

uniqueSimName = createSimName(simName, DUMP, dumpDir, 0);

% Write file for R to read it.
fileID = fopen('simName.R','a');
fprintf(fileID,'SIM <- "%s"\n', uniqueSimName);
fclose(fileID);

folderName = [ dumpDir uniqueSimName '/' ];

if (comp == compLSF)
    % How many sequential simulations in one task.
    SIMS4TASK = 10;
    % How many tasks group in one job.
    TASKS4JOB = 2;
    
    % Task = container of many simulations.
    taskCount = 1;
    % Container = container of many tasks.
    jobCount = 1;
    
    logFolder = ['log/' uniqueSimName];
    mkdir(logFolder); % The name is unique under the dump directory.
    dumpFolder = [ dumpDir uniqueSimName];
    
    % Local
    % sched = parcluster();
    % sched = findResource('scheduler', 'type', 'local');
    
    if (dryLSF ~= 1)
        % Remote.
        parallel.importProfile('/cluster/apps/matlab/support/BrutusLSF8h.settings')
        sched = findResource('scheduler','type','lsf');
        % sched=parcluster('BrutusLSF8h');
        submitArgs = [' -W 1:00 -R "rusage[mem=2000]" -o ' logFolder '/' uniqueSimName '.log'];
        % submitArgs = [' -W 8:00 -R "rusage[mem=2000]" '];
        set(sched, 'SubmitArguments', submitArgs);
        set(sched, 'DataLocation', [logFolder '/']);
        
        j = createJob(sched);
    end
    
elseif (comp == compParallel)

    matlabpool open
end

simCount = 1;


% Nest several loops to simulate parameter sets.
for i1=1:length(CAR_NUMBER)
    my_CAR_NUMBER = CAR_NUMBER(i1);
    
    for i2=1:length(PAYOFF_BUS)
        my_PAYOFF_BUS = PAYOFF_BUS(i2);        
                    
    for i9=1:length(TIME_INTERVAL)
        my_TIME_INTERVAL = TIME_INTERVAL(i9);
                
    for i10=1:length(REWARD_GOT_CAR)
        my_REWARD_GOT_CAR = REWARD_GOT_CAR(i10);
        
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
        
    for i17b=1:length(rho1_relative_to_bus)
        my_rho1_relative_to_bus = rho1_relative_to_bus(i17b);
        
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
                'FIT', FIT, ...
                'INIT_T1', INIT_T1, ...
                'nRuns', nRuns, ...
                'seed', seed, ...
                'CAR_NUMBER', my_CAR_NUMBER, ...
                'PAYOFF_BUS', my_PAYOFF_BUS, ...          
                'TIME_INTERVAL', my_TIME_INTERVAL, ... 
                'REWARD_GOT_CAR', my_REWARD_GOT_CAR, ...
                'INCREASE_SHOCK', my_INCREASE_SHOCK, ...
                'DECREASE_SHOCK', my_DECREASE_SHOCK, ...
                'S1', my_S1, ...
                'epsilon', my_epsilon, ...
                'phi', my_phi, ...                
                'rho1', my_rho1, ...           
                'rho1_relative_to_bus', my_rho1_relative_to_bus, ...
                'wPlus', my_wPlus, ...
                'wMinus', my_wMinus, ...
                'upsilon', my_upsilon ...
            );
            
            fprintf('\nSim count: %i\n',simCount);            
            fprintf('------------------------------------\n');            
            display(paramsObj);        
            fprintf('------------------------------------\n');
            
            
            if (comp == compLocal)        
                simulation_belief(paramsObj);
            else
                % It is convenient to group together more simulations in one
                % task if simulations are short. Matlab overhead to start on
                % each cluster node is about 1 minute.
                taskIdx = mod(simCount, SIMS4TASK);
                
                if (taskIdx == 0)
                    
                    paramObjs{SIMS4TASK} = paramsObj;
                    
                    if (dryLSF ~=1)
                        createTask(j, @wrappersim, 0, {{paramObjs}});
                    end
                    
                    jobIdx = mod(taskCount, TASKS4JOB);
                    % Submit the job to the scheduler in batches
                    if (jobIdx == 0)
                        fprintf('Starting Job %d with %d tasks with %d jobs.\n', ...
                                 jobCount, TASKS4JOB, SIMS4TASK);
                        
                        if (dryLSF ~= 1)
                            submit(j);
                            % TODO: Should create a new job if it is the last one.
                            j = createJob(sched);
                        end
                        
                        % TODO: Should increment only it is the not last one.
                        jobCount = jobCount + 1;
                        
                    end
                    
                    % Update task count after checking to submit job
                    paramObjs = cell(SIMS4TASK, 1);
                    taskCount = taskCount + 1;
                    
                else
                    paramObjs{taskIdx} = paramsObj;
                end
                                
            end
            % Update simulation count.
            simCount = simCount + 1;
            fprintf('\n\n');            
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
end

% Submit the left-over tasks.
if (comp == compLSF)
    remainingSims = mod((simCount - 1), TASKS4JOB*SIMS4TASK);
    if (remainingSims ~= 0)
        if (taskIdx ~= 0 && dryLSF ~= 1)
            createTask(j, @wrappersim, 0, {{paramObjs}});
        end
        
        fprintf('Starting Last Job %d with %d jobs.\n', ...
                                 jobCount, remainingSims);
        if (dryLSF ~= 1)
            submit(j);
        end
    end
    
elseif (comp == compParallel)
     matlabpool close
end



fprintf('Execution completed correctly\n');
% Exit Matlab when invoked from command line with -r option
%exit
