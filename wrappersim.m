function wrappersim( paramArgs )
    
    %% Add other directories to path
    path(path,'util/'); % Help functions

    for i=1:length(paramArgs)
        simulation_belief(paramArgs{i})
    end

end

