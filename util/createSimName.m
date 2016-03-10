function simName = createSimName (simName, DUMP, dumpDir, ~) 
%  CREATESIMNAME Adds identifiers to the name of the simulation

if (DUMP)
    if (nargin > 3)
        simName = [simName '-' time2name()];
    end

    simName = createUniqueDumpDir(dumpDir,simName);    
end

end



function newSimName = createUniqueDumpDir(dumpDir,simName)
%CREATEUNIQUEDUMPDIR Creates a unique name for the dump directory 

variants = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n', ...
            'o','p','q','r','s','t','u','v','w','x','y','z'];

newSimName = simName;
folderName = [ dumpDir simName ];

i=1;
while (exist(folderName,'dir')~=0 )
    newSimName = [simName '-' variants(i)];
    folderName = [ dumpDir newSimName];
    i=i+1;
    
    if (i==length(variants))
        error('Too Many dump directories under the same name');
    end
end

mkdir(folderName);
%fprintf('Created dump dir: %s\n',folderName);

end