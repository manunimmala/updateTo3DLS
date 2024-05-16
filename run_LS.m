%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LSMODEL_RUNSIMULATIONS Run Lagrangian stochastic model simulations for stable 
% and unstable conditions
%   This script loads the necessary variables from 'runFile.mat' and runs
%   the Lagrangian stochastic model simulations for stable and unstable
%   atmospheric conditions. The simulations are run in parallel using the
%   'parfor' loop.
%
%   The script first determines the indices of the simulations that still
%   need to be run by comparing the simulation numbers in the 'runsdone'
%   folder with the total number of simulations specified in 'simNum'.
%
%   For each simulation that needs to be run, the script checks the Obukhov
%   length (L) to determine whether the case is stable (L > 0) or unstable
%   (L < 0). It then calls the appropriate function (LS_stableModel or
%   LS_unstableModel) with the corresponding input parameters.
%
%   The output variables (cgrid, xgrid, zgrid, depgrid, and GLC) are saved
%   as binary files in the 'results' folder within the specified output
%   folder. The binary files are named according to the 'outputFileNames'
%   array.
%
%   Additionally, a binary file is created in the 'runsdone' folder to keep
%   track of the completed simulations. The file contains the simulation
%   number, the elapsed time for the simulation, and the current date and
%   time.
%
%   Note: Make sure that the 'runFile.mat' file contains all the necessary
%   variables for running the simulations, including:
%   - ustar, wstar, L, z_i, z0, xmin, xmax, zmin, zmax, np, vs, x0, h0,
%     nxgrid, nzgrid, C0: Input parameters for the LS_stableModel and
%     LS_unstableModel functions
%   - outputFolder: Path to the output folder where the results will be saved
%   - outputFileNames: Cell array containing the names for the output files
%   - simNum: Total number of simulations to be run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load 'runFile.mat' which contains all the variables necessary to run
% simulations
load('runFile')

% Find the indices of the simulations that still need to be run
doneSims = sort(sscanf(ls(strcat(outputFolder,'/runsdone')),'%d.bin'))';
notDoneSims = setdiff(1:simNum,doneSims);

%%
% Run simulations
parfor index = 1:length(notDoneSims)
    tic

    i = notDoneSims(index);
    % Stable cases    
    if L(i) > 0 
        [xgrid, zgrid, cgrid,depgrid] = LS_stableModel(ustar(i), wstar(i), ...
            L(i), z_i(i), z0(i), xmin(i), xmax(i), zmin(i), zmax(i), np(i),...
            vs(i), x0(i), h0(i), nxgrid(i), nzgrid(i), C0(i));

    % Unstable cases
    else 
        [xgrid, zgrid, cgrid, depgrid] = LS_unstableModel(ustar(i), wstar(i), ...
            L(i), z_i(i), z0(i), xmin(i), xmax(i), zmin(i), zmax(i), np(i),...
            vs(i), x0(i), h0(i), nxgrid(i), nzgrid(i), C0(i));
    end

    % Save output variables into binary files: 
    % cgrid, xgrid, zgrid, depgrid, and GLC;
    cgridFID = fopen(strcat(outputFolder,'/results/cgrid_',outputFileNames(i),'.bin'),'w');
    xgridFID = fopen(strcat(outputFolder,'/results/xgrid_',outputFileNames(i),'.bin'),'w');
    zgridFID = fopen(strcat(outputFolder,'/results/zgrid_',outputFileNames(i),'.bin'),'w');
    glcFID = fopen(strcat(outputFolder,'/results/glc_',outputFileNames(i),'.bin'),'w');    
    depgridFID = fopen(strcat(outputFolder,'/results/depgrid_',outputFileNames(i),'.bin'),'w');
    runNumFID = fopen(strcat(outputFolder,'/runsdone/',num2str(i),'.bin'),'w');
    
    % Write to binary files
    fwrite(cgridFID,cgrid,'double');
    fwrite(xgridFID,xgrid,'double');
    fwrite(zgridFID,zgrid,'double');       
    fwrite(glcFID,cgrid(:,1),'double');
    fwrite(depgridFID,depgrid,'double');
    fwrite(runNumFID,sprintf('%s\t%s',toc,string(datetime)),'char');


    % Close them
    fclose(cgridFID);
    fclose(xgridFID);
    fclose(zgridFID);
    fclose(glcFID);
    fclose(depgridFID);
    fclose(runNumFID);

end



