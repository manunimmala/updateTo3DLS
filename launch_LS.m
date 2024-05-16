%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RUN_LS_BATCH Run Lagrangian stochastic model simulations in batch mode
%   This script runs the Lagrangian stochastic model simulations in batch
%   mode using the 'batch' function. It submits the 'run_LS' script as a
%   batch job with the specified profile, pool size, and diary capture.
%
%   The script performs the following steps:
%   1. Clears the workspace and closes all figures.
%   2. Prints the current date and time to indicate the start of the simulations.
%   3. Creates a batch job using the 'batch' function, specifying the 'run_LS'
%      script, the 'local' profile, a pool size of 70, and diary capture.
%   4. Waits for the batch job to complete using the 'wait' function.
%   5. Retrieves the diary of the batch job using the 'diary' function.
%   6. Loads the completed batch job using the 'load' function.
%   7. Prints the current date and time to indicate the end of the simulations.
%   8. Deletes the batch job and clears the 'my_job' variable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear the workspace and close all figures
clear all
close all

% Print the current date and time to indicate the start of the simulations
fprintf(1,'Running simulations -- %s \n', string(datetime));

% Create a batch job using the 'batch' function
my_job = batch('run_LS','Profile', 'local', 'Pool', 70, 'CaptureDiary',true);

% Wait for the batch job to complete
wait(my_job);

% Retrieve the diary of the batch job
diary(my_job);

% Load the completed batch job
load(my_job);

% Print the current date and time to indicate the end of the simulations
fprintf(1,'Finished running -- %s \n', string(datetime));

% Delete the batch job and clear the 'my_job' variable
delete(my_job)
clear my_job