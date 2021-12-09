function [Par_NextObs,Par_Reward,IsDone,Par_LoggedSignals] = Par_Obj_Fun(Par_Action,Par_LoggedSignals)
% Custom step function to construct cart-pole environment for the function
% name case.
%
 
% This function applies the given action to the environment and evaluates
% the system dynamics for one simulation step.


global Number_SubNets Clusters_APs Clusters_UEs Partitioning_Indeces R_Old h_Equivalent_All Global_Reward;
% This Part is Responsible for Analog Beamsteering at
% each AP.
% This Loop May be deployed through Parallel Programming
 
%% Hybrid Optimization Code

% Analog Beamsteering
for SubNet_n = 1:Number_SubNets
    
    % Calculating Effective Channel gain (Analog Beamsteering)
    Analog_Beamsteering_DRL(Partitioning_Indeces(Par_Action,1),...
        Partitioning_Indeces(Par_Action,2),SubNet_n);
 
end
 
% Digital Beamforming
for SubNet_n = 1:Number_SubNets
    
    R_Beamforming = Digital_Beamforming_DRL(Partitioning_Indeces(Par_Action,1),...
        Partitioning_Indeces(Par_Action,2),SubNet_n);
 
    % Optimized Sum Rate Per SubNet.
    R_SubNets(SubNet_n) = R_Beamforming;
    
end
h_Equivalent_All = [];% Clear Buffer Content of h_Equivalent

R_Reward = sum(R_SubNets)/Number_SubNets;

 


% Unpack the state vector from the logged signals.
% State = Par_LoggedSignals.State;
  
Par_LoggedSignals.State = R_SubNets';

% Transform state to observation.
Par_NextObs = Par_LoggedSignals.State;


% Setting the Reward
Par_Reward = R_Reward;
IsDone = true;

 
 Global_Reward = Global_Reward+ R_Reward;


 

 