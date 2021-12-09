function [Par_InitialObservation, Par_LoggedSignal] = Par_State_Initialization()
% This Function Reset Environement States into Some Random Values
 
global Number_SubNets;


Partitioning_States = 100*rand(Number_SubNets,1);  
 
% % Theta (randomize)
% T0 = 2 * 0.05 * rand() - 0.05;
% % Thetadot
% Td0 = 0;
% % X
% X0 = 0;
% % Xdot
% Xd0 = 0;
 
% Return initial environment state variables as logged signals.
Par_LoggedSignal.State = Partitioning_States;
Par_InitialObservation = Par_LoggedSignal.State;

end