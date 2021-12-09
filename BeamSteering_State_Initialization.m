function [Par_InitialObservation, Par_LoggedSignal] = BeamSteering_State_Initialization()
% This Function Reset Environement States into Some Random Values
 
global Number_SubNets K;
 
Partitioning_States = 100*rand(K-Number_SubNets+1,1);  
 
% Return initial environment state variables as logged signals.
Par_LoggedSignal.State = Partitioning_States;
Par_InitialObservation = Par_LoggedSignal.State;

end