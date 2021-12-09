% This Function Reset Environement States into Some Random Values
function [InitialObservation, LoggedSignal] = State_Initialization(M,a,K,u)



for 

% Theta (randomize)
T0 = 2 * 0.05 * rand() - 0.05;
% Thetadot
Td0 = 0;
% X
X0 = 0;
% Xdot
Xd0 = 0;





% Return initial environment state variables as logged signals.
LoggedSignal.State = [X0;Xd0;T0;Td0];
InitialObservation = LoggedSignal.State;

end