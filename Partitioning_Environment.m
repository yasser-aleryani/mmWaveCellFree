
function Par_ENV = Partitioning_Environment()
  % This Function Defines the Environment of the Problem
global Clusters_APs Clusters_UEs Partitioning_Indeces Number_SubNets M K; 
Clusters_APs         = SetPartition(M,Number_SubNets); % Possible Configurations of APs
Clusters_UEs         = SetPartition(K,Number_SubNets); % Possible Configurations of UEs
APs_Configs          = 1:length(Clusters_APs);
UEs_Configs          = 1:length(Clusters_UEs);
[Par_A,Par_B]        = meshgrid(APs_Configs,UEs_Configs);
partitioning_c       = cat(2,Par_A',Par_B');
Partitioning_Indeces = reshape(partitioning_c,[],2);
Par_Actions_Vector   = 1:length(Partitioning_Indeces(:,1));

%% Defining the Spaces for Partitioning Actions and States
% Defining States Spaces 
Par_StateSpace             = rlNumericSpec([Number_SubNets 1]);
Par_StateSpace.LowerLimit  = 0;
Par_StateSpace.Name        = 'Per-SubNet Sum Rate';
Par_StateSpace.Description = 'f(Per-UE Gammas)';

% Defining Action Space
Par_ActionSpace             = rlFiniteSetSpec(Par_Actions_Vector);
Par_ActionSpace.Name        = 'Partitioning Configuration';

% States Initialization Initialization for Network Partitioning System
[Par_InitialObservation,Par_LoggedSignals] =Par_State_Initialization(); 
  
% Objrctive Function Defined Here (Beamsteering-->Beamforming Functions Inclused Inside the Step Function)
% Initial_Action = 1; 
% Par_Obj_Fun(Initial_Action,Par_LoggedSignals); 

Par_ENV = rlFunctionEnv(Par_StateSpace,Par_ActionSpace,'Par_Obj_Fun','Par_State_Initialization');
end