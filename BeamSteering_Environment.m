
function BeamSteering_ENV = BeamSteering_Environment()
  % This Function Defines the Environment of the Problem
global Clusters_APs Clusters_UEs Partitioning_Indeces Number_SubNets M K u a; 
Clusters_APs         = SetPartition(M,Number_SubNets); % Possible Configurations of APs
Clusters_UEs         = SetPartition(K,Number_SubNets); % Possible Configurations of UEs
APs_Configs          = 1:length(Clusters_APs);
UEs_Configs          = 1:length(Clusters_UEs);
[Par_A,Par_B]        = meshgrid(APs_Configs,UEs_Configs);
partitioning_c       = cat(2,Par_A',Par_B');
Partitioning_Indeces = reshape(partitioning_c,[],2);
Par_Actions_Vector   = 1:length(Partitioning_Indeces(:,1));

BeamSteering_Actions_Vector = a*(M-Number_SubNets+1)+u*(K-Number_SubNets+1);

%% Defining the Spaces for Both Actions and States
% Defining States Spaces 
BeamSteering_StateSpace             = rlNumericSpec([K-Number_SubNets+1 1]);
BeamSteering_StateSpace.LowerLimit  = 0;
BeamSteering_StateSpace.Name        = 'Per-SubNet Sum Rate';
BeamSteering_StateSpace.Description = 'f(Per-UE Gammas of SubNet)';

% Defining Action Space

BeemSteering_ActionSpace             = rlNumericSpec(Par_Actions_Vector);
BeemSteering_ActionSpace.UpperLimit  = 1;
BeemSteering_ActionSpace.LowerLimit  = 0;
BeemSteering_ActionSpace.Name        = 'BeamSteering Vector of SubNet';

% States Initialization Initialization for Network Partitioning System
[BeamSteering_InitialObservation,BeamSteering_LoggedSignals] = BeamSteering_State_Initialization(); 
  
% Objrctive Function Defined Here (Beamsteering-->Beamforming Functions Inclused Inside the Step Function)
% Initial_Action = 1; 
% Par_Obj_Fun(Initial_Action,Par_LoggedSignals); 

BeamSteering_ENV = rlFunctionEnv(BeamSteering_StateSpace,BeemSteering_ActionSpace,'BeamSteering_Obj_Fun','BeamSteering_State_Initialization');
end