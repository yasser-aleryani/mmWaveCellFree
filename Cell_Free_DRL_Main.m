clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Written by: Yasser Al Eryani%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Downlink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global H Number_SubNets M a K u L Count_Iter...
    Global_Reward Clusters_APs Clusters_UEs P a u;

Iterations     = 500;    % Number of Monte-Carlo Iterations
M = 5; a  = 1;  K = 3; u = 1; L = 3; % Network Size (Antenna Numbers should have integer root square)

Number_SubNets = 3; % # of Cell Free Subnetworks

Clusters_APs = SetPartition(M,Number_SubNets); % Possible Configurations of APs
Clusters_UEs = SetPartition(K,Number_SubNets); % Possible Configurations of UEs

% Maximum Allowable Transmission Power Per UE
P_dB = 37; %20:5:50; % In dBm
P = 10.^((P_dB-30)/10);

% Opts = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%     'MaxFunEvals',1000000,'MaxIter',1000,'TolFun',1.e-10,'TolX',1.e-10);
% R_i = zeros(1,Iterations);R=zeros(1,length(P));
% Defining Updated Vectors as
% global Delta_All;

%===========================Code Start Here===============================%
for Count_P=1:length(P)
    Global_Reward = 0;
    for Count_Iter=1:Iterations
        load('H_mmWaves');
        %         H = mmWave_Ch_Generation(M,a,K,u,L);
        
        %% Network Partitioning DRL Model
        
        % Initializing The Network Partitioning Environement
        Par_Env         = Partitioning_Environment();
        
        ObservationInfo_Par = getObservationInfo(Partitioning_Environment);
        
        ActionInfo_Par      = getActionInfo(Partitioning_Environment);
        
        % Defining The Agent (Can be any Agent)
        % Approximating the Critic Function
        Par_DDQN_Agent = rlDQNAgent(ObservationInfo_Par,ActionInfo_Par);
        
        
        % Setting Training Parameters
        Par_TrainOpts = rlTrainingOptions(...
            'MaxEpisodes',1000, ...
            'MaxStepsPerEpisode',100, ...
            'Verbose',false, ...
            'Plots','training-progress',...
            'StopTrainingCriteria','AverageReward',...
            'StopTrainingValue',480);
         
        % Train the DDQN Agent
        doTraining = false;
        if doTraining
            trainingStats = train(Par_DDQN_Agent,Par_Env,Par_TrainOpts);
            save('Trained_Par_Agent','Par_DDQN_Agent');
        else
            % Load A Previously Trained .
            load('Trained_Par_Agent.mat','Par_DDQN_Agent');
        end
        
    end
    Ave_Global_Reward(Count_P) = Global_Reward/Iterations;
end






