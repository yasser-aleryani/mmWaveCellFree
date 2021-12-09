% This function accept the equivalent channel matrices and calculate the
% digital beamforming matices of each AP under network configuration C_j
function [R_Beamforming_DRL] = Digital_Beamforming_DRL(Count_AP_Configs,Count_UE_Configs,SubNet_n)


 
Opts = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    'MaxFunEvals',1000000,'MaxIter',1000,'TolFun',1.e-10,'TolX',1.e-10);

global W_All P Clusters_APs Clusters_UEs Number_SubNets Count_Iter;
% Calculating Network Parameters
 
Indeces_AP_SubNet_Main  = cell2mat(Clusters_APs{Count_AP_Configs,1}(SubNet_n));
Numbers_AP_SubNet_Main  = length(Indeces_AP_SubNet_Main);
Indeces_UE_SubNet_Main  = cell2mat(Clusters_UEs{Count_UE_Configs,1}(SubNet_n));
Numbers_UE_SubNet_Main  = length(Indeces_UE_SubNet_Main);

Indeces_AP_Outer_SubNets           = 1:Number_SubNets;
Indeces_AP_Outer_SubNets(SubNet_n) = [];
Indeces_AP_SubNet_Outer         = cell2mat(Clusters_APs{Count_AP_Configs,1}(Indeces_AP_Outer_SubNets));

Indeces_Outer_SubNets           = 1:Number_SubNets;
Indeces_Outer_SubNets(SubNet_n) = [];
Indeces_UE_SubNet_Outer         = cell2mat(Clusters_UEs{Count_UE_Configs,1}(Indeces_Outer_SubNets));
 
% Beamforming Initialization for Outer Network APs
if Number_SubNets > 1 || Count_Iter == 1 || SubNet_n == 1
    for Count_Net_l = 1:Number_SubNets
        Indeces_AP_SubNet_l  = cell2mat(Clusters_APs{Count_AP_Configs,1}(Count_Net_l));
        Indeces_UE_SubNet_l  = cell2mat(Clusters_UEs{Count_UE_Configs,1}(Count_Net_l));
        Numbers_UE_SubNet_l  = length(Indeces_UE_SubNet_l);
        for AP_ml = Indeces_AP_SubNet_l
            % Each Column of this Matrix Related to One User within the j'th Subnetwork
            W_All{AP_ml}  =  randn(Numbers_UE_SubNet_l,Numbers_UE_SubNet_l);
            % To Ensure That The Sum of Portions of Signals of each User is 1
            Normalizing_Vec  =  vecnorm(W_All{AP_ml});
            W_All{AP_ml} = W_All{AP_ml}*diag(Normalizing_Vec);
        end
    end
end

lb = 0*ones(Numbers_AP_SubNet_Main*(Numbers_UE_SubNet_Main*Numbers_UE_SubNet_Main) ,1);
ub =   ones(Numbers_AP_SubNet_Main*(Numbers_UE_SubNet_Main*Numbers_UE_SubNet_Main) ,1);
Aeq = []; Beq = [];  nonlcon = [];
% Defining The Norm Condition
Count_Norm = 1;
Aineq = zeros(Numbers_AP_SubNet_Main*Numbers_UE_SubNet_Main,Numbers_AP_SubNet_Main*(Numbers_UE_SubNet_Main*Numbers_UE_SubNet_Main));
for User_kn = 1: Numbers_UE_SubNet_Main
    for AP_mn = 1 : Numbers_AP_SubNet_Main
        Aineq(Count_Norm,(Count_Norm-1)*Numbers_UE_SubNet_Main + 1:(Count_Norm-1)*Numbers_UE_SubNet_Main + Numbers_UE_SubNet_Main) = 1;
        Count_Norm = Count_Norm + 1;
    end
end
Bineq = ones(Numbers_AP_SubNet_Main*Numbers_UE_SubNet_Main,1);
F_x = @(X) Obj_Fun_Beamforming_DRL(X,P,Numbers_UE_SubNet_Main,...
    Indeces_AP_SubNet_Main,Indeces_UE_SubNet_Main,SubNet_n,Clusters_APs,Count_AP_Configs,Clusters_UEs,Count_UE_Configs,Number_SubNets);
X0 = randn(Numbers_AP_SubNet_Main*(Numbers_UE_SubNet_Main*Numbers_UE_SubNet_Main) ,1);

X = fmincon(F_x,X0,Aineq,Bineq,Aeq,Beq,lb,ub,nonlcon,Opts);

% Per-UE Sum Rate
R_Beamforming_DRL = -F_x(X)/Numbers_UE_SubNet_Main; 
end