% This function accept CSI matrices and return effective CSI vectors after
% conducting analog beamsteering
function  Analog_Beamsteering_DRL(Count_AP_Configs,Count_UE_Configs,SubNet_n)
global H Delta_All A_All h_Equivalent_All Count_Iter a K u...
       Clusters_APs Clusters_UEs Number_SubNets;


Opts = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
    'MaxFunEvals',1000000,'MaxIter',1000,'TolFun',1.e-10,'TolX',1.e-10);

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
Numbers_UE_SubNet_Outer         = length(Indeces_UE_SubNet_Outer);

% Beamforming Initialization
if Count_Iter == 1 || Count_AP_Configs == 1 || Count_UE_Configs == 1 || SubNet_n == 1
    Delta_All  =  2*pi*rand( K*(u*1)                     ,1);
end
% Optimize and Stor Dependent Variables
% Beamsteering Matrices Initialization
lb = 0*ones( Numbers_AP_SubNet_Main*(a*Numbers_UE_SubNet_Main) + (Numbers_UE_SubNet_Main*u),1);
ub = 2*pi*ones(Numbers_AP_SubNet_Main*(a*Numbers_UE_SubNet_Main) + (Numbers_UE_SubNet_Main*u),1);
Aeq = []; Beq = []; Aineq = []; Bineq = []; NLCON = [];
X0 = 2*pi*rand(Numbers_AP_SubNet_Main*(a*Numbers_UE_SubNet_Main) + (Numbers_UE_SubNet_Main*u),1);

% Implementing Objective Function for the i'th Subnetwork
F_x = @(X) Obj_Fun_Beamsteering_DRL(X,H,Numbers_UE_SubNet_Main,SubNet_n,a,u,...
    Indeces_AP_SubNet_Main,Indeces_UE_SubNet_Main, Indeces_UE_SubNet_Outer,Numbers_UE_SubNet_Outer);
X_Optimal = fmincon(F_x,X0,Aineq,Bineq,Aeq,Beq,lb,ub,NLCON,Opts);
X_A_Opt = X_Optimal(1:Numbers_AP_SubNet_Main*(a*Numbers_UE_SubNet_Main));
X_Delta = X_Optimal(Numbers_AP_SubNet_Main*(a*Numbers_UE_SubNet_Main)+1:end);


Count_AP   = 1;
for AP_mn = Indeces_AP_SubNet_Main
    A_All{AP_mn,SubNet_n} = reshape(X_A_Opt( (Count_AP-1)*(a*Numbers_UE_SubNet_Main) + 1:...
        (Count_AP-1)*(a*Numbers_UE_SubNet_Main) + a*Numbers_UE_SubNet_Main ),a,Numbers_UE_SubNet_Main);
    Count_User = 1;
    
    for User_kn = Indeces_UE_SubNet_Main
        Delta_All( (User_kn-1)*(u*1)+1:((User_kn-1)*(u*1)) + u*1) = ...
            X_Delta( (Count_User  -1)*(u*1)+1: ((Count_User  -1)*(u*1)) + u*1);
        
        % Calculating Equivalent Channel Gain of AP_i-->User_i links
        h_Equivalent_All{AP_mn,User_kn}(1:Numbers_UE_SubNet_Main)        = Delta_All( (User_kn-1)*(u*1)+1:((User_kn-1)*(u*1)) + u*1)'...
            *H{AP_mn,User_kn}*A_All{AP_mn,SubNet_n};
        
        Count_User = Count_User + 1;
    end
    
    for User_i_Outer = Indeces_UE_SubNet_Outer
        
        % Calculating Equivalent Channel Gain of AP_i-->User_i links
        h_Equivalent_All{AP_mn,User_i_Outer}(1:Numbers_UE_SubNet_Main) = Delta_All( (User_i_Outer-1)*(u*1)+1:((User_i_Outer-1)*(u*1)) + u*1)'...
            *H{AP_mn,User_i_Outer}*A_All{AP_mn,SubNet_n};
        
    end
    Count_AP = Count_AP + 1;
end

 
end


