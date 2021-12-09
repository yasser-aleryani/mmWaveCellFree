% This function calculate the objective function related to analog
% beamsteering problem
function F_x = Obj_Fun_Beamsteering_DRL(X,H,Numbers_UE_SubNet_Main,SubNet_n,a,u,...
    Indeces_AP_SubNet_Main,Indeces_UE_SubNet_Main, Indeces_UE_SubNet_Outer,Numbers_UE_SubNet_Outer)
global Delta_All H_SubNet H_SubNet_Interference;

Alpha = 0.3; % Trade-Off Factor Between Beamsteering and Interference Nulling
F_x_Non_Orthogonal = 0;
F_x_Orthogonal = 0;
for AP_mn = 1:length(Indeces_AP_SubNet_Main)

    % Defining Beamsteering Matrices for Obtimization
    A_SubNet_Main{AP_mn} = reshape(X( (AP_mn-1)*(a*Numbers_UE_SubNet_Main) + 1:(AP_mn-1)*(a*Numbers_UE_SubNet_Main)+ a*Numbers_UE_SubNet_Main ),a,Numbers_UE_SubNet_Main);
    A_SubNet_Main{AP_mn} = exp(A_SubNet_Main{AP_mn}*1i);
    
    % Calculating CSI Inside Each Subnetwork
    for User_kn = 1:Numbers_UE_SubNet_Main
        % Defining UEs Beamsteering/Combining Matrix
        Delta_SubNet_Main{User_kn} = reshape(X( (User_kn-1)*(u*1)  + 1:(User_kn-1)*(u*1) + u*1),u,1);
        Delta_SubNet_Main{User_kn} = exp(Delta_SubNet_Main{User_kn}*1i);

        % Channel Matrix at A_i-->U_i Link
        H_SubNet{SubNet_n,AP_mn,User_kn}                     = H{Indeces_AP_SubNet_Main(AP_mn),Indeces_UE_SubNet_Main(User_kn)};
        % SVD of the A_i-->U_i Link Channel Matrix
        [U_i, S_i, V_i] = svd(H_SubNet{SubNet_n,AP_mn,User_kn});
        Eigen_User_i = diag(S_i);
        Index_NZ_Eigen = find(Eigen_User_i >= 10^(-3));   % This Value Defines the Degree of Annihilation for the subspace
        U_1_Matrix{AP_mn,User_kn} = U_i(:,Index_NZ_Eigen); % Left Non-annihilating Matrix
        V_1_Matrix{AP_mn,User_kn} = V_i(:,Index_NZ_Eigen); % Right Non-annihilating Matrix
        % Calculating Non-Orthogonal Projections at Users of the i'th Subnetwork
        Delta_1_Vector        = Delta_SubNet_Main{User_kn}'*U_1_Matrix{AP_mn,User_kn}*U_1_Matrix{AP_mn,User_kn}';
        A_1_Matrix            = V_1_Matrix{AP_mn,User_kn}*V_1_Matrix{AP_mn,User_kn}'*A_SubNet_Main{AP_mn};
        
        % Desired Part of Objective Function
        F_x_Non_Orthogonal = F_x_Non_Orthogonal + norm(Delta_1_Vector*S_i*A_1_Matrix,2); 
         
    end
 
    for User_kl_Outer = 1:Numbers_UE_SubNet_Outer
        % Channel Matrix at A_i-->U_i_Outer Link
        H_SubNet_Interference{SubNet_n,AP_mn,User_kl_Outer}  = H{Indeces_AP_SubNet_Main(AP_mn),Indeces_UE_SubNet_Outer(User_kl_Outer)};
        % SVD of the A_i-->U_i_Outer Link Channel Matrix
%         whos H_SubNet
        [U_Outer_kl, S_Outer_kl, V_Outer_kl] = svd(H_SubNet_Interference{SubNet_n,AP_mn,User_kl_Outer});
        Eigen_Outer_User_kl = diag(S_Outer_kl);
        Index_Z_Eigen = find(Eigen_Outer_User_kl <= 10^(-3));   % This Value Defines the Degree of Annihilation for the subspace
        U_0_Matrix{AP_mn,User_kl_Outer} = U_Outer_kl(:,Index_Z_Eigen); % Right Non-annihilating Matrix
        V_0_Matrix{AP_mn,User_kl_Outer} = V_Outer_kl(:,Index_Z_Eigen); % Right Non-annihilating Matrix
        % Calculating Non-Orthogonal Projections at Users of the i'th Subnetwork
        Delta_0_Vector    = Delta_All( (Indeces_UE_SubNet_Outer(User_kl_Outer)-1)*(u*1)+1: ((Indeces_UE_SubNet_Outer(User_kl_Outer)-1)*(u*1)) + (u*1) )'*...
                                                         U_0_Matrix{AP_mn,User_kl_Outer}*U_0_Matrix{AP_mn,User_kl_Outer}';
        A_0_Matrix        = V_0_Matrix{AP_mn,User_kl_Outer}*V_0_Matrix{AP_mn,User_kl_Outer}'*A_SubNet_Main{AP_mn};
        
        % Undesired Part of Objective Function
        F_x_Orthogonal = F_x_Orthogonal + norm(Delta_0_Vector*S_Outer_kl*A_0_Matrix,2);
  
    end
end
 
F_x = Alpha*F_x_Non_Orthogonal + (1-Alpha)*F_x_Orthogonal;

end



 