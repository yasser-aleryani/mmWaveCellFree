% This function calculate the objective function related to digital
% beamforming problem**************Try Proportional Fairness Function
function [R_Objective] = Obj_Fun_Beamforming_DRL(X,P,Numbers_UE_SubNet_Main,...
    Indeces_AP_SubNet_Main,Indeces_UE_SubNet_Main,SubNet_n,Clusters_APs,Count_AP_Configs,Clusters_UEs,Count_UE_Configs,Number_SubNets)
global h_Equivalent_All W_All;
 
% Ordering Users Within the i'th Subnetwork
if Numbers_UE_SubNet_Main == 1
    Indeces_Ord_UE_SubNet_Main = Indeces_UE_SubNet_Main;
else
    for User_kn = Indeces_UE_SubNet_Main
        Norm_User_kn = 0;
        Count_User_kn = 1;
        for AP_mn = Indeces_AP_SubNet_Main
            Norm_User_kn(Count_User_kn) = Norm_User_kn + norm(h_Equivalent_All{AP_mn,User_kn},2);
        end
        Count_User_kn = Count_User_kn + 1;
    end
    [~,Norm_Ord_Users] = sort(Norm_User_kn,'ascend');
    Indeces_Ord_UE_SubNet_Main = Indeces_UE_SubNet_Main(Norm_Ord_Users);
end
 
Count_kn_mn = 1;
Count_in_mn = 1;
for User_kn = Indeces_Ord_UE_SubNet_Main
    
    Desired = 0; IUI =0; ISNI = 0;
    Indeces_Ord_SIC_kn = Indeces_Ord_UE_SubNet_Main;
    [~,kn] = find(Indeces_Ord_SIC_kn == User_kn);
    for AP_mn = Indeces_AP_SubNet_Main
        
        % Calculating the Numerator of the SINR
        w_kn_mn = reshape(X( (Count_kn_mn-1)*Numbers_UE_SubNet_Main + 1:...
            (Count_kn_mn-1)*Numbers_UE_SubNet_Main + Numbers_UE_SubNet_Main),Numbers_UE_SubNet_Main,1);
        Desired = Desired + norm(h_Equivalent_All{AP_mn,User_kn}.^2*w_kn_mn,1);
        Count_kn_mn = Count_kn_mn +1;
        % Calculating the IUI at the kn'th User
        if kn == Numbers_UE_SubNet_Main
            IUI = 0;
        else
            for User_Ord_kn = Indeces_Ord_SIC_kn(kn+1:end)
                w_in_mn = reshape(X( (Count_in_mn-1)*Numbers_UE_SubNet_Main + 1:...
                    (Count_in_mn-1)*Numbers_UE_SubNet_Main + Numbers_UE_SubNet_Main),Numbers_UE_SubNet_Main,1);
                IUI = IUI + norm(h_Equivalent_All{AP_mn,User_Ord_kn}.^2*w_in_mn,1);
                
                Count_in_mn = Count_in_mn +1;
            end
        end
        
    end
        Indeces_SubNets_Outer           = 1:Number_SubNets;
        Indeces_SubNets_Outer(SubNet_n) = [];
    for SubNet_l           = Indeces_SubNets_Outer
 
        Indeces_AP_SubNet_l         = cell2mat(Clusters_APs{Count_AP_Configs,1}(SubNet_l));
        Indeces_UE_SubNet_l         = cell2mat(Clusters_UEs{Count_UE_Configs,1}(SubNet_l));

        for AP_ml = Indeces_AP_SubNet_l
            Count_kl = 1;
            for User_kl = Indeces_UE_SubNet_l
                ISNI = ISNI + norm(h_Equivalent_All{AP_ml,User_kl}.^2*W_All{AP_ml}(:,Count_kl),1);
                Count_kl = Count_kl + 1;
            end
            
        end
    end
 
    Gamma_kn(User_kn) =  Desired/(IUI+ISNI+1/P);
    R_kn(User_kn)     =  log2(1 + Gamma_kn(User_kn)); 
    end  



R_Objective = -sum(R_kn);

end

