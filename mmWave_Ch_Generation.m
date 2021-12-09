function H = mmWave_Ch_Generation(M,a,K,u,L)
% clc; clear all; close all;
% M = 5; a = 4; K = 3; u =9;  L = 3;

% mmWave Frequence and Wavelength
f = 24*10^(9); Lambda = 3*10^(8)/f;
d = Lambda;
% Array Parameters
L_U_1 = sqrt(u); L_U_2 = L_U_1;
L_A_1 = sqrt(a); L_A_2 = L_A_1;


[L_U_m_1,L_U_m_2] = meshgrid(1:L_U_1,1:L_U_2);
[L_A_m_1,L_A_m_2] = meshgrid(1:L_A_1,1:L_A_2);
c_U=cat(2,L_U_m_1',L_U_m_2');
c_A=cat(2,L_A_m_1',L_A_m_2');
u_v_U =reshape(c_U,[],2);
u_v_A =reshape(c_A,[],2);
U_V_U = u_v_U' - 1; 
U_V_A = u_v_A' - 1;
% index_x = 1:u:u*L; index_y = 1:a:a*L;
for k = 1: K
    for m = 1:M
        % mmWave Channel Parameters
        Theta_A = pi.*rand(L,1);           % Elevation Angle [0   Pi]
        Phi_A = -pi + (pi+pi).*rand(L,M);  % Azimuth Angle   [-Pi Pi]
        Theta_U = pi.*rand(L,1);           % Elevation Angle [0   Pi]
        Phi_U = -pi + (pi+pi).*rand(L,M);  % Azimuth Angle   [-Pi Pi]
        b_A = exp(d*2*pi*1i*(U_V_A(1,:).*sin(Theta_A(:,1)).*cos(Phi_A(:,1))+...
            U_V_A(2,:).*sin(Phi_A(:,1))) /(Lambda));
        b_U = exp(d*2*pi*1i*(U_V_U(1,:).*sin(Theta_U(:,1)).*cos(Phi_U(:,1))+...
            U_V_U(2,:).*sin(Phi_U(:,1))) /(Lambda));
        %     B = b_U * b_A';
        % Extracting Phase Matrices Per Signal Path
        H{m,k} = zeros(u,a);
        for i = 1:L
            B = b_U(i,:)'*b_A(i,:);
            alpha(i)  = randn(1,1)+randn(1,1)*1i;
            H{m,k}  = H{m,k} + alpha(i)*B;
        end
        Kappa  = alpha(1)^2/sum(alpha(2:end).^2);
        H{m,k} = sqrt(1/(Kappa + L - 1))*H{m,k}; 
    end
end


end
