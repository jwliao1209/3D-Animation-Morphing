clear;clc;close all;

addpath('./function');
mkdir data
FileName = 'Horse';
load(fullfile('data', [FileName '.mat']));

figure(1);
Tri.Surf(F, V);

%% Define indices of landmarks
Back = [122;299;647;107;926;2775;6725;9705;10645;11268;11516;11641;11935;12129;12306;12429;12319;12537;12601;12927;12807;12958;12980;12833];
legLF = [5529; 6105; 6572;];
legRF = [18055; 17888; 17569;];
legLB = [5920; 5931; 6275;];
legRB = [19376; 18837; 18564;];

V_Back = V(Back,:);
V_legLF = V(legLF,:);
V_legRF = V(legRF,:);
V_legLB = V(legLB,:);
V_legRB = V(legRB,:);

landmark = [Back; legLF; legRF; legLB; legRB];

BackTheta = (pi/2) * (-1:6);
s_theta = pi/6;
Theta_main_LF = (pi/4) * [0,1,2,3,4,5,6,7]; Theta_small_LF = s_theta * [0,1,2,1,0,0,0,0];
Theta_main_RF = (pi/4) * [4,5,6,7,0,1,2,3]; Theta_small_RF = s_theta * [0,0,0,0,0,1,2,1];
Theta_main_LB = (pi/4) * [2,3,4,5,6,7,0,1]; Theta_small_LB = s_theta * [2,1,0,0,0,0,0,1];
Theta_main_RB = (pi/4) * [6,7,0,1,2,3,4,5]; Theta_small_RB = s_theta * [0,0,0,0,1,2,1,0];

D_Back = -0.03 * cos(BackTheta);

d_main_leg = 0.1;
D_main_LF_x = d_main_leg * cos(Theta_main_LF); 
D_main_LF_y = d_main_leg * max(sin(Theta_main_LF), 0);
D_main_RF_x = d_main_leg * cos(Theta_main_RF);  
D_main_RF_y = d_main_leg * max(sin(Theta_main_RF), 0);
D_main_LB_x = d_main_leg * cos(Theta_main_LB);  
D_main_LB_y = d_main_leg * max(sin(Theta_main_LB), 0);
D_main_RB_x = d_main_leg * cos(Theta_main_RB);  
D_main_RB_y = d_main_leg * max(sin(Theta_main_RB), 0);

d_small_leg = abs(V_legLF(3,3)-V_legLF(1,3));
D_small_LF_x = d_small_leg * ((Theta_small_LF > 0) .* cos(Theta_small_LF)); 
D_small_LF_y = d_small_leg * ((Theta_small_LF > 0) .* sin(Theta_small_LF));
D_small_RF_x = d_small_leg * ((Theta_small_RF > 0) .* cos(Theta_small_RF)); 
D_small_RF_y = d_small_leg * ((Theta_small_RF > 0) .* sin(Theta_small_RF));
D_small_LB_x = d_small_leg * ((Theta_small_LB > 0) .* cos(Theta_small_LB)); 
D_small_LB_y = d_small_leg * ((Theta_small_LB > 0) .* sin(Theta_small_LB));
D_small_RB_x = d_small_leg * ((Theta_small_RB > 0) .* cos(Theta_small_RB)); 
D_small_RB_y = d_small_leg * ((Theta_small_RB > 0) .* sin(Theta_small_RB));

%% Define target points
for N_Frame = 1 : 8
    load(fullfile('data', [FileName '.mat']));
    V_Back = V(Back,:);
    V_legLF = V(legLF,:);
    V_legRF = V(legRF,:);
    V_legLB = V(legLB,:);
    V_legRB = V(legRB,:);

    d_back = D_Back(N_Frame);
    V_Back(1:6,2) = V_Back(1:6,2) + d_back; % Head has the same placement
    V_Back(7:11,2) = V_Back(7:11,2) + (2:-1:-2)' * d_back / 3;
    V_Back(12,2) = V_Back(12,2) - d_back;
    V_Back(13:23,2) = V_Back(13:23,2) + (-5:1:5)' * d_back / 5;
    V_Back(end,2) = V_Back(end,2) + d_back;
    
    V_legLF(:,3) = V_legLF(:,3) + D_main_LF_x(N_Frame);
    V_legLB(:,3) = V_legLB(:,3) + D_main_LB_x(N_Frame);
    V_legRF(:,3) = V_legRF(:,3) + D_main_RF_x(N_Frame);
    V_legRB(:,3) = V_legRB(:,3) + D_main_RB_x(N_Frame);
    
    V_legLF(:,2) = V_legLF(:,2) + D_main_LF_y(N_Frame);
    V_legLB(:,2) = V_legLB(:,2) + D_main_LB_y(N_Frame);
    V_legRF(:,2) = V_legRF(:,2) + D_main_RF_y(N_Frame);
    V_legRB(:,2) = V_legRB(:,2) + D_main_RB_y(N_Frame);
    
    if (Theta_small_LF(N_Frame) > 0)
        V_legLF(1,3) = V_legLF(1,3) - D_small_LF_x(N_Frame);
        V_legLF(2,3) = V_legLF(3,3) + D_small_LF_x(N_Frame);
        V_legLF(1,2) = V_legLF(1,2) - D_small_LF_y(N_Frame);
        V_legLF(2,2) = V_legLF(3,2) + D_small_LF_y(N_Frame);
    end
    if (Theta_small_RF(N_Frame) > 0)
        V_legRF(1,3) = V_legRF(1,3) - D_small_RF_x(N_Frame);
        V_legRF(2,3) = V_legRF(3,3) + D_small_RF_x(N_Frame);
        V_legRF(1,2) = V_legRF(1,2) - D_small_RF_y(N_Frame);
        V_legRF(2,2) = V_legRF(3,2) + D_small_RF_y(N_Frame);
    end
    if (Theta_small_LB(N_Frame) > 0)
        V_legLB(1,3) = V_legLB(1,3) - D_small_LB_x(N_Frame);
        V_legLB(2,3) = V_legLB(3,3) + D_small_LB_x(N_Frame);
        V_legLB(1,2) = V_legLB(1,2) - D_small_LB_y(N_Frame);
        V_legLB(2,2) = V_legLB(3,2) + D_small_LB_y(N_Frame);
    end
    if (Theta_small_RB(N_Frame) > 0)
        V_legRB(1,3) = V_legRB(1,3) - D_small_RB_x(N_Frame);
        V_legRB(2,3) = V_legRB(3,3) + D_small_RB_x(N_Frame);
        V_legRB(1,2) = V_legRB(1,2) - D_small_RB_y(N_Frame);
        V_legRB(2,2) = V_legRB(3,2) + D_small_RB_y(N_Frame);
    end
    
    % Perform Mesh Editing
    V_landmark = [V_Back; V_legLF; V_legRF; V_legLB; V_legRB];
    Vnew = MeshEditing(F, V, landmark, V_landmark);
    
    Tmp = V;
    V = Vnew;
    save(fullfile('data', ['Horse' num2str(N_Frame)]), 'V', 'F');
    V = Tmp;

    % Rotation for a better visualization
    theta = -pi/2;
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    V(:,[1,3]) = V(:,[1,3])*R.';
    Vnew(:,[1,3]) = Vnew(:,[1,3])*R.';

    figure(2)
    subplot(1,2,1)
    Tri.Surf(F, V);
    subplot(1,2,2)
    Tri.Surf(F, Vnew);
end