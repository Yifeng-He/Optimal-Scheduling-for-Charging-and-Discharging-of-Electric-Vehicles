function x=func_group(Current_T, Par_set, Current_EV, Current_slot, P_L_b_mic)

%     save T_Par_set.txt Par_set -ascii;
%     save T_Current_EV.txt Current_EV -ascii;
%     save T_Current_slot.txt Current_slot -ascii;
% clc;
% clear;
% close all;
% % test only
% load T_Par_set.txt;
% load T_Current_EV.txt;
% load T_Current_slot.txt;
% load T_P_L_b_mic.txt;
% Par_set=T_Par_set;
% Current_EV=T_Current_EV;
% Current_slot=T_Current_slot;
% P_L_b_mic=T_P_L_b_mic;
% Current_T=1;

% handle arguments: Par_set=[k_con, tau, Cap_battery_org, gamma, P_max, omega]';
k_con=Par_set(1);
tau=Par_set(2);
Cap_battery_org=Par_set(3);
gamma=Par_set(4);
P_max=Par_set(5);
omega=Par_set(6);
max_CHG_EV=Par_set(7); % the CHG EVs, from No. 1 to No. num_CHG_EV, the other are V2G EVs, BUG found on Dec 16, 2010.

Cap_battery=gamma*Cap_battery_org;

num_slot=length(Current_slot);
num_EV=length(Current_EV(:,1));

% determine the number of CHG EVs in the Current_EV matrix
num_CHG_EV=0;
for i=1:length(Current_EV(:,1))
    if Current_EV(i,7)==0
        num_CHG_EV=num_CHG_EV+1;
    end
end

L_b_mic=P_L_b_mic(Current_slot(1):Current_slot(num_slot));

% construction the F matrix
% Note: the time slots beginns from 1
F=zeros(num_EV, num_slot);
for i=1:num_EV
    for j=max(Current_T,Current_EV(i,4))-Current_T+1:Current_EV(i,5)-Current_T+1
        F(i,j)=1;
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  optimization (EV Charging) %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the equality constraint: Ax=b
% optimization variables x=[z1, z2, ..., z_24, x11, x12, ...., x_100,24]',a colum vextor
num_OptVar=1*num_slot+num_slot*num_EV;
b_a=L_b_mic; % the matrix for the first equality constraint
A1_a=zeros(num_slot, num_OptVar-1*num_slot);
A1=[eye(num_slot) A1_a];

A2_a=zeros(num_slot, num_OptVar-1*num_slot);
s_temp=0;
for i=1:num_slot
    for j=1:num_EV
        A2_a(i, (j-1)*num_slot+i)=F(j,i);
        % fprintf('Assign F(%g,%g)=%g, to A2_a(%g, %g).\n',j,i,F(j,i),i,(j-1)*num_slot+1);
        s_temp=s_temp+F(j,i);
    end
end
A2_b=zeros(num_slot, num_slot);
A2=[A2_b A2_a];
    
A_a=A1-A2;  % the matrix for the first equality constraint
clear A1 A2 A1_a A2_a A2_b;

% the second equality constarint
B_1=zeros(num_EV, num_OptVar-1*num_slot);
for i=1:num_EV
    B_1(i,(i-1)*num_slot+1:(i-1)*num_slot+num_slot)=F(i,:);
end
temp_1=zeros(num_EV, num_slot);
B1=[temp_1 B_1];    % the matrix for the second equality constraint
b_b=(Cap_battery/tau)*ones(num_EV,1)-Current_EV(:,3);% the matrix for the second equality constraint
clear  B_1  temp_1; 
    
% % combine the equality matrix
% Eq_left=[A_a' B1']';
% Eq_right=[b_a' b_b']';


% the equlity constraint ************
Eq_L=A_a;
Eq_R=b_a;

% the inequlity constraint ************
% 1) the first inequality
In_1=zeros(num_EV*num_slot, num_OptVar);
for i=1:num_slot
    for j=1:num_EV
        In_1((i-1)*num_EV+j,num_slot+(j-1)*num_slot+1:num_slot+(j-1)*num_slot+i)=F(j,1:i);
%         fprintf('set row %g, col %g:%g by using F(%g,1:%g).\n',(i-1)*num_EV+j,num_slot+(j-1)*num_slot+1,num_slot+(j-1)*num_slot+i,j,i);
    end
end
In_1=-1*In_1;  % the first inequality, left side
In_b1=zeros(num_EV*num_slot, 1);    % the first inequality, right side, [EV1_slot1, EV2_slot1, ..., EV1_slot2, EV2_slot2,...]'
for i=1:num_slot
    In_b1( (i-1)*num_EV+1:(i-1)*num_EV+num_EV, 1 )= (1/tau)*Current_EV(:,3); 
end

% 2) the second inequality
In_2=-1*In_1; % the second inequality, left side
In_b2=zeros(num_EV*num_slot, 1);    % the second inequality, right side, [EV1_slot1, EV2_slot1, ..., EV1_slot2, EV2_slot2,...]'
temp_b2=Cap_battery_org - Current_EV(:,3);
for i=1:num_slot
    In_b2( (i-1)*num_EV+1:(i-1)*num_EV+num_EV, 1 )= (1/tau)*temp_b2; 
end

% 3) the third inequality
In_3=-1*B1; % the third inequality, left side
In_b3=-1*b_b;    % the third inequality, right side, 

% combine all the inequality constraints
InEq_L=[In_1' In_2' In_3']';
InEq_R=[In_b1' In_b2' In_b3']';

% solve the optimization problem using quadratic programming: [x,f_obj] = quadprog(H,f,A,b,Aeq,beq,lb,ub) 
% the initial value
x_initial=zeros(num_OptVar,1);
x_initial(1:num_slot,1)=L_b_mic;
% lower and upper bound
x_lb=zeros(num_OptVar,1); % the lower bound for EV charging

% the lower bound for V2G
x_lb(num_slot+1:num_OptVar,1)=-1*P_max*ones(num_OptVar-num_slot,1);
% distinguish the CHG Evs and V2G EVs ( changed here, on Dec 16, 2010)
for i=1:length(Current_EV(:,1))
    if Current_EV(i,7)==0
        for j=1:num_slot
           x_lb(i*num_slot+j,1)=0; 
        end
    end
end

x_ub=P_max*ones(num_OptVar,1);
for i=1:num_slot
    x_ub(i,1)=2*omega;
end

% the objective function is changed to f=z_1^2+z_2^2+...
f_1=zeros(num_OptVar,1);
H_1=zeros(num_OptVar,num_OptVar);
for i=1:num_slot
    H_1(i,i)=2;
end
% solve the problem 
% [x,f_obj] = quadprog(H_1,f_1,[],[],Eq_left,Eq_right,x_lb,x_ub);% equality constraint only
[x,f_obj,exitflag,output] = quadprog(H_1,f_1,InEq_L,InEq_R,Eq_L,Eq_R,x_lb,x_ub); % equality + inequality constraint

% % the original objective function value
% obj_value=0;
% for i=1:num_slot
%    obj_value=obj_value + k_con* (x(i)^(theta+1) - L_b_mic(i)^(theta+1));
% end

% organize and verify the results
x_temp=x(num_slot+1:num_OptVar,1);
x_Matrix = reshape(x_temp,num_slot, num_EV);
x_Matrix=x_Matrix';  % charging rate matrix: num_EV * num_slot
% verify if all the EVs have been fully charged
E_Charged=zeros(num_EV,3);  % 1) initial energy, 2) charged energy, 3) final energy
E_Charged(:,1)=Current_EV(:,3);
for i=1:num_EV
    for j=1:num_slot
       E_Charged(i,2)= E_Charged(i,2) + x_Matrix(i,j)*F(i,j);
    end
    E_Charged(i,3)=E_Charged(i,1)+E_Charged(i,2);
    if abs(E_Charged(i,3)-Cap_battery) > 0.1
        fprintf('EV %g, final energy is %g, less than the battery capacity %g.\n',i,E_Charged(i,2),Cap_battery);
    end
end
% verify the charging load at each interval
Charged_Load=zeros(num_slot,6); % 1) the base load, 2) the charged load, 3) the total load (from EV rates), 4) the total load (variables z), 5) difference of 3) and 4)
Charged_Load(:,1)=L_b_mic; % the base load
for i=1:num_slot
    for j=1:num_EV
        Charged_Load(i,2)=Charged_Load(i,2)+x_Matrix(j,i)*F(j,i);
    end
    Charged_Load(i,3)=Charged_Load(i,1)+Charged_Load(i,2); % total load calculated from charged loads of individual EVs
    Charged_Load(i,4)=x(i,1); % total load calculated from optimization variables z_i
    Charged_Load(i,5)=Charged_Load(i,4)-Charged_Load(i,3); % the difference of the results
    Charged_Load(i,6)=Charged_Load(i,5)/Charged_Load(i,1); % the relative error
%     if abs(Total_Load(i,1)-Total_Load(i,2)) > 0.1
%         fprintf('At charging interval %g, the total load is different, %g, vs %g.\n',i,Total_Load(i,1),Total_Load(i,2));
%     end
end

% verify the energy level at each interval for each EV
% the evolution of energy level for each EV
Energy_variation=zeros(num_EV,num_slot+1);
for i=1:num_EV
    Energy_variation(i,1)=E_Charged(i,1);  % the initial energy
    for j=1:num_slot
        Energy_variation(i,j+1)=Energy_variation(i,j)+F(i,j)*x_Matrix(i,j);
    end
end

% returen the results
x=[];
x=x_Matrix; % the charging rate at each interval
