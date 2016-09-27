close all;
clc;
clear;
% cd('C:\Program Files\MATLAB\cvx');
% cvx_setup;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  the parameters      %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the base load vector for Toronto (data: Toronto 21/08/2009, temperature avg/max/min=22/25/20, unit=KW)
L_b=[
1648181.19
1510114.92
1404600.51
1348140.62
1355294.53
1380434.41
1472398.73
1619859.88
1774726.71
1891739.09
2024589.36
2154006.65
2214969.44
2248348.3
2257444.05
2253898.58
2251095.52
2240891.11
2187658.11
2060455.55
2032531.95
1984614.47
1846244.76
1676395.15];

% predicted base load (mean relative error =0.089)
P_L_b1=[
1737223.863
1603267.075
1501043.33
1434323.708
1427034.265
1445808.13
1535646.84
1703557
1878913.403
2041190.428
2189203.243
2300813.958
2376838.56
2431838.998
2456820.033
2472628.88
2497750.498
2502256.753
2487853.418
2333072.535
2313458.498
2278059.543
2121749.573
1908084.005];


% predicted base load (mean relative error =0.0414)
P_L_b2=[
1599948.91
1468345.518
1370502.199
1310946.046
1308525.559
1330298.288
1417044.254
1592177.276
1774327.31
1935791.394
2081450.673
2194127.508
2266347.528
2310264.895
2317806.969
2328853.886
2350763.395
2356239.943
2340220.599
2215920.234
2192252.838
2168526.739
2014122.374
1803843.804];

% second time predicted base load (better prediction)(mean relative error =0.0234)
P_L_b3=[
1648037.519
1458944.669
1366366.043
1310276.476
1309105.99
1329282.055
1411672.493
1587762.671
1728429.463
1888917.875
2032770.183
2133804.711
2199858.463
2240117.111
2241852.165
2245048.574
2259487.69
2260142.395
2243557.341
2145562.018
2122157.145
2104126.554
1948760.596
1740961.875];

% choose which forecased load is used
P_L_b=P_L_b2;

% the base load for a micogrid
Scale_factor=1/1500;
L_b_mic=L_b*Scale_factor;  % the base load
P_L_b_mic=P_L_b*Scale_factor;  % the predicted base load

% the price model
% alfa=1.0;
% theta=1;
omega=1.2*max(L_b_mic);
k_0=0.0001;
k_1=0.00012;
k_2=0;
% the life reduction cost
beta=0.001;
beta=0;

% the constant multiplication
% k_con=alfa/(omega^theta*(theta+1));

% the lengh of charging interval
tau=1; % in hours

% number of charging intervals
num_slot=length(L_b_mic);

% the basic price
price_basic=zeros(num_slot,1); % the price based on basic load
for i=1:num_slot
   price_basic(i)=k_0+k_1* L_b_mic(i);
end
fprintf('The price, min price=%g, max price=%g.\n',min(price_basic), max(price_basic));


% EV battery capacity
Cap_battery_org=16; % in KWh
gamma=0.9; % the percentage of battery when charging is completed
Cap_battery=gamma*Cap_battery_org;

% the max charging rate
P_max=5; % in KW
% number of EVs
num_EV=200;

% the percentage of EVs that charge the battery only
P_Chg=0;
% number of CHG EVs
num_CHG_EV=round(P_Chg*num_EV);  % the CHG EVs will be in the front part of the EV info matrix.
% number of V2G EVs
num_V2G_EV=num_EV-num_CHG_EV;

% the EV charging pattern
%************************************************
% 30% of the EVs connect to charging stations before interval 1, and others uniformly distributed 
%************************************************
% % % % EV info matrix: 1) arrival, 2) departure time, 3) initial energy, 4)
% % % % charging period, 5) min charging time 
% 
% EV_info=zeros(num_EV,3);
% % percentage of EVS connected to station before interval 1
% Per_EV=0.1;
% % other EVs arrival times are uniformly distributed between [1, 20]
% for i=1:num_EV
%     temp_00=rand;
%     if temp_00<=Per_EV
%         T_arrival(i,1)=1;
%     else
%         T_arrival(i,1)=round(1 + (20-1).*rand);
%     end
% end
% 
% % the charging periods are uniformly distributed between [4, 12] hours
% T_charging= round(4 + (12-4).*rand(num_EV,1));
% T_charging=-1*sort(-1*T_charging);
% 
% % the departure time
% for i=1:num_EV
%     T_departure(i,1)=min(24, T_arrival(i,1)+T_charging(i,1));
% end
% % the initial energy is uniformly distributed between [0 0.8] of battery capacity
% Ini_percentage=0+ (0.8-0).*rand(num_EV,1);
% % fill the EV_info
% EV_info(:,1)=T_arrival;
% EV_info(:,2)=T_departure;
% EV_info(:,3)=Cap_battery_org*Ini_percentage;
% 
% for i=1:num_EV
%     EV_info(i,4)=EV_info(i,2)-EV_info(i,1)+1; % charging period
%     EV_info(i,5)=EV_info(i,3)/P_max; % min charging time
%     if EV_info(i,4) < EV_info(i,5)
%         fprintf('EV %g charging time is not valis.\n',i);
%     end
% end
% 
% % save and load EV_info
% save EV_info.txt EV_info -ascii;

load EV_info.txt;
EV_info=EV_info(1:num_EV,:);

% the relationship between EV and charging interval
F=zeros(num_EV, num_slot);
G=ones(num_EV, num_slot);
for i=1:num_EV
    for j=EV_info(i,1):EV_info(i,2)
        F(i,j)=1;
        G(i,j)=0;
    end
end
F1=reshape(F',1,[]);
% F=ones(num_EV, num_slot);

% plot the base load
xx_1=1:num_slot;
figure;
yy_1(:,1)=L_b_mic;
yy_1(:,2)=P_L_b_mic;
plot(xx_1,yy_1);
ylabel('Load [KW]');
xlabel('Hour No.');
legend('Real base load','Predicted base load',1);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Global optimal scheme for V2G using CVX tool
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
b_b=(Cap_battery/tau)*ones(num_EV,1)-EV_info(:,3);% the matrix for the second equality constraint
clear  B_1  temp_1; 
    
% % combine the equality matrix
% Eq_left=[A_a' B1']';
% Eq_right=[b_a' b_b']';


% the equlity constraint ************
Eq_L=A_a;
Eq_R=b_a;
clear  A_a  b_a; 
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
    In_b1( (i-1)*num_EV+1:(i-1)*num_EV+num_EV, 1 )= (1/tau)*EV_info(1:num_EV,3); 
end

% 2) the second inequality
In_2=-1*In_1; % the second inequality, left side
In_b2=zeros(num_EV*num_slot, 1);    % the second inequality, right side, [EV1_slot1, EV2_slot1, ..., EV1_slot2, EV2_slot2,...]'
temp_b2=Cap_battery_org - EV_info(1:num_EV,3);
for i=1:num_slot
    In_b2( (i-1)*num_EV+1:(i-1)*num_EV+num_EV, 1 )= (1/tau)*temp_b2; 
end

% 3) the third inequality
In_3=-1*B1; % the third inequality, left side
In_b3=-1*b_b;    % the third inequality, right side, 

% combine all the inequality constraints
% InEq_L=[In_1' In_2' In_3']';
% InEq_R=[In_b1' In_b2' In_b3']';

% solve the optimization problem using quadratic programming: [x,f_obj] = quadprog(H,f,A,b,Aeq,beq,lb,ub) 
% the initial value
x_initial=zeros(num_OptVar,1);
x_initial(1:num_slot,1)=L_b_mic;
% lower and upper bound
x_lb=zeros(num_OptVar,1); % the lower bound for EV charging

% the lower bound for V2G
x_lb(num_slot+1:num_OptVar,1)=-1*P_max*ones(num_OptVar-num_slot,1);
% distinguish the CHG Evs and V2G EVs
for i=1:num_CHG_EV*num_slot
    x_lb(num_slot+i,1)=0;
end

x_ub=P_max*ones(num_OptVar,1);
for i=1:num_slot
    x_ub(i,1)=2*omega;
end

% the objective function is changed to f=z_1^2+z_2^2+...
% f_1=zeros(num_OptVar,1);
% H_1=zeros(num_OptVar,num_OptVar);
% for i=1:num_slot
%     H_1(i,i)=2;
% end

% *****************************************************************
%******** solve the optimization problem using CVX tool
% *****************************************************************

% quadratic program
disp('Computing the solution of the quadratic program...');
cvx_begin
    variable v_x(num_OptVar);
    minimize(  k_0*sum(pow_p(v_x(1:num_slot),1)) + (k_1/2)*sum(pow_p(v_x(1:num_slot),2)) + beta*sum(square(F1)*square(v_x(num_slot+1:num_OptVar))) ...
        -k_0*sum(pow_p(L_b_mic(1:num_slot),1)) - (k_1/2)*sum(pow_p(L_b_mic(1:num_slot),2)) )
    %     minimize( quad_form(v_x,H_1) )
    Eq_L * v_x == Eq_R;
    In_1 * v_x <= In_b1;
    In_2 * v_x <= In_b2;
    In_3 * v_x <= In_b3;
%     InEq_L * v_x <= InEq_R;
    v_x >= x_lb;
    v_x <= x_ub;
cvx_end

clear  Eq_L  In_1 In_2 In_3 Eq_R In_b1 In_b2 In_b3 x_lb x_ub; 

% [x,f_obj] = quadprog(H_1,f_1,[],[],Eq_left,Eq_right,x_lb,x_ub);% equality constraint only
% [v_x2,f_obj,exitflag,output] = quadprog(H_1,f_1,InEq_L,InEq_R,Eq_L,Eq_R,x_lb,x_ub); % equality + inequality constraint


% organize and verify the results
v_x_temp=v_x(num_slot+1:num_OptVar,1);
v_x_Matrix = reshape(v_x_temp,num_slot, num_EV);
v_x_Matrix=v_x_Matrix';  % charging rate matrix: num_EV * num_slot
% verify if all the EVs have been fully charged
v_E_Charged=zeros(num_EV,3);  % 1) initial energy, 2) charged energy, 3) final energy
v_E_Charged(:,1)=EV_info(1:num_EV,3);
for i=1:num_EV
    for j=1:num_slot
       v_E_Charged(i,2)= v_E_Charged(i,2) + v_x_Matrix(i,j)*F(i,j);
    end
    v_E_Charged(i,3)=v_E_Charged(i,1)+v_E_Charged(i,2);
    if abs(v_E_Charged(i,3)-Cap_battery) > 0.1
        fprintf('EV %g, final energy is %g, less than the battery capacity %g.\n',i,v_E_Charged(i,2),Cap_battery);
    end
end
% verify the charging load at each interval
v_Charged_Load=zeros(num_slot,6); % 1) the base load, 2) the charged load, 3) the total load (from EV rates), 4) the total load (variables z), 5) difference of 3) and 4)
v_Charged_Load(:,1)=L_b_mic; % the base load
for i=1:num_slot
    for j=1:num_EV
        v_Charged_Load(i,2)=v_Charged_Load(i,2)+v_x_Matrix(j,i)*F(j,i);
    end
    v_Charged_Load(i,3)=v_Charged_Load(i,1)+v_Charged_Load(i,2); % total load calculated from charged loads of individual EVs
    v_Charged_Load(i,4)=v_x(i,1); % total load calculated from optimization variables z_i
    v_Charged_Load(i,5)=v_Charged_Load(i,4)-v_Charged_Load(i,3); % the difference of the results
    v_Charged_Load(i,6)=v_Charged_Load(i,5)/v_Charged_Load(i,1); % the relative error
%     if abs(Total_Load(i,1)-Total_Load(i,2)) > 0.1
%         fprintf('At charging interval %g, the total load is different, %g, vs %g.\n',i,Total_Load(i,1),Total_Load(i,2));
%     end
end

% verify the energy level at each interval for each EV
% the evolution of energy level for each EV
v_Energy_variation=zeros(num_EV,num_slot+1);
for i=1:num_EV
    v_Energy_variation(i,1)=v_E_Charged(i,1);  % the initial energy
    for j=1:num_slot
        v_Energy_variation(i,j+1)=v_Energy_variation(i,j)+F(i,j)*v_x_Matrix(i,j);
    end
end

% save the charging rate for V2G
% save h_x_Matrix.txt v_x_Matrix -ascii;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  The discrete-time optimization loop
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_Charged=zeros(num_EV,3);  % 1) initial energy, 2) charged energy, 3) final energy
E_Charged(:,1)=EV_info(1:num_EV,3);

% monitor the status of all EVs
Rate_Matrix=zeros(num_EV, num_slot); % the actual charging rates for each EV

% EV_Matrix: col. 1=EV No., col2=initial energy, col3=current energy,
% col4=starting slot No., col5=ending slot No., col6=status (0-not started,
% 1-ongoing, 2-completed), col.7=CHG or V2G (0=CHG, 1=V2G), col. 8=group ID(1, 2, 3...)
EV_Matrix=zeros(num_EV,6); 
EV_Matrix(:,1)=1:num_EV; %EV No.
EV_Matrix(:,2)=EV_info(:,3); %initial energy
EV_Matrix(:,3)=EV_info(:,3); %current energy
EV_Matrix(:,4)=EV_info(:,1); %starting slot No.
EV_Matrix(:,5)=EV_info(:,2); %ending slot No.
for i=1:length(EV_Matrix(:,1))
    if EV_Matrix(i,1)<=num_CHG_EV
        EV_Matrix(i,7)=0;  % CHG
    else
        EV_Matrix(i,7)=1; % V2G
    end
end
        
%%%%%%%%% The group level **********************
% the group size
group_size=100; % number of EVs in a group
% the group number
group_num=ceil(num_EV/group_size);
% shuffle the EVs
temp_ID=EV_Matrix(:,1);
temp_No=shuffle(temp_ID);
% % save the EV_ID
% % save temp_No.txt temp_No -ascii;
% 
% load temp_No.txt;

for i=1:length(temp_No)
    EV_Matrix_1(i,:)=EV_Matrix(temp_No(i),:);
end
% identify the group no.
EV_Matrix_1(:,8)=group_num*ones(num_EV,1);
for i=1:group_num-1
    EV_Matrix_1((i-1)*group_size+1:i*group_size,8)=i*ones(group_size,1);
end
EV_Matrix=[];
    
% iteration from groups
for gg=1:group_num
    
    EV_Matrix=[];
    % constrcut the matrix EV_Matrix
    temp_ss=0;
    for i=1:num_EV
        if EV_Matrix_1(i,8)==gg
            temp_ss=temp_ss+1;
            EV_Matrix(temp_ss,:)=EV_Matrix_1(i,:);
        end
    end
    EV_Matrix=sortrows(EV_Matrix,7);
    % number of EVs in the current group
    num_EV_g=length(EV_Matrix(:,1));
    % the number of CHG EVs
    num_CHG_g=0;
    for i=1:num_EV_g
        if EV_Matrix(i,7)==0
           num_CHG_g=num_CHG_g+1;
        end
    end
    % EV ID mapping
    EV_Mapping=zeros(num_EV_g,2);
    EV_Mapping(:,1)=1:num_EV_g;
    EV_Mapping(:,2)=EV_Matrix(:,1);
    % change the EV IDs
    EV_Matrix(:,1)=1:num_EV_g;
    % the rate matrix for the group
    Rate_Matrix_g=zeros(num_EV_g, num_slot); 

    
    for iii=1:num_slot

        Current_T=iii;  % the current time instant

        % prepare the parameters which will be input to the optimization
        % function, Par_set=[k_con, tau, Cap_battery_org, gamma, P_max, omega, num_CHG_g, theta]';
        Par_set=[beta, tau, Cap_battery_org, gamma, P_max, omega, num_CHG_g, k_0, k_1]';

        % find out the charging EVs at the current time
        Current_EV=[];
        temp_d1=0;
        for i=1:num_EV_g
            if Current_T >= EV_Matrix(i,4) & EV_Matrix(i,6)<2
                temp_d1=temp_d1+1;
                Current_EV(temp_d1,:)=EV_Matrix(i,:);
                % fprintf('EV %g is put into Current_EV set.\n',i);
            end
        end


        % save data for test
    %     save T_Par_set.txt Par_set -ascii;
    %     save T_Current_EV.txt Current_EV -ascii;
    %     save T_Current_slot.txt Current_slot -ascii;
        if isempty(Current_EV)
            fprintf('At time slot %g, the Current_EV set is empty.\n',Current_T);
            continue;
        else
            % the sliding window
            Current_slot=[];
            Current_slot=Current_T: max(Current_EV(:,5));
            Current_slot=Current_slot';
            % the optimization function, return charging rates during the sliding window
            % % case 1): EV charging only 
    %         x=func_EV_Chg(Current_T, Par_set, Current_EV, Current_slot, P_L_b_mic);
        %     % % case 2): V2G
            x=func_group_cvx(Current_T, Par_set, Current_EV, Current_slot, P_L_b_mic);  
        end

        % update the EV status: Rate_Matrix_g
        for i=1:length(Current_EV(:,1))
            Rate_Matrix_g(Current_EV(i,1),Current_T:Current_T+length(x(1,:))-1)=x(i,:);
        end
        % update the EV status: EV_Matrix
        for i=1:length(Current_EV(:,1))
            EV_Matrix(Current_EV(i,1),3)=EV_Matrix(Current_EV(i,1),3)+Rate_Matrix_g(Current_EV(i,1),Current_T); % update the current energy
            if abs(EV_Matrix(Current_EV(i,1),3)-Cap_battery)<0.01
                EV_Matrix(Current_EV(i,1),6)=2;
                elseif Current_T >= EV_Matrix(Current_EV(i,1),4)
                    EV_Matrix(Current_EV(i,1),6)=1;
            end
            % fprintf('EV %g is updated.\n',EV_Matrix(Current_EV(i,1),1));
        end


    end  % end of iii iteration
 
    % verify the final energy for each EV
    for i=1:num_EV_g
        E_Charged(EV_Mapping(i,2),2)=sum(Rate_Matrix_g(i,:)); % charged energy
        E_Charged(EV_Mapping(i,2),3)=E_Charged(EV_Mapping(i,2),1)+E_Charged(EV_Mapping(i,2),2);
    end

    % update the global rate matrix
    for i=1:num_EV_g
        Rate_Matrix(EV_Mapping(i,2),:)=Rate_Matrix_g(EV_Mapping(i,1),:);
        %fprintf('Group %g: EV %g is completed.\n',gg, EV_Mapping(i,2));
    end

    
end  % end of gg iteration

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Verify the results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[];
x_Matrix=[];
x_Matrix=Rate_Matrix;  % charging rate matrix: num_EV * num_slot
% verify the charging load at each interval
Charged_Load=zeros(num_slot,6); % 1) the base load, 2) the charged load, 3) the total load (from EV rates), 4) the total load (variables z), 5) difference of 3) and 4)
Charged_Load(:,1)=L_b_mic; % the base load
for i=1:num_slot
    for j=1:num_EV
        Charged_Load(i,2)=Charged_Load(i,2)+x_Matrix(j,i)*F(j,i);
    end
    Charged_Load(i,3)=Charged_Load(i,1)+Charged_Load(i,2); % total load calculated from charged loads of individual EVs
%     Charged_Load(i,4)=x(i,1); % total load calculated from optimization variables z_i
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
        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Naive charging scheme without discharging %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_inital_en(:,1)=E_Charged(:,1); % the initial energy
% N_inital_en(:,2)=Cap_battery-N_inital_en(:,1); % the energy requied to be charged
% N_x_Matrix=zeros(num_EV,num_slot);
% % the charging rate allocation
% for i=1:num_EV
%     temp_11(i,1)=N_inital_en(i,1);  % accumulated energy
%     for j=1: num_slot
%         if F(i,j)==1
%             temp_11(i,1)=temp_11(i,1)+P_max;
%             if temp_11(i,1)<=Cap_battery
%                 N_x_Matrix(i,j)=P_max;
%             else
%                 N_x_Matrix(i,j)=Cap_battery-(temp_11(i,1)-P_max);
%                 break;
%             end
%         end
%     end
% end
% % the evolution of energy level for each EV
% N_Energy_variation=zeros(num_EV,num_slot+1);
% for i=1:num_EV
%     N_Energy_variation(i,1)=N_inital_en(i,1);  % the initial energy
%     for j=1:num_slot
%         N_Energy_variation(i,j+1)=N_Energy_variation(i,j)+F(i,j)*N_x_Matrix(i,j);
%     end
% end             
% % the charged load
% N_Charged_Load=zeros(num_slot,3); % 1) the base load, 2) the charged load, 3) the total load (from EV rates), 
% N_Charged_Load(:,1)=L_b_mic; % the base load
% for i=1:num_slot
%     for j=1:num_EV
%         N_Charged_Load(i,2)=N_Charged_Load(i,2)+N_x_Matrix(j,i)*F(j,i);
%     end
%     N_Charged_Load(i,3)=N_Charged_Load(i,1)+N_Charged_Load(i,2); % total load calculated from charged loads of individual EVs
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Naive charging scheme with discharging %%% revised on June 03,
%%%%%%%%%%  2011
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_inital_en(:,1)=E_Charged(:,1); % the initial energy
N_inital_en(:,2)=Cap_battery-N_inital_en(:,1); % the energy requied to be charged
H_x_Matrix=zeros(num_EV,num_slot);

% naive method: constant chg/dichg rate at each interval, the number of
% dicharging intervals are pre-dertermined
% allocate chg/dischg rate
num_dischg_slot=1;
for i=1:num_CHG_EV % for the CHG-only EVs
    % number of parking slots
    num_p_slots(i,1)=sum(F(i,:));
    % the constant rate
    C_rate(i,1)=N_inital_en(i,2)/num_p_slots(i,1);
    if C_rate(i,1)>P_max
        fprintf('At EV %g, the constant rate %g exceeded the max rate %g.\n', i, C_rate(i,1), P_max);
    end
    % rate allocation
    for j=1:num_slot
        if F(i,j)==1
             H_x_Matrix(i,j)=C_rate(i,1);
        end
    end
end
for i=num_CHG_EV+1:num_EV  % for the V2G EVs
    % number of parking slots
    num_p_slots(i,1)=sum(F(i,:));
    % the constant rate
    if num_p_slots(i,1)-2*num_dischg_slot>0
        C_rate(i,1)=N_inital_en(i,2)/(num_p_slots(i,1)-2*num_dischg_slot);
    else
        fprintf('The number of discharging slots %g is too high.\n',num_dischg_slot);
    end
    if C_rate(i,1)>P_max
        fprintf('At EV %g, the constant rate %g exceeded the max rate %g.\n', i, C_rate(i,1), P_max);
    end
    % rate allocation
    price_vector=[];
    i_slot=0;
    for j=1:num_slot
        if F(i,j)==1
            i_slot=i_slot+1;
            price_vector(i_slot,1)=j;
            price_vector(i_slot,2)=price_basic(j);
            H_x_Matrix(i,j)=C_rate(i,1);
        end
    end
    % determine the discharging slot
    [value ind_1]=max(price_vector(:,2)); % max price
    H_x_Matrix(i,price_vector(ind_1,1))=-1*C_rate(i,1);
end
N_x_Matrix=H_x_Matrix;

% accumulative energy
accum_en=N_inital_en(:,1); % the initial energy
for i=1:num_slot
   accum_en(:,i+1)= accum_en(:,i)+H_x_Matrix(:,i);
end
% adjustment to remove negative accumulative energy
for i=1:num_EV
    for j=2:num_slot+1
        if accum_en(i,j)<0
            H_x_Matrix(i,j-1)=-1*H_x_Matrix(i,j-1);
            H_x_Matrix(i,j)=-1*H_x_Matrix(i,j);
        end
    end
end
% adjustment to remove overcharging (exceeding the battery capacity)
for i=1:num_EV
    tem_diff=0;
    for j=2:num_slot+1
        if accum_en(i,j)> Cap_battery_org
            tem_diff=accum_en(i,j)-Cap_battery_org;
            temp_old_chg=H_x_Matrix(i,j-1);
            H_x_Matrix(i,j-1)=H_x_Matrix(i,j-1)-tem_diff; % reduce the charging power
            %fprintf('EV %g at slot %g, accumulative energy %g, changed chg rate from %g to %g, cut %g.\n',i,j,accum_en(i,j),temp_old_chg, H_x_Matrix(i,j-1),tem_diff);
            for k=j:num_slot
                if H_x_Matrix(i,k)<0 % the discharging slot
                    temp_old_dischg=H_x_Matrix(i,k);
                    H_x_Matrix(i,j)=H_x_Matrix(i,k)+tem_diff;
                    %fprintf('EV %g at slot %g, changed chg rate from %g to %g, increased by %g.\n',i,k,temp_old_dischg,H_x_Matrix(i,j),tem_diff);
                end
            end
        end
    end
end

N_x_Matrix= H_x_Matrix;

% the evolution of energy level for each EV
N_Energy_variation=zeros(num_EV,num_slot+1);
for i=1:num_EV
    N_Energy_variation(i,1)=N_inital_en(i,1);  % the initial energy
    for j=1:num_slot
        N_Energy_variation(i,j+1)=N_Energy_variation(i,j)+F(i,j)*N_x_Matrix(i,j);
    end
end             
% the charged load
N_Charged_Load=zeros(num_slot,3); % 1) the base load, 2) the charged load, 3) the total load (from EV rates), 
N_Charged_Load(:,1)=L_b_mic; % the base load
for i=1:num_slot
    for j=1:num_EV
        N_Charged_Load(i,2)=N_Charged_Load(i,2)+N_x_Matrix(j,i)*F(j,i);
    end
    N_Charged_Load(i,3)=N_Charged_Load(i,1)+N_Charged_Load(i,2); % total load calculated from charged loads of individual EVs
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Plot the results %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Plot the results %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot the base load without EV charging
% load Glb_Chg_Load.txt;
figure;
xx=1:num_slot;
yy(:,1)=L_b_mic;  % the base load
yy(:,2)=v_Charged_Load(:,3); % globally optimal scheme
yy(:,3)=Charged_Load(:,3); %locally optimal scheme
yy(:,4)=N_Charged_Load(:,3); % Equal allocation scheme
plot(xx,yy);
ylabel('Load [KW]');
xlabel('Hour No.');
legend('Base load without EV charging','Total load with globally optimal EV charging','Total load with locally optimal EV charging','Total load with naive EV charging',1);

% conversion of charging rate matrix
G_x_Matrix_1D=reshape(v_x_Matrix', [], 1);  % globall optimal
L_x_Matrix_1D=reshape(x_Matrix', [], 1);  % locally optimal
N_x_Matrix_1D=reshape(N_x_Matrix', [], 1);  % naive

% the comparison of objective values    
obj_value=zeros(3,1);
obj_value(1)=k_0*sum(pow_p(yy(:,2),1)) + (k_1/2)*sum(pow_p(yy(:,2),2)) + beta*sum(square(F1)*square(G_x_Matrix_1D)) ...
        -k_0*sum(pow_p(L_b_mic(1:num_slot),1)) - (k_1/2)*sum(pow_p(L_b_mic(1:num_slot),2)); % total cost in globally optimal scheme
obj_value(2)=k_0*sum(pow_p(yy(:,3),1)) + (k_1/2)*sum(pow_p(yy(:,3),2)) + beta*sum(square(F1)*square(L_x_Matrix_1D)) ...
        -k_0*sum(pow_p(L_b_mic(1:num_slot),1)) - (k_1/2)*sum(pow_p(L_b_mic(1:num_slot),2)); % total cost in locally optimal scheme
obj_value(3)=k_0*sum(pow_p(yy(:,4),1)) + (k_1/2)*sum(pow_p(yy(:,4),2)) + beta*sum(square(F1)*square(N_x_Matrix_1D)) ...
        -k_0*sum(pow_p(L_b_mic(1:num_slot),1)) - (k_1/2)*sum(pow_p(L_b_mic(1:num_slot),2)); % total cost in Equal allocation scheme
    
imp_1=(obj_value(3)-obj_value(1))/obj_value(3);
imp_2=(obj_value(3)-obj_value(2))/obj_value(3);
fprintf('The objective value comparison: globally optimal scheme=%g, locally optimal scheme (group size %g) =%g, Equal allocation scheme=%g.\n',obj_value(1),group_size, obj_value(2),obj_value(3) );
fprintf('The improvement of the global optimal scheme is %g, the improvement of locally optimal scheme is %g.\n',imp_1, imp_2);


% plot the charged load in each interval
figure;
xx=1:num_slot;
zz(:,1)=v_Charged_Load(:,2); % globally optimal scheme
zz(:,2)=Charged_Load(:,2); %locally optimal scheme
zz(:,3)=N_Charged_Load(:,2); % Equal allocation scheme
plot(xx,zz);
ylabel('Charged load [KW]');
xlabel('Hour No.');
legend('Globally optimal scheme','Locally optimal scheme','Equal allocation scheme',1);


% the total charged energy
total_charged=0;
for i=1:num_EV
    total_charged=total_charged+ (Cap_battery-E_Charged(i,1));
end
fprintf('The total charged energy should be %g.\n',total_charged);
fprintf('The actual total charged energy:  globally optimal scheme=%g, locally optimal scheme=%g, Equal allocation scheme=%g.\n',...
    sum(yy(:,2)-yy(:,1)), sum(yy(:,3)-yy(:,1)), sum(yy(:,4)-yy(:,1)) );
    


% the load peak
Peak_based=max(L_b_mic); % the base load
Peak_Charged=max(Charged_Load(:,3)); % the load with EV charging
Peak_reduction=(Peak_based-Peak_Charged)/Peak_based;
fprintf('The peak comparison: base load=%g, globally optimal scheme=%g, locally optimal scheme=%g, Equal allocation scheme=%g KW.\n',...
    max(yy(:,1)), max(yy(:,2)), max(yy(:,3)), max(yy(:,4)) );

% the load standard deviation
std_based=std(L_b_mic); % the base load
std_Charged=std(Charged_Load(:,3)); % the load with EV charging
std_reduction=(std_based-std_Charged)/std_based;
fprintf('The standard deviation comparison: base load=%g, globally optimal scheme=%g, locally optimal scheme=%g, Equal allocation scheme=%g.\n',...
    std(yy(:,1)), std(yy(:,2)), std(yy(:,3)), std(yy(:,4)) );


fprintf('The number of the EVs is %g. The percentage of charging-only EVs is %g.\n', num_EV, P_Chg);

% plot the evolution of energy level for each EV (globally optimal scheme)
figure;
xxx=0:num_slot;
plot(xxx,v_Energy_variation(1:40,:));
ylabel('Energy [KWH]');
xlabel('Hour No.');
legend('EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','EV9','EV10',1);
title('The energy evolution in globally optimal scheme');

% plot the evolution of energy level for each EV (locally optimal scheme)
figure;
xxx=0:num_slot;
plot(xxx,Energy_variation);
ylabel('Energy [KWH]');
xlabel('Hour No.');
legend('EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','EV9','EV10',1);
title('The energy evolution in locally optimal scheme');

% plot the evolution of energy level for each EV (Equal allocation scheme)
figure;
plot(xxx,N_Energy_variation);
ylabel('Energy [KWH]');
xlabel('Hour No.');
legend('EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','EV9','EV10',1);
title('The energy evolution in Equal allocation scheme');

% % plot the charging rates for each EV
figure; 
EV_ID=65;
energy_mmm(1,:)=v_Energy_variation(EV_ID,:);
energy_mmm(2,:)=Energy_variation(EV_ID,:);
energy_mmm(3,:)=N_Energy_variation(EV_ID,:);
plot(xxx,energy_mmm);
ylabel('Energy [KWH]');
xlabel('Time (Hours)');
legend('Globally optimal scheme','Locally optimal scheme','Equal allocation scheme',1);
title('The energy evolution for an EV');

figure;
nnn(:,1)=v_x_Matrix(EV_ID,:)';%globally optimal scheme
nnn(:,2)=x_Matrix(EV_ID,:)';%locally optimal scheme
nnn(:,3)=N_x_Matrix(EV_ID,:)';%Equal allocation scheme
h=bar(xx,nnn);
ylabel('Rate [KW]');
xlabel('Time (Hours)');
legend('Globally optimal scheme','Locally optimal scheme','Equal allocation scheme',1);
title('The charging/discharging rate in globally optimal scheme');

% save glabol_x_Matrix.txt v_x_Matrix -ascii;
% save local_x_Matrix.txt x_Matrix -ascii;
% save naive_x_Matrix.txt N_x_Matrix -ascii;

