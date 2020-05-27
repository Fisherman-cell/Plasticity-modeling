% 设置材料参数（copper）
% 参考：http://www.matweb.com/search/datasheet_print.aspx?matguid=9aebe83845c04c1db5126fada6f76f7e
% 参考：https://en.wikipedia.org/wiki/Bulk_modulus
YOUNG = 110*10^3; % E=110*10^3MPa
POISS = 0.343; % nu =0.343
SIGMAY_0 = 33.3; % Tensile Strength, Yield = 33.3 MPa
ULTI_TENSILE = 210; % Tensile Strength, Ultimate =210 MPa
BREAK_LONGATION = 0.6; % Elongation at Break 	

% 计算体积模量和其他必须的参量
% Assume that pure copper is Homogeneous isotropic linear elastic materials
GMODU = YOUNG/(2*(1+POISS)); % shear modulu G,真实为46.0 GPa	
BULK = YOUNG/(3*(1-2*POISS)); % bulk modulu K 真实为140 GPa	
%H = (GMODU/5)*(1-EPSILONP*10^3);  %(ULTI_TENSILE-SIGMAY_0)/(BREAK_LONGATION-SIGMAY_0/YOUNG);% hardening modulus

% 定义单轴拉伸条件
FACTOR = 1*10^(-5); %每一步的增量应变
NUM = 100; %1/4部分的数据点数目
DELTA_EPSILON_1 = ones(1,NUM);
DELTA_EPSILON_2 = FACTOR.*DELTA_EPSILON_1; %设置一个周期中拉伸部分的累积塑性
DELTA_EPSILON_3 = -FACTOR.*DELTA_EPSILON_1; %设置一个周期中压缩部分的累积塑性
DELTA_EPSILON = cat(2,DELTA_EPSILON_2,DELTA_EPSILON_3);
NUM_ALL = 100;
% CYCLE_NUM = 1;
% NUM_ALL = 4*NUM*CYCLE_NUM; %整个循环中的数据点数目
% DELTA_EPSILON = zeros(3);
% DELTA_EPSILON(1,1) = 5*10^(-5);% 增量应变
% DATA_NUM = 30; %提取数据点的数目

%初始化输出数据
EPSILON_CA = zeros(1,NUM_ALL+1); % 单轴拉伸曲线上的应变，标量
SIGMA_CA = zeros(1,NUM_ALL+1); % 单轴拉伸曲线上的应力，标量

% 定义初始的应力及应变
EPSILON = zeros(3); %初始化应变张量
EPSILONP = 0;% accumulated plastic strain
I = eye(3); % 设置单位矩阵
EPSILONV = trace(EPSILON).*I./3; % volumetric strain tensor，张量形式
EPSILOND = EPSILON-EPSILONV; % deviatoric strain tensor
P = BULK*EPSILONV;% hydrostatic pressure
S = 2*GMODU*EPSILOND; % deviatoric stress tensor
SIGMA = P+S; % stress tensor
EPSILON_trial=zeros(3);


for i=1:NUM_ALL
    % 弹性预测，完成弹性测试阶段
    % -------------------------------------------
    % 假设状态为弹性，这累积塑性乘积为0
    DGAMA = 0; % The incremental plastic multiplier
    EPSILONP_trial = EPSILONP;  % trial accumulated plastic strain
    H = (GMODU/5)*(1-EPSILONP*10^3);  %(ULTI_TENSILE-SIGMAY_0)/(BREAK_LONGATION-SIGMAY_0/YOUNG);% hardening modulus
    SIGMAY = SIGMAY_0+H*EPSILONP_trial^1; % 塑性硬化函数：线性硬化
    SIGMAY_CA(1,i) =  SIGMAY;
    
    % evaluate the elastic trial state
    EPSILON_trial = EPSILON;
    EPSILON_trial(1,1) =EPSILON_trial(1,1) +DELTA_EPSILON(1,i); % 设置the elastic strain tensor
    EPSILONV_trial = trace(EPSILON_trial).*I/3; % trial volumetric strain tensor，张量形式
    EPSILOND_trial = EPSILON_trial-EPSILONV_trial; % trial deviatoric strain tensor
    P_trial = BULK.*EPSILONV_trial;% trial hydrostatic pressure
    S_trial = 2*GMODU.*EPSILOND_trial; % trial deviatoric stress tensor
    SIGMA_trial = P_trial+S_trial; % trial stress tensor
    Q_trial = sqrt(3/2*dot(dot(S_trial,S_trial),[1 1 1])); % the elastic trial Von Mises effect stress
    PHI_trial = Q_trial-SIGMAY; % yield function for Von Mises model

    % check plastic admissibility
    if PHI_trial <= 0 | DELTA_EPSILON(1,i)<=0
        % Update stress and strain
        EPSILON = EPSILON_trial;
        EPSILONV = trace(EPSILON).*I/3;
        EPSILOND = EPSILOND_trial;
        EPSILONP = EPSILONP_trial;
        P = P_trial;
        S = S_trial;
        Q = Q_trial;
        SIGMA = P+S; % upadte stress tensor
    else
        % plastic correction
        %----------------------------------------------------------
        % use Newton-Raphson algorithm to solve the return mapping equation
        % calculate DGAMA
        j=0;
        while (abs(PHI_trial) >= 1*10^(-5))
            j=j+1
            %H = (GMODU/5)*(1-EPSILONP);%*EPSILONP;  %(ULTI_TENSILE-SIGMAY_0)/(BREAK_LONGATION-SIGMAY_0/YOUNG);% hardening modulus
            d = -3*GMODU-H; % Compute residual derivative,prefect plasticity
            DGAMA = DGAMA-PHI_trial/d; %计算迭代之后的增量塑性乘积
            SIGMAY = SIGMAY_0+H*(EPSILONP_trial+DGAMA); % 塑性硬化函数：线性硬化
            PHI_trial =  Q_trial-3*GMODU*DGAMA-SIGMAY; %求解屈服函数
            
            FACTOR = (1-(DGAMA*3*GMODU)/Q_trial);
            S_trial = FACTOR.*S_trial; % upadte deviatoric stress tensor
            EPSILONP_trial = EPSILONP_trial+DGAMA; % upadte accumulated plastic strain
            
        end
        Q_trial = sqrt(3/2.*dot(dot(S_trial,S_trial),[1 1 1])); % the elastic trial Von Mises effect stress
        P = P_trial; % update hydrostatic pressure
        S = S_trial; % upadte deviatoric stress tensor
        SIGMA = P+S; % upadte stress tensor
        EPSILOND = S./(2*GMODU); % upadte volumetric strain tensor
        EPSILON = EPSILOND+EPSILONV_trial; % upadte the elastic strain    
        EPSILONP = EPSILONP_trial; % upadte accumulated plastic strain
    end
%         % update the state variables
%         P = P_trial; % update hydrostatic pressure
%         S = (1-(DGAMA*3*GMODU)/Q_trial).*S_trial; % upadte deviatoric stress tensor
%         SIGMA = P+S; % upadte stress tensor
%         EPSILON = S./(2*GMODU)+EPSILONV; % upadte the elastic strain
%         EPSILOND = S./(2*GMODU); % upadte volumetric strain tensor
%         EPSILONP = EPSILONP_trial+DGAMA; % upadte accumulated plastic strain 

    % 将数据储存到单轴拉伸曲线图的数据点（应力+应变）
    EPSILOND_CA(1,i+1) = EPSILOND(1,1);
    EPSILONP_CA(1,i+1) = EPSILONP;
    EPSILON_CA(1,i+1) = EPSILON(1,1);
    SIGMA_CA(1,i+1) = S_trial(1,1);
    a=SIGMA_CA./EPSILON_CA;
end
scatter(EPSILON_CA,SIGMA_CA);
hold on;
plot(EPSILON_CA,SIGMA_CA);
hold off;

    
    
