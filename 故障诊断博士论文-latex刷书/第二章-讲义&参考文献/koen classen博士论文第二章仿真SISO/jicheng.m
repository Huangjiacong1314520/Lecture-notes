% 第二章仿真
clear all; close all; clc;

%% 0.1. 定义系统参数
s = tf('s');

% 这里作为例子可能不太合适，这是example2.24，而2.24是直接给出\tilde{G}_{d}(s,\Delta)
% 修正传递函数定义（根据文献精确表达式）
G1 = (1/pi * s)/(s + 2*pi)^2;      % G1 = (1/π * s)/(s+2π)²
G2 = (1/(10*pi) * s)/(s + 20*pi)^2; % G2 = (1/10π * s)/(s+20π)²  
G3 = (1/(100*pi) * s)/(s + 200*pi)^2; % G3 = (1/100π * s)/(s+200π)²

% bode(G1)
% hold on
% bode(G2)
% hold on
% bode(G3)
% hold on
Gu_s0 = 7*G1*0.25 + 3*G2*1 + 3*G3*0.75;
bode(Gu_s0)

% 动态权重函数
W = 0.2*(s + 10*pi)/(s + 20*pi);

fprintf('========== Example 2.24 & 2.25 仿真 ==========\n');

%% 0.2. 定义频率范围
freq = logspace(-1, 3, 500);  % 0.1 Hz to 1000 Hz
omega = 2*pi*freq;
n_freq = length(omega);

%% 0.3. 计算名义系统 Gu(s,0)
Gu_s0 = 7*G1*0.25 + 3*G2*1 + 3*G3*0.75;
[mag_nominal, phase_nominal] = bode(Gu_s0, omega);
mag_nominal_dB = 20*log10(squeeze(mag_nominal));
bode(Gu_s0)


%% 0.4. 计算不确定性系统 GU(s,Δ)
% 使用不确定性对象定义 δ1 和 Δ1
delta1 = ureal('delta1', 0, 'Range', [-1 1]);   % δ1 ∈ [-1, 1] 不确定的实参数
Delta1 = ultidyn('Delta1', [1 1]);             % Δ1(s) 动态不确定性，||Δ1||∞ ≤ 1

% 构造不确定系统传函: Gu(s, Δ) = [7*G1*(0.25+δ1) + 3*G2*(1-δ1) + 3*G3*(0.75-δ1)] * [1 + W(s)*Δ1(s)]
N = 7*G1*(0.25 + delta1) + 3*G2*(1 - delta1) + 3*G3*(0.75 - delta1);
Gu_sDelta = N * (1 + W * Delta1);
figure
bode(Gu_sDelta)
 
%% 0.5. 计算不确定性系统外界干扰d和故障f对输出y的含不确定性传递函数
% y = Gu(s, ∆)u + Gd(s, ∆)d + Gf(s, ∆)f
% Gd(s, ∆)怎么得到？只有系统不确定性模型的时候？
% Gd_s_Delta=
% Gf(s, ∆)怎么得到？只有系统不确定性模型的时候？
% Gf_s_Delta=

%% 1. 计算标称系统Gu(s,0)的左互质分解LCF，得到\tilde{M}_{u}, \tilde{N}_{u}
% Bop = bodeoptions;
% Bop.Grid = 'on';
% Bop.XLim = {[0.0001 100000]};
% Bop.FreqUnits='Hz';
% %
% Atemp = ss(Gu_s0).A;
% Ctemp = ss(Gu_s0).C;
% Btemp = ss(Gu_s0).B;
% Dtemp = ss(Gu_s0).D;
% % Lpole = [-0.0001 -0.0005 -0.0006 -0.0009 -0.001 -0.0012];
% Lpole = [-0.00001 -0.00005 -0.00006 -0.00009 -0.0001 -0.00012];
% L_lncf = place(Atemp',Ctemp',Lpole);
% 
% eig_value = eig(Atemp - L_lncf'*Ctemp);
% 
% Mu_tilde = ss(Atemp + L_lncf'*Ctemp,L_lncf',Ctemp,eye(3));
% Nu_tilde = ss(Atemp + L_lncf'*Ctemp,Btemp + L_lncf'*Dtemp,Ctemp,Dtemp);
% 
% 
% figure;subplot(221);bode(Mu_tilde,Bop);title('Ml求');
% subplot(222);bode(Nu_tilde,Bop);title('Nl求');
% subplot(223);bode(inv(Mu_tilde)*Nu_tilde,Bop);title('Ml^-1*Nl求');

% 上面是自己拿因子拼接，可能就不是最小状态空间实现
% 返回Gu_s0的最小状态空间实现fact和互质因子
[fact,Mu_tilde,Nu_tilde] = lncf(Gu_s0);
% 作图：用于验证[Ml Nl] 是否是酉的（即 𝑀 𝑙 𝑀 𝑙 ∗ + 𝑁 𝑙 𝑁 𝑙 ∗ = 𝐼 ）
figure
sigma(fact)
% 作图：验证验证是否sys=Ml^(−1)​Nl​，
figure
sigma(Gu_s0,'b-', Mu_tilde\Nu_tilde,'r--')


%% 2. 计算系统的不确定性部分：$\tilde{G}_{u}(s, \Delta):=G_{u}(s, \Delta)-G_{u}(s, 0)$

Gu_s_tilde=Gu_sDelta-Gu_s0;



%% 2. 鲁棒控制器设计





robustcontroller;

%% 2. 计算rdf到预残差$\tilde{\epsilon}$的不确定性传递矩阵
% 鲁棒控制器
K = K_hinf;%H无穷鲁棒控制器
% K = 其他鲁棒控制器
% 这个eye肯定是有问题的，注意修改
S_Delta = (eye(3) + Gu_sDelta*K)^(-1);

M = Mu_tilde;
N = Nu_tilde;

% Tr = M*Gu_s_tilde*K*(eye(3)-S_Delta*Gu*K);

Tr = M*Gu_s_tilde*K*S_Delta;

Td = M*(Gd_s_Delta - Gu_s_tilde*K*S_Delta*Gd_s_Delta);
Tf = M*(Gf_s_Delta - Gu_s_tilde*K*S_Delta*Gf_s_Delta);
% Gd_s_tilde
Gd_s_tilde = [Tr Td];

% 师兄给的滤波器，应该不适用这里
filter2 = tf([1 0 0],[1 0.707*2*0.001*2*pi 0.001*2*pi*0.001*2*pi]);














wc = wcgain(Gu_sDelta);
    wcGain = wc.UpperBound;
    wcFreq = wc.CriticalFrequency;
wcGain_dB = 20*log10(wcGain);

fprintf('最坏情况增益 = %.4f (%.2f dB)，发生在频率 ω = %.4f rad/s。\n', wcGain, wcGain_dB, wcFreq);
% 注：wcgain 自动搜索整个频率范围的最坏增益及不确定值。