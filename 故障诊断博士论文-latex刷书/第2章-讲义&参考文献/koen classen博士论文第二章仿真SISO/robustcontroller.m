%% 3. 鲁棒控制器设计（H∞ 控制器设计）
% 参考：Kemin Zhou, Essentials of Robust Control, Ch.14
% 思路：构造广义植物 P(s)，选择性能权重函数 W1, W2, W3，进行 H∞ 综合。
% 目标：在满足鲁棒稳定性的前提下，最小化闭环传递函数 ||Tzw||∞。

fprintf('\n========== H∞ 鲁棒控制器设计 ==========\n');

%% 3.1 构造加权系统（混合灵敏度设计）
% -------------------------------------------
% 性能指标：||Tzw||∞ = || [ W1*S; W2*K*S; W3*T ] ||∞ < γ
% S = (I+PK)^(-1)，T = PK(I+PK)^(-1)
%
% 权重选择：
%   W1: 低频增强跟踪性能（大增益，抑制稳态误差）
%   W2: 高频限制控制输入（防止过大控制能量）
%   W3: 建模不确定性加权（确保鲁棒稳定性）

W1 = makeweight(10, 100, 0.01);  % 低频10，高频0.01，交叉频率≈100 rad/s
W2 = makeweight(0.1, 50, 10);    % 限制控制能量
% 根据鲁棒控制的需求重新选择参数
s = tf('s');
W3 = 0.5*(s/200 + 1)/(s/2000 + 1);


% 名义被控对象
P_nom = Gu_s0;

% 使用 connect 或 sysic 构造广义plant
systemnames = 'P_nom W1 W2 W3';
inputvar = '[r; n; d; u]';          % r:参考输入, n:测量噪声, d:扰动
outputvar = '[W1; W2; W3; r-P_nom]'; % 控制性能输出 + 误差信号
input_to_P_nom = '[n + u]';      % 控制器输入为测量信号+控制
input_to_W1 = '[r-P_nom]';
input_to_W2 = '[u]';
input_to_W3 = '[P_nom]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;

%% 3.2 设计 H∞ 控制器
% -------------------------------------------
% hinfsyn(P, nmeas, ncont)
% nmeas = 1：测量输出维度（y）
% ncont = 1：控制输入维度（u）
nmeas = 1;
ncont = 1;
[K_hinf, CL_hinf, gamma_opt] = hinfsyn(P, nmeas, ncont);

fprintf('最优 H∞ 控制器求解完成：γ* = %.4f\n', gamma_opt);

%% 3.3 检查闭环鲁棒性能
% -------------------------------------------
% 检查闭环稳定性与鲁棒裕度
figure;
sigma(CL_hinf);
title('闭环奇异值图 (H∞ 设计结果)');

% 计算最坏情况增益
wc = wcgain(CL_hinf);
fprintf('闭环系统最坏增益上界 = %.3f (%.2f dB)\n', wc.UpperBound, 20*log10(wc.UpperBound));

%% 3.4 时域响应验证
% -------------------------------------------
% 对单位阶跃输入进行时域仿真
T_closed = feedback(P_nom*K_hinf, 1);  % 闭环传递函数
figure;
step(T_closed);
title('闭环阶跃响应 (H∞ 控制器)');

%% 3.5 鲁棒稳定性与鲁棒性能分析
% -------------------------------------------
% 使用 Robust Control Toolbox 的 µ-分析函数：
[stabmarg, destabunc, report, info] = robuststab(CL_hinf);
disp(report);   % 显示鲁棒稳定性报告

[perfmarg, pertu, reportp, infop] = robustperf(CL_hinf);
disp(reportp);  % 显示鲁棒性能报告

fprintf('\n鲁棒稳定性裕度：%.3f\n鲁棒性能裕度：%.3f\n', ...
    stabmarg.LowerBound, perfmarg.LowerBound);

%% 3.6 结果可视化总结
figure;
subplot(2,1,1);
bode(P_nom, K_hinf);
legend('被控对象 G(s)','H∞ 控制器 K(s)');
title('H∞ 控制器与被控对象频率响应');

subplot(2,1,2);
margin(P_nom*K_hinf);
title('开环增益与相位裕度分析');
grid on;

% %% 4. µ-综合（D–K 迭代）
% % 目标：在存在结构化不确定性 Δ = diag(delta1, Delta1) 的情况下，
% %       通过 D–K 迭代最小化 sup_ω μ( M(jω) )，获得更强鲁棒性的控制器。
% % 说明：若有 musyn（R2016a+），将使用其内置 D–K 迭代；否则退回 dksyn。
% 
% fprintf('\n========== µ-综合（D–K 迭代） ==========\n');
% 
% 
% % 复用上面定义的权重和不确定系统
% % Gu_sDelta: 不确定对象
% % W1, W2, W3: 权重函数
% systemnames = 'Gu_sDelta W1 W2 W3';
% inputvar    = '[r; n; u]';   % 三个输入：r参考、n噪声、u控制信号
% 
% % 定义中间信号:
% % e1 = r - y （跟踪误差信号）
% % y = Gu_sDelta * (n + u)
% % z1 = W1 * e1
% % z2 = W2 * u
% % z3 = W3 * y
% %
% % 输出顺序：性能输出[z1; z2; z3] 和 测量输出 y（供控制器使用）
% outputvar   = '[ W1*(r - Gu_sDelta*(n + u)); W2*u; W3*(Gu_sDelta*(n + u)); Gu_sDelta*(n + u) ]';
% 
% % 各系统输入
% input_to_Gu_sDelta = '[n + u]';
% input_to_W1 = '[r - Gu_sDelta*(n + u)]';
% input_to_W2 = '[u]';
% input_to_W3 = '[Gu_sDelta*(n + u)]';
% 
% sysoutname  = 'P_mu';
% cleanupsysic = 'yes';
% sysic;
% 
% 
% % 计数：测量输出 y=Gu_sDelta → nmeas = 1；控制输入 u → ncont = 1
% nmeas = 1; 
% ncont = 1;
% 
% % --- 4.2 执行 µ-综合（优先 musyn，备用 dksyn） ---
% use_musyn = exist('musyn','file') == 2;
% if use_musyn
%     musyn_opts = musynOptions('Display','brief');    % 简洁输出
%     [K_mu, CL_mu, muinfo] = musyn(P_mu, nmeas, ncont, musyn_opts);
%     fprintf('musyn 已完成 D–K 迭代。近似最优性能（upper bound）：%.4f\n', muinfo.UpperBound);
% else
%     % dksyn 会自动在内部做 D–K 迭代，可能比 musyn 慢，报告在 info 中
%     [K_mu, CL_mu, muinfo] = dksyn(P_mu, nmeas, ncont);
%     if isfield(muinfo,'gamma')
%         fprintf('dksyn 已完成 D–K 迭代。近似最优 γ：%.4f\n', muinfo.gamma);
%     else
%         fprintf('dksyn 已完成 D–K 迭代。\n');
%     end
% end
% 
% % --- 4.3 闭环验证（不确定闭环） ---
% % CL_mu 是从 [r; n; d] → [z; y] 的闭环不确定系统
% % 进行鲁棒稳定性与鲁棒性能分析
% fprintf('\n--- µ-综合：鲁棒性报告 ---\n');
% [stabmarg_mu, destabunc_mu, report_mu_stab] = robuststab(CL_mu);
% disp(report_mu_stab);
% [perfmarg_mu, pert_mu, report_mu_perf] = robustperf(CL_mu);
% disp(report_mu_perf);
% 
% % 最坏情况增益（含不确定性）与关键频率
% wc_mu = wcgain(CL_mu);
% fprintf('µ-综合闭环最坏增益上界 = %.3f (%.2f dB) ，关键频率 ω = %.4f rad/s\n', ...
%     wc_mu.UpperBound, 20*log10(wc_mu.UpperBound), wc_mu.CriticalFrequency);
% 
% % --- 4.4 可视化对比：H∞ vs µ ---
% % 若你前面已得到 H∞ 控制器 K_hinf 与闭环 CL_hinf，这里给出对比图；
% % 若没有 H∞ 部分，这段对比可以注释掉。
% try
%     figure('Name','闭环奇异值对比 H∞ vs µ'); 
%     sigma(CL_hinf,'b-', CL_mu,'r--'); grid on;
%     legend('H_\infty 闭环','\mu-综合 闭环','Location','Best');
%     title('闭环奇异值对比（H_\infty vs \mu-综合）');
% 
%     % 阶跃响应对比（名义模型下）
%     % 以名义对象 Gu_s0 形成两者的名义闭环（便于时域对比）
%     T_hinf_nom = feedback(Gu_s0*K_hinf, 1);
%     T_mu_nom   = feedback(Gu_s0*K_mu,   1);
%     figure('Name','阶跃响应对比 H∞ vs µ（名义模型）');
%     step(T_hinf_nom,'b-', T_mu_nom,'r--', 0.5); grid on;  % 视系统快慢可调时间轴
%     legend('H_\infty 名义闭环','\mu-综合 名义闭环','Location','Best');
%     title('阶跃响应对比（名义模型下）');
% catch
%     % 若 H∞ 变量不存在，则仅绘制 µ-综合结果
%     figure('Name','µ-综合：闭环奇异值');
%     sigma(CL_mu); grid on; title('\mu-综合 闭环奇异值');
%     figure('Name','µ-综合：名义阶跃响应');
%     step(feedback(Gu_s0*K_mu,1)); grid on; title('\mu-综合 名义阶跃响应');
% end
% 
% % --- 4.5 结果小结 ---
% % 一般经验：
% % - 若模型存在显著的结构化不确定性（如你的 Delta1: ultidyn；delta1: ureal），
% %   则 µ-综合较 H∞ 能获得更大的鲁棒裕度（robuststab LowerBound 更大），
% %   同时 robustperf 报告中的 LowerBound 往往也更接近或超过 1。
% % - 代价是控制器复杂度与带宽可能上升。可用 reduce/balred 做阶次约简：
% %   K_mu_red = balred(K_mu, 5);   % 例如约到 5 阶（按需要调整）



%% 4. µ-综合（D–K 迭代）
fprintf('\n========== µ-综合（D–K 迭代） ==========\n');

% 使用之前定义的权重和不确定系统
% Gu_sDelta, W1, W2, W3

% ---- 4.1 定义信号接口 ----
% 输入端口：[r; n; u]  (参考输入、噪声、控制)
% 输出端口：[z1; z2; z3; y]
%
% 其中：
%   e = r - y
%   y = Gu_sDelta*(n + u)
%   z1 = W1*e
%   z2 = W2*u
%   z3 = W3*y
%
% 我们用 connect() 手动搭建信号路径

% 定义信号节点名
Gu_sDelta.InputName  = {'nin','uin'};   % 两个输入：噪声n，控制u
Gu_sDelta.OutputName = 'y';

W1.InputName = 'e';
W1.OutputName = 'z1';

W2.InputName = 'u';
W2.OutputName = 'z2';

W3.InputName = 'y';
W3.OutputName = 'z3';

% 定义误差信号 e = r - y
Sum1 = sumblk('e = r - y');
% 定义对象输入 y = Gu_sDelta*(n + u)
Sum2 = sumblk('uin = u');
Sum3 = sumblk('nin = n');

% 构造广义植物
P_mu = connect(Gu_sDelta, W1, W2, W3, Sum1, Sum2, Sum3, {'r','n','u'}, {'z1','z2','z3','y'});

% ---- 4.2 µ-综合（D–K 迭代） ----
nmeas = 1;  % 测量输出 y
ncont = 1;  % 控制输入 u

if exist('musyn','file') == 2
    [K_mu, CL_mu, muinfo] = musyn(P_mu, nmeas, ncont);
    fprintf('musyn 已完成 D–K 迭代，UpperBound = %.4f\n', muinfo.UpperBound);
else
    [K_mu, CL_mu, muinfo] = dksyn(P_mu, nmeas, ncont);
    fprintf('dksyn 已完成 D–K 迭代。\n');
end

% ---- 4.3 µ-鲁棒性能验证 ----
fprintf('\n--- µ-综合鲁棒性报告 ---\n');
[stabmarg_mu,~,report_mu_stab] = robuststab(CL_mu);
disp(report_mu_stab);
[perfmarg_mu,~,report_mu_perf] = robustperf(CL_mu);
disp(report_mu_perf);

wc_mu = wcgain(CL_mu);
fprintf('µ-综合闭环最坏增益上界 = %.3f (%.2f dB) ，关键频率 ω = %.4f rad/s\n', ...
    wc_mu.UpperBound, 20*log10(wc_mu.UpperBound), wc_mu.CriticalFrequency);

% ---- 4.4 可视化对比 ----
try
    figure('Name','闭环奇异值对比 H∞ vs µ');
    sigma(CL_hinf,'b-', CL_mu,'r--'); grid on;
    legend('H_\infty 闭环','\mu-综合 闭环','Location','Best');
    title('闭环奇异值对比（H_\infty vs \mu-综合）');

    figure('Name','阶跃响应对比 H∞ vs µ（名义模型）');
    T_hinf_nom = feedback(Gu_s0*K_hinf,1);
    T_mu_nom   = feedback(Gu_s0*K_mu,1);
    step(T_hinf_nom,'b-',T_mu_nom,'r--',0.5);
    legend('H_\infty 名义闭环','\mu-综合 名义闭环');
    grid on; title('阶跃响应对比（名义模型）');
catch
    figure; sigma(CL_mu); title('\mu-综合 闭环奇异值'); grid on;
end
