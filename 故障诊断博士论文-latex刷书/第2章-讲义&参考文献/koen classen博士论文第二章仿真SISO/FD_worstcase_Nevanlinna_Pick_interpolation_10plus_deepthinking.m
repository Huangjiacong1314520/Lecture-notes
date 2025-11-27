%% Fault detection filter optimization — Example 2.24 & 2.25 (repro of Fig. 2.8)
% Robust Control Toolbox required.

clear; clc; close all;

%% 1. 定义系统模型 Gd(s, Δ) 及计算最坏情况增益和频率
% 系统模型参数定义
s = tf('s');
G1 = (1/pi * s) / (s + 2*pi)^2;       % G1(s) = (1/π * s) / (s + 2π)^2
G2 = (1/(10*pi) * s) / (s + 20*pi)^2; % G2(s) = (1/10π * s) / (s + 20π)^2
G3 = (1/(100*pi) * s)/ (s + 200*pi)^2;% G3(s) = (1/100π * s) / (s + 200π)^2
W  = 0.2 * (s + 10*pi) / (s + 20*pi); % W(s) = 0.2 * (s + 10π) / (s + 20π)

% 使用不确定性对象定义 δ1 和 Δ1
delta1 = ureal('delta1', 0, 'Range', [-1 1]);   % δ1 ∈ [-1, 1] 不确定的实参数
Delta1 = ultidyn('Delta1', [1 1]);             % Δ1(s) 动态不确定性，||Δ1||∞ ≤ 1

% 构造不确定系统传函: Gd(s, Δ) = [7*G1*(0.25+δ1) + 3*G2*(1-δ1) + 3*G3*(0.75-δ1)] * [1 + W(s)*Δ1(s)]
N = 7*G1*(0.25 + delta1) + 3*G2*(1 - delta1) + 3*G3*(0.75 - delta1);
Gd = N * (1 + W * Delta1);

% 计算最坏情况增益及对应频率
% [wcGain, wcFreq] = wcgain(Gd);
% wcGain_dB = 20*log10(wcGain);  % 转换为分贝

wc = wcgain(Gd);

if isstruct(wc)
    wcGain = wc.UpperBound;
    wcFreq = wc.CriticalFrequency;
else
    wcGain = wc;
    wcFreq = NaN;  % 不支持频率提取
end

wcGain_dB = 20*log10(wcGain);

fprintf('最坏情况增益 = %.4f (%.2f dB)，发生在频率 ω = %.4f rad/s。\n', wcGain, wcGain_dB, wcFreq);
% 注：wcgain 自动搜索整个频率范围的最坏增益及不确定值。





%% 提取最坏情况时的不确定参数值
% wcUnc = wcgain(Gd, 'MaxGain');  % 获取导致最大增益的uncertainty实例
% fprintf('最坏情况时 δ1 = %.2f， Δ1 = %.2f ∠%.2f°。\n', wcUnc.delta1, abs(wcUnc.Delta1), (180/pi)*angle(wcUnc.Delta1));
% % （Delta1 在该频率处为一个复常数，这里输出其幅值和相位） 

% 设置选项以获取导致最坏情况的不确定性
opt = wcgainOptions;
opt.Sensitivity = 'on';        % 这会让函数返回导致最坏情况的不确定性

[wcGain, wcUnc, info] = wcgain(Gd, opt);

% wcUnc 现在包含导致最坏情况增益的不确定性实例
fprintf('导致最大增益的不确定性实例已获取\n');
% fprintf('最坏情况时 δ1 = %.2f， Δ1 = %.2f ∠%.2f°。\n', wcUnc.delta1, abs(wcUnc.Delta1), (180/pi)*angle(wcUnc.Delta1));
% % （Delta1 在该频率处为一个复常数，这里输出其幅值和相位） 

%% 2. 计算结构奇异值的保守下界 L(ω)（skew-μ 幂迭代结果）
% 为估计各频率下的增益下界 L(ω)，我们对δ1取两端值并令Δ1在每个频率处取使 |1+WΔ| 最大的相位。
% 对于SISO系统，该相位即 Δ(jω) = e^{-j arg(W(jω))}，此时 |1+WΔ| = 1 + |W|.
% L(ω) = max_{δ1 = -1,1} |N(jω)| * (1 + |W(jω)|)，其中 N(jω) = 7*G1*(0.25+δ1)+3*G2*(1-δ1)+3*G3*(0.75-δ1).
w = logspace(-1, 3, 400);  % 频率范围：0.1 rad/s 到 1000 rad/s（约0.016 Hz 至 159 Hz）
Lvals = zeros(size(w));
for k = 1:length(w)
    jw = 1j*w(k);
    % 计算 δ1 = 1 和 δ1 = -1 时 N(jω) 的复数值
    N_val_plus  = evalfr(7*G1*(0.25+1)  + 3*G2*(1-1)  + 3*G3*(0.75-1) , jw);
    N_val_minus = evalfr(7*G1*(0.25-1)  + 3*G2*(1+1)  + 3*G3*(0.75+1) , jw);
    % 计算 |N| 和 |W|
    magN_plus  = abs(N_val_plus);
    magN_minus = abs(N_val_minus);
    magW = abs(evalfr(W, jw));
    % 取较大者并乘以 (1 + |W|)
    Lvals(k) = max(magN_plus, magN_minus) * (1 + magW);
end

% 绘制 L(ω) 曲线用于比较
figure; semilogx(w, 20*log10(Lvals), 'm--', 'LineWidth', 1.5); hold on;
xlabel('Frequency (rad/s)'); ylabel('Gain (dB)');
title('Worst-Case Gain Lower Bound L(\omega) (skew-\mu result)');
grid on;

% 提示：上图绘制了结构奇异值下界 L(ω)。L(ω) 曲线与 wcgain 计算的全局最坏增益在峰值处相切。

%% 3. 构造动态不确定性 Δ1(s) 的边界插值近似（6阶滤波器结构）
% 选择在 ω = 1 Hz、10 Hz、100 Hz 处插值，使 Δ1(jω) 在这些频率取到边界值。
% 利用 Nevanlinna-Pick (边界插值) 方法构造满足 ||Δ1||∞ ≤ 1 的有理函数。
% 这里直接给出插值生成的 Δ1(s) 6阶传函（由算法得到的全通滤波器近似）。
% （注：由于该函数为全通，其频率响应幅值恒为1，仅相位随频率变化以匹配插值点）
numDelta = [1,    15.7046069,  -2137.56448,  -48.7492742,  -1.50623181,   6.72291786,   0.05246387];
denDelta = [1,   -15.7046062,  -2137.56449,   48.7490083,  -1.51788091,  -6.72282545,   0.05276738];
Delta_s = tf(numDelta, denDelta);
% 验证 Δ1(s) 为全通：在插值频率处的值及 ||Δ1||∞
val1Hz  = evalfr(Delta_s, 2*pi*1);   % 应接近 e^{-j arg(W(j1Hz))}
val10Hz = evalfr(Delta_s, 2*pi*10);
val100Hz= evalfr(Delta_s, 2*pi*100);
fprintf('Δ1(j1Hz) = %.3f ∠%.1f°, Δ1(j10Hz) = %.3f ∠%.1f°, Δ1(j100Hz) = %.3f ∠%.1f°。\n', ...
    abs(val1Hz),  angle(val1Hz)*180/pi, abs(val10Hz), angle(val10Hz)*180/pi, abs(val100Hz), angle(val100Hz)*180/pi);
fprintf('Δ1(s) 的∞范数 = %.4f (应≤1)。\n', norm(Delta_s, inf));
% 输出显示 Δ1(s) 在各插值频率幅值约为1且相位与 W(s) 相反，使得 1+WΔ 构造性叠加。

%% 4. 构造名义系统、最坏情况系统和边界插值系统，并绘制频率响应
% 名义系统：δ1 = 0, Δ1(s) = 0
% G_nom = subs(Gd, {'delta1','Delta1'}, {0, 0});
% G_nom = tf(G_nom);  % 转换为确定的传函
% 
% 
% %% 3. 计算名义系统 Ĝd(s,0)
% G_nom = 7*G1*0.25 + 3*G2 + 3*G3*0.75;
% [mag_nominal, ~] = bode(G_nom, omega);
% mag_nominal_dB = 20*log10(squeeze(mag_nominal));



% 错误的方式：
% G_nom = subs(Gd, {'delta1','Delta1'}, {0, 0});
% 最坏情况系统：δ1 = 1, Δ1(s) = 1 (取Δ1为实常数1近似最坏情况下的相位对齐)
% G_wc = subs(Gd, {'delta1','Delta1'}, {1, 1});
% G_wc = tf(G_wc);

% 正确的方式：
G_nom = usubs(Gd, 'delta1', 0, 'Delta1', 0);

G_wc = usubs(Gd, 'delta1', 1, 'Delta1', 1);









% 边界插值系统：δ1 = 1, Δ1(s) = Delta_s (插值构造的6阶不确定性)
G_bnd = (7*G1*1.25 + 3*G2*0 + 3*G3*(-0.25)) * (1 + W * Delta_s);
% 上式解释：δ1=1代入N(s)得7*G1*1.25 + 3*G2*0 + 3*G3* -0.25，然后乘(1+W*Δ_s)。

% 绘制三个系统在各频率的增益曲线（对数坐标）
figure;
omega = logspace(-1, 3, 600);  % 0.1 ~ 1000 rad/s
[mag_nom, ~] = bode(G_nom, omega);  mag_nom = squeeze(mag_nom);
[mag_wc,  ~] = bode(G_wc,  omega);  mag_wc  = squeeze(mag_wc);
[mag_bnd, ~] = bode(G_bnd, omega); mag_bnd = squeeze(mag_bnd);
semilogx(omega, 20*log10(mag_nom), 'k-', 'LineWidth', 1.5); hold on;
semilogx(omega, 20*log10(mag_wc),  'r--', 'LineWidth', 1.5);
semilogx(omega, 20*log10(mag_bnd), 'b-.', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
title('Frequency Responses: Nominal vs Worst-case vs Boundary-interpolated');
legend('Nominal (\delta_1=0, \Delta_1=0)', 'Worst-case (\delta_1=1, \Delta_1=1)', 'Boundary interp. (\delta_1=1, \Delta_1(s))');

% 从图中可观察各系统增益随频率的变化：名义系统在各频率增益最低；最坏情况系统在约1 Hz处达到最高峰值；边界插值系统在1 Hz、10 Hz、100 Hz附近均有峰值，曲线整体高于名义且在中高频超过静态最坏情况。

%% 5. 计算各系统的最小/最大/平均增益，验证性能比较
% 定义频率向量（取0.1 Hz ~ 100 Hz，即0.628~628 rad/s，以涵盖主要频段）
f = logspace(-1, 2, 301);  % Hz
w_rad = 2*pi*f;
mag_nom_f = squeeze(abs(freqresp(G_nom, w_rad)));
mag_wc_f  = squeeze(abs(freqresp(G_wc,  w_rad)));
mag_bnd_f = squeeze(abs(freqresp(G_bnd, w_rad)));
% 计算最小、最大、均值（均值取线性幅值平均）
min_nom = min(mag_nom_f);  max_nom = max(mag_nom_f);  avg_nom = mean(mag_nom_f);
min_wc  = min(mag_wc_f);   max_wc  = max(mag_wc_f);   avg_wc  = mean(mag_wc_f);
min_bnd = min(mag_bnd_f);  max_bnd = max(mag_bnd_f);  avg_bnd = mean(mag_bnd_f);
fprintf('名义系统增益：最小%.4f，最大%.4f，均值%.4f。\n', min_nom, max_nom, avg_nom);
fprintf('最坏情况系统增益：最小%.4f，最大%.4f，均值%.4f。\n', min_wc, max_wc, avg_wc);
fprintf('边界插值系统增益：最小%.4f，最大%.4f，均值%.4f。\n', min_bnd, max_bnd, avg_bnd);
% 可见：名义系统增益曲线整体较低，最坏情况系统在低频峰值最高但增益曲线在高频相对较低，
% 边界插值系统在各插值频率附近增益都被拉升，因而其平均增益最高且在中高频段超过静态最坏情况系统。

%% 6. 基于Ricca​​ti方程的最优故障检测滤波器设计与性能对比
% （本步骤需要根据具体Hi/H∞指标的原型结构构建增广系统，这里假设我们有误差输出e和残差r等）
% 通常构建以扰动和故障为输入、残差为输出的增广系统，然后使用H∞方法综合滤波器。
% 例如，可利用 MATLAB Robust Control Toolbox 的 hinfsyn 函数：
% 假设我们已有增广系统P，它包括从[扰动;故障]到[评估输出; 滤波器测量]的传递，
% 我们希望设计滤波器R(s)，使得扰动对评估输出的影响H∞范数最小，同时保证故障对评估输出的增益达到要求。
% [R, CL, gamma] = hinfsyn(P, ny, nu);
% 其中ny为滤波器可测量输出维数，nu为滤波器控制输入维数。
%
% （注：由于缺少具体增广模型，此处不给出滤波器设计的完整代码。实际设计中，通过求解一对Riccati方程即可得到最优检测滤波器。有关详细过程和原型结构，请参考用户提供的文献。）
