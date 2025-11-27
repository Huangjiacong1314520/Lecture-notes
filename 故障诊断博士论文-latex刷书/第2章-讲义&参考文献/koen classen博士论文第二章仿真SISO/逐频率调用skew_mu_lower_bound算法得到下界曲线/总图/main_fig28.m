%% ===============================================================
%  Figure 2.8 style replication (nominal / lower bound / worst-case init)
%  System: \tilde{G}_d(s,Δ) = (7G1(0.25+δ1)+3G2(1-δ1)+3G3(0.75-δ1)) * (1 + W Δ1(s))
%  Δ = diag(δ1, Δ1),  δ1 ∈ R, |δ1|≤1;  Δ1 dynamic complex, ||Δ1||∞≤1
% ===============================================================
clear; close all; clc;

%% 1) Build uncertain plant (uss)
s = tf('s');

% G1 = (1/pi * s) / (s + 2*pi)^2;
% G2 = (1/(10*pi) * s) / (s + 20*pi)^2;
% G3 = (1/(100*pi) * s) / (s + 200*pi)^2;
G1 = 0.36*100*(1/(1*pi )* s)/(s + 2*pi)^2;  
G2 = 0.36*100*(10/(1*pi) * s)/(s + 20*pi)^2; 
G3 = 0.36*100*(100/(1*pi) * s)/(s + 200*pi)^2; 
W  = 0.2*(s + 10*pi)/(s + 20*pi);

delta1 = ureal('delta1', 0, 'Range', [-1, 1]);     % param block
Delta1 = ultidyn('Delta1', [1 1], 'Bound', 1);     % dynamic complex

polyA  = 7*G1*(0.25 + delta1) + 3*G2*(1 - delta1) + 3*G3*(0.75 - delta1);
Gtilde = polyA * (1 + W*Delta1);                   % uss

%% 2) Extract M(s) seen by Δ via lftdata
[Mlin, ~, BlkStruct] = lftdata(Gtilde);   % Gtilde = lft(Δ, Mlin)

% ---- auto-detect block sizes (fixed=ureal, varying=ultidyn/ucomplex) ----
idx_fixed = []; idx_vary = [];
nf = 0; nv = 0; off = 0;
% 对应论文对块结构的定义：Kf​(mrf​,mcf​,mCf​),Kv​(mrv​,mcv​,mCv​)，这里为了对不同的传递函数都适用，所以用了一个循环自动计算块结构
for i = 1:numel(BlkStruct)
    dim = BlkStruct(i).Size(1);        % repeated count
    typ = lower(BlkStruct(i).Type);
    ch  = off + (1:dim);
    switch typ
        case 'ureal'
            idx_fixed = [idx_fixed, ch]; nf = nf + dim;
        case {'ultidyn','ucomplex','ucomplexm'}
            idx_vary  = [idx_vary,  ch]; nv = nv + dim;
        otherwise
            error('Unsupported uncertainty type: %s', BlkStruct(i).Type);
    end
    off = off + dim;
end
nM = size(Mlin.B,2);
if nf+nv ~= nM
    % 补齐遗漏通道（极少见），统归 varying
    extra = setdiff(1:nM, [idx_fixed idx_vary], 'stable');
    idx_vary = [idx_vary extra];
    nv = numel(idx_vary);
end
perm = [idx_fixed, idx_vary];  % reorder to [fixed, varying]

% ---- block vectors for power iteration (we treat all fixed as real-repeated; all varying as complex-repeated) ----
k_f = [nf 0 0];
k_v = [0  nv 0];

%% 3) Frequency sweep: lower bound L(ω) via skew-μ power iteration
%  论文算法是静态的（针对单个M矩阵），而实际系统M(jω)随频率变化，因此需在每个 ω 点调用一次幂迭代。
fHz = logspace(-1, 3, 160);   % 0.1~1e3 Hz
w   = 2*pi*fHz; %角频率点
opts = struct('maxIter', 200, 'tol', 1e-3, 'verbose', 0); %设置幂迭代选项

Lw       = zeros(size(w));    % 初始化数组
iterOuts = cell(size(w));     % 存储各频点迭代结果 store iteration info (to recover phase/sign)。Lw(k) 保存每个 ω 的 μₛ 下界 𝐿 ( 𝜔 ) L(ω)，iterOuts{k} 保存收敛时的 a,w 向量（论文式 (21) 的 a,w 对应右、左特征向量，用于恢复 Δ₍wc₎ 相位）。
% 逐频率点调用skew_mu_lower_bound_Holland函数进行频率扫描
for k = 1:numel(w)
    Mjw = squeeze(freqresp(Mlin, w(k)));   % n×n 复矩阵 M(jω) 返回在由向量 w(k) 指定的实频率网格上的频率响应。论文算法输入即为某频率下的复矩阵𝑀∈𝐶𝑛×𝑛
    % channel reorder to [fixed, varying]
    % 对Mjw的行列重新排序，使其分块结构为[M11 M12; M21 M22]      %关于为什么要重新排序：在 MATLAB 中，用
    % lftdata 抽取 M(s) 后，得到的矩阵通道顺序往往是按照不确定块出现顺序自动生成的。而Δ = diag(D₍f₎, D₍v₎),
    % where D₍f₎ is the fixed-range part and D₍v₎ is the variable part.”要求把把矩阵 M 相应分成块
    P = speye(nM); P = P(:, perm);
    Mjw = P' * Mjw * P;
    % 调用实现幂迭代算法
    out = skew_mu_lower_bound_Holland(Mjw, k_f, k_v, opts);  % returns out.mu_lb and out.last vectors
    Lw(k)        = out.mu_lb;    % 保存该频率的下界 L(ω)
    iterOuts{k}  = out;          % 保存最后的 a,w 向量用于 Δwc 重构                 % save a,w for wc reconstruction
end

% worst-case frequency by lower bound (for demo)
[~, idx_wc] = max(Lw);
f_wc = fHz(idx_wc);
fprintf('Estimated worst-case frequency: %.4g Hz\n', f_wc);

%% 4) Reconstruct Δ_wc from iteration info at ω_wc
out_wc = iterOuts{idx_wc};
nf_loc = sum(k_f); nv_loc = sum(k_v);
a_wc = out_wc.last.a;     % last right stage vector
w_wc = out_wc.last.w;

% -- δ1_wc (real repeated): take global sign of sum(Re(a1 .* conj(w1)))
a1 = a_wc(1:nf_loc);
w1 = w_wc(1:nf_loc);
sig = sign(real(sum(a1 .* conj(w1))));
if sig==0, sig=1; end
delta1_wc = sig;             % ±1

% -- Δ1_wc (complex repeated): take phase of inner product a2' * w2
a2 = a_wc(nf_loc+1:end);
w2 = w_wc(nf_loc+1:end);
th = angle(a2' * w2);
Delta1_wc = exp(1j*th);      % static complex gain, |.|=1

fprintf('Recovered delta1_wc = %+d,  angle(Delta1_wc) = %.2f deg\n', delta1_wc, rad2deg(th));

%% 5) Build curves: nominal, lower bound, and worst-case-init candidate
% Nominal (Δ=0)
G_nom = (7*G1*(0.25 + 0) + 3*G2*(1 - 0) + 3*G3*(0.75 - 0)) * (1 + W*0);

% Worst-case-init candidate using (delta1_wc, Delta1_wc) as static choices
G_init = (7*G1*(0.25 + delta1_wc) + 3*G2*(1 - delta1_wc) + 3*G3*(0.75 - delta1_wc)) ...
         * (1 + W * tf(Delta1_wc));   % tf of a complex static gain

% Magnitudes
[mag_nom, ~]  = bode(G_nom,  w);  mag_nom  = squeeze(mag_nom);
[mag_init, ~] = bode(G_init, w);  mag_init = squeeze(mag_init);

% % dB conversion
% L_dB    = 20*log10(Lw);
% nom_dB  = 20*log10(mag_nom);
% init_dB = 20*log10(mag_init);
% 
% %% 6) Plot (Fig. 2.8 style)
% figure('Color','w'); hold on; grid on; box on;
% % lower bound as small squares
% plot(fHz, L_dB, 'ks', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'DisplayName','Lower bound L(\omega)');


% % nominal
% semilogx(fHz, nom_dB,  'b-', 'LineWidth', 1.6, 'DisplayName','\tilde{G}_d(s,0)');
% % worst-case-init candidate
% semilogx(fHz, init_dB, 'Color',[1 0.5 0], 'LineWidth', 2.0, 'DisplayName','\tilde{G}_{d,init}(s)');
% % mark ω_wc
% plot(f_wc, L_dB(idx_wc), 'kx', 'MarkerSize', 10, 'LineWidth', 1.4, 'DisplayName','\omega_{wc}');
% 
% xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
% title('Nominal / Lower bound / Worst-case-init (Holland-style)');
% legend('Location','best');
% ylim([min([nom_dB; init_dB; L_dB])-5, max([nom_dB; init_dB; L_dB])+5]);
% --- convert to dB and ensure column vectors ---
L_dB    = 20*log10(abs(Lw(:)));        % column
nom_dB  = 20*log10(abs(squeeze(mag_nom(:))));
init_dB = 20*log10(abs(squeeze(mag_init(:))));

% --- ensure same length (truncate/pad if needed) ---
Nmin = min([numel(L_dB), numel(nom_dB), numel(init_dB)]);
L_dB    = L_dB(1:Nmin);
nom_dB  = nom_dB(1:Nmin);
init_dB = init_dB(1:Nmin);
fHz     = fHz(1:Nmin);

%% 6) Plot (Fig. 2.8 style)
figure('Color','w'); hold on; grid on; box on;

% lower bound markers
plot(fHz, L_dB, 'ks', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'DisplayName','Lower bound L(\omega)');
% nominal
semilogx(fHz, nom_dB,  'b-', 'LineWidth', 1.6, 'DisplayName','\tilde{G}_d(s,0)');
% worst-case-init
semilogx(fHz, init_dB, 'Color',[1 0.5 0], 'LineWidth', 2.0, 'DisplayName','\tilde{G}_{d,init}(s)');
% mark ω_wc
plot(f_wc, L_dB(idx_wc), 'kx', 'MarkerSize', 10, 'LineWidth', 1.4, 'DisplayName','\omega_{wc}');

xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Nominal / Lower bound / Worst-case-init (Holland-style)');
legend('Location','best');

ylim([min([L_dB; nom_dB; init_dB])-5, max([L_dB; nom_dB; init_dB])+5]);
