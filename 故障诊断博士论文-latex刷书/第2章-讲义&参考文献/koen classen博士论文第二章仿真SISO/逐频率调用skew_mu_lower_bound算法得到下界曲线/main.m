%% ===============================================================
%  Skew-μ lower bound demo for the system in the screenshots
%  Δ = diag(δ1, Δ1(s)),  δ1 ∈ R, |δ1|≤1;  Δ1(s) dynamic, ||Δ1||∞ ≤ 1
%  \tilde{G}_d(s,Δ) = (7 G1(s)(0.25+δ1) + 3 G2(s)(1-δ1) + 3 G3(s)(0.75-δ1)) * (1 + W(s) Δ1(s))
%  This script:
%    1) builds the uncertain system (uss)
%    2) extracts M(s) via lftdata (the Δ–Δ interconnection)
%    3) runs Holland-style skew-μ lower bound per frequency
% ===============================================================

clear; close all; clc;

% -------- transfer functions --------
s = tf('s');

G1 = (1/pi * s) / (s + 2*pi)^2;        % corresponds to text in your notes
G2 = (1/(10*pi) * s) / (s + 20*pi)^2;
G3 = (1/(100*pi) * s) / (s + 200*pi)^2;
% 换成下面这个传递函数出的图很奇怪
% G1 = 0.36*100*(1/(1*pi )* s)/(s + 2*pi)^2;  
% G2 = 0.36*100*(10/(1*pi) * s)/(s + 20*pi)^2; 
% G3 = 0.36*100*(100/(1*pi) * s)/(s + 200*pi)^2; 

W  = 0.2*(s + 10*pi)/(s + 20*pi);      % dynamic weight，作用是在不同频率点给动态不确定性加权，这样才能显示出所谓的动态来

% -------- uncertainties: δ1 real scalar; Δ1 dynamic complex with ||.||∞≤1 --------
delta1 = ureal('delta1', 0, 'Range', [-1, 1]);    % parametric (real) block
Delta1 = ultidyn('Delta1', [1 1], 'Bound', 1);    % dynamic complex block (scalar)

% -------- uncertain plant (your screenshot’s structure) --------
% \tilde G_d(s,Δ) = (...δ1...) * (1 + W(s) Δ1(s))
polyA = 7*G1*(0.25 + delta1) + 3*G2*(1 - delta1) + 3*G3*(0.75 - delta1);
Gt    = polyA * (1 + W*Delta1);                   % uss

% % -------- extract M(s) mapping for Δ-channels via lftdata --------
% % If T(s,Δ) = lft(Δ, M(s)), lftdata(T) returns M(s) & the Δ block info.
% [Mlin, DeltaBlk, BlkStruct] = lftdata(Gt);
% % Here, M(s) is a 2x2 LTI (ss) because we have 2 scalar Δ blocks (δ1, Δ1)
% 
% % -------- partition sizes for skew-μ: fixed = δ1(real), varying = Δ1(complex) ----
% k_f = [1 0 0];   % fixed side:   1 real, 0 complex, 0 full
% k_v = [0 1 0];   % varying side: 0 real, 1 complex, 0 full
% 
% % -------- frequency grid (rad/s) --------
% % Screenshot uses Hz on x-axis; pick 0.1–1e3 Hz => multiply by 2π
% fHz = logspace(-1, 3, 120);
% w   = 2*pi*fHz;
% 
% % -------- run skew-μ lower bound per frequency --------
% opts = struct('maxIter', 200, 'tol', 1e-3, 'verbose', 0);
% Lw = zeros(size(w));
% for k = 1:numel(w)
%     Mjw = freqresp(Mlin, w(k));      % complex 2x2 at this freq
%     Mjw = squeeze(Mjw);
%     out = skew_mu_lb_onefreq(Mjw, k_f, k_v, opts);
%     Lw(k) = out.mu_lb;               % lower bound at this ω
% end
% 
% % -------- find worst-case frequency by lower bound --------
% [~, idx_wc] = max(Lw);
% f_wc = fHz(idx_wc);
% 
% % -------- plots --------
% figure; 
% semilogx(fHz, 20*log10(Lw), 'LineWidth', 1.8); grid on;
% xlabel('Frequency (Hz)');
% ylabel('Lower bound L(\omega)  [dB]');
% title('Skew-\mu Lower Bound vs Frequency (Holland-style)');
% hold on; 
% plot(f_wc, 20*log10(Lw(idx_wc)), 'kx', 'MarkerSize', 9, 'LineWidth', 1.5);
% text(f_wc, 20*log10(Lw(idx_wc))+0.8, sprintf('\\leftarrow \\omega_{wc} \\approx %.3g Hz', f_wc));


% -------- extract M(s) mapping for Δ-channels via lftdata --------
[Mlin, DeltaBlk, BlkStruct] = lftdata(Gt);

% ===== 自动识别并统计块维度（关键修正）=====
% 将 ureal 视为 fixed（实重复），ultidyn 视为 varying（复重复）
idx_fixed = []; idx_vary = [];
nf = 0; nv = 0; off = 0;

for i = 1:numel(BlkStruct)
    bs  = BlkStruct(i).Size;        % e.g. [m n]，标量块就是 [1 1]
    dim = bs(1);                         % 标量块取 1，重复实标量会被展开成 dim 次
    typ = lower(BlkStruct(i).Type);      % 'ureal' 或 'ultidyn'（也可能是 'ucomplex' 等）

    % 这一步构造“原始通道索引区间”
    ch_range = off + (1:dim);

    switch typ
        case 'ureal'
            idx_fixed = [idx_fixed, ch_range]; %#ok<AGROW>
            nf = nf + dim;                     % 记入 fixed 维数
        case {'ultidyn','ucomplex','ucomplexm'}  % 统统当作“复重复”varying
            idx_vary  = [idx_vary,  ch_range]; %#ok<AGROW>
            nv = nv + dim;
        otherwise
            error('Unsupported uncertainty type: %s', BlkStruct(i).Type);
    end
    off = off + dim;
end


% % 现在总通道数应该等于 size(Mjw,1)（也即 size(Mlin, 'Frequency response')）
% % 我们需要把 M 的通道顺序“变换”为 [fixed, varying]
% perm = [idx_fixed, idx_vary];            % 目标顺序：先 fixed 再 varying
% 
% % -------- frequency grid --------
% fHz = logspace(-1, 3, 120);
% w   = 2*pi*fHz;
% 
% opts = struct('maxIter', 200, 'tol', 1e-3, 'verbose', 0);
% Lw = zeros(size(w));
% 
% % 幂迭代的 k_f/k_v：把所有 fixed 统为“实重复”，varying 统为“复重复”
% k_f = [nf, 0, 0];        % nf 个 real-repeated
% k_v = [0,  nv, 0];       % nv 个 complex-repeated
% 
% for k = 1:numel(w)
%     Mjw = freqresp(Mlin, w(k));      % 得到  (n×n×1×1)
%     Mjw = squeeze(Mjw);              % 变成 n×n 复矩阵
% 
%     % ===== 通道重排到 [fixed, varying]（关键修正）=====
%     P   = eye(nf+nv);
%     P   = P(:, perm);                % 把旧顺序的通道列重排成 [fixed, varying]
%     Mjw = P' * Mjw * P;
% 
%     out = skew_mu_lb_onefreq(Mjw, k_f, k_v, opts);
%     Lw(k) = out.mu_lb;
% end
% 
% [~, idx_wc] = max(Lw);
% f_wc = fHz(idx_wc);
% 
% figure;
% semilogx(fHz, 20*log10(Lw), 'LineWidth', 1.8); grid on;
% xlabel('Frequency (Hz)'); ylabel('Lower bound L(\omega) [dB]');
% title('Skew-\mu Lower Bound vs Frequency (Auto-detected block sizes)');
% hold on;
% plot(f_wc, 20*log10(Lw(idx_wc)), 'kx', 'MarkerSize', 9, 'LineWidth', 1.5);
% text(f_wc, 20*log10(Lw(idx_wc))+0.8, sprintf('\\leftarrow \\omega_{wc} \\approx %.3g Hz', f_wc));

% % ===== 自动识别块维度 =====
% idx_fixed = []; idx_vary = [];
% nf = 0; nv = 0; off = 0;
% 
% for i = 1:numel(BlkStruct)
%     bs  = BlkStruct(i).Size;         % 每个块的尺寸
%     dim = bs(1);
%     typ = lower(BlkStruct(i).Type);
% 
%     ch_range = off + (1:dim);             % 通道编号区间
%     switch typ
%         case 'ureal'
%             idx_fixed = [idx_fixed, ch_range]; %#ok<AGROW>
%             nf = nf + dim;
%         case {'ultidyn','ucomplex','ucomplexm'}
%             idx_vary  = [idx_vary,  ch_range]; %#ok<AGROW>
%             nv = nv + dim;
%         otherwise
%             error('Unsupported block type: %s', BlkStruct(i).Type);
%     end
%     off = off + dim;
% end
% 
% % 总通道数
% nM = size(Mlin.B,2);     % 直接取 LTI 系统的 Δ 输入数
% nDelta = nf + nv;
% if nM ~= nDelta
%     warning('检测到 M 矩阵维度 %d 与 nf+nv=%d 不符，自动以 nM 为准', nM, nDelta);
%     nDelta = nM;
% end
% 
% perm = [idx_fixed, idx_vary];             % 通道重排索引
% perm = perm(:)';                          % 强制成行向量
% perm = perm(1:nDelta);                    % 截取到正确长度

% ===== 自动识别块维度 (修正版：考虑重复块) =====
idx_fixed = []; idx_vary = [];
nf = 0; nv = 0; off = 0;

for i = 1:numel(BlkStruct)
    bs  = BlkStruct(i).Size;         % [m n]，例如 [2 2] 表示重复两次
    dim = bs(1);                          % 通道数量
    typ = lower(BlkStruct(i).Type);
    ch_range = off + (1:dim);             % 通道区间
    
    switch typ
        case 'ureal'   % 实重复块
            idx_fixed = [idx_fixed, ch_range]; %#ok<AGROW>
            nf = nf + dim;
        case {'ultidyn','ucomplex','ucomplexm'} % 复或动态块
            idx_vary = [idx_vary, ch_range]; %#ok<AGROW>
            nv = nv + dim;
        otherwise
            warning('忽略未知类型块 %s', typ);
    end
    off = off + dim;
end

% 实际 M 的总通道数
nM = size(Mlin.B,2);
nDelta = nf + nv;

if nM ~= nDelta
    warning('检测到 M 矩阵维度 %d 与 nf+nv=%d 不符，自动调整为 nM=%d', nM, nDelta, nM);
    nDelta = nM;
    % 如果 idx_fixed+idx_vary 通道总数 < nM，则补齐剩余索引
    all_idx = 1:nM;
    known_idx = [idx_fixed idx_vary];
    extra_idx = setdiff(all_idx, known_idx, 'stable');
    idx_vary = [idx_vary extra_idx];  % 把多余通道作为 varying
    nv = numel(idx_vary);
end

% 通道重排索引：[fixed, varying]
perm = [idx_fixed, idx_vary];
perm = perm(:)';                      % 转行向量
perm = perm(1:nDelta);                % 防止越界
fprintf('检测结果: nf=%d (fixed), nv=%d (varying), M size=%dx%d\n', nf, nv, size(Mlin.B,2), size(Mlin.B,2));
disp('perm ='); disp(perm);           % fixed varying 通道总数


% ===== 扫频 =====
fHz = logspace(-1, 3, 120);
w   = 2*pi*fHz;
opts = struct('maxIter', 200, 'tol', 1e-3, 'verbose', 0);

Lw = zeros(size(w));
k_f = [nf 0 0];
k_v = [0 nv 0];

for k = 1:numel(w)
    Mjw = squeeze(freqresp(Mlin, w(k)));  % n×n 复矩阵
    % 保证 Mjw 是二维方阵
    sz = size(Mjw);
    if numel(sz) > 2
        Mjw = reshape(Mjw, sz(1), sz(2));
    end
    % 如果 perm 与 Mjw 维度不匹配，则只截取交集
    nM = size(Mjw,1);
    perm_safe = perm(perm <= nM);

    % 构造置换矩阵 P（稀疏安全方式）
    P = speye(nM);
    P = P(:,perm_safe);

    % 重排 M 的通道顺序：[fixed, varying]
    Mjw = P' * Mjw * P;

    % 调用 skew-μ 幂迭代
    out = skew_mu_lb_onefreq(Mjw, k_f, k_v, opts);
    Lw(k) = out.mu_lb;
end

[~, idx_wc] = max(Lw);
f_wc = fHz(idx_wc);

figure;
semilogx(fHz, 20*log10(Lw), 'LineWidth', 1.8);
grid on;
xlabel('Frequency (Hz)');
ylabel('Lower bound L(\omega) [dB]');
title('Skew-\mu Lower Bound vs Frequency (auto partition)');
hold on;
plot(f_wc, 20*log10(Lw(idx_wc)), 'kx', 'MarkerSize', 9, 'LineWidth', 1.5);
text(f_wc, 20*log10(Lw(idx_wc))+1, sprintf('\\leftarrow \\omega_{wc}\\approx %.2f Hz',f_wc));
