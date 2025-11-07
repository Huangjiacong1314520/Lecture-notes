% 第三版，修改了第二版的部分错误
function [mu_s_lower, v_total, history, b, a, z, w] = skew_mu_lower_bound(M, k_f, k_v, max_iter, tol)
% SKEW_MU_LOWER_BOUND - 严格按照 Holland (2005) 论文实现的 Skew μ 下界算法
%
% 本实现修正了所有关键问题：
%   1. v 公式：使用 ||g2||^2 而非 ||g2|| (论文方程15/16)
%   2. 缩放方向：M <- S_L * M * S_R (论文方程23/24)
%   3. q 参数：使用 sign(Re(a*w)) (论文方程14)
%   4. 算法顺序：右迭代 -> 左迭代 -> 缩放 (论文4.2节)
%   5. 终止条件：检测 v < 1 和 ||g1|| >= 1 (论文4.3节)
%   6. 数值稳定：使用 \ 而非 inv()
%   7. 归一化：每轮保持 |a|=|w|=|b|=|z|=1
%   8. 历史记录：保存 v, g1, c1 用于调试
%
% Inputs:
%   M        - 系统矩阵 (n x n 复数矩阵)
%   k_f      - 固定扰动块结构 [k1_f, k2_f, k3_f]
%              k1_f: 实重复标量块大小
%              k2_f: 复重复标量块大小
%              k3_f: 全复数块大小
%   k_v      - 变化扰动块结构 [k1_v, k2_v, k3_v]
%   max_iter - 最大迭代次数 (默认: 100)
%   tol      - 收敛容差 (默认: 1e-6)
%
% Outputs:
%   mu_s_lower - skew μ 的下界
%   v_total    - 累积缩放因子
%   history    - 迭代历史 (v, g1_norm, c1_norm, v_star, v_sharp)
%   b, a, z, w - 收敛的特征向量

%% ========================================================================
%  初始化与参数检查
%  ========================================================================
if nargin < 4, max_iter = 100; end
if nargin < 5, tol = 1e-6; end

% 提取块尺寸
k1_f = k_f(1); k2_f = k_f(2); k3_f = k_f(3);
k1_v = k_v(1); k2_v = k_v(2); k3_v = k_v(3);

n_f = k1_f + k2_f + k3_f;  % 固定扰动总维数
n_v = k1_v + k2_v + k3_v;  % 变化扰动总维数
n = n_f + n_v;

if size(M,1) ~= n || size(M,2) ~= n
    error('Matrix M dimension (%d×%d) does not match structure n_f+n_v=%d', ...
          size(M,1), size(M,2), n);
end

% 初始化历史记录结构
history.v = [];
history.g1_norm = [];
history.c1_norm = [];
history.v_star = [];
history.v_sharp = [];
history.converged = false;

% 初始化特征向量（随机化并归一化）
b = randn(n,1) + 1i*randn(n,1); b = b / norm(b);
z = randn(n,1) + 1i*randn(n,1); z = z / norm(z);
w = randn(n,1) + 1i*randn(n,1); w = w / norm(w);

% 初始化缩放
v_total = 1;
M_scaled = M;

fprintf('=== Skew μ Lower Bound Algorithm (Holland 2005) ===\n');
fprintf('Fixed blocks: n_f=%d [real:%d, complex:%d, full:%d]\n', ...
        n_f, k1_f, k2_f, k3_f);
fprintf('Varying blocks: n_v=%d [real:%d, complex:%d, full:%d]\n', ...
        n_v, k1_v, k2_v, k3_v);
fprintf('Starting iterations...\n\n');

%% ========================================================================
%  主循环：Power Algorithm (论文 Algorithm 4.2)
%  ========================================================================
for iter = 1:max_iter
    
    %% ====================================================================
    %  阶段 I: 右特征向量迭代
    %  目标: 从 b_k 计算 a_{k+1}, 并得到 v_star
    %  ====================================================================
    
    % 步骤 1: 计算 g = M * b (论文方程 S_{k+1}*a_{k+1} = M*b_k)
    g = M_scaled * b;
    
    % 步骤 2: 分割为固定/变化部分 (论文方程18)
    g1 = g(1:n_f);       % 固定扰动部分
    g2 = g(n_f+1:end);   % 变化扰动部分
    
    g1_norm = norm(g1);
    g2_norm = norm(g2);
    
    % 步骤 3: 检查 ||g1||^2 >= 1 (论文4.3节，固定扰动能单独失稳)
    if g1_norm >= 1
        fprintf('[Iter %d] ||g1|| = %.6f >= 1: Fixed perturbation destabilizes!\n', ...
                iter, g1_norm);
        % 论文方程(23): 应用缩放 M <- S_L * M * S_R
        % 这里 S_L = diag(I_f, 1/√2*I_v), S_R = diag(I_f, √2*I_v)
        v_total = 2 * v_total;
        S_L = blkdiag(eye(n_f), 1/sqrt(2)*eye(n_v));
        S_R = blkdiag(eye(n_f), sqrt(2)*eye(n_v));
        M_scaled = S_L * M_scaled * S_R;
        
        % 重新归一化向量（因为坐标系变了）
        b = b / norm(b);
        z = z / norm(z);
        w = w / norm(w);
        continue;
    end
    
    % 步骤 4: 计算 v_star (论文定理3，方程15) - 修正：分子是 ||g2||^2
    if g1_norm < 1
        v_star = (g2_norm^2) / sqrt(1 - g1_norm^2);
    else
        v_star = inf;
    end
    
    % 步骤 5: 构造 S_star 并计算 a (论文方程11)
    S_star = blkdiag(eye(n_f), v_star*eye(n_v));
    a = S_star \ g;  % 使用 \ 而非 inv() 提高数值稳定性
    a = a / norm(a); % 归一化
    
    % 步骤 6: 分块向量 (论文方程10)
    [a1_f, a2_f, a3_f] = partition_vector(a, 1:n_f, k1_f, k2_f, k3_f);
    [a1_v, a2_v, a3_v] = partition_vector(a, n_f+1:n, k1_v, k2_v, k3_v);
    [w1_f, w2_f, w3_f] = partition_vector(w, 1:n_f, k1_f, k2_f, k3_f);
    [w1_v, w2_v, w3_v] = partition_vector(w, n_f+1:n, k1_v, k2_v, k3_v);
    
    % 步骤 7: 计算 q_star (论文方程14) - 修正：使用 sign(Re(a*w))
    q_star_f = compute_q_correct(a1_f, w1_f);
    q_star_v = compute_q_correct(a1_v, w1_v);
    
    % 步骤 8: 更新 z (论文方程13，右侧)
    z = zeros(n, 1);
    z = assemble_vector(z, 1:n_f, ...
        update_real_block(q_star_f, w1_f), ...
        update_complex_block_correct(w2_f, a2_f), ...
        update_full_block_correct(w3_f, a3_f));
    z = assemble_vector(z, n_f+1:n, ...
        update_real_block(q_star_v, w1_v), ...
        update_complex_block_correct(w2_v, a2_v), ...
        update_full_block_correct(w3_v, a3_v));
    z = z / norm(z); % 归一化
    
    %% ====================================================================
    %  阶段 II: 左特征向量迭代
    %  目标: 从 z_{k+1} 计算 w_{k+1}, 并得到 v_sharp
    %  ====================================================================
    
    % 步骤 1: 计算 c = M^* * z (论文方程 S_{k+1}*w_{k+1} = M^**z_{k+1})
    c = M_scaled' * z;
    
    % 步骤 2: 分割
    c1 = c(1:n_f);
    c2 = c(n_f+1:end);
    
    c1_norm = norm(c1);
    c2_norm = norm(c2);
    
    % 步骤 3: 计算 v_sharp (论文定理3，方程16) - 修正：分子是 ||c2||^2
    if c1_norm < 1
        v_sharp = (c2_norm^2) / sqrt(1 - c1_norm^2);
    else
        v_sharp = inf;
    end
    
    % 步骤 4: 构造 S_sharp 并计算 w
    S_sharp = blkdiag(eye(n_f), v_sharp*eye(n_v));
    w_new = S_sharp \ c;
    w_new = w_new / norm(w_new);
    
    % 步骤 5: 分块新的 w
    [w1_f, w2_f, w3_f] = partition_vector(w_new, 1:n_f, k1_f, k2_f, k3_f);
    [w1_v, w2_v, w3_v] = partition_vector(w_new, n_f+1:n, k1_v, k2_v, k3_v);
    
    % 步骤 6: 计算 q_sharp
    q_sharp_f = compute_q_correct(a1_f, w1_f);
    q_sharp_v = compute_q_correct(a1_v, w1_v);
    
    % 步骤 7: 更新 b (论文方程13，左侧)
    b_new = zeros(n, 1);
    b_new = assemble_vector(b_new, 1:n_f, ...
        update_real_block(q_sharp_f, a1_f), ...
        update_complex_block_correct(a2_f, w2_f), ...
        update_full_block_correct(w3_f, a3_f));
    b_new = assemble_vector(b_new, n_f+1:n, ...
        update_real_block(q_sharp_v, a1_v), ...
        update_complex_block_correct(a2_v, w2_v), ...
        update_full_block_correct(w3_v, a3_v));
    b_new = b_new / norm(b_new);
    
    %% ====================================================================
    %  阶段 III: 收敛判定与缩放更新
    %  ====================================================================
    
    % 计算 v 的平均值和差异
    v_avg = (v_star + v_sharp) / 2;
    v_diff = abs(v_star - v_sharp);
    
    % 记录历史
    history.v(iter) = v_avg;
    history.g1_norm(iter) = g1_norm;
    history.c1_norm(iter) = c1_norm;
    history.v_star(iter) = v_star;
    history.v_sharp(iter) = v_sharp;
    
    % 打印当前状态
    fprintf('[Iter %3d] v*=%.6f, v#=%.6f, v_avg=%.6f, ||g1||=%.4f, ||c1||=%.4f, v_total=%.6f\n', ...
            iter, v_star, v_sharp, v_avg, g1_norm, c1_norm, v_total);
    
    % 【修正5】终止条件1: v < 1 立即停止 (论文4.3节)
    if v_avg < 1
        fprintf('\n>>> Converged: v_avg < 1, returning mu_s = v_total = %.6f\n', v_total);
        history.converged = true;
        mu_s_lower = v_total;
        b = b_new;
        w = w_new;
        return;
    end
    
    % 终止条件2: v_star 与 v_sharp 收敛 (论文备注)
    if v_diff < tol
        fprintf('\n>>> Converged: |v* - v#| < tol\n');
        
        % 当 v >= 1 时，需要继续缩放 (论文方程24)
        if v_avg >= 1
            % 更新累积因子
            v_total = v_avg * v_total;
            
            % 论文方程(24): M <- S_L * M * S_R
            % S_L = diag(I_f, 1/√v*I_v), S_R = diag(I_f, √v*I_v)
            S_L = blkdiag(eye(n_f), 1/sqrt(v_avg)*eye(n_v));
            S_R = blkdiag(eye(n_f), sqrt(v_avg)*eye(n_v));
            M_scaled = S_L * M_scaled * S_R;
            
            fprintf('    Scaling M with v=%.6f, new v_total=%.6f\n', v_avg, v_total);
            
            % 检查是否需要继续迭代
            if v_total > 1e6
                fprintf('\n>>> v_total >> 1: Fixed perturbation destabilizes (mu_s = inf)\n');
                mu_s_lower = inf;
                history.converged = true;
                return;
            end
            
            % 重置收敛标志，继续迭代
            % (因为缩放后需要重新寻找新的特征向量)
        else
            % v < 1 且已收敛，返回结果
            mu_s_lower = v_total;
            history.converged = true;
            b = b_new;
            w = w_new;
            return;
        end
    end
    
    % 更新向量进入下一轮
    b = b_new;
    w = w_new;
    
end

%% ========================================================================
%  达到最大迭代次数
%  ========================================================================
warning('Maximum iterations (%d) reached without full convergence', max_iter);
% 论文备注: 即使未完全收敛，max(v*, v#) 仍是有效下界
mu_s_lower = max(v_star, v_sharp) * v_total;
fprintf('\n>>> Max iter reached: returning mu_s = %.6f (may be conservative)\n', mu_s_lower);

end

%% ========================================================================
%  辅助函数：向量分块
%  ========================================================================
function [v1, v2, v3] = partition_vector(v, range, k1, k2, k3)
% 将向量 v(range) 按块大小 k1, k2, k3 分割
v_sub = v(range);
v1 = v_sub(1:k1);
v2 = v_sub(k1+1:k1+k2);
v3 = v_sub(k1+k2+1:end);
end

function v = assemble_vector(v, range, v1, v2, v3)
% 将子块 v1, v2, v3 组装回 v(range)
idx = range(1) - 1;
if ~isempty(v1), v(idx+1:idx+length(v1)) = v1; idx = idx+length(v1); end
if ~isempty(v2), v(idx+1:idx+length(v2)) = v2; idx = idx+length(v2); end
if ~isempty(v3), v(idx+1:idx+length(v3)) = v3; end
end

%% ========================================================================
%  辅助函数：q 参数计算 (论文方程14) - 修正版
%  ========================================================================
function q = compute_q_correct(a1, w1)
% 论文公式: q = sign(Re(a^* w))
% 条件 (14): Re(a^* w) ≤ 0 for q=-1; ≥0 for q=+1; =0 for |q|<1
if isempty(a1) || isempty(w1)
    q = 1;
    return;
end

% 计算 Re(a^* * w)
s = real(a1' * w1);

% 取符号
if abs(s) < 1e-12
    q = 1;  % 默认 +1
else
    q = sign(s);
end

% 确保 q 在 [-1, 1]
q = max(-1, min(1, q));
end

%% ========================================================================
%  辅助函数：块更新 (论文方程13) - 修正版
%  ========================================================================
function z_out = update_real_block(q, w)
% 实标量块: z = q * w
if isempty(w)
    z_out = [];
else
    z_out = q * w;
end
end

function z_out = update_complex_block_correct(w, a)
% 复标量块: z = (w^* a / |w^* a|) * w
% 修正：内积方向应为 w^* * a (论文方程13)
if isempty(w) || isempty(a)
    z_out = [];
    return;
end

inner = w' * a;  % w^* * a
if abs(inner) > 1e-12
    z_out = (inner / abs(inner)) * w;
else
    z_out = w;
end
end

function z_out = update_full_block_correct(w, a)
% 全复数块: z = (|w| / |a|) .* a
% 论文方程13: 逐元素比例缩放
if isempty(w) || isempty(a)
    z_out = [];
    return;
end

% 防止除零
a_abs = abs(a);
a_abs(a_abs < 1e-12) = 1e-12;

z_out = (abs(w) ./ a_abs) .* a;
end

%% ========================================================================
%  测试函数
%  ========================================================================
function test_skew_mu_corrected()
% 测试用例：4×4 系统，2固定 + 2变化

fprintf('\n========== Test Case 1: Simple 4x4 system ==========\n');

% 构造测试矩阵
M = [0.5+0.1i, 0.2,      0.1,      0.05;
     0.2,      0.6-0.1i, 0.15,     0.1;
     0.1,      0.15,     0.4+0.2i, 0.2;
     0.05,     0.1,      0.2,      0.5-0.1i];

% 块结构: 各1个实标量块 + 1个复标量块 (固定和变化各一组)
k_f = [1, 1, 0];  % 固定: 1实 + 1复 + 0满
k_v = [1, 1, 0];  % 变化: 1实 + 1复 + 0满

[mu_s, v_total, history] = skew_mu_lower_bound(M, k_f, k_v, 50, 1e-6);

fprintf('\n=== Final Results ===\n');
fprintf('Skew μ lower bound: %.6f\n', mu_s);
fprintf('Total scaling factor: %.6f\n', v_total);
fprintf('Converged: %s\n', string(history.converged));
fprintf('Total iterations: %d\n', length(history.v));

% 绘制收敛曲线
if length(history.v) > 1
    figure('Name', 'Skew μ Convergence History');
    subplot(2,1,1);
    semilogy(1:length(history.v), history.v_star, 'b-o', ...
             1:length(history.v), history.v_sharp, 'r-x', ...
             1:length(history.v), history.v, 'k--');
    xlabel('Iteration'); ylabel('v value');
    legend('v^*', 'v^\#', 'v_{avg}');
    title('v convergence');
    grid on;
    
    subplot(2,1,2);
    plot(1:length(history.g1_norm), history.g1_norm, 'b-o', ...
         1:length(history.c1_norm), history.c1_norm, 'r-x');
    yline(1, 'k--', 'LineWidth', 1.5);
    xlabel('Iteration'); ylabel('Norm');
    legend('||g_1||', '||c_1||', '||·|| = 1');
    title('Fixed perturbation norms');
    grid on;
end

fprintf('\n========== Test Case 2: Larger system (6×6) ==========\n');
M2 = randn(6,6) + 0.3i*randn(6,6);
M2 = 0.3 * M2 / norm(M2);  % 缩放到合理范围

k_f2 = [2, 1, 0];  % 固定: 2实 + 1复
k_v2 = [1, 1, 1];  % 变化: 1实 + 1复 + 1满

[mu_s2, v_total2, history2] = skew_mu_lower_bound(M2, k_f2, k_v2, 50, 1e-6);

fprintf('\n=== Final Results ===\n');
fprintf('Skew μ lower bound: %.6f\n', mu_s2);
fprintf('Total scaling factor: %.6f\n', v_total2);
fprintf('Converged: %s\n', string(history2.converged));

end















































% % 第二版，有详细注释
% function [mu_s_lower, v_total, b, a, z, w] = skew_mu_lower_bound(M, k_f, k_v, max_iter, tol)
% % SKEW_MU_LOWER_BOUND - Compute lower bound for skewed structured singular value
% %
% % 实现论文 "Development of a skew μ lower bound" (Holland et al., 2005) 中的算法
% %
% % Inputs:
% %   M        - 系统矩阵 (n x n 复数矩阵)
% %   k_f      - 固定扰动的块结构 [k1_f, k2_f, k3_f]
% %              k1_f: 重复实标量块的大小
% %              k2_f: 重复复标量块的大小
% %              k3_f: 全复数块的大小
% %   k_v      - 变化扰动的块结构 [k1_v, k2_v, k3_v]
% %              (结构同上)
% %   max_iter - 最大迭代次数 (默认: 100)
% %   tol      - 收敛容差 (默认: 1e-6)
% %
% % Outputs:
% %   mu_s_lower - skew μ 的下界 (论文中的主要结果)
% %   v_total    - 总缩放因子 (累积的缩放)
% %   b, a, z, w - 收敛后的特征向量
% %
% % 理论基础：
% %   - Theorem 1: power算法的存在性和结构
% %   - Theorem 2: v是μ_s(M)的下界
% %   - Theorem 3: v的直接计算公式 (方程15和16)
% 
% if nargin < 4
%     max_iter = 100;
% end
% if nargin < 5
%     tol = 1e-6;
% end
% 
% % ========================================================================
% % 第1步: 提取块结构并初始化
% % ========================================================================
% % 从输入参数中提取各个块的大小
% k1_f = k_f(1); k2_f = k_f(2); k3_f = k_f(3);  % 固定扰动块
% k1_v = k_v(1); k2_v = k_v(2); k3_v = k_v(3);  % 变化扰动块
% 
% n_f = k1_f + k2_f + k3_f;  % 固定扰动的总维数
% n_v = k1_v + k2_v + k3_v;  % 变化扰动的总维数
% n = n_f + n_v;              % 系统总维数
% 
% % 验证矩阵维度是否匹配
% if size(M,1) ~= n || size(M,2) ~= n
%     error('Matrix M dimension does not match perturbation structure');
% end
% 
% % ========================================================================
% % 第2步: 初始化特征向量
% % ========================================================================
% % 论文备注: "一个可以用任意正实标量ζ乘以b和a，用任意正实标量ψ乘以z和w"
% % 因此可以施加额外约束 |a| = |w| = 1
% b = randn(n,1) + 1i*randn(n,1);  % 随机初始化复向量
% b = b / norm(b);                  % 归一化
% w = randn(n,1) + 1i*randn(n,1);
% w = w / norm(w);
% 
% % ========================================================================
% % 第3步: 初始化缩放
% % ========================================================================
% % v_total: 累积的缩放因子，用于跟踪总的μ_s值
% % M_scaled: 经过缩放的矩阵，用于避免数值问题
% v_total = 1;
% M_scaled = M;
% 
% % ========================================================================
% % 主循环: Power算法迭代
% % ========================================================================
% for iter = 1:max_iter
% 
%     % ====================================================================
%     % 阶段 I: 右特征向量迭代 (论文 4.2节)
%     % ====================================================================
%     % 求解: S_{k+1} * a_{k+1} = M * b_k (论文方程12的第一个等式)
% 
%     % 步骤 1: 计算 g = Mb
%     g = M_scaled * b;
% 
%     % 步骤 2: 将 g 分割为固定部分和变化部分
%     % 这对应于论文方程(18): g = [g1; g2]
%     g1 = g(1:n_f);      % 固定扰动对应的部分
%     g2 = g(n_f+1:end);  % 变化扰动对应的部分
% 
%     % ====================================================================
%     % 步骤 3: 检查固定扰动是否能使系统失稳 (论文 4.3节)
%     % ====================================================================
%     % 如果 ||g1||² > 1，说明固定扰动可以单独使系统失稳
%     % 这时需要缩放变化部分以避免数值问题
%     g1_norm = norm(g1);
% 
%     if g1_norm > 1
%         % 论文方程(23): 应用缩放矩阵来处理这种情况
%         % v_total *= 2: 记录缩放历史
%         v_total = 2 * v_total;
%         % 构造缩放矩阵: S = diag(I_f, sqrt(2)*I_v)
%         S_scale = blkdiag(eye(n_f), sqrt(2)*eye(n_v));
%         % M_scaled = S^(-1) * M * S^(-1)
%         M_scaled = inv(S_scale) * M * inv(S_scale);
%         continue;  % 重新开始迭代
%     end
% 
%     % ====================================================================
%     % 步骤 4: 计算 v_star (论文定理3，方程15)
%     % ====================================================================
%     % v = ||g2||₂ / sqrt(1 - ||g1||₂²)
%     g2_norm = norm(g2);
%     if g1_norm < 1
%         v_star = g2_norm / sqrt(1 - g1_norm^2);
%     else
%         v_star = inf;  % 边界情况
%     end
% 
%     % ====================================================================
%     % 步骤 5: 构造 S_star 并计算 a
%     % ====================================================================
%     % 论文方程(11): S = diag(I_f, v*I_v)
%     S_star = blkdiag(eye(n_f), v_star*eye(n_v));
%     % 从 Mb = Sa 得到 a = S^(-1)g
%     a = inv(S_star) * g;
%     a = a / norm(a);  % 归一化使 |a| = 1
% 
%     % ====================================================================
%     % 步骤 6: 向量分块 (论文方程10)
%     % ====================================================================
%     % 根据块结构将向量分割成对应的子块
%     % 分块顺序: [实标量块; 复标量块; 全复数块] 对于固定和变化部分
% 
%     % 固定扰动部分的分块
%     a1_f = a(1:k1_f);                           % 实标量块
%     a2_f = a(k1_f+1:k1_f+k2_f);                % 复标量块
%     a3_f = a(k1_f+k2_f+1:n_f);                 % 全复数块
% 
%     % 变化扰动部分的分块
%     a1_v = a(n_f+1:n_f+k1_v);                  % 实标量块
%     a2_v = a(n_f+k1_v+1:n_f+k1_v+k2_v);       % 复标量块
%     a3_v = a(n_f+k1_v+k2_v+1:end);            % 全复数块
% 
%     % 同样分块 w 向量
%     w1_f = w(1:k1_f);
%     w2_f = w(k1_f+1:k1_f+k2_f);
%     w3_f = w(k1_f+k2_f+1:n_f);
%     w1_v = w(n_f+1:n_f+k1_v);
%     w2_v = w(n_f+k1_v+1:n_f+k1_v+k2_v);
%     w3_v = w(n_f+k1_v+k2_v+1:end);
% 
%     % ====================================================================
%     % 步骤 7: 计算 q 参数 (论文方程14)
%     % ====================================================================
%     % q 是实标量块的相位参数，取值范围 [-1, 1]
%     % 条件: Re(a₁*w₁) ≤ 0 when q = -1
%     %       Re(a₁*w₁) ≥ 0 when q = +1
%     %       Re(a₁*w₁) = 0 when |q| < 1
%     q_star_f = compute_q(a1_f, w1_f, b(1:k1_f));
%     q_star_v = compute_q(a1_v, w1_v, b(n_f+1:n_f+k1_v));
% 
%     % ====================================================================
%     % 步骤 8: 更新 z 向量 (论文方程13，右侧部分)
%     % ====================================================================
%     % 根据不同的块类型应用不同的更新规则
%     z = zeros(n,1);
% 
%     % 固定扰动块的更新
%     if k1_f > 0
%         % 实标量块: z₁ = q*w₁
%         z(1:k1_f) = q_star_f * w1_f;
%     end
%     if k2_f > 0
%         % 复标量块: z₂ = (w₂*a₂)/(|w₂*a₂|) * w₂
%         z(k1_f+1:k1_f+k2_f) = update_complex_block(w2_f, a2_f);
%     end
%     if k3_f > 0
%         % 全复数块: z₃ = (|w₃|/|a₃|) * a₃
%         z(k1_f+k2_f+1:n_f) = update_full_block(w3_f, a3_f);
%     end
% 
%     % 变化扰动块的更新（规则相同）
%     if k1_v > 0
%         z(n_f+1:n_f+k1_v) = q_star_v * w1_v;
%     end
%     if k2_v > 0
%         z(n_f+k1_v+1:n_f+k1_v+k2_v) = update_complex_block(w2_v, a2_v);
%     end
%     if k3_v > 0
%         z(n_f+k1_v+k2_v+1:end) = update_full_block(w3_v, a3_v);
%     end
% 
%     % ====================================================================
%     % 阶段 II: 左特征向量迭代 (论文 4.2节)
%     % ====================================================================
%     % 求解: S_{k+1}* w_{k+1} = M* * z_{k+1} (论文方程12的第二个等式)
% 
%     % 步骤 1: 计算 c = M* * z
%     c = M_scaled' * z;
% 
%     % 步骤 2: 分割 c
%     c1 = c(1:n_f);
%     c2 = c(n_f+1:end);
% 
%     % ====================================================================
%     % 步骤 3: 计算 v_sharp (论文定理3，方程16)
%     % ====================================================================
%     % 这是从左特征向量得到的 v 值
%     c1_norm = norm(c1);
%     c2_norm = norm(c2);
% 
%     if c1_norm < 1
%         v_sharp = c2_norm / sqrt(1 - c1_norm^2);
%     else
%         v_sharp = inf;
%     end
% 
%     % ====================================================================
%     % 步骤 4: 构造 S_sharp 并计算新的 w
%     % ====================================================================
%     S_sharp = blkdiag(eye(n_f), v_sharp*eye(n_v));
%     w_new = inv(S_sharp) * c;
%     w_new = w_new / norm(w_new);  % 归一化使 |w| = 1
% 
%     % ====================================================================
%     % 步骤 5: 分块新的 w 向量
%     % ====================================================================
%     w1_f = w_new(1:k1_f);
%     w2_f = w_new(k1_f+1:k1_f+k2_f);
%     w3_f = w_new(k1_f+k2_f+1:n_f);
%     w1_v = w_new(n_f+1:n_f+k1_v);
%     w2_v = w_new(n_f+k1_v+1:n_f+k1_v+k2_v);
%     w3_v = w_new(n_f+k1_v+k2_v+1:end);
% 
%     % ====================================================================
%     % 步骤 6: 计算新的 q_sharp
%     % ====================================================================
%     q_sharp_f = compute_q(a1_f, w1_f, b(1:k1_f));
%     q_sharp_v = compute_q(a1_v, w1_v, b(n_f+1:n_f+k1_v));
% 
%     % ====================================================================
%     % 步骤 7: 更新 b 向量 (论文方程13，左侧部分)
%     % ====================================================================
%     b_new = zeros(n,1);
% 
%     % 固定扰动块的更新
%     if k1_f > 0
%         % 实标量块: b₁ = q*a₁
%         b_new(1:k1_f) = q_sharp_f * a1_f;
%     end
%     if k2_f > 0
%         % 复标量块: b₂ = (a₂*w₂)/(|a₂*w₂|) * a₂
%         b_new(k1_f+1:k1_f+k2_f) = update_complex_block(a2_f, w2_f);
%     end
%     if k3_f > 0
%         % 全复数块: b₃ = (|a₃|/|w₃|) * w₃
%         b_new(k1_f+k2_f+1:n_f) = update_full_block(w3_f, a3_f);
%     end
% 
%     % 变化扰动块的更新
%     if k1_v > 0
%         b_new(n_f+1:n_f+k1_v) = q_sharp_v * a1_v;
%     end
%     if k2_v > 0
%         b_new(n_f+k1_v+1:n_f+k1_v+k2_v) = update_complex_block(a2_v, w2_v);
%     end
%     if k3_v > 0
%         b_new(n_f+k1_v+k2_v+1:end) = update_full_block(w3_v, a3_v);
%     end
% 
%     % ====================================================================
%     % 步骤 8: 检查收敛性
%     % ====================================================================
%     % 论文备注: "如果 v_star ≠ v_sharp，则 max(v_star, v_sharp) 仍然给出下界"
%     v_avg = (v_star + v_sharp) / 2;
% 
%     if abs(v_star - v_sharp) < tol
%         % 收敛！v_star 和 v_sharp 已经足够接近
%         if v_avg >= 1
%             % 情况1: v ≥ 1，已找到最终答案
%             % μ_s = v * v_total (论文定理2)
%             mu_s_lower = v_avg * v_total;
%             b = b_new;
%             w = w_new;
%             return;
%         else
%             % 情况2: v < 1，需要继续缩放
%             % 论文 4.3节: 更新 v_total 并缩放矩阵
%             v_total = v_avg * v_total;
%             S_scale = blkdiag(eye(n_f), (1/sqrt(v_avg))*eye(n_v));
%             M_scaled = inv(S_scale) * M * inv(S_scale);
%         end
%     end
% 
%     % ====================================================================
%     % 步骤 9: 更新迭代变量
%     % ====================================================================
%     b = b_new;
%     w = w_new;
% 
%     % ====================================================================
%     % 步骤 10: 发散检测
%     % ====================================================================
%     % 如果 v_total 趋向无穷，说明固定扰动可以使系统失稳
%     % 论文 4.3节: "如果 ||g1||² ≥ 1 或 ||c1||² ≥ 1，则 μ_s(M) = ∞"
%     if v_total > 1e10
%         warning('v_total diverging - fixed perturbation can destabilize');
%         mu_s_lower = inf;
%         return;
%     end
% end
% 
% % ========================================================================
% % 达到最大迭代次数但未收敛
% % ========================================================================
% warning('Maximum iterations reached without convergence');
% % 即使未收敛，max(v_star, v_sharp) 仍是有效的下界
% mu_s_lower = max(v_star, v_sharp) * v_total;
%     if v_total > 1e10
%         warning('v_total diverging - fixed perturbation can destabilize');
%         mu_s_lower = inf;
%         return;
%     end
% % end
% 
% % If max iterations reached
% warning('Maximum iterations reached without convergence');
% mu_s_lower = max(v_star, v_sharp) * v_total;
% 
% end
% 
% %% Helper functions
% function q = compute_q(a1, w1, b1)
%     % Compute q parameter for real scalar blocks (Equation 14)
%     if isempty(a1) || isempty(w1)
%         q = 1;
%         return;
%     end
% 
%     alpha = abs(b1) / abs(a1) + real(a1' * w1);
% 
%     if abs(alpha) < 1
%         q = alpha / abs(alpha);
%     else
%         q = alpha;
%     end
% 
%     % Ensure q is in [-1, 1]
%     q = max(-1, min(1, real(q)));
% end
% 
% function z_out = update_complex_block(w, a)
%     % Update for repeated complex scalar blocks
%     if isempty(w) || isempty(a)
%         z_out = [];
%         return;
%     end
% 
%     inner_prod = w' * a;
%     if abs(inner_prod) > 0
%         z_out = (inner_prod / abs(inner_prod)) * w;
%     else
%         z_out = w;
%     end
% end
% 
% function z_out = update_full_block(w, a)
%     % Update for full complex blocks
%     if isempty(w) || isempty(a)
%         z_out = [];
%         return;
%     end
% 
%     z_out = (abs(w) ./ abs(a)) .* a;
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 第一版，没有详细注释
% % function [mu_s_lower, v_total, b, a, z, w] = skew_mu_lower_bound(M, k_f, k_v, max_iter, tol)
% % % SKEW_MU_LOWER_BOUND - Compute lower bound for skewed structured singular value
% % %
% % % Inputs:
% % %   M        - System matrix (n x n complex matrix)
% % %   k_f      - Structure of fixed perturbations [k1_f, k2_f, k3_f]
% % %              k1_f: size of repeated real scalar block
% % %              k2_f: size of repeated complex scalar block  
% % %              k3_f: size of full complex block
% % %   k_v      - Structure of varying perturbations [k1_v, k2_v, k3_v]
% % %   max_iter - Maximum number of iterations (default: 100)
% % %   tol      - Convergence tolerance (default: 1e-6)
% % %
% % % Outputs:
% % %   mu_s_lower - Lower bound for skew mu
% % %   v_total    - Total scaling factor
% % %   b, a, z, w - Converged eigenvectors
% % 
% % if nargin < 4
% %     max_iter = 100;
% % end
% % if nargin < 5
% %     tol = 1e-6;
% % end
% % 
% % % Extract block sizes
% % k1_f = k_f(1); k2_f = k_f(2); k3_f = k_f(3);
% % k1_v = k_v(1); k2_v = k_v(2); k3_v = k_v(3);
% % 
% % n_f = k1_f + k2_f + k3_f;  % Total size of fixed perturbations
% % n_v = k1_v + k2_v + k3_v;  % Total size of varying perturbations
% % n = n_f + n_v;
% % 
% % % Verify matrix dimensions
% % if size(M,1) ~= n || size(M,2) ~= n
% %     error('Matrix M dimension does not match perturbation structure');
% % end
% % 
% % % Initialize eigenvectors with random values
% % b = randn(n,1) + 1i*randn(n,1);
% % b = b / norm(b);
% % w = randn(n,1) + 1i*randn(n,1);
% % w = w / norm(w);
% % 
% % % Initialize scaling
% % v_total = 1;
% % M_scaled = M;
% % 
% % % Power algorithm iteration
% % for iter = 1:max_iter
% % 
% %     %% Right eigenvector iteration
% %     % Compute S_k+1 * a_k+1 = M * b_k
% %     g = M_scaled * b;
% % 
% %     % Partition g into fixed and varying parts
% %     g1 = g(1:n_f);
% %     g2 = g(n_f+1:end);
% % 
% %     % Check if fixed perturbation can destabilize
% %     g1_norm = norm(g1);
% % 
% %     if g1_norm > 1
% %         % Scale down varying part
% %         v_total = 2 * v_total;
% %         S_scale = blkdiag(eye(n_f), sqrt(2)*eye(n_v));
% %         M_scaled = inv(S_scale) * M * inv(S_scale);
% %         continue;
% %     end
% % 
% %     % Compute v_star
% %     g2_norm = norm(g2);
% %     if g1_norm < 1
% %         v_star = g2_norm / sqrt(1 - g1_norm^2);
% %     else
% %         v_star = inf;
% %     end
% % 
% %     % Form S_star
% %     S_star = blkdiag(eye(n_f), v_star*eye(n_v));
% %     a = inv(S_star) * g;
% %     a = a / norm(a);  % Normalize
% % 
% %     % Partition vectors according to block structure
% %     a1_f = a(1:k1_f);
% %     a2_f = a(k1_f+1:k1_f+k2_f);
% %     a3_f = a(k1_f+k2_f+1:n_f);
% %     a1_v = a(n_f+1:n_f+k1_v);
% %     a2_v = a(n_f+k1_v+1:n_f+k1_v+k2_v);
% %     a3_v = a(n_f+k1_v+k2_v+1:end);
% % 
% %     w1_f = w(1:k1_f);
% %     w2_f = w(k1_f+1:k1_f+k2_f);
% %     w3_f = w(k1_f+k2_f+1:n_f);
% %     w1_v = w(n_f+1:n_f+k1_v);
% %     w2_v = w(n_f+k1_v+1:n_f+k1_v+k2_v);
% %     w3_v = w(n_f+k1_v+k2_v+1:end);
% % 
% %     % Compute q_star for real blocks
% %     q_star_f = compute_q(a1_f, w1_f, b(1:k1_f));
% %     q_star_v = compute_q(a1_v, w1_v, b(n_f+1:n_f+k1_v));
% % 
% %     % Update z according to Equation (13)
% %     z = zeros(n,1);
% % 
% %     % Fixed blocks
% %     if k1_f > 0
% %         z(1:k1_f) = q_star_f * w1_f;
% %     end
% %     if k2_f > 0
% %         z(k1_f+1:k1_f+k2_f) = update_complex_block(w2_f, a2_f);
% %     end
% %     if k3_f > 0
% %         z(k1_f+k2_f+1:n_f) = update_full_block(w3_f, a3_f);
% %     end
% % 
% %     % Varying blocks
% %     if k1_v > 0
% %         z(n_f+1:n_f+k1_v) = q_star_v * w1_v;
% %     end
% %     if k2_v > 0
% %         z(n_f+k1_v+1:n_f+k1_v+k2_v) = update_complex_block(w2_v, a2_v);
% %     end
% %     if k3_v > 0
% %         z(n_f+k1_v+k2_v+1:end) = update_full_block(w3_v, a3_v);
% %     end
% % 
% %     %% Left eigenvector iteration
% %     c = M_scaled' * z;
% % 
% %     % Partition c
% %     c1 = c(1:n_f);
% %     c2 = c(n_f+1:end);
% % 
% %     % Compute v_sharp
% %     c1_norm = norm(c1);
% %     c2_norm = norm(c2);
% % 
% %     if c1_norm < 1
% %         v_sharp = c2_norm / sqrt(1 - c1_norm^2);
% %     else
% %         v_sharp = inf;
% %     end
% % 
% %     % Form S_sharp
% %     S_sharp = blkdiag(eye(n_f), v_sharp*eye(n_v));
% %     w_new = inv(S_sharp) * c;
% %     w_new = w_new / norm(w_new);
% % 
% %     % Partition new w
% %     w1_f = w_new(1:k1_f);
% %     w2_f = w_new(k1_f+1:k1_f+k2_f);
% %     w3_f = w_new(k1_f+k2_f+1:n_f);
% %     w1_v = w_new(n_f+1:n_f+k1_v);
% %     w2_v = w_new(n_f+k1_v+1:n_f+k1_v+k2_v);
% %     w3_v = w_new(n_f+k1_v+k2_v+1:end);
% % 
% %     % Compute q_sharp
% %     q_sharp_f = compute_q(a1_f, w1_f, b(1:k1_f));
% %     q_sharp_v = compute_q(a1_v, w1_v, b(n_f+1:n_f+k1_v));
% % 
% %     % Update b according to Equation (13)
% %     b_new = zeros(n,1);
% % 
% %     % Fixed blocks
% %     if k1_f > 0
% %         b_new(1:k1_f) = q_sharp_f * a1_f;
% %     end
% %     if k2_f > 0
% %         b_new(k1_f+1:k1_f+k2_f) = update_complex_block(a2_f, w2_f);
% %     end
% %     if k3_f > 0
% %         b_new(k1_f+k2_f+1:n_f) = update_full_block(w3_f, a3_f);
% %     end
% % 
% %     % Varying blocks
% %     if k1_v > 0
% %         b_new(n_f+1:n_f+k1_v) = q_sharp_v * a1_v;
% %     end
% %     if k2_v > 0
% %         b_new(n_f+k1_v+1:n_f+k1_v+k2_v) = update_complex_block(a2_v, w2_v);
% %     end
% %     if k3_v > 0
% %         b_new(n_f+k1_v+k2_v+1:end) = update_full_block(w3_v, a3_v);
% %     end
% % 
% %     % Check convergence
% %     v_avg = (v_star + v_sharp) / 2;
% % 
% %     if abs(v_star - v_sharp) < tol
% %         % Converged
% %         if v_avg >= 1
% %             mu_s_lower = v_avg * v_total;
% %             b = b_new;
% %             w = w_new;
% %             return;
% %         else
% %             % Scale and continue
% %             v_total = v_avg * v_total;
% %             S_scale = blkdiag(eye(n_f), (1/sqrt(v_avg))*eye(n_v));
% %             M_scaled = inv(S_scale) * M * inv(S_scale);
% %         end
% %     end
% % 
% %     % Update for next iteration
% %     b = b_new;
% %     w = w_new;
% % 
% %     % Check for divergence
% %     if v_total > 1e10
% %         warning('v_total diverging - fixed perturbation can destabilize');
% %         mu_s_lower = inf;
% %         return;
% %     end
% % end
% % 
% % % If max iterations reached
% % warning('Maximum iterations reached without convergence');
% % mu_s_lower = max(v_star, v_sharp) * v_total;
% % 
% % end
% % 
% % %% Helper functions
% % function q = compute_q(a1, w1, b1)
% %     % Compute q parameter for real scalar blocks (Equation 14)
% %     if isempty(a1) || isempty(w1)
% %         q = 1;
% %         return;
% %     end
% % 
% %     alpha = abs(b1) / abs(a1) + real(a1' * w1);
% % 
% %     if abs(alpha) < 1
% %         q = alpha / abs(alpha);
% %     else
% %         q = alpha;
% %     end
% % 
% %     % Ensure q is in [-1, 1]
% %     q = max(-1, min(1, real(q)));
% % end
% % 
% % function z_out = update_complex_block(w, a)
% %     % Update for repeated complex scalar blocks
% %     if isempty(w) || isempty(a)
% %         z_out = [];
% %         return;
% %     end
% % 
% %     inner_prod = w' * a;
% %     if abs(inner_prod) > 0
% %         z_out = (inner_prod / abs(inner_prod)) * w;
% %     else
% %         z_out = w;
% %     end
% % end
% % 
% % function z_out = update_full_block(w, a)
% %     % Update for full complex blocks
% %     if isempty(w) || isempty(a)
% %         z_out = [];
% %         return;
% %     end
% % 
% %     z_out = (abs(w) ./ abs(a)) .* a;
% % end
% % 
% 
% 
% 
% 
