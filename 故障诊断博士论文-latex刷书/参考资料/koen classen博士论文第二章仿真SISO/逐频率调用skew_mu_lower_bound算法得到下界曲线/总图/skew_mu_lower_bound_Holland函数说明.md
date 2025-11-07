下面把这段 `skew_mu_lower_bound_Holland.m` 逐行（按功能块）对照 Holland (2005) 的推导解释清楚；我会在关键位置标注论文里的对应式子/段落。

------

# 入口与参数检查

```matlab
function out = skew_mu_lower_bound_Holland(M, k_f, k_v, opts)
% Skew-μ lower bound (Holland 2005) with last-iteration vectors
```

实现 Holland 提出的 **skew-μ 下界幂迭代**。论文目标：把 skew-μ 的下界转化为一次在给定频率下矩阵 (M) 上的幂法迭代计算【】。

```matlab
if nargin < 4, opts = struct; end
if ~isfield(opts,'maxIter'), opts.maxIter = 200; end
if ~isfield(opts,'tol'),     opts.tol     = 1e-6; end
if ~isfield(opts,'verbose'), opts.verbose = 0;    end
```

算法停止准则与最大迭代次数（论文 §4.2 给出 “Detailed steps for a skew-μ lower bound power algorithm”，需要迭代到收敛；代码用 `tol` 控制收敛精度）【】。

```matlab
nf = sum(k_f); nv = sum(k_v); n = nf+nv;
[nr,nc] = size(M);
if nr~=n || nc~=n, error('M size %dx%d inconsistent with nf+nv=%d.',nr,nc,n); end
```

按照 skew 结构把不确定性通道分成 **固定块**（size 合成 `nf`）与 **可变块**（size 合成 `nv`），从而将 (M) 写成
\[
M=\begin{bmatrix}M_{11}&M_{12}\\M_{21}&M_{22}\end{bmatrix},
 \quad g=\begin{bmatrix}g_1\\ g_2\end{bmatrix}=Mb
\]
 这是论文 §3 的块结构前提，也是 §4.2 幂迭代中计算 \(||g_1||\)、\(||g_2||\) 的必要分区【】【】。

------

# 初始向量与缩放量

```matlab
rng(1);
b = randn(n,1)+1j*randn(n,1); b = b/norm(b);
w = randn(n,1)+1j*randn(n,1); w = w/norm(w);
```

初始化右向量 (b) 与左向量 (w)，并单位范数化。论文在幂迭代讨论中允许对向量做正实尺度缩放而不影响平衡点（用以便利收敛），即可以规定 (|a|=|w|=1)【】【】。

```matlab
v_total = 1; flag = 0; M_scaled = M;
```

`v_total` 是论文算法里累计缩放量的实现（“wrap into vtotal”）；当系统可被固定块单独致不稳时，`v_total → ∞`，即 (\mu_s(M)\to\infty) 的情形【】。

```matlab
hist_v=[]; hist_g1=[]; hist_c1=[];
a_last = []; w_last = [];
```

用于记录轨迹与最后迭代的向量（后续重构 (\Delta_{wc}) 会用到）。

------

# 主迭代：右步（b→a）

```matlab
for iter=1:opts.maxIter
    % ---- RIGHT: b -> a ----
    g = M_scaled*b; g1=g(1:nf); g2=g(nf+1:end);
```

对应论文 (17)(18)：令 \(g=Mb=[g_1^T\ g_2^T]^T\)，其中 \(g_1\) 与 \(g_2\) 分别是固定/可变子空间对应分量【】【】。

```matlab
    n_g1 = norm(g1); n_g2 = norm(g2);
    if n_g1 >= 1-1e-9, flag=1; break; end
```

当 \(||g_1||_2 \ge 1\) 时，表示**仅凭固定块**就能使 (\det(I-M\Delta)=0)（系统已越界），这是论文定理 3 后的讨论结论之一：此时 (\mu_s(M)=1)（或趋于无穷大下界情况视定义归一），算法按“固定块致不稳”退出【】【】。

```matlab
    v_star = (n_g2^2)/sqrt(max(1e-16,1-n_g1^2));
```

**核心公式（右步）**：论文定理 3 的式 (15)
\[
v_*=\frac{||g_2||_2^2}{\sqrt{1-||g_1||_2^2}}
\]
 这是把 (a=S^{-1}g)（见下）与 (|a|=1) 条件结合推得的显式下界计算式【】【】。

```matlab
    S_star = blkdiag(eye(nf), v_star*eye(nv));
    a = S_star \ g; a = a/norm(a);
```

对应论文 (17)：(S^{-1}=\mathrm{diag}(I_f,,\tfrac{1}{v}I_v))；于是
\[
a=S^{-1}g=\begin{bmatrix}g_1\\ \tfrac{1}{v}g_2\end{bmatrix},
 \quad |a|=1
\]
 代码用\(S=\mathrm{diag}(I_f,v I_v)\) 写法，求解 `a = S\g` 等价【】【】。

------

# 主迭代：左步（z→w）

```matlab
    % ---- LEFT: z -> w ----
    z = local_update_z(a,w,k_f,k_v,nf,nv);
```

`z` 的构造满足论文 (13)(14) 的 **块对齐（alignment / scaling）规则**：对不同块类型（实重复/复重复/全复）采用不同的相位/符号/幅值对齐（Q、D 集的选择）以保证块约束与幂法的对偶关系，这些规则在定理叙述与式 (21) 的关系中起到“把 a、w 与 Q、D 对齐”的作用【】【】【】。

```matlab
    c = M_scaled' * z; c1=c(1:nf); c2=c(nf+1:end);
    n_c1 = norm(c1); n_c2 = norm(c2);
    if n_c1 >= 1-1e-9, flag=1; break; end
```

左步完全对称：令 \(c=M^H z=[c_1^T\ c_2^T]^T\)。若 \(||c_1||\ge 1\) 亦意味着固定块即可致不稳（同上论断）【】。

```matlab
    v_sharp = (n_c2^2)/sqrt(max(1e-16,1-n_c1^2));
    S_sharp = blkdiag(eye(nf), v_sharp*eye(nv));
    w_new = S_sharp \ c; w_new = w_new/norm(w_new);
```

**核心公式（左步）**：论文定理 3 的式 (16)
\[
v_c=\frac{|c_2|_2^2}{\sqrt{1-|c_1|_2^2}}
\]
 同理由 (w=S^{-1}c)、(|w|=1) 推得【】。

------

# 合并两侧下界并缩放

```matlab
    v_local = max(v_star, v_sharp);
    hist_v(end+1)=v_local; hist_g1(end+1)=n_g1; hist_c1(end+1)=n_c1;

    % save for reconstruction if we stop next
    a_last = a; w_last = w_new;
```

论文对齐性质表明：若 (v_*=v_#) 则已达到方程 (12) 的分解；若二者不等，**取 (\max(v_*,v_#)) 仍是有效下界**（Young–Doyle 混合-μ 下界中的标准结论，论文式 (21) 后给出结论）【】。

```matlab
    if v_local < 1+opts.tol
        break;
    end
```

若当前迭代的 (v<1)（容差内），说明在当前缩放下已“入界”，可停止；累计缩放量 `v_total` 即可作为下界。对应 Holland 的实现细节与 4.3 的“避免数值问题/停止条件”描述【】。

```matlab
    v_total = v_total * v_local;
```

把这一轮的 (v) 乘到累计量（论文把“不断把 2 或 (v) 包装进 (v_{total})”的口径写得很直白），最后输出的就是 **skew-μ 下界**【】。

```matlab
    S_L = blkdiag(eye(nf), 1/sqrt(v_local)*eye(nv));
    S_R = blkdiag(eye(nf), sqrt(v_local)*eye(nv));
    M_scaled = S_L * M_scaled * S_R;
```

这就是论文 (23)(24) 的**配平/缩放步骤**。当 (|g_1|>1) 时按 (23) 用 (\sqrt{2}) 缩放；当 (v\ge 1) 继续迭代时按 (24) 用 (\sqrt{v}) 缩放，使数值保持良性并维持分块结构的相对尺度【】【】。
 本实现用统一的 (v_{local}) 来做 (24) 形态的缩放（当检测到 (|g_1|) 越界时在外层亦会触发“固定块致不稳”的退出）。

```matlab
    % update right vector using block rules
    b = local_update_b(a, w_new, k_f, k_v, nf, nv); b = b/norm(b);
    w = w_new;
```

更新右向量 (b)。这一步与 `local_update_z` 一样，是论文 (13)(14) 的**块对齐规则**（把 a、w 的相位/符号对准，使下一次 (g=Mb) 能继续朝满足 (21) 的方向靠拢）【】。

```matlab
    if v_total>1e10, flag=1; break; end
end
```

按照论文 §4.3 的“三种可能结果”，若 `v_total` 爆长（趋向 ∞），就意味着“可由固定块致不稳”的情形，直接退出【】。

------

# 输出

```matlab
out.mu_lb   = v_total;
out.v_total = v_total;
out.iter    = iter;
out.flag    = flag;
out.history = struct('v',hist_v,'g1',hist_g1,'c1',hist_c1);
out.last    = struct('a',a_last,'w',w_last);   % <—— for Δ_wc reconstruction
```

* `mu_lb=v_total`：最后累计缩放量就是 **skew-μ 下界**（论文结论）【】。
* `last.a/w` 用来按式 (21) 的对齐关系恢复 (\Delta_{wc}) 的**相位/符号方向**（实块取 (\operatorname{sign}(\Re(a_1\bar w_1)))，复块取 (\angle(a_2^H w_2))），这对应“当达到平衡点时存在 (Q,D) 使得 (Mb=!*Sa,, M^H z=#S w)”的关系【】。

------

# 两个 helper：块对齐 = 论文 (13)(14) 的实现

```matlab
function z = local_update_z(a,w,kf,kv,nf,nv)
```

对 **固定/可变**、**实/复/全复** 三类块分开更新（注释里 `kf(1:3)`、`kv(1:3)`）。规则要点：

* **实重复块**：

  ```matlab
  z = sign(real(a.*conj(w))).*w;
  ```

  对齐符号，使 (\Re(a^\ast \circ w)) 非负；对应论文里实块的 alignment 条件（Q 集合）【】【】。

* **复重复块**：

  ```matlab
  inner = a'*w; z = (inner/|inner|) * w;
  ```

  用内积相位对齐，使相对相位一致（D/Q 集的复相位对齐）。

* **全复块**：

  ```matlab
  z = (|w|./|a|).*a;
  ```

  幅值配平，使 (|a|=|w|)（论文 remark 允许把 a、w 缩放到单位范数；全复块用幅值配平保持块约束）【】。

```matlab
function b = local_update_b(a,w,kf,kv,nf,nv)
```

与上面对称：把 **b 的各块**沿着 a、w 的对齐方向更新，以便下一轮右步计算 (g=Mb) 时满足式 (21) 那组“Q、D 与向量”的一致性关系（论文式 (21)）【】。

------

# 小结（论文—代码映射）

* **定理 3 / 式 (15)(16)**：`v_star` 与 `v_sharp` 的两种计算（右/左步），取 (\max) 仍为有效下界【】【】。
* **式 (17)**：`a = S^{-1} g`，代码用 `S\g` 实现【】。
* **(13)(14) 的块约束/对齐**：`local_update_z` / `local_update_b`；保证实/复/全复块的相位和幅值满足集合 (Q_{K_f,K_v})、(D_{K_f,K_v}) 的构造要求【】。
* **式 (23)(24)**：`M_scaled = S_L * M_scaled * S_R` 的缩放，防止数值发散并实现 “wrap into (v_{total})”【】【】。
* **三种结局**：收敛到 (v\to 1)（找到 `v_total`）、固定块越界（`v_total→∞`）、一直低于 1（`v_total→0`）【】。

这样，这个函数就完整地把 Holland 的 skew-μ 下界幂迭代按公式逐步落到了可运行的 MATLAB 实现上。