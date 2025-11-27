## 对什么样的系统进行辨识

线性时不变（LTI）、因果（Causal）、连续时间的 MIMO 系统，并且是加性系统

### **输入与输出**：

系统有 $p$ 个输入（$u(t) \in \mathbb{R}^p$）和 $q$ 个输出（$y(t) \in \mathbb{R}^q$）。

### **加性系统**：

系统被建模为 $K$ 个子系统的并联叠加。直接对应于物理世界中的“模态叠加原理”。

* [ ] 总系统表达式

### **子系统参数化**：

每一个子系统 $G^*_i(p)$ 都被表示为一个分式：
$$
G^*_i(p, \boldsymbol{\theta}) = \frac{\textbf{B}^*_i(p)}{A^*_i(p)}
$$

* **分母 $A^*_i(p)$**：这是一个**标量（Scalar）**多项式，代表该子系统的极点（Pole）。
  注意它是“反首一”（anti-monic）的（常数项为1），并且假设是稳定的（根在左半平面）。
  
  * **物理含义**：每一个模态的**极点**（自然频率和阻尼比）对于所有输入输出通道都是**共有**的。比如，你敲击一个桌子的不同位置，桌子发出的声音频率（固有频率）是一样的。因此，分母必须是标量多项式。
  
  * **Anti-monic 约束**：常数项设为1。这是连续时间辨识中的标准做法，为了保证参数的唯一性并避免平凡解。
  
  \[
  A^*_i(p) = 1 + a^*_{i,1}p + \dots + a^*_{i,n_i}p^{n_i}
  \]
  
* **分子 $B^*_i(p)$**：这是一个**矩阵（Matrix）**多项式，维度为 $q \times p$，代表该子系统的零点和增益结构。

  * **物理含义**：虽然频率一样，但在不同位置激励和测量时，**振型（Mode Shapes）**的幅值和相位是不同的。这个矩阵捕捉了空间上的输入输出增益关系（即留数矩阵）。它是 $q \times p$ 维的。


\[
  \textbf{B}^*_i(p) = \textbf{B}^*_{i,0} + \textbf{B}^*_{i,1}p + \dots + \textbf{B}^*_{i,m_i}p^{m_i}
\]

* **互质性 (Coprimeness)**：假设 $A^*_i$ 和 $\textbf{B}^*_i$ 是互质的，且不同的子系统之间也是互质的，以保证模型表示的唯一性。

**物理意义**：这种结构允许我们将一个复杂的机械系统拆解为一个个独立的“模态”。例如，一个子系统可能代表“刚体模态”（二重积分），另一个代表“一阶弯曲模态”（二阶振荡）。

**每个子系统的参数**：每个子系统的参数被收集在向量 $\boldsymbol{\theta}_i^*$ 中：
$$
\boldsymbol{\theta}_i^* = [a_{i,1}^* \dots a_{i,n_i}^* \quad b_{i,0}^* \dots b_{i,m_i}^*]^\top \tag{}
$$
参数合并为：
$$
\beta^* := \begin{bmatrix} \boldsymbol{\theta}_1^{*\top} & \boldsymbol{\theta}_2^{*\top} & \dots & \boldsymbol{\theta}_K^{*\top} \end{bmatrix}^\top
$$



## 连续无噪声输出（理想输出）与采样有噪声输出（实际输出）

### **理想输出**：连续且无噪声

连续且无噪声

无噪声的系统输出 $x(t)$ 是 $K$ 个子系统输出的总和：
$$
 x(t) = \sum_{i=1}^{K} G_i^*(p)u(t) \tag{9.1}
$$
**算子 $p$**：这里用 $p$ 表示微分算子 $\frac{d}{dt}$，即 $p u(t)$ 代表对 $u(t)$ 求导。在频域分析时，你可以把它看作拉普拉斯算子 $s$。



### **实际输出**：离散的采样输出，且含有噪声

离散的采样输出，且含有噪声

* [ ] 论文在这里处理了一个连续时间辨识中的经典难题：**连续时间白噪声在物理上是不存在的（功率无穷大）**。


公式 (9.4) 描述了我们实际测量到的输出 $y(t_k)$：
$$
 \mathbf{y}(t_k) = \mathbf{x}(t_k) + \mathbf{v}(t_k) 
$$

*   $t_k$：采样时刻。
*   $\mathbf{x}(t_k)$：真实的连续系统在采样时刻的输出。
*   $\mathbf{v}(t_k)$：离散的测量噪声。这里假设它是**平稳随机过程**（统计特性不随时间变化）。
*   $\mathbf{y}(t_k)$：我们在 $t_k$ 时刻实际采样得到的带噪测量值。



### **采样信号重构**：计算机VS实际

这是一个很多学生容易忽略的细节。计算机只能存储离散的 $u(t_k)$，但连续系统 $G(p)$ 需要连续的输入 $u(t)$ 才能演化。
文中明确假设：**输入信号 $u(t)$ 在采样点之间保持不变**，即 **零阶保持器 (ZOH)** 行为。

系统是连续的（微分方程），数据是离散的（采样点），需要处理采样带来的影响（ZOH假设）



## （对整体论文理解无影响，可略过）对带噪数据直接差分？引入状态变量滤波器（数据预处理）

在连续时间系统辨识（Continuous-Time System Identification）中，SVF 是连接“理论上的微分方程”与“现实中的采样数据”的**桥梁**。没有它，直接辨识微分方程系数几乎是不可能的。

### 差分

我们的目标模型是微分方程（公式 2.1）：
$$ \mathbf{x}(t) = \sum \mathbf{G}_i^*(p)\mathbf{u}(t) $$
这其中包含算子 $p^k$，意味着我们需要求 $k$ 阶导数。

**困境**：

1.  我们只有离散采样数据 $\mathbf{u}(t_k)$ 和 $\mathbf{y}(t_k)$。
2.  数据中必然包含测量噪声 $\mathbf{v}(t_k)$。

**灾难性的后果**：
如果你尝试直接用差分来近似微分（例如 $\dot{y} \approx \frac{y_k - y_{k-1}}{T}$）：

*   **有用信号**：低频为主，差分后幅值变化符合物理规律。
*   **无用噪声**：通常是宽带（白）噪声，高频能量很大。
*   **数学事实**：微分操作在频域等价于乘以 $j\omega$。频率 $\omega$ 越高，增益越大。
*   **结论**：直接对带噪数据求导，会**极度放大高频噪声**，导致信噪比（SNR）急剧下降，数据完全不可用。



我们想要辨识的物理模型通常是微分方程的形式，例如：
$$ a_2 \frac{d^2y}{dt^2} + a_1 \frac{dy}{dt} + y = b_0 u $$
这就要求我们需要知道输入 $u(t)$ 和输出 $y(t)$ 的**时间导数**（$\dot{u}, \ddot{u}, \dot{y}, \ddot{y}$ 等）。

然而，在现实中，我们只有离散的采样数据 $u(t_k), y(t_k)$，且不可避免地含有**噪声**。

*   **直接差分的灾难**：如果你尝试用差分（例如 $\frac{y_k - y_{k-1}}{\Delta t}$）来近似导数，高频噪声会被极度放大（因为微分在频域相当于乘以 $j\omega$，频率越高增益越大）。这会导致信噪比极低，辨识出的参数完全错误。

**SVF 的解决思路**：
与其直接求导，不如先对信号进行**低通滤波**以压制高频噪声，然后再（或者同时）计算导数。数学上，我们在原微分方程两边同时除以一个多项式 $F(p)$，方程的等式关系依然成立，但每一项都变成了“滤波后的导数”。



### 如果对带噪数据直接差分会怎么样？高频噪声超级放大

假设 $y(t) = \sin(t)$，有一个频率很高但幅值很小的噪声 $n(t) = 0.01 \sin(1000t)$。

$$ y_{meas}(t) = \sin(t) + 0.01 \sin(1000t) $$

如果我们尝试计算二阶导数 $\ddot{y}$（加速度）：

1.  **真实信号部分**：$\frac{d^2}{dt^2}\sin(t) = -\sin(t)$ （幅值是 1）。
2.  **噪声部分**：$\frac{d^2}{dt^2} 0.01 \sin(1000t) = -0.01 \times 1000^2 \times \sin(1000t) = -10000 \sin(1000t)$。

**结果**：噪声被放大了 10000 倍！你的真实信号会被淹没在噪声的海洋里，导致辨识出的 $a_1, a_0$ 完全错误。

### 引入SVF：本质是带通/低通滤波器

为了解决这个问题，我们设计一个**二阶低通滤波器** $F(p)$。
假设我们要辨识的频率范围是 0-5Hz，我们选一个截止频率稍高一点的滤波器。
比如选 $F(p) = (p + \lambda)^2$，设 $\lambda = 5$。

$$ F(p) = p^2 + 10p + 25 $$
对应论文公式 (2.13)，我们的滤波器系数是：

*   $c_1 = 10$
*   $c_2 = 25$

我们的目标不再是求 $\ddot{y}$，而是求**滤波后的导数**：

*   $y_f^{(0)}(t) = \frac{1}{F(p)} y(t)$ （滤波后的位移）
*   $y_f^{(1)}(t) = \frac{p}{F(p)} y(t)$ （滤波后的速度）
*   $y_f^{(2)}(t) = \frac{p^2}{F(p)} y(t)$ （滤波后的加速度）



### 公式 (2.9) 和 (2.10)：滤波后的导数重建

$$ \mathbf{u}^{(l)}(t_k) = \frac{p^l}{F(p)} \mathbf{u}(t_k), \quad l = 0, \dots, m_i $$
$$ \mathbf{y}^{(i)}(t_k) = \frac{p^i}{F(p)} \mathbf{y}(t_k), \quad i = 0, \dots, n_i $$

*   **$p$ (Heaviside Operator)**：这里的 $p$ 代表微分算子 $\frac{d}{dt}$。因此 $p^l$ 代表 $l$ 阶导数。
*   **$\frac{1}{F(p)}$**：这是一个**低通滤波器**。
*   **$\frac{p^l}{F(p)}$**：这是一个组合操作。它表示“对信号进行 $l$ 阶微分”**并且**“经过滤波器 $1/F(p)$”。
    *   在频域上看，$\frac{(j\omega)^l}{F(j\omega)}$。
    *   当频率较低（系统带宽内）时，它是微分器。
    *   当频率很高（噪声区域）时，由于 $F(p)$ 的阶数够高，幅值会衰减，从而抑制噪声。
*   **$\mathbf{u}^{(l)}(t_k)$**：这是计算出来的、在 $t_k$ 时刻的“滤波后输入信号的 $l$ 阶导数”。注意，这**不是**原始信号的导数，而是经过预处理后的信号的导数。由于方程两边都做了同样的处理，参数估计依然是无偏的（在理想情况下）。

### 公式 (2.11)：滤波器的设计

$$ F(p) = p^n + c_1 p^{n-1} + \dots + c_{n-1}p + c_n $$

*   **用户自定义**：这个多项式 $F(p)$ 是由作为工程师的你来设计的，而不是辨识出来的。
*   **Monic（首一多项式）**：最高次项 $p^n$ 的系数固定为 1。
*   **稳定性**：$F(p)$ 的根必须都在左半平面（Stable），否则滤波过程会产生发散的数值。
*   **带宽选择**：$F(p)$ 的截止频率通常选择在比系统期望的带宽稍高一点的地方。如果太低，会滤掉有用的系统动态；如果太高，抑制噪声的效果就不好。

### 关键约束：相对阶 (Relative Degree)

文中提到：“The relative degree of the filter $1/F(p)$ must have a relative degree greater or equal to $n_i$”。

*   **数学含义**：$F(p)$ 的阶数 $n$ 必须大于或等于系统模型中最高的微分阶数 $n_i$。
*   **直观理解**：
    *   如果要计算 $n_i$ 阶导数，即算子是 $p^{n_i}$。
    *   如果 $F(p)$ 的阶数小于 $n_i$，那么 $\frac{p^{n_i}}{F(p)}$ 在高频下仍然是一个增益随频率增加而增加的滤波器（非真分式），这会放大噪声。
    *   只有当 $\text{deg}(F) \ge n_i$ 时，$\frac{p^{n_i}}{F(p)}$ 才是**真分式（Proper）**，在高频下的增益才是常数或衰减的，才能物理实现并抑制噪声。

### 3. Remark 1 的深度解读：采样与重构的哲学

这一段 Remark 非常精彩，它触及了数字信号处理中“连续”与“离散”的边界问题。

**混合记号的含义**：
文中使用了混合记号 $\frac{1}{F(p)} z(t_k)$。一个连续时间的滤波器怎么能作用在一个离散的采样序列 $z(t_k)$ 上呢？

**两种解释的对比**：

1.  **解释 A（实际操作）：先插值，再滤波**
    *   **符号**：$(1/F(p))z(t_k)$
    *   **流程**：
        1.  拿到采样点 $z(t_1), z(t_2), \dots$。
        2.  假设采样点之间的行为（**Intersample Behavior**）。比如假设两个点之间是直线连接（FOH）或者保持不变（ZOH）。
        3.  基于这个假设，重构出连续信号 $\hat{z}(t)$。
        4.  将 $\hat{z}(t)$ 送入连续滤波器 $1/F(p)$ 进行模拟。
        5.  对滤波后的输出在 $t_k$ 时刻进行重采样。
    *   **论文观点**：这是实际工程中必须采用的方法，因为我们不知道采样点之间发生了什么。通常对于控制系统的输入信号，我们确切知道它是 ZOH（零阶保持，由DAC产生）；对于输出信号，FOH（一阶保持）通常是更好的近似。

2.  **解释 B（理论理想）：先滤波，再采样**
    *   **符号**：$\{1/F(p)z(t)\}_{t=t_k}$
    *   **流程**：假设存在一个真实的连续信号 $z(t)$，它先经过物理滤波器 $1/F(p)$，然后我们在 $t_k$ 时刻去采样这个滤波后的结果。
    *   **论文观点**：这在数学推导上很方便，但在实际辨识中是**不现实的**。因为我们没法在计算机里拿到真实的 $z(t)$，我们拿到手的时候它已经是离散的 $z(t_k)$ 了。因此，我们必须依赖解释 A 中的假设。

### 教授总结

这一节不仅是定义公式，更是在解决**连续时间系统辨识的三大挑战**：

1.  **导数不可测** $\rightarrow$ 使用 SVF 进行**带通滤波**（低频近似微分，高频截止）。
2.  **因果性与实现** $\rightarrow$ 确保滤波器是**真分式**（Proper），且通过状态空间方程在计算机中数值积分实现。
3.  **离散与连续的鸿沟** $\rightarrow$ 明确指出了必须对**采样点之间的行为**做出假设（如 ZOH 或 FOH），并通过数值积分来精确计算滤波后的值，而不是简单地做离散变换。

这为后续构建线性回归方程 $\mathbf{y}_f = \mathbf{\Phi}_f \boldsymbol{\theta}$ 打下了坚实的基础，其中 $\mathbf{y}_f$ 和 $\mathbf{\Phi}_f$ 都是由上述 SVF 方法处理得到的“干净”数据。





## 对输入的假设

**注记 9.2.** 在本章中，假设输入（对于开环情况）和参考信号（对于闭环情况）是拟平稳的。因此，考虑拟平稳信号的标准期望定义[161, p. 34]
$$
\mathbb{E}\{g(t_k)\} := \lim_{N \to \infty} \frac{1}{N} \sum_{k=1}^N \mathbb{E}\{g(t_k)\}.
$$

## 辨识准则的基础（目标）：最小化OE

辨识效果好=得到了好的模型。好的模型的标准就是**输出误差 (Output Error, OE) 最小化**。



* 首先定义**实际观测输出**与**参数化模型输出**（待估计）之间的残差 $\varepsilon(t_k, \beta)$：

$$
\boldsymbol{\varepsilon}(t_{k},\boldsymbol{\beta})=\mathbf{y}(t_{k})-\sum_{i=1}^ {K}\mathbf{P}_{i}(p,\boldsymbol{\beta})\mathbf{u}(t_{k}), \tag{2.15}
$$

 其中 $P_i$ 是第 $i$ 个子系统的模型。

即：在开环系统中，最自然的想法是让**模型预测的输出**尽可能接近**实际测量的输出**。

定义代价函数，方便得到辨识准则

### 代价函数

\[
V(\boldsymbol{\beta}) = \frac{1}{N} \sum_{k=1}^{N} \| \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta}) \|^2_{\boldsymbol{\Lambda}_0^{-1}} \tag{2.19}
\]

其中：

* **$\boldsymbol{\beta}$**：模型参数（包括所有模态的分子、分母系数）。

* **$\boldsymbol{\varepsilon}$**：**预测误差（残差）**。即：真实测量值 $y$ - 模型计算值 $\hat{y}$。

* 写在右下角的下标 $\boldsymbol{\Lambda}_0^{-1}$，代表这不仅仅是一个普通的平方和，而是一个**加权范数（Weighted Norm）**

  * 在数学上，符号 $\| \mathbf{x} \|_{\mathbf{W}}^2$ 被定义为：
    $$ \| \mathbf{x} \|_{\mathbf{W}}^2 = \mathbf{x}^\top \mathbf{W} \mathbf{x} $$

    具体到你提到的公式 (2.19)，它的完整展开形式是：
    $$ \| \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta}) \|^2_{\boldsymbol{\Lambda}_0^{-1}} = \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta})^\top \boldsymbol{\Lambda}_0^{-1} \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta}) $$

* **$\boldsymbol{\Lambda}_0^{-1}$**：这是一个权重矩阵（噪声协方差的逆）。

  * **直观理解**：如果某个传感器的噪声很大（方差大），我们在计算总误差时就应该少听它的（权重小）；反之亦然。这叫做加权最小二乘，目的是让估计结果符合**最大似然（Maximum Likelihood）**原则。

  * $\boldsymbol{\Lambda}_{0}\in\mathbb{R}^{q\times q}$ 表示真实噪声协方差矩阵
    \[
    \boldsymbol{\Lambda}_{0}=\mathbb{E}\{\mathbf{v}(t_{k})\mathbf{v}^{\top}(t_{k}) \}. \tag{2.17}
    \]

我们想要找到使代价函数最小的\(\boldsymbol{\beta}\)，即
$$
\hat{\boldsymbol{\beta}} = \arg\min_{\boldsymbol{\beta}} V(\boldsymbol{\beta}), \tag{2.18}
$$
为了使得其最小，我们将代价函数的优化器作为一阶最优条件的解得到
为了找到让 $V(\boldsymbol{\beta})$ 最小的参数，我们需要对它求导并令其为 0。
$$
\mathbf{0} = \frac{\partial}{\partial \boldsymbol{\beta}} V(\boldsymbol{\beta}) = \frac{1}{N} \sum_{k=1}^{N} \left( \frac{\partial \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta})}{\partial \boldsymbol{\beta}} \right)^\top \boldsymbol{\Lambda}_0^{-1} \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta}). \tag{2.20}
$$
请仔细看这个式子的结构：它是一个**加权内积**。
它告诉我们：当**残差 $\boldsymbol{\varepsilon}$** 与 **残差的梯度（灵敏度函数）** 正交（内积为0）时，我们达到了最优点。



# 非线性问题详解

通常我们认为：**开环系统 $\rightarrow$ 输入与噪声无关 $\rightarrow$ 最小二乘法（LS）就够了**。那为什么这篇论文（以及著名的 SRIVC 算法）在开环下还要大费周章地使用工具变量（IV）呢？

答案并非在于**物理系统**本身（物理上输入确实是干净的），而在于为了求解问题所采用的**数学变换**。

认知：非线性问题没有解析解只有数值解？线性问题直接解析解。

简单来说：**为了把难解的非线性问题转化为好解的线性问题，我们引入了“含噪回归向量”，这迫使我们必须使用 IV 来修正偏差。**

## 输出误差对参数是非线性的，无法写成 $\mathbf{Y} = \mathbf{\Phi} \boldsymbol{\theta}$形式，不能用解析解求解\(\boldsymbol{\theta}\)

我们的终极目标是最小化输出误差（Output Error, OE）。
模型结构是：
$$ y(t_k) = \frac{B(p)}{A(p)} u(t_k) + v(t_k) $$

我们要最小化的残差是：
$$ \varepsilon(t_k) = y(t_k) - \frac{B(p)}{A(p)} u(t_k) $$

请注意分母上的 $A(p)$。参数（即 $A(p)$ 的系数）出现在分母里，这意味着**残差 $\varepsilon$ 对参数是非线性的**。

其中：

*   $A(p) = 1 + a_1 p + a_2 p^2 + \dots + a_n p^n$ 
*   $B(p) = b_0 + b_1 p + \dots + b_m p^m$

**痛点**：参数 $a_i$ 在分母 $A(p)$ 里。如果你对 $a_i$ 求导，你会得到非线性的结果，这就没法用线性回归 $Ax=b$ 来解了。

*   你不能写成 $\mathbf{Y} = \mathbf{\Phi} \boldsymbol{\theta}$ 的形式。
*   你无法直接用解析解 $\boldsymbol{\theta} = (\mathbf{\Phi}^\top\mathbf{\Phi})^{-1}\mathbf{\Phi}^\top\mathbf{Y}$ 一步算出答案。

通常解决非线性优化需要用梯度下降法（如 Gauss-Newton），但这些方法对初始值非常敏感，容易陷入局部最优。

因此，我们是想把优化问题写成类似线性的形式的

## 解决方法：通过变换写成伪线性形式（How）

## 如何把残差写成伪线性形式

### .1 提出分母 $A(p)$ 

SRIVC 算法做了一个极其聪明的操作：**把分母 $A(p)$ 提出来**。

我们要最小化的残差是：
$$ \varepsilon(t_k) = y(t_k) - \frac{B(p)}{A(p)} u(t_k) $$

我们把残差公式通分一下（或者说提取公因式 $\frac{1}{A(p)}$）：

$$ \begin{aligned} \varepsilon(t) &= \frac{1}{A(p)} \left[ A(p) y(t) - B(p) u(t) \right] \\ &= \frac{1}{A(p)} \Big[ \underbrace{(1 + a_1 p + \dots + a_n p^n)}_{\text{展开 } A(p)} y(t) - \underbrace{(b_0 + b_1 p + \dots + b_m p^m)}_{\text{展开 } B(p)} u(t) \Big] \end{aligned} $$

注意，这里 $A(p)$ 的第一项是 **1**。我们把这一项单独拆出来：

$$ \varepsilon(t) = \frac{1}{A(p)} \left[ y(t) + (a_1 p + \dots) y(t) - (b_0 + \dots) u(t) \right] $$

现在，利用线性系统的性质，把外面的滤波器 $\frac{1}{A(p)}$ 乘进去：

$$ \varepsilon(t) = \underbrace{\frac{1}{A(p)} y(t)}_{\text{滤波后输出 } y_f} + \underbrace{\frac{a_1 p}{A(p)} y(t) + \dots}_{\text{含参数 } a_i \text{ 的项}} - \underbrace{\frac{b_0}{A(p)} u(t) - \dots}_{\text{含参数 } b_i \text{ 的项}} $$

### .2 整理成线性回归的形式

我们把不含未知参数的项 $y_f$ 留在左边（或者作为基准），把含有未知参数 $a_i, b_i$ 的项整理成向量内积的形式。

\[
 \varepsilon(t) = y_f(t) - \left[ \underbrace{- \frac{p}{A(p)}y(t)}_{\text{对应 } a_1} \cdot a_1 - \dots - \underbrace{\frac{p^n}{A(p)}y(t)}_{\text{对应 } a_n} \cdot a_n + \underbrace{\frac{1}{A(p)}u(t)}_{\text{对应 } b_0} \cdot b_0 + \dots \right] 
\]

我们可以定义一个**参数向量** $\boldsymbol{\theta}$：
$$ \boldsymbol{\theta} = [a_1, a_2, \dots, a_n, b_0, b_1, \dots, b_m]^\top $$

再定义一个**回归向量** $\boldsymbol{\Phi}(t)$：
$$ \boldsymbol{\Phi}(t) = \left[ \underbrace{\frac{-p}{A(p)}y(t), \frac{-p^2}{A(p)}y(t), \dots}_{\text{对应分母参数}}, \underbrace{\frac{1}{A(p)}u(t), \frac{p}{A(p)}u(t), \dots}_{\text{对应分子参数}} \right]^\top $$

于是，残差就神奇地变成了：
$$ \varepsilon(t) = y_f(t) - \boldsymbol{\Phi}^\top(t) \boldsymbol{\theta} $$

这就是论文中的 **公式 (2.32)**。



## 这种伪线性回归存在的问题1：在回归矩阵中引入噪声→参数估计有偏

残差重写成的**伪线性形式**（公式 2.32）：
$$ \varepsilon(t_k) \approx y_f(t_k) - \boldsymbol{\Phi}^\top(t_k) \boldsymbol{\theta} $$

**回归向量 $\boldsymbol{\Phi}(t_k)$** （公式 2.34）：
$$ \boldsymbol{\Phi}(t_k) = \left[ \underbrace{\frac{-p}{A(p)}\mathbf{y}(t_k), \dots, \frac{-p^n}{A(p)}\mathbf{y}(t_k)}_{\text{注意这里！}}, \quad \frac{1}{A(p)}\mathbf{u}(t_k), \dots \right]^\top $$

**陷阱出现了！**
为了让方程看起来像线性的，我们把**输出信号 $\mathbf{y}(t_k)$** 放进了回归向量 $\boldsymbol{\Phi}$ 中。

**由于事实1和2，所以事实3和4**

*   **事实 1**：回归向量 $\boldsymbol{\Phi}$ 中包含测量输出 $\mathbf{y}(t_k)$ 。
*   **事实 2**：测量输出 $\mathbf{y}(t_k)$ 必然包含噪声 $\mathbf{v}(t_k)$。
*   **事实 3**：所以，回归向量 $\boldsymbol{\Phi}$ 包含噪声 $\mathbf{v}(t_k)$。
*   **事实 4**：这就导致回归矩阵与误差项**高度相关**（Correlation）。



#### **如果用普通最小二乘法解**

如果你在这一步直接使用**普通的最小二乘法**（LS）的解析解去解\(\hat{\boldsymbol{\theta}}_{LS}\)：
\[
\hat{\boldsymbol{\theta}}_{LS} = (\mathbf{\Phi}^\top \mathbf{\Phi})^{-1} \mathbf{\Phi}^\top \mathbf{Y}
\]
其中**真实模型输出=估计模型输出+噪声**，可以表示为
\[
\mathbf{Y} = \mathbf{\Phi} \boldsymbol{\theta} + \mathbf{V}\\
 y(t_k) = x(t_k) + v(t_k)\\
 其中模型估计x(t_k)为\boldsymbol{\Phi}^\top(t_k) \boldsymbol{\theta}
\]
将\(\mathbf{Y}\)代入上式，得到\(\hat{\boldsymbol{\theta}}_{LS}\)的估计值
\[
\begin{aligned}
\hat{\boldsymbol{\theta}}_{LS} &= (\mathbf{\Phi}^\top \mathbf{\Phi})^{-1} \mathbf{\Phi}^\top (\mathbf{\Phi} \boldsymbol{\theta} + \mathbf{V}) \\
&= \underbrace{(\mathbf{\Phi}^\top \mathbf{\Phi})^{-1} \mathbf{\Phi}^\top \mathbf{\Phi}}_{\text{变成单位矩阵 } \mathbf{I}} \boldsymbol{\theta} + (\mathbf{\Phi}^\top \mathbf{\Phi})^{-1} \mathbf{\Phi}^\top \mathbf{V} \\
&= \boldsymbol{\theta} + \underbrace{(\mathbf{\Phi}^\top \mathbf{\Phi})^{-1} \mathbf{\Phi}^\top \mathbf{V}}_{\text{偏差项 (Bias Term)}}
\end{aligned}
\]
**结论来了：**
你的估计值 $\hat{\boldsymbol{\theta}}_{LS}$ 等于 **真值** 加上 **偏差项**。

我们当然是想要这个偏差项是为0的

#### **最小二乘估计的误差公式：**

\[
\hat{\boldsymbol{\theta}}_{LS} - \boldsymbol{\theta} = (\mathbf{\Phi}^\top \mathbf{\Phi})^{-1} \mathbf{\Phi}^\top \mathbf{V}
\]

为了分析 $N \to \infty$ 时的性质，我们在分子分母同时除以 $N$（数据量），这样可以将“求和”转化为“平均”：
\[
\hat{\boldsymbol{\theta}}_{LS} - \boldsymbol{\theta} = \left( \frac{1}{N} \mathbf{\Phi}^\top \mathbf{\Phi} \right)^{-1} \left( \frac{1}{N} \mathbf{\Phi}^\top \mathbf{V} \right)
\]
根据**大数定律**，当数据量 $N$ 趋于无穷大时，**样本的平均值**会收敛到其**数学期望**（Expectation）。我们将极限符号写为 $\text{plim}$ (probability limit)：
\[
\text{plim}_{N \to \infty} (\hat{\boldsymbol{\theta}}_{LS} - \boldsymbol{\theta}) = \underbrace{\left( \mathbb{E}[\boldsymbol{\phi}_k \boldsymbol{\phi}_k^\top] \right)^{-1}}_{\text{矩阵 A}} \cdot \underbrace{\mathbb{E}[\boldsymbol{\phi}_k v_k]}_{\text{向量 B}}
\]
这里：

*   $\boldsymbol{\phi}_k$ 是第 $k$ 个时刻的回归向量（即 $\mathbf{\Phi}$ 的第 $k$ 行）。
*   $v_k$ 是第 $k$ 个时刻的噪声。

#### 分析结果：

1.  **矩阵 A**：$\mathbb{E}[\boldsymbol{\phi}_k \boldsymbol{\phi}_k^\top]$ 是回归向量的自协方差矩阵。只要你的输入信号是“持续激励”的（Rich enough），这个矩阵就是**正定且可逆的**。它是一个常数矩阵，**不等于 0**。
2.  **向量 B**：$\mathbb{E}[\boldsymbol{\phi}_k v_k]$ 是回归向量与噪声的互协方差向量。

**结论：**
要让整体偏差趋于 0，**必须且只需**满足：
\[
 \mathbb{E}[\boldsymbol{\phi}_k v_k] = \mathbf{0} 
\]

这就是“**回归向量与噪声正交（不相关）**”的数学定义。

但是由于**回归向量受到噪声污染**，实际上变成了两部分的叠加：
\[
\phi_k = \underbrace{\frac{1}{A(p)}x_k}_{\text{信号部分}} + \underbrace{\frac{1}{A(p)}v_k}_{\text{噪声部分}}
\]
**计算期望 $\mathbb{E}[\phi_k v_k]$：**

我们将 $\phi_k$ 代入：
$$ \begin{aligned} \mathbb{E}[\phi_k v_k] &= \mathbb{E}\left[ \left( \frac{1}{A(p)}x_k + \frac{1}{A(p)}v_k \right) v_k \right] \\ &= \underbrace{\mathbb{E}\left[ \frac{1}{A(p)}x_k \cdot v_k \right]}_{\text{第一项}} + \underbrace{\mathbb{E}\left[ \left( \frac{1}{A(p)}v_k \right) \cdot v_k \right]}_{\text{第二项}} \end{aligned} $$

* **第一项**：$x_k$ 是真实输出，只由 $u_k$ 决定。开环下 $u_k$ 与 $v_k$ 无关，所以这一项为 **0**。

* **第二项（关键！）**：
  令 $\tilde{v}_k = \frac{1}{A(p)}v_k$（滤波后的噪声）。这一项实际上是计算**滤波后的噪声 $\tilde{v}$ 与原始噪声 $v$ 的相关性**。
  $$ \mathbb{E}[\tilde{v}_k \cdot v_k] $$

  设滤波器 $\frac{1}{A(p)}$ 的单位脉冲响应是 $h_0, h_1, h_2, \dots$，那么：
  $$ \tilde{v}_k = h_0 v_k + h_1 v_{k-1} + h_2 v_{k-2} + \dots $$

  代入计算：
  $$ \mathbb{E}[(h_0 v_k + h_1 v_{k-1} + \dots) \cdot v_k] = h_0 \underbrace{\mathbb{E}[v_k^2]}_{\sigma_v^2} + h_1 \underbrace{\mathbb{E}[v_{k-1}v_k]}_{0} + \dots $$

  结果等于 $h_0 \sigma_v^2$。

  *   $\sigma_v^2$ 是噪声方差，恒大于 0。
  *   $h_0$ 是滤波器的直通项（第一项系数）。在连续系统离散化后或者特定结构下，通常 $h_0 \neq 0$。

**结果**：$\mathbb{E}[\phi_k v_k] = h_0 \sigma_v^2 \neq 0$。
**结论**：偏差存在，且无法随数据量增加而消除。


**即使是开环系统，因为你的算法人为地把含噪输出引入了回归矩阵，最小二乘法也会产生严重的偏差。** 这种偏差通常被称为“方程误差偏差”（Equation Error Bias）。

## 这种伪线性回归存在的问题2：依赖上次迭代参数

虽然上面的公式看起来像是完美的线性回归 $Y = \Phi \theta$，但这里有一个巨大的**逻辑漏洞**：

你在构建回归向量 $\boldsymbol{\Phi}(t)$ 的时候，用到了滤波器 $\frac{1}{A(p)}$。
但是，**$A(p)$ 本身就包含了我们要估计的未知参数 $a_i$ 啊！**

你还没解出 $a_i$，怎么能用它来构造滤波器去解 $a_i$ 呢？这不就是“我揪着自己的头发想把自己提离地面”吗？

**解决方案：迭代 (Iteration)**
这就是 SRIVC 算法的核心机制：

1.  第 $j$ 次迭代时，我们手里有一组旧的参数估计 $\hat{A}^{(j)}(p)$。
2.  我们用这个**旧的** $A^{(j)}$ 来构造滤波器和回归向量 $\boldsymbol{\Phi}$。
3.  此时 $\boldsymbol{\Phi}$ 就变成了已知数（常数矩阵）。
4.  然后我们解线性方程，求出新的参数 $\hat{A}^{(j+1)}$。
5.  用新的 $\hat{A}^{(j+1)}$ 更新滤波器，重复上述过程。

所以它叫“伪”线性回归，因为它依赖于上一次迭代的参数。

## “为了线性化而进行的数学变换”也强行引入了偏差

现在回头看回归向量的定义（公式 2.34）：
$$ \boldsymbol{\Phi}(t_k) = \left[ \frac{-p}{A(p)}\mathbf{y}(t_k), \dots \right]^\top $$

为了把 $A(p)$ 从分母里“拿出来”变成线性的 $a_i$，我们不得不把它乘到 $y(t)$ 这一边。
这就导致**回归向量里直接包含了输出信号 $y(t)$**。

*   **物理现实**：测量到的输出 $y(t)$ = 真实输出 $x(t)$ + 噪声 $v(t)$。
*   **数学后果**：回归向量 $\boldsymbol{\Phi}$ 里混入了噪声 $v(t)$。
*   **致命一击**：在最小二乘法中，$\boldsymbol{\Phi}$ 必须独立于误差。但现在 $\boldsymbol{\Phi}$ 里有 $v$，残差里也有 $v$，两者**自相关**了！

这就是为什么即便是在开环系统中，这种“为了线性化而进行的数学变换”也强行引入了偏差，迫使我们必须祭出**工具变量（IV）**这个大杀器来修正它。







## 问题 1：如果不用线性回归，是否就不需要工具变量了？

**答案是：原则上是的。**

工具变量（IV）是专门为了拯救**线性回归（或伪线性回归）**而诞生的“补丁”。

*   **线性回归的困境**：为了把模型凑成 $Y = \Phi \theta$ 这种好解的线性形式，我们被迫把含噪的输出信号 $y$ 放进了回归矩阵 $\Phi$ 里。这就导致了 $\Phi$ 与噪声相关，从而产生偏差。IV 就是用来修正这个偏差的。
*   **非线性优化的路径**：如果你不强求线性形式，而是直接去解原始的**输出误差（Output Error, OE）**最小化问题：
    $$ \min_{\theta} \sum (y - \hat{y}(\theta))^2 $$
    这里的预测输出 $\hat{y}(\theta)$ 是参数 $\theta$ 的非线性函数（因为有分母 $A(p)$）。
    在这种方法中，我们通常使用梯度下降法或高斯-牛顿法直接搜索最优解。
    **在这种框架下，不需要构造工具变量矩阵。** 只要优化算法收敛到全局最优，得到的估计值通常就是**无偏（Unbiased）**且**一致（Consistent）**的。

### **那为什么还要用线性回归+IV（如 SRIVC）呢？**

因为**非线性优化太难了**（见问题 4）。SRIVC 的精妙之处在于：它用“线性回归+IV”的方式，巧妙地避开了非线性优化的陷阱，却达到了和非线性优化一样好的统计效果。



# 工具变量的引入：由于闭环噪声与非线性问题

由于闭环噪声问题和非线性问题，我们引入工具变量来进行解决

#### 闭环噪声问题

*   **开环 (Open-loop)**：输入 $\mathbf{u}$ 与噪声 $\mathbf{v}$ 不相关。这是最简单的情况。
*   **闭环 (Closed-loop)**：如图2.2所示，由于控制器的存在，输出的噪声 $\mathbf{v}$ 会通过反馈回路影响输入 $\mathbf{u}$。
    *   **严谨推导**：公式 (2.7) $\mathbf{u}(t_k) = \mathbf{S}^*_{uo}(q)\mathbf{r}(t_k) - \mathbf{S}^*_{uo}(q)\mathbf{v}(t_k)$ 展示了这一点。输入 $\mathbf{u}$ 包含了一项与噪声 $\mathbf{v}$ 相关的项。
    *   **后果**：如果在闭环下直接使用普通的最小二乘法，会导致**偏差（Bias）**。这也是后续章节引入“工具变量法（Instrumental Variable）”的主要原因之一。

#### 非线性问题

现在我们已经知道了，为了把难解的**非线性**问题转化为好解的**线性**问题，我们引入了“**含噪回归向量**”，这迫使我们必须使用**工具变量**来修正偏差。

### 工具变量的本质是为了解决噪声问题

对于开环，我们需要用工具变量处理非线性转伪线性回归中的“含噪声回归向量”

对于闭环，我们既需要处理”含噪回归向量”，又需要处理反馈导致的非纯净控制量输入（受噪声影响）

因此，需要引入工具变量

# 工具变量要满足的前提条件

既然是我们为了算法便利（线性化）而引入了相关性偏差，我们就必须用工具变量把这个偏差消除掉。

消除偏差的本质是找到一个新的回归矩阵，它既能是代价函数的最优解，而且与残差的相关性为0？

#### 工具变量法（IV）的**核心估计方程**：

\[
\frac{1}{N} \sum_{k=1}^{N} \hat{\boldsymbol{\Phi}}(t_k) \boldsymbol{\Lambda}_0^{-1} \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta}) = \mathbf{0} \tag{2.16}
\]

#### .1条件：和噪声正交

公式 (2.16) 本质上是在说：**工具变量 $\hat{\boldsymbol{\Phi}}$ 和残差 $\boldsymbol{\varepsilon}$ 的加权相关性（Correlation）必须为零。**

在统计学中，如果我们认为参数 $\boldsymbol{\beta}$ 已经收敛到了真值 $\boldsymbol{\beta}^*$，那么残差 $\boldsymbol{\varepsilon}$ 就只剩下纯粹的噪声 $\mathbf{v}$。
此时，公式 (2.16) 等价于要求：
$$ \mathbb{E} [ \hat{\boldsymbol{\Phi}}^\top \mathbf{v} ] = \mathbf{0} $$

##### **为什么要满足它？**

*   **回顾 LS 的失败**：在普通最小二乘法中，回归矩阵 $\boldsymbol{\Phi}$ 含有噪声，导致 $\mathbb{E}[\boldsymbol{\Phi}^\top \mathbf{v}] \neq \mathbf{0}$。这意味着噪声被误认为是信号的一部分，从而拉偏了参数估计（偏差）。
*   **IV 的修正**：我们强行寻找一组变量 $\hat{\boldsymbol{\Phi}}$，它们是从纯净的输入信号推导出来的，与物理噪声 $\mathbf{v}$ **毫无瓜葛**（统计独立）。
*   **操作层面**：我们在求解参数时，强迫“工具变量”和“剩余误差”正交（垂直）。既然工具变量是干净的，那么这就保证了剩余下来的误差里不包含任何与系统动态有关的信息，只剩下纯噪声。

**本质**：满足这个条件，就是**切断噪声对参数估计的干扰通道**。

##### 为什么公式（2.16）不是直接和噪声正交？

因为噪声是不可测的，我们不知道噪声实际上是多少，我们知道的只有y，真实系统输出我们也不知道。

如果都知道噪声了，那直接用测量系统输出减去噪声就是真实系统输出了，那我们就不用考虑噪声问题了。



虽然公式 (2.16) 没有直接写出来，但要解出这个方程，隐含着**第二个条件**：
**工具矩阵 $\hat{\boldsymbol{\Phi}}$ 必须与原回归矩阵 $\boldsymbol{\Phi}$ 强相关。**

#### .2隐藏条件（工具矩阵 $\hat{\boldsymbol{\Phi}}$ 必须与原回归矩阵 $\boldsymbol{\Phi}$ 强相关）

让我们试着解一下公式 (2.16) 的线性化版本（假设 $\boldsymbol{\varepsilon} \approx \mathbf{y} - \boldsymbol{\Phi}\boldsymbol{\beta}$）：
$$ \hat{\boldsymbol{\Phi}}^\top (\mathbf{y} - \boldsymbol{\Phi}\hat{\boldsymbol{\beta}}) = 0 $$
移项得到：
$$ \hat{\boldsymbol{\Phi}}^\top \mathbf{y} = \hat{\boldsymbol{\Phi}}^\top \boldsymbol{\Phi} \hat{\boldsymbol{\beta}} $$
$$ \hat{\boldsymbol{\beta}} = (\hat{\boldsymbol{\Phi}}^\top \boldsymbol{\Phi})^{-1} \hat{\boldsymbol{\Phi}}^\top \mathbf{y} $$

注意那个求逆的部分 $(\hat{\boldsymbol{\Phi}}^\top \boldsymbol{\Phi})^{-1}$。如果 $\hat{\boldsymbol{\Phi}}$ 和 $\boldsymbol{\Phi}$ 不相关（比如 $\hat{\boldsymbol{\Phi}}$ 全是 0，或者全是随机乱数），这个矩阵就会奇异（不可逆），方程这就没法解了。

##### 论文中的设计

文中提到 $\hat{\boldsymbol{\Phi}}$ 是由 "noise-free estimates of the input and output signals" 构建的。

*   **本质**：它不仅要干净（满足条件一），还必须**长得像**真实的回归矩阵（满足条件二）。它必须包含系统真实的动力学特征（极点、零点信息）。
*   这就是为什么要用**辅助模型（Auxiliary Model）**来生成工具变量的原因：为了确保它足够“像”真数据，从而能捕捉到系统的参数信息。

#### .3为什么要乘 $\boldsymbol{\Lambda}_0^{-1}$？（最优性条件）

你可能注意到了公式中间夹着一个 $\boldsymbol{\Lambda}_0^{-1}$（噪声协方差的逆）。这是为了满足**第三个高阶需求：统计效率（Efficiency）**。

*   **问题**：能消除偏差的工具变量有很多种。只要和噪声无关就行。哪一种最好？
*   **答案**：让估计方差最小的那一种。
*   **原理**：根据 **高斯-马尔可夫定理（Gauss-Markov Theorem）** 的推广，当我们在方程中加入 $\boldsymbol{\Lambda}_0^{-1}$ 进行加权时，我们实际上是在对数据进行“白化”（Whitening）或者说“标准化”。
    *   对于信噪比高的通道，赋予高权重。
    *   对于信噪比低的通道，赋予低权重。

**本质**：这个权重项确保了我们不仅算得**对**（无偏），而且算得**准**（方差最小，利用了所有有效信息）。

#### 总结：三大条件——正交性、相关性、加权

让我们把视角拉高，公式 (2.16) $\sum \hat{\boldsymbol{\Phi}} \boldsymbol{\Lambda}^{-1} \boldsymbol{\varepsilon} = 0$ 的本质是什么？

它是在寻找一个参数 $\hat{\boldsymbol{\beta}}$，使得：
**“模型剩下的残差 $\boldsymbol{\varepsilon}$，在工具变量 $\hat{\boldsymbol{\Phi}}$ 所张成的空间里，没有任何投影。”**

1.  **纯净性（正交性）**：因为 $\hat{\boldsymbol{\Phi}}$ 是纯净空间，如果残差在这个空间里没投影，说明残差里不含任何系统动态信息，剩下的全是纯噪声。
2.  **有效性（相关性）**：因为 $\hat{\boldsymbol{\Phi}}$ 是由模型生成的，它代表了我们对系统动态的“最佳猜测方向”。
3.  **最优性（加权）**：$\boldsymbol{\Lambda}_0^{-1}$ 调整了投影的角度，使其符合噪声的统计分布特性。

**一句话概括**：
这个方程是**用“纯净的模拟信号”作为筛子，在“加权统计空间”中，把测量数据里的“系统动力学”完全过滤出来，只留下“纯噪声”在筛子里。** 此时筛出来的参数，就是真实参数。





# 如何得到工具变量

## 比较双解得到构造依据：选择工具矩阵等于残差梯度的负数

通过比较代价函数最优解（2.20）
\[
\frac{\partial V}{\partial \boldsymbol{\beta}} = \frac{1}{N} \sum_{k=1}^{N} \left( \frac{\partial \boldsymbol{\varepsilon}}{\partial \boldsymbol{\beta}} \right)^\top \boldsymbol{\Lambda}_0^{-1} \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta})=\mathbf{0} \tag{2.20}
\]
与工具变量方法的解（2.16）
\[
\frac{1}{N} \sum_{k=1}^{N} \hat{\boldsymbol{\Phi}}(t_k) \boldsymbol{\Lambda}_0^{-1} \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta}) = \mathbf{0} \tag{2.16}
\]
这两个解我们需要同时满足，因此，也即让 IV 方法的解**恰好**也是上面那个代价函数的**最优解**

为了同时满足，**我们需要选择工具矩阵等于残差梯度的负数（即灵敏度矩阵）！**
\[
 \hat{\boldsymbol{\Phi}} \approx -\frac{\partial \boldsymbol{\varepsilon}}{\partial \boldsymbol{\beta}} 
\]
因此下一步就是：推导残差对参数的梯度，直接作为工具矩阵的构造依据。

## 推导残差对参数的梯度，直接作为工具矩阵的构造依据（可以略）

现在我们的任务变成了单纯的数学演算：对加法模型求导。

#### 3.1 模型的结构

$$ \text{模型} \mathbf{P} = \sum \frac{1}{A_i(p)} \mathbf{B}_i(p) $$
$$ \text{残差} \boldsymbol{\varepsilon} = \mathbf{y} - \mathbf{P}\mathbf{u} $$

因为是加法结构，对参数 $\boldsymbol{\theta}_i$ 求导时，其他子系统 $\mathbf{P}_{j \neq i}$ 的导数都是 0。我们只需要关注第 $i$ 个子系统。

#### 3.2 对分母系数 $a_{i,j}$ 求导 (公式 2.23)

第 $i$ 个子系统的项是 $\frac{\mathbf{B}_i}{A_i} \mathbf{u}$。
我们利用除法求导法则（Quotient Rule）：$\left(\frac{1}{A}\right)' = -\frac{A'}{A^2}$。

$$ \frac{\partial}{\partial a_{i,j}} \left( - \frac{\mathbf{B}_i(p)}{A_i(p)} \mathbf{u} \right) = - \mathbf{B}_i(p) \mathbf{u} \cdot \left( -\frac{p^j}{A_i(p)^2} \right) = \frac{p^j \mathbf{B}_i(p)}{A_i(p)^2} \mathbf{u} $$

*   **直观理解**：这是“灵敏度函数”。它表示如果你稍微改变一下极点（分母），输出会发生多大变化。注意分母变成了平方 $A^2$，这是典型的极点灵敏度特征。

#### 3.3 对分子系数 $\mathbf{B}_{i,j}$ 求导 (公式 2.24)

这一步比较tricky，因为分子是一个**矩阵**。论文在这里使用了向量化算子 ($\text{vec}$) 和克罗内克积 ($\otimes$)。

我们要对 $\frac{1}{A_i(p)} \mathbf{B}_i(p) \mathbf{u}(t_k)$ 中的 $\mathbf{B}_i$ 求导。
利用恒等式 (1.3) $\mathbf{X}\mathbf{b} = (\mathbf{b}^\top \otimes \mathbf{I})\text{vec}(\mathbf{X})$，我们可以把 $\mathbf{B}_i \mathbf{u}$ 写成：
$$ \mathbf{B}_i(p)\mathbf{u} = (\mathbf{u}^\top \otimes \mathbf{I}_q) \text{vec}(\mathbf{B}_i) $$

令 $\mathbf{U}^\top = \mathbf{u}^\top \otimes \mathbf{I}_q$（这是一个扩展的输入矩阵），那么模型项就是：
$$ \frac{1}{A_i(p)} \mathbf{U}^\top \text{vec}(\mathbf{B}_i) $$

这一项对 $\text{vec}(\mathbf{B}_i)$ 的系数 $\mathbf{B}_{i,j}$ 求导就非常简单了，就是把系数拿掉，剩下：
$$ \frac{p^j}{A_i(p)} \mathbf{U}^\top $$

*   **直观理解**：分子参数主要影响增益，所以它的灵敏度就是经过滤波（$1/A$）后的输入信号。

---

### 4. 组装得到最终结果：最优工具矩阵 (The Optimal Instrument)

将上述两个导数结果拼起来，就得到了公式 (2.27)，即对于第 $i$ 个模态的最优工具块 $\hat{\boldsymbol{\Phi}}_i$：

$$ \hat{\boldsymbol{\Phi}}_i(t_k, \boldsymbol{\theta}_i) = \left[ \underbrace{\frac{-p \mathbf{B}_i}{A_i^2}\mathbf{u}, \dots}_{\text{分母部分}}, \underbrace{\frac{1}{A_i}\mathbf{U}^\top, \dots}_{\text{分子部分}} \right]^\top $$

**这个公式的物理意义非常明确：**

1.  **分母部分（左边）**：这是通过模拟得到的信号。它相当于把输入 $\mathbf{u}$ 通过一个传递函数 $\frac{-p^j \mathbf{B}_i}{A_i^2}$。这意味着我们在用模型当前的估计值来生成“无噪声的输出梯度”。
2.  **分子部分（右边）**：这是通过模拟得到的输入。它相当于把输入 $\mathbf{u}$ 通过预滤波器 $\frac{p^j}{A_i}$。

**教授的特别提醒**：
细心的你可能会发现一个**悖论**：
要计算这个最优工具矩阵 $\hat{\boldsymbol{\Phi}}$，我们需要知道 $\mathbf{B}_i$ 和 $A_i$。
但这正是我们要去估计的参数啊！如果我知道了，还辨识什么？

**解决办法**：这就是为什么 SRIVC 是一个**迭代算法 (Iterative Algorithm)**。

1.  先猜一组参数 $\boldsymbol{\beta}^{(0)}$。
2.  用这就组参数构建工具矩阵 $\hat{\boldsymbol{\Phi}}(\boldsymbol{\beta}^{(0)})$。
3.  解线性方程得到新的参数 $\boldsymbol{\beta}^{(1)}$。
4.  重复上述过程，直到参数不再变化。





# 辨识准则：工具变量准则

由于（2.20）和（2.16）一样了，所以OE准则也就是工具变量准则了

所考虑的辨识准则是辅助变量准则[15, 50]。不是直接最小化成本函数，而是通过求解辅助变量方程来获得参数估计。该方程通过将辅助信号矩阵与残差函数(2.15)相关联而得到

$$\hat{\boldsymbol{\beta}}\in\operatorname*{sol}_{\boldsymbol{\beta}}\left\{\frac {1}{N}\sum_{k=1}^{N}\hat{\boldsymbol{\Phi}}\left(t_{k}\right)\boldsymbol{\Lambda }^{-1}_{0}\boldsymbol{\varepsilon}(t_{k},\boldsymbol{\beta})=\mathbf{0}\right\}, \tag{2.16}$$

其中 $\hat{\boldsymbol{\Phi}}\left(t_{k}\right)\in\mathbb{R}^{\sum_{i}n_{i}+pq\left(m_ {i}+K\right)\times q}$ 是辅助矩阵，$\boldsymbol{\Lambda}_{0}\in\mathbb{R}^{q\times q}$ 表示真实噪声协方差矩阵

$$\boldsymbol{\Lambda}_{0}=\mathbb{E}\{\mathbf{v}(t_{k})\mathbf{v}^{\top}(t_{k}) \}. \tag{2.17}$$

辅助矩阵 $\hat{\boldsymbol{\Phi}}(t_{k})$ 旨在消除渐近偏差。这通常通过使用输入和输出信号的无噪声估计来构建矩阵来实现，从而避免与噪声扰动 $\mathbf{v}(t_{k})$ 的相关性。辅助矩阵的最优表达式可以根据系统运行的具体环境推导。在接下来的小节中，将推导在开环或闭环情况下运行的加性模型结构(2.1)的最优性准则和由此产生的最优辅助表达式。



# 如何利用工具变量进行辨识？**SRIVC 算法**

你已经手里握着一把“屠龙刀”——**最优工具矩阵 $\hat{\boldsymbol{\Phi}}$**。你知道它能消除偏差，也知道它长什么样。

但现在你面临一个尴尬的局面：**这把刀还在石头里拔不出来。**

为什么？因为 $\hat{\boldsymbol{\Phi}}$ 的定义里包含了系统参数 $\boldsymbol{\beta}$（你需要用参数来构建滤波器 $1/A(p)$ 和模拟输出 $\hat{\mathbf{x}}$），而 $\boldsymbol{\beta}$ 正是你想要计算的东西。

这就构成了“先有鸡还是先有蛋”的死循环。打破这个循环，并最终算出数值的步骤，就是论文 **2.4节 "Proposed solution procedure"** 的内容。

这一步的核心是：**SRIVC 算法（迭代精炼工具变量法）**。

下面我分三步详细解释这“下一步”具体是怎么走的。

---

### 第一步：伪线性化 (Pseudo-Linear Form) —— 把靶子竖起来

我们在之前讨论“陷阱”时提到过，为了求解方便，我们需要把非线性的残差强行写成线性的样子。这是算法运行的基础。

在 **2.4.1节** 中，论文通过代数变换，把残差写成了：
$$ \boldsymbol{\varepsilon}(t_k, \boldsymbol{\beta}) = \mathbf{y}_f(t_k, \boldsymbol{\beta}) - \boldsymbol{\Phi}^\top(t_k, \boldsymbol{\beta}) \boldsymbol{\theta} $$

*   **$\mathbf{y}_f$**：滤波后的输出（作为因变量）。
*   **$\boldsymbol{\Phi}$**：回归矩阵（包含滤波后的输入和输出）。
*   **$\boldsymbol{\theta}$**：我们要解的线性参数。

**关键点**：这叫做“伪”线性，因为等式右边的 $\mathbf{y}_f$ 和 $\boldsymbol{\Phi}$ 本身也都依赖于参数 $\boldsymbol{\beta}$（藏在滤波器 $1/A(p)$ 里）。

但我们可以耍一个赖皮：**“如果我们暂时假定滤波器里的参数是已知的，那么这就变成了一个标准的线性方程组 $Ax=b$。”**

---

### 第二步：迭代求解 (Algorithm 1) —— 磨刀不误砍柴工

我们采用**迭代（Iterative）**的策略。既然不知道真实参数，我们就先猜一个，然后一轮一轮地修正。

这个过程就像**对焦**：一开始画面模糊（参数不准），我们根据模糊的画面调整镜头，画面变清晰一点；再根据新画面继续微调，直到画面最清晰。

**具体的“下一步”操作流程（对应论文 Algorithm 1）：**

#### 1. 初始化 (Initialization)
*   **动作**：随便给一组初始参数 $\boldsymbol{\beta}^{(0)}$。
*   **通常做法**：假设滤波器 $A(p)=1$（即不滤波），用普通的最小二乘法算一个初始解。虽然这个解有偏差，但大概方向是对的。

#### 2. 第 $j$ 轮迭代 (Iteration $j$)
假设我们手里已经有了第 $j$ 轮的参数估计 $\boldsymbol{\beta}^{(j)}$。

*   **构造滤波器**：用 $\boldsymbol{\beta}^{(j)}$ 中的分母系数构建滤波器 $F(p) = A^{(j)}(p)$。
*   **构造数据**：
    *   **回归矩阵 $\boldsymbol{\Phi}$**：把测量数据 $u, y$ 通过滤波器，填进去。
    *   **因变量 $\mathbf{y}_f$**：把测量输出 $y$ 通过滤波器。
    *   **关键一步：构造工具矩阵 $\hat{\boldsymbol{\Phi}}$**：用 $\boldsymbol{\beta}^{(j)}$ 进行系统仿真，得到纯净输出 $\hat{x}$，再通过滤波器，填入工具矩阵。

#### 3. 求解更新 (Update)
现在，$\boldsymbol{\Phi}$ 和 $\hat{\boldsymbol{\Phi}}$ 都变成了**已知的数字矩阵**。我们利用 IV 的解析解公式（公式 2.52）算出第 $j+1$ 轮的参数：

$$ \hat{\boldsymbol{\beta}}^{(j+1)} = \left[ \sum_{k=1}^N \hat{\boldsymbol{\Phi}}(t_k) \boldsymbol{\Lambda}^{-1} \boldsymbol{\Phi}^\top(t_k) \right]^{-1} \left[ \sum_{k=1}^N \hat{\boldsymbol{\Phi}}(t_k) \boldsymbol{\Lambda}^{-1} \mathbf{y}_f(t_k) \right] $$

*   这是一个标准的线性代数运算（形如 $x = (A^T A)^{-1} A^T b$）。
*   这一步算出来的是**当前条件下**的最优解。

#### 4. 稳不稳定？(Stability Check)
*   因为是连续系统，随机生成的参数可能会导致 $A(p)$ 不稳定（极点在右半平面）。
*   如果 $1/A(p)$ 不稳定，滤波后的信号会爆炸（无穷大）。
*   **操作**：论文提到（Algorithm 1 第12行），如果发现极点不稳定，就把它的实部翻转一下（$Re(s) \to -Re(s)$），强制把极点推回左半平面。这只改变稳定性，不改变幅频特性。

#### 5. 循环判断
*   比较 $\boldsymbol{\beta}^{(j+1)}$ 和 $\boldsymbol{\beta}^{(j)}$。如果变化很小（比如小于 $10^{-5}$），说明收敛了，停止。否则，$j = j+1$，回到第2步继续。

---

### 第三步：处理 MIMO 的复杂性 (Stacking)

如果你的系统是多输入多输出（MIMO）的，还有一个额外的步骤：**堆叠（Stacking）**。

论文在 **2.4.1 节后半部分** 提到，因为有 $K$ 个模态（子系统），我们不能一个一个单独算，而要一起算。

*   **子问题**：先把每个模态 $i$ 看作独立的，写成伪线性形式。
*   **大矩阵**：把所有模态的方程堆在一起，形成一个巨大的矩阵方程（公式 2.48）：
    $$ \mathcal{E}^\top = \boldsymbol{\Upsilon}^\top - \boldsymbol{\Phi}^\top \mathcal{B} $$
    其中 $\mathcal{B}$ 是包含所有模态参数的大块对角矩阵。

这一步在编程实现时很重要，但在概念上，它依然遵循上面的迭代逻辑









# 疑问

## 线性vs非线性，工具变量

如果不用线性回归，是否就不需要工具变量了？

除了线性回归还有哪些线性优化方法？

非线性优化方法有哪些？

为什么非线性优化方法对初始值非常敏感，容易陷入局部最优？

线性回归方法不会对初始值敏感，容易陷入局部最优吗？

你好！我是你的教授。

这五个问题问得非常棒，逻辑环环相扣，直击系统辨识的核心痛点——**计算复杂性与统计准确性之间的权衡（Trade-off）**。

让我们逐层剖析，彻底弄清楚这些概念之间的纠葛。

---

### 问题 1：如果不用线性回归，是否就不需要工具变量了？

**答案是：原则上是的。**

工具变量（IV）是专门为了拯救**线性回归（或伪线性回归）**而诞生的“补丁”。

*   **线性回归的困境**：为了把模型凑成 $Y = \Phi \theta$ 这种好解的线性形式，我们被迫把含噪的输出信号 $y$ 放进了回归矩阵 $\Phi$ 里。这就导致了 $\Phi$ 与噪声相关，从而产生偏差。IV 就是用来修正这个偏差的。
*   **非线性优化的路径**：如果你不强求线性形式，而是直接去解原始的**输出误差（Output Error, OE）**最小化问题：
    $$ \min_{\theta} \sum (y - \hat{y}(\theta))^2 $$
    这里的预测输出 $\hat{y}(\theta)$ 是参数 $\theta$ 的非线性函数（因为有分母 $A(p)$）。
    在这种方法中，我们通常使用梯度下降法或高斯-牛顿法直接搜索最优解。
    **在这种框架下，不需要构造工具变量矩阵。** 只要优化算法收敛到全局最优，得到的估计值通常就是**无偏（Unbiased）**且**一致（Consistent）**的。

**那为什么还要用线性回归+IV（如 SRIVC）呢？**
因为**非线性优化太难了**（见问题 4）。SRIVC 的精妙之处在于：它用“线性回归+IV”的方式，巧妙地避开了非线性优化的陷阱，却达到了和非线性优化一样好的统计效果。

---

### 问题 2：除了线性回归还有哪些线性优化方法？

在系统辨识和广义优化领域，除了最小二乘（LS）这种标准的线性回归，还有以下几类“线性”或“凸”方法：

1.  **递归最小二乘 (Recursive Least Squares, RLS)**：
    *   这就是线性回归的“在线版”。数据来一个处理一个，参数实时更新。常用于自适应控制。
2.  **总体最小二乘 (Total Least Squares, TLS)**：
    *   普通线性回归假设只有 $Y$ 有噪声，$\Phi$ 是准的。TLS 假设 $Y$ 和 $\Phi$ 都有噪声（Errors-in-variables 模型）。它通过奇异值分解（SVD）求解。
3.  **子空间辨识法 (Subspace Identification, e.g., N4SID, MOESP)**：
    *   **这是MIMO系统的神器**。它不显式地解回归方程，而是利用线性代数中的投影和 SVD 技术，直接从数据矩阵的行空间或列空间中提取状态空间矩阵 $(A, B, C, D)$。它全是线性运算，非常稳健。
4.  **凸优化 (Convex Optimization)**：
    *   **线性规划 (LP)**：带约束的参数估计。
    *   **二次规划 (QP)**：带约束的最小二乘。
    *   **半定规划 (SDP)**：利用线性矩阵不等式（LMI）求解。这在现代鲁棒控制和稀疏辨识中用得很多。

---

### 问题 3：非线性优化方法有哪些？

当代价函数 $V(\theta)$ 对参数是非线性的时候（比如有很多局部坑洼），我们需要用迭代搜索算法：

1.  **梯度下降法 (Gradient Descent)**：
    *   最基础的方法。沿着梯度的反方向走。缺点：收敛速度慢，容易走出“之”字形。
2.  **牛顿法 (Newton-Raphson)**：
    *   利用二阶导数（海森矩阵 Hessian）。收敛极快，但计算海森矩阵很贵，且如果海森矩阵非正定会失效。
3.  **高斯-牛顿法 (Gauss-Newton, GN)**：
    *   **系统辨识最常用**。它利用雅可比矩阵（梯度）近似海森矩阵。专门用于解决非线性最小二乘问题。
4.  **列文伯格-马夸特法 (Levenberg-Marquardt, LM)**：
    *   **工业界标准**。它是梯度下降和高斯-牛顿的结合体。离最优值远时像梯度下降（稳），近时像高斯-牛顿（快）。MATLAB 的 `lsqnonlin` 默认就用这个。
5.  **信赖域方法 (Trust Region)**：
    *   限制每一步走的距离，保证算法稳定性。
6.  **全局优化算法**：
    *   遗传算法 (GA)、粒子群算法 (PSO)、模拟退火 (SA)。这些属于“乱枪打鸟”，计算量极大，通常只用于极其复杂的黑箱模型。

---

### 问题 4：为什么非线性优化方法对初始值非常敏感，容易陷入局部最优？

想象你在爬山（或者下山），目标是找到**地球上的最低点（马里亚纳海沟）**。

1.  **非凸性 (Non-convexity)**：
    *   线性回归的代价函数像一个**碗**（凸函数）。无论你从哪扔个球，它最终都会滚到碗底（全局最优）。
    *   非线性优化（如输出误差 OE）的代价函数像**连绵起伏的山脉**。有无数个小山谷（局部最优）。
    *   **数学原因**：对于传递函数 $\frac{B(p)}{A(p)}$，参数 $A(p)$ 在分母里。这就构成了有理分式，导致误差曲面非常复杂。

2.  **梯度的短视**：
    *   高斯-牛顿等方法是基于**梯度（坡度）**的。
    *   如果你如果不小心被空投到了喜马拉雅山的一个小坑里（初始值选错了），由于四周都是上坡，算法会以为这里就是最低点，于是停在这里不动了。它看不见隔壁就是马里亚纳海沟。

3.  **参数的强耦合**：
    *   在极点和零点靠得很近，或者模型阶次选高了的时候，参数之间会“打架”。误差曲面会出现很平坦的区域（Flat Valley），算法在上面走不动，不知道该往哪走。

所以，非线性优化**必须**要有一个非常好的**初始猜测值（Initial Guess）**，让算法一开始就落在全局最优解附近的“坑”里。

---

### 问题 5：线性回归方法不会对初始值敏感，容易陷入局部最优吗？

**对于标准的线性回归（OLS）：**
*   **完全不敏感**。
*   **完全没有局部最优**。
*   因为它的解是**解析解（Closed-form solution）**：$\hat{\theta} = (\Phi^\top \Phi)^{-1} \Phi^\top Y$。
*   不管你心里想的初始值是什么，公式一步算出来的结果永远是同一个（那个碗底）。这是线性代数保证的。

**但是（注意这里有个转折）：**

这篇论文里用的 **SRIVC 算法**，虽然名字里带“回归”，但它属于**迭代伪线性回归**。
1.  它依然是迭代的（先猜参数 $\to$ 建滤波器 $\to$ 算回归 $\to$ 得到新参数 $\to$ 重复）。
2.  **它也需要初始值**。
3.  **但是**，SRIVC 的收敛区域（Basin of Attraction）通常比纯粹的非线性优化（如高斯-牛顿）要**大得多**，且地形更平滑。

**结论**：
*   **标准线性回归**：一步到位，无视初始值，无局部最优。
*   **非线性优化**：步步惊心，极度依赖初始值，到处是坑。
*   **SRIVC (论文方法)**：介于两者之间。它利用了线性回归的稳定性，但也保留了迭代的特性。它比非线性优化鲁棒得多，但在极高阶或极复杂的MIMO系统中，依然需要像论文第4章讨论的那样，进行精细的**初始化（Initialization）**。

## 克罗内克积推导

你好！我是你的教授。

这一步确实是数学推导中的“拦路虎”，因为它涉及到**矩阵微积分（Matrix Calculus）**。很多同学看到克罗内克积（Kronecker Product, $\otimes$）和向量化（Vectorization, $\text{vec}$）头就开始痛了。

但请放心，这其实只是一个**“整理抽屉”**的过程。

我们的目标非常单纯：**把混在一起的“已知数据”和“未知参数”彻底分离开。**

让我们像剥洋葱一样，把这个公式 (2.24) 也就是对分子系数 $\mathbf{B}_{i,j}$ 的求导过程，拆解得清清楚楚。

---

### 1. 为什么要这么麻烦？（维度困境）

首先，我们要明白我们在对什么求导。

*   **模型项**：$\mathbf{y}_{part} = \mathbf{B}_i(p) \mathbf{u}(t)$。
*   **输入** $\mathbf{u}$ 是 $p \times 1$ 的向量。
*   **输出** $\mathbf{y}_{part}$ 是 $q \times 1$ 的向量。
*   **分子参数** $\mathbf{B}_{i,j}$ 是一个 $q \times p$ 的**矩阵**。

**问题来了**：我们要计算一个向量 $\mathbf{y}$ 对一个矩阵 $\mathbf{B}$ 的导数。这在标准微积分里是没法直接定义的（结果会是一个三维张量）。

**解决办法**：把矩阵 $\mathbf{B}$ 拉直（Flatten），变成一个长向量 $\text{vec}(\mathbf{B})$。这样问题就变成了“向量对向量求导”，这就好办了，结果是一个雅可比矩阵（Jacobian Matrix）。

---

### 2. 核心工具：恒等式 (1.3) 的魔力

为了把参数 $\mathbf{B}$ 从乘法 $\mathbf{B}\mathbf{u}$ 中“抠”出来，我们需要用到论文中提到的恒等式：

$$ \mathbf{X}\mathbf{b} = (\mathbf{b}^\top \otimes \mathbf{I}_q) \text{vec}(\mathbf{X}) $$

这个公式看起来很吓人，但它的物理意义非常直观。让我举个最简单的 $2 \times 2$ 例子给你看。

假设：
$$ \mathbf{X} = \begin{bmatrix} x_{11} & x_{12} \\ x_{21} & x_{22} \end{bmatrix}, \quad \mathbf{b} = \begin{bmatrix} b_1 \\ b_2 \end{bmatrix} $$

**左边（正常乘法）**：
$$ \mathbf{X}\mathbf{b} = \begin{bmatrix} x_{11}b_1 + x_{12}b_2 \\ x_{21}b_1 + x_{22}b_2 \end{bmatrix} $$

**右边（分离参数）**：
我们把 $\mathbf{X}$ 拉直：$\text{vec}(\mathbf{X}) = [x_{11}, x_{21}, x_{12}, x_{22}]^\top$（注意是按列堆叠）。
构建克罗内克积 $\mathbf{b}^\top \otimes \mathbf{I}_2$：
$$ [b_1, b_2] \otimes \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} = \begin{bmatrix} b_1 \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} & b_2 \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} \end{bmatrix} = \begin{bmatrix} b_1 & 0 & b_2 & 0 \\ 0 & b_1 & 0 & b_2 \end{bmatrix} $$

现在让右边的矩阵乘以拉直的向量：
$$ \begin{bmatrix} b_1 & 0 & b_2 & 0 \\ 0 & b_1 & 0 & b_2 \end{bmatrix} \begin{bmatrix} x_{11} \\ x_{21} \\ x_{12} \\ x_{22} \end{bmatrix} = \begin{bmatrix} b_1 x_{11} + b_2 x_{12} \\ b_1 x_{21} + b_2 x_{22} \end{bmatrix} $$

**见证奇迹**：结果和左边一模一样！

**本质**：这个恒等式的作用就是把**参数矩阵 $\mathbf{X}$** 变成**参数向量 $\text{vec}(\mathbf{X})$**，同时把**数据向量 $\mathbf{b}$** 变成一个**巨大的稀疏数据矩阵**。这样就实现了**参数与数据的线性分离**。

---

### 3. 回到论文推导过程

现在我们把这个工具用到公式 (2.24) 上。

#### 第一步：分离参数
我们要处理的项是第 $j$ 个微分项：
$$ \text{Term} = \frac{p^j}{A_i(p)} \mathbf{B}_{i,j} \mathbf{u}(t_k) $$

这里 $\mathbf{B}_{i,j}$ 是常数矩阵，$\mathbf{u}(t_k)$ 是随时间变化的信号。
利用恒等式，把 $\mathbf{B}_{i,j} \mathbf{u}$ 拆开：
$$ \mathbf{B}_{i,j} \mathbf{u}(t_k) = (\mathbf{u}(t_k)^\top \otimes \mathbf{I}_q) \text{vec}(\mathbf{B}_{i,j}) $$

论文定义了一个符号 $\mathbf{U}(t_k) = \mathbf{u}(t_k) \otimes \mathbf{I}_q$（公式 2.25），注意这里是列向量 $\mathbf{u}$。
根据转置规则 $(\mathbf{A} \otimes \mathbf{B})^\top = \mathbf{A}^\top \otimes \mathbf{B}^\top$，我们可以写成：
$$ \mathbf{u}(t_k)^\top \otimes \mathbf{I}_q = (\mathbf{u}(t_k) \otimes \mathbf{I}_q)^\top = \mathbf{U}(t_k)^\top $$

所以：
$$ \mathbf{B}_{i,j} \mathbf{u}(t_k) = \mathbf{U}(t_k)^\top \text{vec}(\mathbf{B}_{i,j}) $$

#### 第二步：加入滤波器
现在把前面的系数 $\frac{p^j}{A_i(p)}$ 加回来。
注意：$\text{vec}(\mathbf{B}_{i,j})$ 是待辨识的**常数**参数，它不受微分算子 $p$ 的影响（常数的导数是0，但在乘法里可以提出来）。
受 $p$ 影响的是输入信号 $\mathbf{U}(t_k)$。

$$ \text{Term} = \left[ \frac{p^j}{A_i(p)} \mathbf{U}(t_k)^\top \right] \text{vec}(\mathbf{B}_{i,j}) $$

#### 第三步：求导
现在形式已经变成了标准的线性方程：
$$ \mathbf{y} = \mathbf{H} \cdot \mathbf{x} $$
其中：
*   $\mathbf{y}$ 是模型输出项。
*   $\mathbf{x} = \text{vec}(\mathbf{B}_{i,j})$ 是未知参数。
*   $\mathbf{H} = \frac{p^j}{A_i(p)} \mathbf{U}(t_k)^\top$ 是系数矩阵。

线性代数告诉我们，$\frac{\partial (\mathbf{H}\mathbf{x})}{\partial \mathbf{x}} = \mathbf{H}$。

所以，导数就是：
$$ \frac{\partial \text{Term}}{\partial \text{vec}(\mathbf{B}_{i,j})} = \frac{p^j}{A_i(p)} \mathbf{U}(t_k)^\top $$

这就是公式 (2.24) 的最终结果。

---

### 4. 物理直观：这到底是个啥？

别被 $\otimes$ 吓到了，我们来看看最后这个梯度的物理意义。

$\mathbf{U}(t_k)^\top$ 其实就是一个**数据分发器**。
假设系统有 2 个输入，2 个输出。
$$ \mathbf{U}^\top = \begin{bmatrix} u_1 & 0 & u_2 & 0 \\ 0 & u_1 & 0 & u_2 \end{bmatrix} $$

这个梯度矩阵告诉我们：
1.  **参数 $\mathbf{B}_{i,j}$ 的每个元素**（比如第 1 行第 1 列的元素，代表输入 1 对输出 1 的增益），它的变化如何影响最终输出？
2.  答案是：它影响的大小取决于**输入信号 $u_1$** 的大小。
3.  再考虑前面的 $\frac{p^j}{A_i(p)}$：因为这是第 $j$ 阶微分项，所以还要加上微分和滤波。

**总结**：
分子参数 $\mathbf{B}_{i,j}$ 实际上就是**增益**。
在控制理论中，**增益的灵敏度（梯度） = 输入信号**。
因为我们有滤波器 $1/A(p)$ 和微分算子 $p^j$，所以：
**分子参数的灵敏度 = 经过滤波和微分处理后的输入信号**。

这与我们的直观完全一致！复杂的 Kronecker 积推导，只是为了在 MIMO 系统的多维空间中，严谨地把这个直观结论表达出来而已。
