# 参数辨识 vs 非参数辨识

“参数辨识”（Parametric Identification）与“非参数辨识”（Non-parametric Identification）是系统辨识（System Identification）领域的两大核心方法。

简而言之，**参数辨识**是先假设系统符合某种数学公式（结构），然后算出公式里的具体系数（参数）；而**非参数辨识**则是不预设具体的数学公式，直接用曲线或图表（如频率响应曲线）来描述系统的行为。

以下结合您提供的文件《Modal Identification of Multivariable Systems》以及控制理论知识，为您详细解释并对比这两者。

## 1. 参数辨识 (Parametric Identification)

**定义：**
 参数辨识是指在已知或假设模型结构（如传递函数、状态空间方程）的前提下，利用输入输出数据来估计模型中有限个未知参数（如增益、时间常数、阻尼比等）的过程。[journals.sagepub+1](https://journals.sagepub.com/doi/pdf/10.1260/026309206778494274)​

* **核心特征：** 模型结构固定，参数数量有限。

* **数学形式：** 例如，假设系统是如下传递函数：
  \[
  G(s) = \frac{b_0}{s^2 + a_1 s + a_0}
  \]
   参数辨识的任务就是找到 $a_1, a_0, b_0$ 这三个具体的数值。[control.dii.unisi](https://control.dii.unisi.it/sysid/pdf/notes_system_identification.pdf)​

* **文件中的应用：** 您上传的 thesis 正是致力于解决**MIMO系统的参数辨识**问题。文中提到，通过参数辨识得到的“模态模型”（Modal Model），可以将系统描述为物理意义明确的质量、刚度、阻尼矩阵或模态参数（频率 $\omega_i$、阻尼比 $\xi_i$）。

* **常见方法：** 最小二乘法 (Least Squares)、极大似然法 (Maximum Likelihood)、预测误差法 (PEM)、以及文件中提到的 SRIVC (Simplified Refined Instrumental Variable for Continuous-time systems) 算法。

## 2. 非参数辨识 (Non-parametric Identification)

**定义：**
 非参数辨识是指在不假设系统具体数学模型结构的情况下，直接确定系统的动态特性。它通常用函数曲线（如脉冲响应、频率响应）来描述系统，而不涉及具体的方程系数。[sciencedirect+1](https://www.sciencedirect.com/science/article/pii/0005109881900844)​

* **核心特征：** 没有固定的方程形式，模型本质上是无限维的（是一条曲线或一组数据点）。
* **数学形式：** 例如，直接测量出系统的**频率响应函数 (FRF)** $G(j\omega_k)$，它是一组随频率变化的复数数据点，而不是一个紧凑的公式。[journals.sagepub](https://journals.sagepub.com/doi/pdf/10.1260/026309206778494274)
* **文件中的应用：** 文中提到，非参数模型（如 FRF）通常是辨识的第一步。通过正弦扫描或多正弦激励实验，先得到非参数的 FRF 数据（Bode 图），然后以此为基础去拟合参数模型。
* **常见方法：** 谱分析 (Spectral Analysis)、相关分析 (Correlation Analysis)、傅里叶变换 (FFT)、经验传递函数估计 (ETFE)。

------

## 3. 参数辨识 vs. 非参数辨识

| 比较维度       | 参数辨识 (Parametric)                                        | 非参数辨识 (Non-parametric)                            |
| :------------- | :----------------------------------------------------------- | :----------------------------------------------------- |
| **模型结构**   | **固定结构**（需预知阶数、物理结构）                         | **灵活/无结构**（仅依赖数据）                          |
| **输出形式**   | 具体的**方程**与**参数值**（如 $\omega_n = 10$ rad/s）       | **曲线**或**图表**（如 Bode 图、阶跃响应曲线）         |
| **物理意义**   | **强**。参数通常对应物理量（质量、摩擦系数等），如文件中的模态参数。 | **弱**。难以直接对应具体物理部件，仅反映输入输出关系。 |
| **先验知识**   | 需要。必须先知道系统大致是几阶的、线性的还是非线性的。       | 不需要。不需要预先知道系统的物理机理。                 |
| **计算复杂度** | **高**。通常涉及非线性优化、迭代求解（如 SRIVC 算法）。      | **低**。通常基于统计或变换（如 FFT），计算快。         |
| **精度与效率** | 如果模型结构假设正确，**精度极高**且数据利用率高（方差小）。 | 对噪声敏感，为了平滑曲线通常需要大量数据平均。         |
| **主要用途**   | 控制器设计（如 PID 参数整定）、故障诊断、仿真预测。          | 初步分析、验证参数模型、直接频域分析。                 |

## 4. 两者的关系与实际工作流

在实际工程（如您文件中的光刻机运动台辨识）中，这两者通常是**互补**且**串行**使用的，而不是对立的。

**典型工作流：**

1. **非参数辨识（第一步）：** 对系统施加激励（如 chirp 信号），通过 FFT 分析得到**非参数的 FRF 曲线**。这能让你直观地看到系统有哪些共振峰，大约在什么频率。
2. **模型结构确定：** 观察 FRF 曲线，发现有 3 个共振峰，于是决定采用一个 6 阶（3个模态）的参数模型结构。
3. **参数辨识（第二步）：** 使用算法（如文件中提出的频域 SRIVC 方法），让参数模型的曲线去逼近第一步得到的非参数 FRF 数据，从而解算出具体的 $M, C, K$ 矩阵或模态参数。1486187-ModalIdentification_of_MultivariableSystems.pdf

**总结：** 非参数辨识是“看清现象”（画出曲线），参数辨识是“挖掘本质”（算出方程）。您上传的文件核心工作在于**如何更精确、更符合物理意义地完成第二步（参数辨识）**，特别是针对复杂的多输入多输出（MIMO）机械系统。

# 连续时间辨识 vs 离散时间辨识

根据您提供的文件《Modal Identification of Multivariable Systems》以及控制理论的核心概念，连续时间辨识（Continuous-Time Identification）与离散时间辨识（Discrete-Time Identification）的主要区别在于**底层数学模型的形式**以及**物理意义的保留程度**。

简而言之：**连续时间辨识**直接恢复具有物理意义的微分方程参数（如质量、阻尼），适合物理建模；**离散时间辨识**通过差分方程逼近系统，适合计算机直接控制和预测。

------

## 1. 核心定义

### 什么是连续时间辨识 (Continuous-Time Identification)?

连续时间辨识是指利用采样数据（输入 $u(t_k)$ 和输出 $y(t_k)$），直接估计**连续时间微分方程**模型参数的过程。

* **数学模型**：基于微分算子 $p = \frac{d}{dt}$ 或拉普拉斯算子 $s$。
  * 形式示例：$a_n \frac{d^n y}{dt^n} + \dots + a_0 y(t) = b_m \frac{d^m u}{dt^m} + \dots + b_0 u(t)$。
* **关键技术**：由于计算机无法直接处理“微分”，通常使用**状态变量滤波器 (SVF)** 或积分方法来重构信号的导数，从而避免直接微分带来的噪声放大问题。
* **物理意义**：参数直接对应物理系统的特性（如机械系统的刚度 $k$、阻尼 $c$、质量 $m$）。

### 什么是离散时间辨识 (Discrete-Time Identification)?

离散时间辨识是指利用采样数据，估计**离散时间差分方程**模型参数的过程。这是目前最主流的工程辨识方法。

* **数学模型**：基于移位算子 $q$ ($q x(k) = x(k+1)$) 或 Z 变换算子 $z$。
  * 形式示例：$y(k) + a_1 y(k-1) + \dots = b_1 u(k-1) + \dots$。
* **关键技术**：最小二乘法 (Least Squares)、极大似然法等。由于数据本身就是离散的，这类方法在数值计算上非常自然且稳定。
* **物理意义**：参数是采样时间 $T_s$ 的函数，**没有直接的物理直观意义**。例如，离散极点 $z = e^{sT_s}$，这使得参数难以直接解读为物理量。

------

## 2. 详细对比表

| 特性           | 连续时间辨识 (CT Identification)                             | 离散时间辨识 (DT Identification)                             |
| :------------- | :----------------------------------------------------------- | :----------------------------------------------------------- |
| **核心模型**   | 微分方程 ($dx/dt$)                                           | 差分方程 ($x_{k+1} - x_k$)                                   |
| **数学算子**   | $s$ (Laplace) 或 $p$ (微分算子)                              | $z$ (Z-transform) 或 $q$ (移位算子)                          |
| **物理意义**   | **强**。参数直接对应物理量（频率、阻尼、时间常数）。         | **弱**。参数是物理参数与采样时间的复杂组合。                 |
| **采样率依赖** | **不依赖**。模型参数不仅与采样率无关，而且在高采样率下更稳定。 | **强依赖**。采样率改变，模型参数全变。                       |
| **主要难点**   | 需要处理“导数”计算（使用SVF滤波器）；随机噪声在连续域难以建模（通常采用混合模型）。 | 极高采样率下会出现“病态问题”（极点聚集在 $z=1$ 附近），导致数值不稳定。 |
| **辨识方法**   | 直接法 (SRIVC) 或 间接法 (先辨识离散再转换)。                | 预测误差法 (PEM)、子空间法 (N4SID)。                         |
| **适用场景**   | **物理参数估计**（如模态分析、灰箱建模）、非均匀采样数据。   | **数字控制**、纯预测模型、计算机实现。                       |

------

## 3. 深度解析（基于您的附件内容）

在您的文件背景（**多变量模态辨识**）下，连续时间辨识具有特殊的优势和处理方式：

### A. 为什么在模态分析中偏向连续时间？

1. **直接获取物理参数**：模态分析的核心是获取固有频率 ($\omega_n$) 和阻尼比 ($\zeta$)。
   * **CT 方法**：直接估计出的系数矩阵可以直接分解为 $\omega_n$ 和 $\zeta$。
   * **DT 方法**：必须先估计离散参数，再通过 $s = \frac{1}{T_s} \ln(z)$ 转换。这种转换在快速采样或高频模态下容易产生偏差（Bias）和数值误差。
2. **处理刚性系统 (Stiff Systems)**：对于同时包含快慢动态的系统，离散时间模型往往需要极高的采样率，导致离散极点全部挤在单位圆 $1$ 附近，计算机难以区分。连续时间模型不存在此数值敏感性问题。

### B. "直接法" (Direct Approach) vs "间接法" (Indirect Approach)

您的文件特别强调了**直接连续时间辨识 (Direct CT Identification)**：

* **间接法**：先辨识离散模型 $\rightarrow$ 再转换为连续模型。
  * *缺点*：如果离散模型估计有微小误差，转换后的连续参数可能会有巨大偏差（特别是在高频段）。
* **直接法**（推荐）：直接构建连续时间模型的误差函数。
  * *核心技术*：使用 **SRIVC (Simplified Refined Instrumental Variable Continuous-time)** 算法。
  * *混合策略*：文件中提到一种**混合模型 (Hybrid Model)** 策略——**“对象是连续的，噪声是离散的”**。这是因为在连续时间下对“白噪声”建模在数学上非常困难（功率无穷大），所以用离散的 ARMA 模型来描述噪声，用连续的微分方程描述系统动态。

### C. 状态变量滤波器 (SVF) 的作用

在连续时间辨识中，不能直接对采样数据求导（$dy/dt$）来获得回归向量，因为这会极度放大高频噪声。

* **解决方法**：在方程两边同时通过一个低通滤波器 $F(s)$（即 SVF）。
* **效果**：$s^n Y(s)$ 变成了 $\frac{s^n}{F(s)} Y(s)$。这是一个带通滤波操作，既提取了动态信息，又压制了高频噪声，使得构建线性回归方程成为可能。

## 总结建议

* 如果您做的是**纯数字控制**（如设计一个运行在单片机上的 PID），**离散时间辨识**更方便。
* 如果您做的是**设备健康监测、模态分析、或者需要理解系统的物理本质**（如您的附件所示），**连续时间辨识**是更优的选择，特别是采用 SRIVC 这类能直接处理采样数据的直接法。

1. https://ieeexplore.ieee.org/document/8411997/
2. https://www.nature.com/articles/s41598-021-85320-4
3. https://www.sciencedirect.com/science/article/pii/S0165188922002263
4. https://ntrs.nasa.gov/citations/19780034746
5. https://forum.allaboutcircuits.com/threads/relationship-between-continuous-and-discrete-time-models-of-a-lti-system.158497/
6. https://portal.research.lu.se/files/62894316/8146063.pdf
7. https://www.sciencedirect.com/topics/mathematics/discrete-time-system
8. https://www.diva-portal.org/smash/get/diva2:1429291/FULLTEXT01.pdf
9. https://skoge.folk.ntnu.no/prost/proceedings/acc04/Papers/0114_WeM01.3.pdf
10. https://www.mathworks.com/matlabcentral/answers/256722-continuous-time-vs-discrete-time-identification
11. https://arxiv.org/abs/2308.11933
12. https://fiveable.me/adaptive-and-self-tuning-control/unit-10/discrete-time-system-models-identification/study-guide/zsmC3UQmWySV3dWd
13. https://math.berkeley.edu/~chorin/LLC16.pdf
14. http://publications.pvandenhof.nl/Paperfiles/Dankers&etal_CDC2014.pdf
15. https://www.reddit.com/r/ControlTheory/comments/x2whey/discrete_controller_vs_continuous_time_controller/
16. https://www.sciencedirect.com/science/article/pii/S1474667017429703
17. https://engineering.purdue.edu/VISE/ee438L/lab2/pdf/lab2.pdf
18. https://jblevins.org/research/ctgames.pdf
19. https://www.sciencedirect.com/science/article/pii/S1474667017479328
20. https://arxiv.org/html/2511.02701v1
21. https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/147441806/05f43905-bd31-4304-be5d-a63140827e87/1486187-ModalIdentification_of_MultivariableSystems.pdf

# 差分方程和移位算子

为了更深入地理解离散时间辨识，我们需要从数学本质上弄清楚其描述工具：**差分方程**和**移位算子**。它们就像是连续时间系统中的“微分方程”与“拉普拉斯算子 $s$”的关系。

------

## 1. 什么是差分方程 (Difference Equation)?

**差分方程**是描述离散时间系统输入输出动态关系的代数方程。它是连续时间系统中“微分方程”的离散化对应物。

* **定义**：一个方程，它将系统在当前时刻的输出 $y(k)$ 与过去的输出 $y(k-1), \dots$ 以及当前和过去的输入 $u(k), u(k-1), \dots$ 联系起来。
* **物理直观**：它描述了系统如何“一步一步”地演变。
  * 例如：$y(k) = 0.9 y(k-1) + 0.1 u(k-1)$。
  * **解读**：当前的输出是上一步输出的 90% 加上上一步输入的 10%。这就是一个典型的一阶惯性环节的离散化形式。
* **一般形式**（$n$ 阶线性常系数差分方程）：
   $y(k) + a_1 y(k-1) + \dots + a_n y(k-n) = b_0 u(k) + \dots + b_m u(k-m)$
* **在辨识中的作用**：系统辨识的核心任务就是利用观测数据 $u(k), y(k)$，解算出上述方程中的系数 $a_i$ 和 $b_i$。

------

## 2. 什么是移位算子 (Shift Operator)?

**移位算子**是为了简化差分方程的书写和推导而引入的一种**符号运算工具**。它将复杂的递推关系转化为简单的代数多项式运算。

### A. 前向移位算子 $q$ (Forward Shift Operator)

这是现代控制理论中最常用的定义。

* **定义**：$q$ 作用在信号上，表示“**时间向前推进一步**”。
   $q y(k) = y(k+1)$
   $q^{-1} y(k) = y(k-1) \quad (\text{后向移位 / 延迟})$
* **作用**：它把差分方程变成了**传递函数**形式。
  * 差分方程：$y(k+1) + a y(k) = b u(k)$
  * 算子形式：$(q + a) y(k) = b u(k)$
  * 传递函数：$G(q) = \frac{y(k)}{u(k)} = \frac{b}{q + a}$

### B. Z 变换算子 $z$

在频域分析中，我们通常把 $q$ 替换为复变量 $z$。

* **本质**：$z$ 是 $q$ 的频域对等物。
* **关系**：对于线性时不变系统，算子 $q$ 的运算规则与复变量 $z$ 的代数运算规则完全一致。因此工程上常混用 $G(q)$ 和 $G(z)$。

------

## 3. 核心对比：差分方程 vs 移位算子

| 维度         | 差分方程 (Difference Equation)       | 移位算子 (Shift Operator $q$)              |
| :----------- | :----------------------------------- | :----------------------------------------- |
| **本质**     | **时域**的递推公式。                 | **代数域**（算子域）的符号工具。           |
| **表现形式** | $y(k) = -a_1 y(k-1) + b_1 u(k-1)$    | $A(q) y(k) = B(q) u(k)$                    |
| **计算特点** | 适合计算机**编程实现**（迭代循环）。 | 适合数学**推导和分析**（极点、零点分析）。 |
| **辨识视角** | 辨识就是做“线性回归”求系数。         | 辨识就是拟合“有理分式”传递函数。           |
| **物理对应** | 对应微分方程 $\frac{dy}{dt}$。       | 对应拉普拉斯算子 $s$ 或微分算子 $p$。      |

------

## 4. 进阶：为什么模态辨识害怕 $q$ 算子？

在您的模态辨识场景中，**移位算子 $q$ 有一个致命弱点**：

* **极点聚集问题**：
   当采样很快（$T_s \to 0$）时，所有稳定的物理极点 $s_i$ 都会映射到 $z$ 平面上的 $1$ 附近（因为 $z = e^{s T_s} \approx 1 + s T_s$）。
  * 例如：$0.99991$ 和 $0.99992$ 可能分别对应完全不同的物理阻尼。
  * **差分方程表现**：$y(k) \approx y(k-1)$。这意味着当前的输出几乎等于上一步的输出，差异极小，被噪声淹没。
  * **结果**：使用 $q$ 算子（或标准差分方程）进行辨识时，计算机因为精度问题（Numerical Ill-conditioning）无法算出准确的参数。
* **解决办法：$\delta$ 算子 (Delta Operator)**
   为了解决这个问题，通常引入 $\delta$ 算子：
   $\delta = \frac{q - 1}{T_s}$
   这个算子在快速采样时趋近于微分算子 $s$（$\lim_{T_s \to 0} \delta = \frac{d}{dt}$），从而恢复了数值稳定性。这也是为什么您的附件推荐直接使用连续时间方法（类似 $s$ 或 $\delta$）而不是直接辨识离散 $q$ 模型的原因。

1. https://blog.csdn.net/weixin_44378835/article/details/109705137
2. https://web.xidian.edu.cn/junli/files/20140304_203333.pdf
3. https://word.baidu.com/view/797a7ce7ce7931b765ce0508763231126edb77a8.html
4. [https://raw.githubusercontent.com/swjtuhub/SWJTU-Courses/main/%E8%AF%BE%E7%A8%8B%E8%B5%84%E6%96%99/%E4%BF%A1%E5%8F%B7%E4%B8%8E%E7%B3%BB%E7%BB%9F%EF%BC%88%E4%BF%A1%E6%81%AF%EF%BC%89/%E8%AF%BE%E4%BB%B6/%E7%AC%AC5%E7%AB%A0%20%E7%A6%BB%E6%95%A3%E6%97%B6%E9%97%B4%E7%B3%BB%E7%BB%9F%E7%9A%84%E5%88%86%E6%9E%90-1.pdf](https://raw.githubusercontent.com/swjtuhub/SWJTU-Courses/main/课程资料/信号与系统（信息）/课件/第5章 离散时间系统的分析-1.pdf)
5. https://www.jove.com/cn/science-education/v/16046/classification-of-systems-ii
6. https://blog.csdn.net/qrsysterm/article/details/115149763
7. http://zhangchuheng123.github.io/assets/files/2016-10-17-Signal-and-Systems-6.pdf
8. https://www.cnblogs.com/TaigaCon/p/8094422.html
9. [http://www.hcii-lab.net/lianwen/Course/Digital%20Signal%20Processing/main/notes/Lecture%203%20%E7%A6%BB%E6%95%A3%E6%97%B6%E9%97%B4%E7%B3%BB%E7%BB%9F%E7%9A%84%E6%97%B6%E5%9F%9F%E5%88%86%E6%9E%90.pdf](http://www.hcii-lab.net/lianwen/Course/Digital Signal Processing/main/notes/Lecture 3 离散时间系统的时域分析.pdf)

# 连续时间 s 平面到离散时间 z 平面



# 极点聚集问题：离散时间系统辨识快速采样

在**快速采样**（Sampling time $T_s \to 0$）的离散时间系统辨识中，最臭名昭著的现象就是**极点聚集问题**。这是数学变换在数值计算上产生的“副作用”，是离散时间模型（特别是使用 $q$ 算子）在物理建模中的致命弱点。

------

## 1. 为什么会发生？数学根源

这个问题的根源在于**连续时间 $s$ 平面**到**离散时间 $z$ 平面**的映射关系：
 $z = e^{s T_s}$

* **连续极点 $s_i$**：物理系统的极点，决定了系统的固有频率和阻尼。对于稳定系统，$s_i$ 分布在左半平面（实部 $<0$）。
* **离散极点 $z_i$**：数字系统的极点，分布在单位圆内（模 $<1$）。

### 泰勒展开看本质

当采样很快（$T_s$ 很小）时，我们可以对指数函数进行一阶近似：
 $z = e^{s T_s} \approx 1 + s T_s$
 这意味着：
 $z - 1 \approx s T_s$

由于 $T_s$ 是一个非常小的数（例如 $10^{-4}$ 秒），任何有限的物理极点 $s$（例如 $-10, -20, -100$），乘以 $T_s$ 后都会变成一个**极其接近 0 的微小负数**。

* $s = -10 \implies z \approx 1 - 0.001 = 0.999$
* $s = -20 \implies z \approx 1 - 0.002 = 0.998$

**结果**：原本在 $s$ 平面上相距甚远的物理极点（$-10$ 和 $-20$），在 $z$ 平面上全部被“压缩”并“聚集”到了 $(1, 0)$ 这个点附近。

------

## 2. 为什么是“致命”弱点？数值计算灾难

极点聚集在 $z=1$ 附近，对计算机辨识参数意味着什么？

### A. 差分方程的退化

对于一阶系统 $y(k+1) + a y(k) = b u(k)$，极点是 $-a$。如果极点在 $0.9999$：
 $y(k+1) \approx 0.9999 y(k) + \dots$
 这说明 $y(k+1)$ 和 $y(k)$ 几乎相等。
 在计算机看来：
 $y(k+1) - y(k) \approx \text{噪声}$
 计算机试图从这些被噪声淹没的微小差值中恢复出系统参数，这在数值上是**病态的 (Ill-conditioned)**。

### B. 精度丢失 (Loss of Significance)

计算机浮点数（double precision）只有约 15 位有效数字。

* 假设真实极点是 $z_1 = 0.99991$ 和 $z_2 = 0.99992$。
* 在辨识过程中，算法需要计算特征多项式 $(z - z_1)(z - z_2)$。
* 系数会变成 $z^2 - 1.99983 z + 0.9998300072$。
* 注意这些系数小数点后好多位都一样！计算机在矩阵求逆运算中，微小的舍入误差（Round-off error）就会导致解出的 $z$ 跑到 $0.99$ 或者 $1.01$ 去。
* **物理后果**：你算出的阻尼比 $\zeta$ 可能是负的（不稳定），或者误差高达 200%。

------

## 3. 您的模态辨识场景中的具体表现

在模态分析中，系统往往具有**宽频响**特性（既有低频刚体模态，又有高频柔性模态）。

1. **为了捕捉高频模态** $\to$ 你必须使用**高采样率**（$T_s$ 很小）。
2. **导致低频模态极点聚集** $\to$ 低频模态（通常也是最重要的主模态）的 $s_i$ 比较小，乘以微小的 $T_s$ 后，$z_i$ 极度接近 1。
3. **辨识失效** $\to$ 使用 $q$ 算子的算法（如离散 ARX/ARMAX）算出的低频模态参数极其不准，经常出现虚假的“负阻尼”或“重根”。

这就是为什么您的附件强烈推荐使用 **$\delta$ 算子** 或 **直接连续时间 (Continuous-Time) 辨识方法 (SRIVC)**。

* **$\delta$ 算子** 定义为 $\delta = \frac{z-1}{T_s}$。
* 当 $z \to 1$ 时，分子趋近于 0，但分母 $T_s$ 也趋近于 0，两者相除刚好把那个微小的差异“放大”回了物理极点 $s$ 的数量级。
* 这在数值上重新拉开了极点间距，解决了病态问题。

## 总结

**极点聚集**是因为 $z=e^{sT_s}$ 映射把广阔的左半平面极点压缩到了 $z=1$ 的微小邻域内。这导致差分方程系数对物理参数变化极不敏感，计算机无法在噪声存在的情况下精确求解。这是高频采样下离散辨识算法（基于 $q$ 算子）的固有缺陷。

1. https://jeit.ac.cn/en/article/doi/10.11999/JEIT240528
2. https://blog.csdn.net/m0_56997192/article/details/146431928
3. [http://www.stat.ucla.edu/~sczhu/Books/Book_4_MCMC_%E4%B8%AD%E6%96%87%E7%89%88.pdf](http://www.stat.ucla.edu/~sczhu/Books/Book_4_MCMC_中文版.pdf)
4. http://amt.amss.cas.cn/kyjz/kycg/
5. https://wulixb.iphy.ac.cn/article/2007/4
6. https://patents.google.com/patent/CN114127843B/zh
7. https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.1770-0-200607-S!!MSW-C.doc

# 第二章 加性-连续时间系统-时域辨识



# 第三章 加性-连续时间系统-频域辨识