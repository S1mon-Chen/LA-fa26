# Assignment 8 学习笔记：分块消元、可逆性与乘积的秩

课程：MATH1205 线性代数  
原题：[assignment8.pdf](assignment8.pdf)  
提交版：[中文解答 PDF](../../output/pdf/assignment8-sol.pdf)  
源文件：[中文 LaTeX](assignment8-sol.tex)  
相关讲义：[§2.6 分块矩阵的初等变换](../../lectures/LA2.6.pdf)

本次范围：**51(1)(2)、68、69、70**。以下用 $E_k$ 表示 $k$ 阶单位矩阵，用 $r(A)$ 表示秩。$A^*$ 是伴随矩阵，不是共轭转置。

## 1. 四道题的一条主线

这些题都在利用同一个想法：通过合法的初等变换，把一个不容易直接判断的矩阵改写成分块三角形或等价标准形。

| 题号 | 要解决的问题 | 关键动作 |
| --- | --- | --- |
| 51 | 分块矩阵 $Q$ 何时可逆 | 消去左下块，留下一个标量 $b-\alpha^TA^{-1}\alpha$ |
| 68 | 从 $E-AB$ 的逆构造 $E-BA$ 的逆 | 使用 $(E-BA)B=B(E-AB)$，保持乘法次序 |
| 69 | 比较 $\lambda E-AB$ 与 $\lambda E-BA$ 的行列式 | 将它们放进同一个大矩阵的两种分块三角形中 |
| 70 | 控制乘积的秩 | 上界用消元；下界把 $A$ 化成标准形，数最多损失多少行 |

## 2. 分块运算先检查什么

### 2.1 尺寸对齐，乘法顺序固定

在分块乘积中，仍然按“左边的块行乘右边的块列”计算。例如

$$
\begin{pmatrix}A&B\\C&D\end{pmatrix}
\begin{pmatrix}X&Y\\Z&W\end{pmatrix}
=
\begin{pmatrix}AX+BZ&AY+BW\\CX+DZ&CY+DW\end{pmatrix}.
$$

每一项乘法都必须有意义，同一位置相加的两项还必须有相同尺寸。块通常是矩阵，不能把 $CX$ 随意改成 $XC$。

### 2.2 “左行右列”是什么意思

- 块行变换：$R_i\leftarrow R_i+KR_j$，系数矩阵 $K$ 在左边。
- 块列变换：$C_i\leftarrow C_i+C_jK$，系数矩阵 $K$ 在右边。

例如，对第 69 题的矩阵做 $C_1\leftarrow C_1+C_2B$，左上块增加的是 $AB$；再做 $R_2\leftarrow R_2-BR_1$，右下块减去的是 $BA$。这两个不同的乘积正是证明的核心。

第三类分块初等变换保持秩；对方阵还保持行列式。它们可分解成普通第三类初等变换，也可以看成乘上一个对角块为单位矩阵、行列式为 $1$ 的分块三角矩阵。

不同块行或块列的尺寸也可以不同，只要系数的尺寸合适。例如第 70 题中 $C_1$ 有 $n$ 列、$C_2$ 有 $s$ 列，$B$ 是 $n\times s$，所以 $C_1B$ 恰好能加到 $C_2$ 上。

### 2.3 不要套用错误的“二阶行列式公式”

当对角块是方阵时，分块三角矩阵满足

$$
\det\begin{pmatrix}A&C\\O&D\end{pmatrix}=\det A\det D.
$$

但一般的分块矩阵不能直接写成
$\det\begin{pmatrix}A&B\\C&D\end{pmatrix}=\det(AD-BC)$。
各块尺寸可能不同，矩阵也通常不可交换。先消去一个非对角块，再使用三角矩阵公式，才是这里可靠的做法。

## 3. 第 51 题：可逆性归结为一个标量

### 3.1 题目与尺寸

$$
P=\begin{pmatrix}E_n&O\\-\alpha^TA^*&|A|\end{pmatrix},
\qquad
Q=\begin{pmatrix}A&\alpha\\\alpha^T&b\end{pmatrix},
$$

其中 $A$ 为 $n$ 阶可逆矩阵，$\alpha$ 是 $n\times1$ 列矩阵，$b$ 是标量。因此 $P,Q$ 都是 $(n+1)$ 阶矩阵。

特别注意：$\alpha^TA^{-1}\alpha$ 的尺寸是

$$
(1\times n)(n\times n)(n\times1)=1\times1,
$$

所以它是一个数，可以与 $b$ 比较。

### 3.2 为什么题目这样设计 $P$

伴随矩阵满足

$$
A^*A=AA^*=|A|E_n.
$$

因此乘积的左下块恰好抵消：

$$
-\alpha^TA^*A+|A|\alpha^T
=-|A|\alpha^T+|A|\alpha^T=O.
$$

完整结果为

$$
\boxed{
PQ=\begin{pmatrix}
A&\alpha\\
O&b|A|-\alpha^TA^*\alpha
\end{pmatrix}
=\begin{pmatrix}
A&\alpha\\
O&|A|\bigl(b-\alpha^TA^{-1}\alpha\bigr)
\end{pmatrix}.}
$$

第二个等号使用了 $A$ 可逆时的公式 $A^*=|A|A^{-1}$。

### 3.3 从乘积得到行列式与充要条件

因为 $P$ 为分块下三角矩阵，

$$
|P|=|A|\ne0.
$$

而 $PQ$ 为分块上三角矩阵，故

$$
|P||Q|=|PQ|=|A|^2\bigl(b-\alpha^TA^{-1}\alpha\bigr).
$$

约去非零的 $|A|$，得到

$$
\boxed{|Q|=|A|\bigl(b-\alpha^TA^{-1}\alpha\bigr)}.
$$

于是

$$
\boxed{Q\text{ 可逆}\iff b-\alpha^TA^{-1}\alpha\ne0.}
$$

这里必须用到 $A$ 可逆。没有这一条件，既不能写 $A^{-1}$，也不能约去 $|A|$。

### 3.4 消元以后留下的数是什么

直接对 $Q$ 做块行变换也能看出结构：

$$
\begin{pmatrix}E_n&O\\-\alpha^TA^{-1}&1\end{pmatrix}
\begin{pmatrix}A&\alpha\\\alpha^T&b\end{pmatrix}
=\begin{pmatrix}A&\alpha\\O&b-\alpha^TA^{-1}\alpha\end{pmatrix}.
$$

右下角的 $b-\alpha^TA^{-1}\alpha$ 称为关于 $A$ 的 **Schur 补**。不必记术语也能理解它：用前 $n$ 行消去最后一行左侧的部分后，最后一个主元就是这个数。

当它为零时，可以直接给出 $Q$ 的非零零向量：

$$
Q\begin{pmatrix}-A^{-1}\alpha\\1\end{pmatrix}
=\begin{pmatrix}0\\b-\alpha^TA^{-1}\alpha\end{pmatrix}
=0.
$$

这个向量最后一项为 $1$，所以确实非零。注意题目没有假设 $A$ 对称或正定，不能擅自断言 $\alpha^TA^{-1}\alpha\ge0$。

**快速检验：** 当 $n=1$，$A=(a)$、$\alpha=(t)$、$a\ne0$ 时，公式退化为

$$
\det\begin{pmatrix}a&t\\t&b\end{pmatrix}
=ab-t^2=a\left(b-\frac{t^2}{a}\right),
$$

正好是熟悉的二阶行列式。

## 4. 第 68 题：逆矩阵怎样从 $AB$ 转到 $BA$

已知 $E-AB$ 可逆。设

$$
C=(E-AB)^{-1}.
$$

关键恒等式是

$$
(E-BA)B=B(E-AB).
$$

两边展开都是 $B-BAB$。这个等式只用了结合律，没有使用交换律。

尝试构造 $D=E+BCA$，则

$$
\begin{aligned}
(E-BA)D
&=E-BA+(E-BA)BCA\\
&=E-BA+B(E-AB)CA\\
&=E-BA+BA=E.
\end{aligned}
$$

因为这里是同阶方阵，乘积为单位矩阵便说明两者都可逆。例如取行列式可得 $|E-BA||D|=1$。因此

$$
\boxed{(E-BA)^{-1}=E+B(E-AB)^{-1}A.}
$$

也可以反向验证。由 $A(E-BA)=(E-AB)A$，

$$
D(E-BA)=E-BA+BC(E-AB)A=E.
$$

### 为什么能想到这个形式

如果逆矩阵写成 $E+X$，乘开以后要消掉原来的 $-BA$。于是让 $X$ 的两端保持为 $B$ 和 $A$，中间放入已知的 $(E-AB)^{-1}$，就有机会利用上面的恒等式抵消。

顺序是 **$B$、$(E-AB)^{-1}$、$A$**，不能随意换成 $A$、$(E-AB)^{-1}$、$B$。

第 69 题令 $\lambda=1$ 也能得到更短的证明：

$$
|E-BA|=|E-AB|\ne0.
$$

提交版使用构造逆矩阵的方法，使本题的证明独立于后一题，同时得到一个可复用的公式。

## 5. 第 69 题：同一个大行列式，消元出两个答案

目标是证明

$$
|\lambda E-AB|=|\lambda E-BA|,\qquad \lambda\ne0.
$$

不需要把两个 $n$ 阶行列式展开。构造

$$
M_0=\begin{pmatrix}\lambda E-AB&A\\O&\lambda E\end{pmatrix}.
$$

第一步做块列变换 $C_1\leftarrow C_1+C_2B$：

$$
M_0\longrightarrow
M_1=\begin{pmatrix}\lambda E&A\\\lambda B&\lambda E\end{pmatrix}.
$$

- 左上块：$\lambda E-AB+AB=\lambda E$。
- 左下块：$O+(\lambda E)B=\lambda B$。
- 第二块列保持不变。

第二步做块行变换 $R_2\leftarrow R_2-BR_1$：

$$
M_1\longrightarrow
M_2=\begin{pmatrix}\lambda E&A\\O&\lambda E-BA\end{pmatrix}.
$$

- 左下块：$\lambda B-B(\lambda E)=O$。
- 右下块：$\lambda E-BA$。

两步都保持行列式，且首尾都是分块上三角矩阵，所以

$$
\lambda^n|\lambda E-AB|
=\lambda^n|\lambda E-BA|.
$$

由于 $\lambda\ne0$，约去 $\lambda^n$，即得

$$
\boxed{|\lambda E-AB|=|\lambda E-BA|.}
$$

### 5.1 $\lambda\ne0$ 用在了哪里

变换本身没有除以 $\lambda$；只有最后约去 $\lambda^n$ 时需要它非零。

事实上，$\lambda=0$ 时结论也成立，因为

$$
|-AB|=(-1)^n|A||B|=(-1)^n|B||A|=|-BA|.
$$

这只是补充理解，原题要求证明 $\lambda\ne0$ 的情形。

### 5.2 不能因此说 $AB=BA$

例如

$$
A=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
B=\begin{pmatrix}0&1\\0&0\end{pmatrix},
$$

则 $AB=B\ne0$，而 $BA=O$，但两者都满足

$$
|\lambda E-AB|=|\lambda E-BA|=\lambda^2.
$$

这个例子还说明 $r(AB)$ 与 $r(BA)$ 不一定相等：这里分别是 $1$ 和 $0$。
若以后学习特征多项式，本题说明同阶方阵的 $AB$ 与 $BA$ 有相同特征多项式；它并不说明两个矩阵相等或一定相似。

## 6. 第 70 题：乘积的秩最多、最少是多少

设 $A$ 为 $m\times n$，$B$ 为 $n\times s$，则 $AB$ 为 $m\times s$。目标是

$$
\boxed{r(A)+r(B)-n\le r(AB)\le\min\{r(A),r(B)\}.}
$$

下界中减去的是 **内侧的公共尺寸 $n$**，不是 $m$ 或 $s$。这个下界称为 Sylvester 秩不等式。

### 6.1 上界：把乘积放进一个更大的矩阵

先横向拼接 $A$ 与 $AB$，再做块列变换：

$$
\begin{pmatrix}A&AB\end{pmatrix}
\xrightarrow{C_2\leftarrow C_2-C_1B}
\begin{pmatrix}A&O\end{pmatrix}.
$$

初等变换保持秩，而 $AB$ 是原来拼接矩阵的子矩阵，所以

$$
r(AB)\le r\begin{pmatrix}A&AB\end{pmatrix}=r(A).
$$

再纵向拼接 $B$ 与 $AB$，做块行变换：

$$
\begin{pmatrix}B\\AB\end{pmatrix}
\xrightarrow{R_2\leftarrow R_2-AR_1}
\begin{pmatrix}B\\O\end{pmatrix}.
$$

同理

$$
r(AB)\le r\begin{pmatrix}B\\AB\end{pmatrix}=r(B).
$$

于是得到上界 $r(AB)\le\min\{r(A),r(B)\}$。

子矩阵的秩为什么不超过原矩阵？子矩阵中的每个非零子式，也是原矩阵中的一个非零子式。

### 6.2 下界：先把 $A$ 化成最容易相乘的形状

记 $r=r(A)$。用初等行、列变换将 $A$ 化成等价标准形：存在可逆矩阵 $U,V$，分别为 $m$ 阶、$n$ 阶，使

$$
UAV=\begin{pmatrix}E_r&O\\O&O\end{pmatrix}.
$$

这里的行分块为 $r+(m-r)$，列分块为 $r+(n-r)$，所以它仍是 $m\times n$ 矩阵。

**关键补偿：** 右乘 $V$ 改变了 $A$，所以必须把 $B$ 换成 $V^{-1}B$，才能保持乘积关系：

$$
UAB=(UAV)(V^{-1}B).
$$

不能只写成 $(UAV)B$，因为那通常不等于 $UAB$。

按行分块

$$
V^{-1}B=\begin{pmatrix}B_1\\B_2\end{pmatrix},
\qquad B_1\in\mathbb F^{r\times s},\quad B_2\in\mathbb F^{(n-r)\times s}.
$$

相乘得到

$$
UAB=
\begin{pmatrix}E_r&O\\O&O\end{pmatrix}
\begin{pmatrix}B_1\\B_2\end{pmatrix}
=\begin{pmatrix}B_1\\O\end{pmatrix}.
$$

左乘可逆矩阵不改变秩，因此

$$
r(AB)=r(B_1).
$$

### 6.3 数出最多可能损失的秩

只对 $B_1$ 内部的行做初等行变换，把它化为阶梯形，设其秩为 $k=r(B_1)$。它只有 $k$ 个非零行；下面的 $B_2$ 一共只有 $n-r$ 行，因此整个矩阵最多有 $k+n-r$ 个非零行。

任意矩阵的秩都不超过非零行数，所以

$$
\begin{aligned}
r(B)
&=r(V^{-1}B)\\
&\le r(B_1)+(n-r)\\
&=r(AB)+n-r(A).
\end{aligned}
$$

移项即得下界

$$
\boxed{r(AB)\ge r(A)+r(B)-n.}
$$

这里没有声称“整个矩阵已经是阶梯形”，也没有把它的秩等同于非零行数；只是在用非零行数给秩一个上界。

**直观解释：** 标准形中的 $A$ 保留 $B_1$，消去 $B_2$。被消去的部分只有 $n-r(A)$ 行，因此秩至多损失 $n-r(A)$。

### 6.4 边界情况与两个实用推论

如果 $r(A)+r(B)-n<0$，下界仍然正确，只是不如 $r(AB)\ge0$ 有信息。可以合写为

$$
\max\{0,r(A)+r(B)-n\}\le r(AB)\le\min\{r(A),r(B)\}.
$$

当 $r=0$，$B_1$ 没有行，$A=O$，$AB=O$；当 $r=n$，$B_2$ 没有行。证明仍成立。

两个直接推论：

1. 如果 $A$ 满列秩，即 $r(A)=n$，上下界相等，得到 $r(AB)=r(B)$。
2. 如果 $B$ 满行秩，即 $r(B)=n$，同理得到 $r(AB)=r(A)$。

这里不要求 $A$ 或 $B$ 都是方阵。“左乘可逆矩阵保持秩”是第一个推论在方阵情形下的特例。

**下界可以取到：**

$$
A=\operatorname{diag}(1,1,0),\qquad
B=\operatorname{diag}(0,1,1).
$$

这时 $r(A)=r(B)=2$、$n=3$，且 $AB=\operatorname{diag}(0,1,0)$，所以

$$
r(AB)=1=r(A)+r(B)-n.
$$

**上界也可以取到：** 若 $A=E_n$，则 $AB=B$，从而 $r(AB)=r(B)$。

## 7. 结果速查与常见错误

| 题号 | 结果或关键公式 | 最容易漏掉的条件 |
| --- | --- | --- |
| 51(1) | $PQ$ 的左下块为 $O$，右下块为 $b\lvert A\rvert-\alpha^TA^*\alpha$ | 使用 $A^*A=\lvert A\rvert E$ |
| 51(2) | $Q$ 可逆当且仅当 $b\ne\alpha^TA^{-1}\alpha$ | $A$ 可逆，才可约去 $\lvert A\rvert$ |
| 68 | $(E-BA)^{-1}=E+B(E-AB)^{-1}A$ | 不能改变三个因子的次序 |
| 69 | $\lvert\lambda E-AB\rvert=\lvert\lambda E-BA\rvert$ | 约去 $\lambda^n$ 时使用 $\lambda\ne0$ |
| 70 | $r(A)+r(B)-n\le r(AB)\le\min\{r(A),r(B)\}$ | 减去的是公共尺寸 $n$ |

复习时重点检查：

- 是否把伴随矩阵 $A^*$ 误当成转置或共轭转置。
- 是否把块行操作中的系数写到右边，或把块列操作中的系数写到左边。
- 是否将分块矩阵直接套成数的二阶行列式公式。
- 是否在未说明非零时就约去行列式或参数因子。
- 是否把第 70 题的 $V^{-1}B$ 漏写成 $B$。
- 是否误认为 $r(AB)=\min\{r(A),r(B)\}$ 总成立；上述对角矩阵例子已说明可能严格小于。

## 8. 自测：能否不看解答说出理由

**问题 1。** 第 51 题若 $\alpha=0$，可逆条件是什么？  
**答案：** $b\ne0$。此时 $Q=\operatorname{diag}(A,b)$，与一般公式一致。

**问题 2。** 若第 51 题 $b=\alpha^TA^{-1}\alpha$，怎样直接证明 $Q$ 不可逆？  
**答案：** 非零向量 $(-A^{-1}\alpha,1)^T$ 被 $Q$ 映为零。

**问题 3。** 第 68 题是否要求 $A$、$B$ 各自可逆？  
**答案：** 不要求。只需要 $E-AB$ 可逆；整个证明没有使用 $A^{-1}$ 或 $B^{-1}$。

**问题 4。** 第 69 题能否推出 $AB$ 与 $BA$ 的秩相同？  
**答案：** 不能。第 5.2 节给出了秩分别为 $1$ 和 $0$ 的例子。

**问题 5。** 若 $A$ 为 $5\times4$ 矩阵，$B$ 为 $4\times6$ 矩阵，且两者的秩都是 $3$，则 $r(AB)$ 有哪些可能值？  
**答案：** 不等式给出 $2\le r(AB)\le3$。两值均可实现：取 $A$ 的前四行为 $\operatorname{diag}(1,1,1,0)$、第五行为零；若 $B$ 的前四列为 $\operatorname{diag}(1,1,1,0)$、末两列为零，则乘积秩为 $3$；若将该对角块改为 $\operatorname{diag}(0,1,1,1)$，则乘积秩为 $2$。

**问题 6。** 若 $A$ 为 $m\times n$ 且满列秩，对任意 $n\times s$ 矩阵 $B$，为何 $r(AB)=r(B)$？  
**答案：** 下界给出 $r(AB)\ge n+r(B)-n=r(B)$，上界给出反向不等式，因此相等。

## 9. 重新编译提交版

在仓库根目录运行以下命令；需要已安装包含 `ctex` 和 Fandol 字体的 TeX Live / MacTeX：

```bash
mkdir -p output/pdf
xelatex -halt-on-error -output-directory=output/pdf assignments/assignment8/assignment8-sol.tex
xelatex -halt-on-error -output-directory=output/pdf assignments/assignment8/assignment8-sol.tex
```

连续编译两遍可更新 PDF 书签等信息。提交文件为 `output/pdf/assignment8-sol.pdf`。
