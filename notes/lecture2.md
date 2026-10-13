# Lecture 2：矩阵、行列式与秩

> 课程：MATH1205 · 线性代数
>
> 范围：按仓库中第二章的六份课件整理，覆盖 §2.1–§2.6。
>
> 依据：[LA2.1.pdf](../lectures/LA2.1.pdf)、[LA2.2.pdf](../lectures/LA2.2.pdf)、[LA2.3.pdf](../lectures/LA2.3.pdf)、[LA2.4.pdf](../lectures/LA2.4.pdf)、[LA2.5.pdf](../lectures/LA2.5.pdf)、[LA2.6.pdf](../lectures/LA2.6.pdf)。下文页码均为各文件的 PDF 页码。
>
> 整理方式：按知识主线重组课件；标为「补充」或「条件说明」的内容用于解释推导、澄清适用范围。
>
> 记号：$F$ 表示本课程中的数域；$E_n$ 表示 $n$ 阶单位矩阵，$O$ 表示相应大小的零矩阵；$r(A)=\operatorname{rank}(A)$。课件中的 $A^*$ 表示**伴随矩阵**，不是共轭转置。

## 1. 本章主线与阅读导航

第一章学习用行变换解方程。本章把这些操作放进矩阵运算中，并引入两个重要量：

- **行列式**：定义在方阵上，判断矩阵是否可逆，并给出一类方程组的解公式。
- **秩**：适用于任意矩阵，等于消元后的主元数，用来统一描述方程组的解的情况。

本章知识之间的关系为

$$
\text{矩阵运算}
\longrightarrow \text{行列式与逆矩阵}
\longrightarrow \text{初等矩阵与秩}
\longrightarrow \text{方程组的判别与求解}.
$$

分块矩阵则把上述运算应用到较大的矩阵上，让结构清楚、计算简便。

| 课件 | 内容 | 对应笔记 |
| --- | --- | --- |
| §2.1 | 矩阵运算、幂、转置、迹 | 第 2–4 节 |
| §2.2 | 行列式、Laplace 展开、Cramer 法则 | 第 5–7 节 |
| §2.3 | 逆矩阵与伴随矩阵 | 第 8 节 |
| §2.4 | 分块矩阵 | 第 9 节 |
| §2.5 | 初等矩阵、秩、方程组与求逆 | 第 10–12 节 |
| §2.6 | 分块初等变换、降阶公式、秩不等式 | 第 13–14 节 |
| 综合复习 | 方法选择、易错点与自测 | 第 15–16 节 |

## 2. 矩阵的基本运算：先看维度，再看顺序

来源：LA2.1，第 2–18、24–28 页。

### 2.1 加减法与数乘

只有**同类型**矩阵才能相加减。若 $A,B\in F^{m\times n}$，则

$$
(A\pm B)_{ij}=a_{ij}\pm b_{ij},\qquad (cA)_{ij}=ca_{ij}.
$$

加法满足交换律和结合律；数乘满足

$$
c(A+B)=cA+cB,\qquad (c+d)A=cA+dA,\qquad c(dA)=(cd)A.
$$

例如，$3(A+C)=2(B-C)$ 可整理为

$$
5C=2B-3A,\qquad C=\frac15(2B-3A).
$$

这里除以的是非零**数** $5$，不是除以一个矩阵。

### 2.2 矩阵乘法的定义

若 $A\in F^{m\times p}$、$B\in F^{p\times n}$，则 $AB\in F^{m\times n}$，且

$$
\boxed{(AB)_{ij}=\sum_{k=1}^{p}a_{ik}b_{kj}}.
$$

即：左矩阵的第 $i$ 行，与右矩阵的第 $j$ 列，对应相乘再相加。

$$
\underbrace{A}_{m\times p}\underbrace{B}_{p\times n}
=\underbrace{AB}_{m\times n}.
$$

**内侧维度必须相等，外侧维度决定结果大小。** 逐项相乘得到的是另一种运算（Hadamard 积），不是这里的矩阵乘法。

课件例题：

$$
A=\begin{bmatrix}3&0&4\\-1&5&2\end{bmatrix},\qquad
B=\begin{bmatrix}1&0\\0&-1\\1&1\end{bmatrix}.
$$

则

$$
AB=\begin{bmatrix}7&4\\1&-3\end{bmatrix},\qquad
BA=\begin{bmatrix}3&0&4\\1&-5&-2\\2&5&6\end{bmatrix}.
$$

两个乘积都存在，但大小不同；一般还可能只有其中一个乘积有意义。

### 2.3 怎样理解 $Ax$ 与 $AB$？

把 $A$ 按列写为 $A=[\alpha_1\ \cdots\ \alpha_n]$，则

$$
\boxed{Ax=x_1\alpha_1+\cdots+x_n\alpha_n}.
$$

所以 $Ax=b$ 就是在问：能否用 $A$ 的各列按某些系数组合出 $b$？

若 $y=Ax$、$x=Bz$，则

$$
y=A(Bz)=(AB)z.
$$

因此乘积 $AB$ 表示**先作用 $B$，再作用 $A$**。这也是矩阵乘法采用“行乘列”的原因之一。

### 2.4 哪些数的运算法则仍然成立？

在维度允许的前提下：

$$
\begin{aligned}
(AB)C&=A(BC),\\
A(B+C)&=AB+AC,\\
(A+B)C&=AC+BC,\\
c(AB)&=(cA)B=A(cB),\\
E_mA&=AE_n=A\quad(A\in F^{m\times n}).
\end{aligned}
$$

结合律允许改变括号位置，**不允许任意交换因子的顺序**。证明矩阵恒等式的一般方法是：先检查两边类型，再比较任意位置的元素。

例如，结合律两边的 $(i,j)$ 元素都等于

$$
\sum_s\sum_t a_{it}b_{ts}c_{sj};
$$

有限求和可以交换求和次序，因此两边相等。

### 2.5 三个不能照搬的性质

**一、乘法一般不交换。** 例如

$$
A=\begin{bmatrix}1&1\\1&0\end{bmatrix},\quad
B=\begin{bmatrix}1&0\\1&1\end{bmatrix}
\quad\Longrightarrow\quad
AB=\begin{bmatrix}2&1\\1&0\end{bmatrix}
\ne
BA=\begin{bmatrix}1&1\\2&1\end{bmatrix}.
$$

只有当 $AB=BA$ 时，才称 $A,B$ **可交换**。

**二、非零矩阵相乘可能为零。**

$$
\begin{bmatrix}1&1\\2&2\end{bmatrix}
\begin{bmatrix}1&1\\-1&-1\end{bmatrix}=O.
$$

所以 $AB=O$ 不能推出 $A=O$ 或 $B=O$。

**三、非零因子不能直接约掉。** $AB=AC$ 且 $A\ne O$，不一定有 $B=C$。第 8 节会说明：若 $A$ 可逆，就可以在两边左乘 $A^{-1}$ 来消去它。

### 2.6 左乘与右乘对角矩阵

若 $D=\operatorname{diag}(d_1,\ldots,d_m)$，则 $DA$ 把 $A$ 的第 $i$ 行乘以 $d_i$。

若 $D'=\operatorname{diag}(c_1,\ldots,c_n)$，则 $AD'$ 把 $A$ 的第 $j$ 列乘以 $c_j$。

记忆：**左乘作用于行，右乘作用于列。** 第 10 节的初等矩阵会推广这一观察。

## 3. 方阵的幂、多项式与递推

来源：LA2.1，第 19–23、29–34 页。

### 3.1 定义与运算条件

对 $n$ 阶方阵 $A$，规定

$$
A^0=E_n,\qquad A^k=\underbrace{A\cdots A}_{k\text{ 个}},\quad k\ge1.
$$

对非负整数 $p,q$，有

$$
A^pA^q=A^{p+q},\qquad (A^p)^q=A^{pq}.
$$

若 $f(t)=c_0+c_1t+\cdots+c_kt^k$，则

$$
f(A)=c_0E_n+c_1A+\cdots+c_kA^k.
$$

**常数项应写成 $c_0E_n$，不能把一个数直接加到矩阵上。**

对两个同阶方阵，通常只有

$$
(A+B)^2=A^2+AB+BA+B^2.
$$

若已知 $AB=BA$，才可以合并中间两项，并使用通常的二项式定理、$(AB)^k=A^kB^k$ 等公式。同一个矩阵的多项式 $f(A),g(A)$ 总是可交换（补充）。

### 3.2 技巧一：把矩阵拆成列向量乘行向量

设

$$
A=\begin{bmatrix}1&-1&2\\-2&2&-4\\1&-1&2\end{bmatrix}
=uv^T,\qquad
u=\begin{bmatrix}1\\-2\\1\end{bmatrix},\quad
v=\begin{bmatrix}1\\-1\\2\end{bmatrix}.
$$

因为 $v^Tu=5$ 是一个数，所以

$$
A^2=u(v^Tu)v^T=5A,\qquad
\boxed{A^k=5^{k-1}A\quad(k\ge1)}.
$$

一般地，若 $A=uv^T$ 是方阵，则 $A^2=(v^Tu)A$。先算中间的 $1\times1$ 乘积，往往比反复计算大矩阵更简单。

### 3.3 技巧二：拆成数量矩阵与幂零矩阵

设

$$
A=\begin{bmatrix}\lambda&1&0\\0&\lambda&1\\0&0&\lambda\end{bmatrix}
=\lambda E+N,\qquad
N=\begin{bmatrix}0&1&0\\0&0&1\\0&0&0\end{bmatrix}.
$$

由于 $N^3=O$，且 $\lambda E$ 与 $N$ 可交换，二项式展开只留下前三项。对 $k\ge2$，

$$
\boxed{
A^k=\begin{bmatrix}
\lambda^k&k\lambda^{k-1}&\binom{k}{2}\lambda^{k-2}\\
0&\lambda^k&k\lambda^{k-1}\\
0&0&\lambda^k
\end{bmatrix}}.
$$

**条件说明：**$k=0,1$ 直接使用 $A^0=E,A^1=A$。上式中零次幂项按常数多项式 $1$ 理解；特别地，$\lambda=0$ 时直接由 $A=N$ 得到 $A^2=N^2$、$A^k=O$（$k\ge3$），无需引入负指数。

### 3.4 应用：奶蛙数量的递推

最初只有一只黄色奶蛙。每天每只奶蛙消失，变成另外两种颜色的奶蛙各一只。记第 $n$ 天黄、绿、红的数量为 $x_n,y_n,z_n$，则

$$
U_{n+1}=AU_n,\quad
A=\begin{bmatrix}0&1&1\\1&0&1\\1&1&0\end{bmatrix},\quad
U_0=\begin{bmatrix}1\\0\\0\end{bmatrix},\quad U_n=A^nU_0.
$$

第 3 天的数量为 $(2,3,3)$。求一般项时，先观察总数与差值：

$$
S_n=x_n+y_n+z_n=2^n,\qquad
y_n=z_n,\qquad x_n-y_n=(-1)^n.
$$

联立得

$$
\boxed{
x_n=\frac{2^n+2(-1)^n}{3},\qquad
y_n=z_n=\frac{2^n-(-1)^n}{3}}.
$$

因为 $2^n$ 不能被 $3$ 整除，三种颜色不可能数量完全相等。第 $2026$ 天只需在上述公式中代入偶数 $n=2026$。

**补充：矩阵幂形式。** 令 $J$ 为三阶全 $1$ 矩阵，利用 $A=J-E$、$J^2=3J$，可验证

$$
A^n=(-1)^nE+\frac{2^n-(-1)^n}{3}J.
$$

解递推时，既可以研究 $A^n$，也可以先找总和、差值等容易递推的量。

## 4. 转置、对称矩阵与迹

来源：LA2.1，第 16–17、35–44 页。

### 4.1 转置

若 $A\in F^{m\times n}$，则 $A^T\in F^{n\times m}$，且 $(A^T)_{ij}=a_{ji}$。转置把行与列互换：

$$
\begin{aligned}
(A^T)^T&=A,&(A+B)^T&=A^T+B^T,\\
(cA)^T&=cA^T,&\boxed{(AB)^T}&=\boxed{B^TA^T}.
\end{aligned}
$$

多个矩阵乘积转置时，因子顺序全部倒过来；$(A^k)^T=(A^T)^k$。

### 4.2 对称与反对称

- $A^T=A$：对称矩阵，即 $a_{ij}=a_{ji}$。
- $A^T=-A$：反对称矩阵，即 $a_{ij}=-a_{ji}$，其对角元全为 $0$。

任意方阵都有分解

$$
\boxed{A=\frac{A+A^T}{2}+\frac{A-A^T}{2}}.
$$

前一项对称，后一项反对称。**补充：**这一分解唯一，因为同时对称又反对称的矩阵只能是零矩阵。

若 $A$ 是三阶对称矩阵、$x=(x_1,x_2,x_3)^T$，则

$$
x^TAx=a_{11}x_1^2+a_{22}x_2^2+a_{33}x_3^2
+2a_{12}x_1x_2+2a_{13}x_1x_3+2a_{23}x_2x_3.
$$

这说明对称矩阵可以紧凑地表示二次式，交叉项系数是相应非对角元的两倍。

### 4.3 例：$B=E-2uu^T$

若列向量 $u$ 满足 $u^Tu=1$，则

$$
B^T=B,\qquad
B^2=E-4uu^T+4u(u^Tu)u^T=E.
$$

这里关键的一步是先把中间的 $u^Tu$ 化成数 $1$。由 $B^2=E$ 还可知 $B^{-1}=B$（补充）。

### 4.4 迹

方阵的迹是对角元之和：

$$
\operatorname{tr}(A)=\sum_{i=1}^n a_{ii}.
$$

它满足线性性、转置不变性，以及

$$
\boxed{\operatorname{tr}(AB)=\operatorname{tr}(BA)}.
$$

**补充：**这一等式也适用于 $A\in F^{m\times n},B\in F^{n\times m}$，因为两边都等于 $\sum_{i,j}a_{ij}b_{ji}$；但不能由迹相等推出 $AB=BA$。

对实矩阵 $A\in\mathbb R^{m\times n}$，

$$
\operatorname{tr}(AA^T)=\sum_{i=1}^m\sum_{j=1}^n a_{ij}^2\ge0.
$$

因此

$$
\boxed{\operatorname{tr}(AA^T)=0\iff A=O\quad(A\text{ 为实矩阵})}.
$$

**条件说明：**“实矩阵”不能省略。复数平方之和可能为零，例如 $A=[1\ i]$ 满足 $AA^T=0$，但 $A\ne O$。

## 5. 行列式：定义、展开与基本性质

来源：LA2.2，第 2、4–22、33–41 页。

### 5.1 行列式是一个数

对 $n$ 阶方阵 $A$，其行列式记作 $\det A$ 或 $|A|$。**矩阵是一张数表，行列式是由方阵算出的一个数。** 本课程不给非方阵定义行列式。

一阶时 $\det[a]=a$。二阶与三阶的公式为

$$
\begin{vmatrix}a&b\\c&d\end{vmatrix}=ad-bc,
$$

$$
\begin{aligned}
\det A={}&a_{11}a_{22}a_{33}+a_{12}a_{23}a_{31}+a_{13}a_{21}a_{32}\\
&-a_{13}a_{22}a_{31}-a_{11}a_{23}a_{32}-a_{12}a_{21}a_{33}.
\end{aligned}
$$

**条件说明：**二、三阶的对角线计算方法不能直接推广到四阶及以上。

### 5.2 余子式与代数余子式

对 $A$ 中的位置 $(i,j)$：

- 删去第 $i$ 行、第 $j$ 列，剩下的 $(n-1)$ 阶**行列式**叫余子式 $M_{ij}$。
- 加上符号得到代数余子式 $A_{ij}=(-1)^{i+j}M_{ij}$。

符号按棋盘交替：

$$
\begin{bmatrix}+&-&+&\cdots\\-&+&-&\cdots\\+&-&+&\cdots\\\vdots&\vdots&\vdots&\ddots\end{bmatrix}.
$$

余子式是一个数，不是删去行列后留下的子矩阵本身。$A_{ij}$ 也不是元素 $a_{ij}$。

**关键观察：**$A_{ij}$ 不依赖于第 $i$ 行和第 $j$ 列的元素。特别地，改变某一整行，不会改变这一行各位置对应的代数余子式。

### 5.3 按任意一行或一列展开

课件先按第一行递归定义行列式，再证明任意一行、任意一列均可展开：

$$
\boxed{\det A=\sum_{j=1}^n a_{ij}A_{ij}\quad\text{（固定行 }i\text{）}},
$$

$$
\boxed{\det A=\sum_{i=1}^n a_{ij}A_{ij}\quad\text{（固定列 }j\text{）}}.
$$

计算时优先选零最多的行或列。例如

$$
\begin{vmatrix}3&0&4\\1&1&2\\2&1&0\end{vmatrix}
=3\begin{vmatrix}1&2\\1&0\end{vmatrix}
+4\begin{vmatrix}1&1\\2&1\end{vmatrix}
=-6-4=-10.
$$

### 5.4 行列式的基本性质

下面所有关于行的性质，对列也成立，因为 $\det(A^T)=\det A$。

| 操作或结构 | 行列式的变化 |
| --- | --- |
| 交换两行 | 变号 |
| 一行乘数 $c$ | 整个行列式乘 $c$ |
| 一行加上另一行的 $c$ 倍 | 不变 |
| 有一行全零 | 行列式为 $0$ |
| 两行相同或成比例 | 行列式为 $0$ |
| 上三角、下三角或对角矩阵 | 等于主对角元的乘积 |

行列式对**单独一行**具有线性性：保持其他行不动时，

$$
\det(\ldots,u+v,\ldots)
=\det(\ldots,u,\ldots)+\det(\ldots,v,\ldots).
$$

它不意味着 $\det(A+B)=\det A+\det B$。例如 $A=B=E_2$ 时，左边为 $4$，右边为 $2$。

由于 $cA$ 把 $n$ 行都乘了 $c$，

$$
\boxed{\det(cA)=c^n\det A}.
$$

**与解方程的区别：**初等行变换保持方程组同解，但不一定保持行列式的值。计算行列式时，换行的负号、提取的公因子必须记录。

### 5.5 “本行配本行”与“本行配别行”

将第 $i$ 行的元素与第 $j$ 行各位置的代数余子式相乘后求和，有

$$
\boxed{\sum_{k=1}^n a_{ik}A_{jk}=\delta_{ij}\det A},\qquad
\delta_{ij}=\begin{cases}1,&i=j,\\0,&i\ne j.\end{cases}
$$

当 $i=j$，就是按行展开；当 $i\ne j$，把第 $j$ 行替换为第 $i$ 行后展开，所得行列式有两行相同，因而为零。

对应的列公式为

$$
\sum_{k=1}^n a_{ki}A_{kj}=\delta_{ij}\det A.
$$

这两个公式是 Cramer 法则和伴随矩阵恒等式的核心。

### 5.6 Laplace 定理：一次选多行展开

从 $A$ 中选出 $k$ 行、$k$ 列，保持原相对顺序，交叉处形成的 $k$ 阶行列式称为 **$k$ 阶子式**。

记选中的行标为 $I=\{i_1<\cdots<i_k\}$、列标为 $J=\{j_1<\cdots<j_k\}$，子矩阵记为 $A[I,J]$。对应的余子式是 $\det A[I^c,J^c]$，代数余子式还要乘

$$
(-1)^{\sum_{i\in I}i+\sum_{j\in J}j}.
$$

**Laplace 定理：**固定 $k$ 行，遍历所有 $k$ 列的选法，则

$$
\boxed{
\det A=\sum_{|J|=k}
(-1)^{\sum I+\sum J}\det A[I,J]\det A[I^c,J^c]}.
$$

$k=1$ 就是普通的按行展开。课件不要求此定理的证明，重点是利用大片零元素，使大部分子式为零。

例如，若 $B$ 为 $p$ 阶方阵、$D$ 为 $q$ 阶方阵，则

$$
\boxed{
\det\begin{bmatrix}B&O\\C&D\end{bmatrix}
=\det B\det D},\qquad
\boxed{
\det\begin{bmatrix}O&B\\D&C\end{bmatrix}
=(-1)^{pq}\det B\det D}.
$$

第二式交换两组列需跨过 $pq$ 对列，因此带有 $(-1)^{pq}$。这些公式**不要求对角块可逆**。

### 5.7 乘积的行列式

对同阶方阵 $A,B$，

$$
\boxed{\det(AB)=\det A\det B},\qquad
\det(A^k)=(\det A)^k.
$$

因此虽然 $AB$ 与 $BA$ 一般不同，但它们的行列式相同。

**课件证明思路：**对同一个分块行列式作两种计算：

$$
\det\begin{bmatrix}A&O\\-E&B\end{bmatrix}
=\det A\det B.
$$

把下方各行按 $A$ 的系数加到上方，使它变为

$$
\det\begin{bmatrix}O&AB\\-E&B\end{bmatrix}
=(-1)^{n^2}\det(AB)\det(-E)=\det(AB).
$$

因为 $n^2+n$ 为偶数，两个符号恰好抵消。

## 6. 行列式的典型计算与证明

来源：LA2.2，第 11–12、23–32、42–45 页。

### 6.1 化成三角形

$$
\begin{vmatrix}3&1&5\\3&2&5\\0&1&5\end{vmatrix}
\xrightarrow{R_2\leftarrow R_2-R_1}
\begin{vmatrix}3&1&5\\0&1&0\\0&1&5\end{vmatrix}
\xrightarrow{R_3\leftarrow R_3-R_2}
\begin{vmatrix}3&1&5\\0&1&0\\0&0&5\end{vmatrix}
=15.
$$

这里两次都是倍加，行列式不变。一般数字行列式优先考虑消元，避免机械地反复展开。

若反对角矩阵 $A$ 的第 $(i,n+1-i)$ 个元素为 $\lambda_i$，其余元素为零，则

$$
\det A=(-1)^{n(n-1)/2}\lambda_1\cdots\lambda_n.
$$

符号来自把列顺序完全颠倒时产生的交换次数。

### 6.2 对角元全为 $x$、其余元素全为 $a$

对 $n\ge2$，

$$
D_n=\begin{vmatrix}
x&a&\cdots&a\\a&x&\cdots&a\\
\vdots&\vdots&\ddots&\vdots\\a&a&\cdots&x
\end{vmatrix}
=\boxed{\bigl(x+(n-1)a\bigr)(x-a)^{n-1}}.
$$

推导步骤：

1. 把第 $2,\ldots,n$ 列全部加到第一列，第一列都成为 $x+(n-1)a$。
2. 从第一列提取该公因子，使第一列成为全 $1$。
3. 对 $j=2,\ldots,n$ 作 $C_j\leftarrow C_j-aC_1$，得到下三角行列式，其对角元为 $1,x-a,\ldots,x-a$。

**条件说明：**提取公因子是行列式的线性性，不是除以参数；即使公因子为零，恒等式仍成立。

### 6.3 Vandermonde（范德蒙）行列式

$$
V_n=\begin{vmatrix}
1&1&\cdots&1\\
x_1&x_2&\cdots&x_n\\
x_1^2&x_2^2&\cdots&x_n^2\\
\vdots&\vdots&\ddots&\vdots\\
x_1^{n-1}&x_2^{n-1}&\cdots&x_n^{n-1}
\end{vmatrix}
=\boxed{\prod_{1\le i<j\le n}(x_j-x_i)}.
$$

**记住差的方向：后面的 $x_j$ 减前面的 $x_i$。**

证明思路：从下往上依次作 $R_k\leftarrow R_k-x_1R_{k-1}$，使第一列除首项外都为零；沿第一列展开，再从各列提取 $x_j-x_1$，得到

$$
V_n(x_1,\ldots,x_n)=\prod_{j=2}^n(x_j-x_1)V_{n-1}(x_2,\ldots,x_n).
$$

归纳即得公式。必须从下往上做，才能保证使用的上一行还没被改变。

直接推论：$V_n\ne0$ 当且仅当 $x_1,\ldots,x_n$ 两两不同。

### 6.4 用行列式表示多项式

课件中的行列式

$$
F_n(\lambda)=\begin{vmatrix}
\lambda&0&\cdots&0&a_n\\
-1&\lambda&\cdots&0&a_{n-1}\\
0&-1&\ddots&0&a_{n-2}\\
\vdots&\vdots&\ddots&\lambda&a_2\\
0&0&\cdots&-1&\lambda+a_1
\end{vmatrix}
$$

按第一行展开满足 $F_n=\lambda F_{n-1}+a_n$，其中 $F_1=\lambda+a_1$，故

$$
\boxed{F_n(\lambda)=\lambda^n+a_1\lambda^{n-1}+\cdots+a_n}.
$$

### 6.5 证明题：不交换时保留乘法顺序

设 $A,B$ 都是三阶实方阵，且 $A^2=B^2$。不能直接写 $(A+B)(A-B)=A^2-B^2$。正确计算是

$$
(A+B)(A-B)=BA-AB=-(A-B)(A+B).
$$

令 $d=\det(A+B)\det(A-B)$，取行列式得到

$$
d=(-1)^3d=-d,
$$

因此 $d=0$。这个论证同样适用于任意**奇数阶**方阵。

另一类常见思路是先构造乘积。例如 $AA^T=E$ 且 $\det A<0$ 时，

$$
\det(A+E)=\det\bigl(A(E+A^T)\bigr)
=\det A\det(A+E),
$$

所以 $(1-\det A)\det(A+E)=0$，从而 $\det(A+E)=0$。

## 7. Cramer 法则与插值

来源：LA2.2，第 47–56 页。

### 7.1 Cramer（克拉默）法则

对 $n$ 个方程、$n$ 个未知数的方程组 $Ax=b$，若

$$
D=\det A\ne0,
$$

则方程组有唯一解

$$
\boxed{x_j=\frac{D_j}{D},\qquad j=1,\ldots,n},
$$

其中 $D_j$ 是把 $A$ 的**第 $j$ 列**替换成常数列 $b$ 后的行列式。

证明的关键是将第 $k$ 个方程乘 $A_{kj}$ 后求和，利用第 5.5 节的列公式，得到

$$
Dx_j=\sum_{k=1}^n b_kA_{kj}=D_j.
$$

这说明任何解都只能取上述值；再代回并使用行公式，可验证上述值确实组成一个解。因此存在性与唯一性都得到证明。

**使用前提：**系数矩阵必须是方阵，且 $D\ne0$。$D=0$ 时不能使用比值公式，更不能仅由此判为无解；应改用消元或秩。

### 7.2 例：多项式插值的存在唯一性

求次数不超过 $3$ 的多项式，使

$$
f(-1)=0,\quad f(1)=4,\quad f(2)=3,\quad f(3)=16.
$$

设 $f(x)=a_0+a_1x+a_2x^2+a_3x^3$，则系数满足

$$
\begin{bmatrix}
1&-1&1&-1\\1&1&1&1\\1&2&4&8\\1&3&9&27
\end{bmatrix}
\begin{bmatrix}a_0\\a_1\\a_2\\a_3\end{bmatrix}
=\begin{bmatrix}0\\4\\3\\16\end{bmatrix}.
$$

系数矩阵是 Vandermonde 矩阵的转置，其行列式为

$$
(1+1)(2+1)(3+1)(2-1)(3-1)(3-2)=48\ne0.
$$

因此所求多项式存在且唯一。**补充计算：**解得

$$
\boxed{f(x)=2x^3-5x^2+7}.
$$

一般地，在 $n$ 个互不相同的点上指定函数值，会唯一确定一个次数不超过 $n-1$ 的多项式。

### 7.3 齐次方程组

对方阵系数的齐次方程组 $Ax=0$：

- $\det A\ne0$，则只有零解。
- 有非零解，则 $\det A=0$。

第 11 节的秩判别还会给出逆向结论，从而得到

$$
Ax=0\text{ 有非零解}\iff\det A=0.
$$

**课件例题的技巧：**若能证明某个实矩阵满足 $AA^T=sE$，且 $s>0$，那么

$$
(\det A)^2=\det(AA^T)=s^n>0,
$$

因此 $A$ 可逆，$Ax=0$ 只有零解。这个方法可避免直接展开复杂的行列式。

## 8. 可逆矩阵与伴随矩阵

来源：LA2.3，第 2–21 页。

### 8.1 定义、唯一性与单边检验

对 $n$ 阶方阵 $A$，若存在同阶方阵 $B$ 使

$$
AB=BA=E,
$$

则称 $A$ 可逆，$B$ 是其逆矩阵，记为 $A^{-1}$。

逆矩阵唯一：若 $B,C$ 都是逆，则 $B=B(AC)=(BA)C=C$。

本节的核心判别为

$$
\boxed{A\text{ 可逆}\iff\det A\ne0}.
$$

所以“可逆矩阵”与“非奇异矩阵”是同一个概念。非零矩阵未必可逆，例如 $\operatorname{diag}(1,0)$。

**单边检验：**对同阶方阵，只需验证 $AB=E$ 或 $BA=E$ 中的一个，就可推出 $B=A^{-1}$。例如 $AB=E$ 给出 $\det A\det B=1$，故 $A$ 可逆，再左乘 $A^{-1}$ 即得结论。

**条件说明：**不能把单边检验直接套到矩形矩阵上。例如 $A=[1\ 0]$、$B=[1\ 0]^T$ 满足 $AB=E_1$，但 $BA\ne E_2$。

### 8.2 伴随矩阵：代数余子式排好后再转置

设 $A_{ij}$ 为 $a_{ij}$ 的代数余子式，定义

$$
A^*=\begin{bmatrix}
A_{11}&A_{21}&\cdots&A_{n1}\\
A_{12}&A_{22}&\cdots&A_{n2}\\
\vdots&\vdots&\ddots&\vdots\\
A_{1n}&A_{2n}&\cdots&A_{nn}
\end{bmatrix}.
$$

注意 $A_{ij}$ 放在伴随矩阵的 **$(j,i)$ 位置**。

利用第 5.5 节的行、列公式，立刻得到

$$
\boxed{AA^*=A^*A=(\det A)E}.
$$

因此当 $\det A\ne0$ 时，

$$
\boxed{A^{-1}=\frac{1}{\det A}A^*}.
$$

这同时证明了“行列式非零”足以保证可逆。**伴随矩阵在 $A$ 不可逆时也有定义；逆矩阵则不一定存在。**

二阶公式特别常用：

$$
\boxed{
\begin{bmatrix}a&b\\c&d\end{bmatrix}^{-1}
=\frac1{ad-bc}\begin{bmatrix}d&-b\\-c&a\end{bmatrix}},\qquad ad-bc\ne0.
$$

### 8.3 逆矩阵的运算性质

以下 $A,B$ 为同阶可逆方阵，$c\ne0$：

$$
\begin{aligned}
(A^{-1})^{-1}&=A,\\
(cA)^{-1}&=c^{-1}A^{-1},\\
\boxed{(AB)^{-1}}&=\boxed{B^{-1}A^{-1}},\\
(A^T)^{-1}&=(A^{-1})^T,\\
\det(A^{-1})&=(\det A)^{-1},\\
(A^k)^{-1}&=(A^{-1})^k.
\end{aligned}
$$

乘积求逆与乘积转置都要**倒序**。逆矩阵一般不能逐项求倒数，也没有 $(A+B)^{-1}=A^{-1}+B^{-1}$ 这样的加法规则。

### 8.4 矩阵方程：从正确的一侧乘逆

若 $A\in F^{m\times m}$、$B\in F^{n\times n}$ 均可逆，$C\in F^{m\times n}$，则

$$
\boxed{AXB=C\iff X=A^{-1}CB^{-1}}.
$$

对应地，$AX=C$ 得 $X=A^{-1}C$；$XB=C$ 得 $X=CB^{-1}$。

课件例题：

$$
\begin{bmatrix}1&-5\\-1&4\end{bmatrix}X
=\begin{bmatrix}3&2\\1&4\end{bmatrix}
\quad\Longrightarrow\quad
X=\begin{bmatrix}-4&-5\\-1&-1\end{bmatrix}
\begin{bmatrix}3&2\\1&4\end{bmatrix}
=\begin{bmatrix}-17&-28\\-4&-6\end{bmatrix}.
$$

### 8.5 从已知恒等式构造逆矩阵

**例 1：**若 $A^3+A-2E=O$，则

$$
A\left(\frac12(A^2+E)\right)=E,
\qquad\boxed{A^{-1}=\frac12(A^2+E)}.
$$

思路是把已知式改写成“$A$ 乘某矩阵等于 $E$”。

**例 2：**若 $A+B+AB=O$，则

$$
(A+E)(B+E)=E.
$$

单边检验说明 $(B+E)(A+E)=E$ 也成立，于是 $BA=-A-B=AB$，即两矩阵可交换。

**例 3：**若 $A,B,A+B$ 均可逆，则

$$
A^{-1}+B^{-1}=A^{-1}(A+B)B^{-1}=B^{-1}(A+B)A^{-1},
$$

所以

$$
\boxed{(A^{-1}+B^{-1})^{-1}
=B(A+B)^{-1}A=A(A+B)^{-1}B}.
$$

这里的相等来自两种合法的分解，不能解释成任意交换因子。

### 8.6 伴随矩阵的性质与计算

对同阶方阵，

$$
(AB)^*=B^*A^*,\qquad
\boxed{\det(A^*)=(\det A)^{n-1}\quad(n\ge2)}.
$$

第二式的证明要分类：

- 若 $\det A\ne0$，由 $A^*=(\det A)A^{-1}$ 取行列式即可。
- 若 $\det A=0$，假设 $A^*$ 可逆，则由 $AA^*=O$ 得 $A=O$。当 $n\ge2$ 时，零矩阵的全部 $(n-1)$ 阶余子式均为零，于是 $A^*=O$，矛盾。

例：$A$ 为五阶方阵，$\det A=3$，则

$$
\det(2A^*)=2^5\cdot3^4,
$$

$$
\det(2A^*-7A^{-1})
=\det\bigl((6-7)A^{-1}\bigr)
=(-1)^5\frac13=-\frac13.
$$

遇到 $A^*$ 与 $A^{-1}$ 混合的题，优先用 $A^*=(\det A)A^{-1}$ 统一表示。

## 9. 分块矩阵：按结构计算

来源：LA2.4，第 2–18 页。

### 9.1 分块只是换一种记法

用横线和竖线把一个矩阵分成若干子矩阵，每个子矩阵称为一个**块**。例如

$$
A=\begin{bmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{bmatrix}.
$$

本节的 $A_{ij}$ 表示子矩阵，**不是代数余子式**，含义由上下文决定。

分块的目的在于利用零块、单位块、对角块等结构。它不改变原矩阵，也不改变运算结果。

### 9.2 加法、乘法与转置

加法要求对应块大小相同；数乘逐块进行。乘法仍是“行乘列”，但要求左边的列分块尺寸与右边的行分块尺寸相容：

$$
\begin{bmatrix}A&B\\C&D\end{bmatrix}
\begin{bmatrix}P&Q\\R&S\end{bmatrix}
=\begin{bmatrix}AP+BR&AQ+BS\\CP+DR&CQ+DS\end{bmatrix}.
$$

块内部仍是矩阵乘法，**不可随意交换 $BR$ 与 $RB$**。

转置既要交换块的位置，又要转置每个块：

$$
\boxed{
\begin{bmatrix}A&B\\C&D\end{bmatrix}^{T}
=\begin{bmatrix}A^T&C^T\\B^T&D^T\end{bmatrix}}.
$$

常用的按列分块公式是

$$
A[\beta_1\ \cdots\ \beta_k]
=[A\beta_1\ \cdots\ A\beta_k].
$$

所以求解 $AX=B$ 可以理解为同时求解 $B$ 的各列对应的多个方程组。

### 9.3 分块对角矩阵

若 $A=\operatorname{diag}(A_1,\ldots,A_s)$，各 $A_i$ 都为方阵，则

$$
\det A=\prod_{i=1}^s\det A_i,\qquad
A^k=\operatorname{diag}(A_1^k,\ldots,A_s^k).
$$

当且仅当各块都可逆时，$A$ 可逆，且

$$
\boxed{A^{-1}=\operatorname{diag}(A_1^{-1},\ldots,A_s^{-1})}.
$$

例如

$$
\begin{bmatrix}1&0&0&0\\2&-1&0&0\\0&0&1&0\\0&0&0&2\end{bmatrix}^{-1}
=\begin{bmatrix}1&0&0&0\\2&-1&0&0\\0&0&1&0\\0&0&0&\tfrac12\end{bmatrix},
$$

因为左上 $2\times2$ 块的平方为 $E_2$。

### 9.4 分块反对角矩阵与分块三角矩阵

若 $B,C$ 都可逆，则

$$
\boxed{
\begin{bmatrix}O&B\\C&O\end{bmatrix}^{-1}
=\begin{bmatrix}O&C^{-1}\\B^{-1}&O\end{bmatrix}}.
$$

**条件说明：**$B,C$ 可以有不同的阶数。此时逆矩阵的行、列分块尺寸随之交换，不能把各块机械地放回原来的位置。

对分块上三角矩阵，若 $B,C$ 为可逆方阵，则

$$
\boxed{
\begin{bmatrix}B&D\\O&C\end{bmatrix}^{-1}
=\begin{bmatrix}B^{-1}&-B^{-1}DC^{-1}\\O&C^{-1}\end{bmatrix}}.
$$

推导：设逆矩阵为 $\begin{bmatrix}X&Z\\W&Y\end{bmatrix}$，由乘积为单位阵得

$$
CW=O,\quad CY=E,\quad BX+DW=E,\quad BZ+DY=O.
$$

依次解得 $W=O,Y=C^{-1},X=B^{-1},Z=-B^{-1}DC^{-1}$。

**注意：**一般分块矩阵的行列式不能直接套“$AD-BC$”，一般分块矩阵的逆也不能逐块求逆。第 13 节会给出带条件的正确降阶方法。

## 10. 初等矩阵与相抵标准形

来源：LA2.5，第 2–12、19–20、24 页。

### 10.1 初等变换与初等矩阵

矩阵的初等变换包括初等**行**变换和初等**列**变换，均有三类：交换、乘非零数、倍加。

由单位矩阵经过**一次**初等变换得到的矩阵，称为初等矩阵。

| 操作 | 初等矩阵的构造 | 逆操作 |
| --- | --- | --- |
| $R_i\leftrightarrow R_j$ | 交换 $E$ 的第 $i,j$ 行 | 再交换一次 |
| $R_i\leftarrow cR_i$，$c\ne0$ | 将 $E$ 的第 $i$ 个对角元换成 $c$ | 乘 $1/c$ |
| $R_i\leftarrow R_i+cR_j$，$i\ne j$ | 在 $E$ 的 $(i,j)$ 位置填入 $c$ | 加上 $-c$ 倍的第 $j$ 行 |

因此初等矩阵都可逆，逆矩阵仍是同类初等矩阵。

### 10.2 左乘作行变换，右乘作列变换

若 $A$ 为 $m\times n$ 矩阵：

- 左乘一个 $m$ 阶初等矩阵，相当于对 $A$ 作一次初等行变换。
- 右乘一个 $n$ 阶初等矩阵，相当于对 $A$ 作一次初等列变换。

**倍加时要核对下标方向。** 记 $e_i$ 为第 $i$ 个标准列向量，令

$$
P=E+ce_ie_j^T\quad(i\ne j).
$$

则

$$
PA:\ R_i\leftarrow R_i+cR_j,\qquad
AP:\ C_j\leftarrow C_j+cC_i.
$$

同一个 $P$ 左乘与右乘时，被修改的行、列下标不同。

连续先作 $P_1$、再作 $P_2$、最后作 $P_s$ 对应的行变换，结果为

$$
P_s\cdots P_2P_1A.
$$

越后执行的左乘矩阵，写得越靠左。

### 10.3 为什么行变换保持方程组同解？

行变换把 $[A\mid b]$ 变为 $[PA\mid Pb]$，其中 $P$ 可逆。因此

$$
Ax=b\iff PAx=Pb.
$$

逆向只需左乘 $P^{-1}$。这给出了第一章“可逆变换保持同解”的矩阵证明。

**补充：列变换的含义。** 对系数矩阵作列变换得到 $AQ$ 时，$AQy=b$ 对应换元 $x=Qy$。若不记录换元关系，就不能把新未知量的取值直接当成原方程的解。求原通解和用 $[A\mid E]$ 求逆时，应按规定做行变换。

### 10.4 相抵与标准形

如果 $A$ 经过有限次初等行、列变换可变成 $B$，就称 $A,B$ **相抵**。相抵具有反身性、对称性与传递性。

任意 $m\times n$ 矩阵都能通过行、列变换化为

$$
\boxed{J_r=\begin{bmatrix}E_r&O\\O&O\end{bmatrix}_{m\times n}},
$$

称为**相抵标准形**。其中某些零块可以没有，$r=0$ 时就是零矩阵。

构造思路：先用行变换得到规范阶梯形；交换列，把主元列移到前面；再用列倍加消去剩余元素。

等价地，存在可逆矩阵 $P\in F^{m\times m},Q\in F^{n\times n}$，使

$$
PAQ=J_r.
$$

**不要混淆：**规范阶梯形只需要行变换；相抵标准形允许同时使用行、列变换。第 11 节会证明其中的 $r$ 就是矩阵的秩。

## 11. 矩阵的秩与线性方程组

来源：LA2.5，第 13–30 页。

### 11.1 秩的定义

矩阵 $A$ 的秩，是其**非零子式的最高阶数**，记为 $r(A)$。规定 $r(O)=0$。

因此证明 $r(A)=r$，需要说明：

1. 至少有一个 $r$ 阶子式非零；
2. 所有更高阶的子式都为零；实际只需检查所有 $r+1$ 阶子式为零，更高阶结论由展开得到。

当 $r=\min(m,n)$ 时，没有更高阶子式需要检查。

课件例子：

$$
A=\begin{bmatrix}3&2&1&-1\\0&2&3&0\\-3&4&8&1\end{bmatrix}.
$$

左上角二阶子式为 $6\ne0$；又有 $R_3=-R_1+3R_2$，所以所有三阶子式为零，故 $r(A)=2$。

### 11.2 秩的基本性质与计算方法

$$
\begin{gathered}
0\le r(A_{m\times n})\le\min(m,n),\qquad r(A)=0\iff A=O,\\
r(A^T)=r(A),\qquad r(cA)=r(A)\quad(c\ne0).
\end{gathered}
$$

若 $B$ 是 $A$ 的子矩阵，则 $r(B)\le r(A)$。

核心定理是：**初等行、列变换不改变秩。** 因而

$$
\boxed{r(A)=\text{行阶梯形的非零行数}=\text{主元数}}.
$$

求秩通常只需消元到阶梯形，没必要继续化成规范阶梯形，更没必要枚举所有子式。

**证明思路：**利用子式的线性性，说明初等变换后不会使最大非零子式的阶数下降；再用逆变换得到反向不等式，于是秩相等。特别要注意：保持秩不等于保持每一个子式的值。

### 11.3 满秩、可逆与相抵

- $r(A_{m\times n})=m$：行满秩，要求 $m\le n$。
- $r(A_{m\times n})=n$：列满秩，要求 $n\le m$。
- 对 $n$ 阶方阵，$r(A)=n$ 称为满秩。

将本章结论串联起来，对 $n$ 阶方阵有

$$
\boxed{
A\text{ 可逆}
\iff\det A\ne0
\iff r(A)=n
\iff A\text{ 可用初等行变换化为 }E_n}.
$$

此外，可逆矩阵恰好是有限个初等矩阵的乘积；可逆矩阵左乘或右乘都不改变秩：

$$
\boxed{r(PAQ)=r(A)\quad(P,Q\text{ 可逆且维度相容})}.
$$

**补充：**同类型矩阵相抵，当且仅当秩相同，因为它们可化为同一个 $J_r$。

### 11.4 秩为 $1$ 与满秩分解

$$
\boxed{r(A)=1\iff A=\alpha\beta^T,\quad\alpha\ne0,\ \beta\ne0}.
$$

右边说明每一列都是同一个非零列向量 $\alpha$ 的倍数，并且至少一列非零。反方向可由相抵标准形的分解得到。

更一般地，若 $r(A)=r>0$，则存在

$$
A=BC,\qquad B\in F^{m\times r},\quad C\in F^{r\times n},\qquad
r(B)=r(C)=r.
$$

这称为**满秩分解**，一般不唯一。它把一个矩阵写成列满秩矩阵与行满秩矩阵的乘积。

### 11.5 用秩判别 $Ax=b$ 的解

设 $A\in F^{m\times n}$，增广矩阵为 $\widetilde A=[A\mid b]$，则

$$
\boxed{Ax=b\text{ 有解}\iff r(A)=r(\widetilde A)}.
$$

| 秩的关系 | 解的情况 |
| --- | --- |
| $r(A)<r(\widetilde A)$ | 无解 |
| $r(A)=r(\widetilde A)=n$ | 唯一解 |
| $r(A)=r(\widetilde A)=r<n$ | 无穷多解，自由变量有 $n-r$ 个 |

**一定先检查相容性，再数自由变量。** 这里的 $n$ 是未知数个数，即系数矩阵的列数，不是方程数。

对齐次方程组，零解总存在，故

$$
\boxed{Ax=0\text{ 只有零解}\iff r(A)=n},
\qquad
\boxed{Ax=0\text{ 有非零解}\iff r(A)<n}.
$$

**条件说明：**课件第 26 页末尾“方程个数小于未知数个数，必有非零解”的推论应在**齐次方程组**下使用。非齐次方程组在 $m<n$ 时可能无解；只能说若有解，就有无穷多解。

### 11.6 例：用 Vandermonde 直接判无解

考虑

$$
\begin{cases}
x_1+x_2+x_3=1,\\
ax_1+bx_2+cx_3=d,\\
a^2x_1+b^2x_2+c^2x_3=d^2,\\
a^3x_1+b^3x_2+c^3x_3=d^3,
\end{cases}
$$

其中 $a,b,c,d$ 两两不同。增广矩阵是四阶 Vandermonde 矩阵，故

$$
r(\widetilde A)=4,\qquad r(A)\le3.
$$

两秩不同，所以无解，无需逐项消元。

### 11.7 含参数方程组：避免提前除以参数

讨论

$$
\begin{cases}
x_1+x_2-x_3=1,\\
2x_1+(a+2)x_2-(b+2)x_3=3,\\
-3ax_2+(a+2b)x_3=-3.
\end{cases}
$$

先作 $R_2\leftarrow R_2-2R_1$，再作 $R_3\leftarrow R_3+3R_2$，得到

$$
\left[\begin{array}{ccc|c}
1&1&-1&1\\0&a&-b&1\\0&0&a-b&0
\end{array}\right].
$$

现在才对可能的零主元 $a$ 与 $a-b$ 分类：

| 条件 | 结论 | 解 |
| --- | --- | --- |
| $a\ne0$ 且 $a\ne b$ | 唯一解 | $(x_1,x_2,x_3)=(1-1/a,1/a,0)$ |
| $a=b\ne0$ | 无穷多解 | $(x_1,x_2,x_3)=(1-1/a,1/a+t,t)$，$t\in F$ |
| $a=0,b\ne0$ | 无解 | 第二行要求 $-bx_3=1$，第三行要求 $-bx_3=0$ |
| $a=b=0$ | 无解 | 第二行直接为 $0=1$ |

于是也可概括为：$a=0$ 总无解；$a\ne0$ 时，再看 $a=b$ 是否成立。

## 12. 用初等行变换求逆与解矩阵方程

来源：LA2.5，第 31–33 页。

### 12.1 增广矩阵求逆

对方阵 $A$，构造 $[A\mid E]$，对整个增广矩阵同步作初等行变换：

$$
\boxed{[A\mid E]\ \xrightarrow{\text{初等行变换}}\ [E\mid A^{-1}]}.
$$

原理：若全部行变换的总效果是左乘 $P$，且 $PA=E$，则 $P=A^{-1}$，右半边恰好变成 $PE=P$。

若左半边无法化成 $E$，即 $r(A)<n$，则 $A$ 不可逆。

### 12.2 完整例题

求

$$
A=\begin{bmatrix}4&2&3\\3&1&2\\2&1&1\end{bmatrix}
$$

的逆。先作 $R_1\leftarrow R_1-2R_3$、$R_2\leftarrow R_2-R_3$：

$$
\left[\begin{array}{ccc|ccc}
4&2&3&1&0&0\\3&1&2&0&1&0\\2&1&1&0&0&1
\end{array}\right]
\longrightarrow
\left[\begin{array}{ccc|ccc}
0&0&1&1&0&-2\\1&0&1&0&1&-1\\2&1&1&0&0&1
\end{array}\right].
$$

再作 $R_2\leftarrow R_2-R_1$、$R_3\leftarrow R_3-R_1$，然后作 $R_3\leftarrow R_3-2R_2$：

$$
\longrightarrow
\left[\begin{array}{ccc|ccc}
0&0&1&1&0&-2\\1&0&0&-1&1&1\\0&1&0&1&-2&1
\end{array}\right].
$$

最后按原第 $2,3,1$ 行的顺序重新排列，左侧成为 $E_3$，故

$$
\boxed{A^{-1}=\begin{bmatrix}-1&1&1\\1&-2&1\\1&0&-2\end{bmatrix}}.
$$

可通过 $AA^{-1}=E_3$ 检查结果。

### 12.3 直接解 $AX=B$

若 $A$ 可逆，可以直接化简

$$
\boxed{[A\mid B]\ \xrightarrow{\text{初等行变换}}\ [E\mid A^{-1}B]}.
$$

右半边就是 $X$，不必先单独求出整个 $A^{-1}$。

例如课件给出

$$
A=\begin{bmatrix}5&1&-5\\3&-3&2\\1&-2&1\end{bmatrix},\qquad
B=\begin{bmatrix}-8&-5\\3&9\\0&0\end{bmatrix},
$$

化简得到

$$
X=\begin{bmatrix}1&4\\2&5\\3&6\end{bmatrix}.
$$

**补充：**$A$ 不可逆时，矩阵方程未必无解；需要分别检查 $B$ 的各列是否使方程相容，此时不能套 $A^{-1}B$。

## 13. 分块初等变换与行列式降阶

来源：LA2.6，第 2–10 页。

### 13.1 块行、块列的操作规则

分块初等变换有三类：

1. 交换两个块行或两个块列。
2. 某块行左乘一个可逆矩阵，或某块列右乘一个可逆矩阵。
3. 一个块行左乘某矩阵后加到另一块行；或一个块列右乘某矩阵后加到另一块列。

所有操作都必须满足尺寸要求。仍然要记住：**块行的系数在左边，块列的系数在右边。**

这些操作可由普通初等变换实现，所以保持秩。第三类块倍加还保持行列式不变。

**补充：**对相邻、分别含 $p,q$ 行的两个块行，交换后行列式乘 $(-1)^{pq}$，不能一律写“变号”；某块行左乘可逆矩阵 $P$ 时，整个方阵的行列式乘 $\det P$。

### 13.2 行列式的降阶公式

设

$$
M=\begin{bmatrix}A&B\\C&D\end{bmatrix},
\qquad A\in F^{m\times m},\quad D\in F^{n\times n}.
$$

**若 $A$ 可逆**，作块行变换

$$
\mathcal R_2\leftarrow\mathcal R_2-CA^{-1}\mathcal R_1,
$$

得到

$$
\begin{bmatrix}A&B\\C&D\end{bmatrix}
\longrightarrow
\begin{bmatrix}A&B\\O&D-CA^{-1}B\end{bmatrix}.
$$

于是

$$
\boxed{\det M=\det A\cdot\det(D-CA^{-1}B)}.
$$

**若 $D$ 可逆**，类似地，

$$
\boxed{\det M=\det D\cdot\det(A-BD^{-1}C)}.
$$

这两个公式分别只要求所求逆的那个对角块可逆。若 $A,D$ 都可逆，则

$$
\det A\det(D-CA^{-1}B)
=\det D\det(A-BD^{-1}C).
$$

**补充术语与结论：**$D-CA^{-1}B$ 称为关于 $A$ 的 Schur 补。继续用块列变换消去 $B$ 后，可得

$$
r(M)=m+r(D-CA^{-1}B)\quad(A\text{ 可逆}).
$$

因此在 $A$ 可逆的前提下，$M$ 可逆当且仅当该 Schur 补可逆。

### 13.3 将大行列式换成小行列式

若 $A\in F^{m\times n},B\in F^{n\times m}$，则

$$
\boxed{\lambda^n\det(\lambda E_m-AB)
=\lambda^m\det(\lambda E_n-BA)}.
$$

课件用分块初等变换证明：从

$$
K=\begin{bmatrix}\lambda E_m&A\\\lambda B&\lambda E_n\end{bmatrix}
$$

出发，有两种化简：

$$
\mathcal C_1\leftarrow\mathcal C_1-\mathcal C_2B
\quad\Longrightarrow\quad
\begin{bmatrix}\lambda E_m-AB&A\\O&\lambda E_n\end{bmatrix},
$$

$$
\mathcal R_2\leftarrow\mathcal R_2-B\mathcal R_1
\quad\Longrightarrow\quad
\begin{bmatrix}\lambda E_m&A\\O&\lambda E_n-BA\end{bmatrix}.
$$

两者行列式相等，得到上述恒等式。证明没有除以 $\lambda$，因此也覆盖 $\lambda=0$。

取 $\lambda=1$，得到特别常用的

$$
\boxed{\det(E_m-AB)=\det(E_n-BA)}.
$$

把 $B$ 换成 $-B$，也有 $\det(E_m+AB)=\det(E_n+BA)$。两边矩阵阶数可以不同；选阶数较小的一边计算。

### 13.4 例：秩一结构的行列式

令 $u=(a_1,\ldots,a_n)^T$，则

$$
\begin{vmatrix}
1+a_1^2&a_1a_2&\cdots&a_1a_n\\
a_2a_1&1+a_2^2&\cdots&a_2a_n\\
\vdots&\vdots&\ddots&\vdots\\
a_na_1&a_na_2&\cdots&1+a_n^2
\end{vmatrix}
=\det(E_n+uu^T)
=\boxed{1+u^Tu=1+\sum_{i=1}^n a_i^2}.
$$

关键是识别 $uu^T$，把 $n$ 阶行列式降为一阶。

### 13.5 例：降到二阶

若 $G$ 的对角元为 $a_i^2$，非对角元为 $a_ia_j+1$，令

$$
U=\begin{bmatrix}a_1&1\\a_2&1\\\vdots&\vdots\\a_n&1\end{bmatrix},
\qquad G=-E_n+UU^T.
$$

设 $s_1=\sum_i a_i$、$s_2=\sum_i a_i^2$，则

$$
U^TU=\begin{bmatrix}s_2&s_1\\s_1&n\end{bmatrix}.
$$

于是

$$
\begin{aligned}
\det G
&=(-1)^n\det(E_n-UU^T)\\
&=(-1)^n\det(E_2-U^TU)\\
&=\boxed{(-1)^n\bigl((1-s_2)(1-n)-s_1^2\bigr)}.
\end{aligned}
$$

不要漏掉从 $G=-(E_n-UU^T)$ 提出的 $(-1)^n$。

## 14. 秩不等式与典型应用

来源：LA2.6，第 11–19 页。

### 14.1 乘积的秩：上界与下界

对 $A\in F^{m\times n},B\in F^{n\times p}$，

$$
\boxed{r(A)+r(B)-n\le r(AB)\le\min\{r(A),r(B)\}}.
$$

下界中的 $n$ 是**乘积的中间维度**，即 $A$ 的列数、$B$ 的行数。若该下界为负，还可与 $r(AB)\ge0$ 合用。

上界的理解：先把 $A$ 行化简成只有 $r(A)$ 个非零行的阶梯形，再右乘 $B$，零行仍为零，故 $r(AB)\le r(A)$。对转置使用同一结论，得到 $r(AB)\le r(B)$。

下界的证明用到分块矩阵。下节的分块三角秩不等式给出

$$
r(A)+r(B)\le
r\begin{bmatrix}A&O\\-E_n&B\end{bmatrix}.
$$

作块行、块列倍加可化为

$$
\begin{bmatrix}A&O\\-E_n&B\end{bmatrix}
\longrightarrow
\begin{bmatrix}O&AB\\-E_n&B\end{bmatrix}
\longrightarrow
\begin{bmatrix}O&AB\\-E_n&O\end{bmatrix},
$$

其秩为 $n+r(AB)$，从而得到下界。

**常用推论：**若 $AB=O$，则

$$
r(A)+r(B)\le n.
$$

### 14.2 和、拼接与分块矩阵的秩

若 $A,B$ 同类型，则

$$
\boxed{r(A+B)\le r([A\ B])\le r(A)+r(B)}.
$$

第一步可由

$$
A+B=[A\ B]\begin{bmatrix}E\\E\end{bmatrix}
$$

和乘积秩上界得到。第二步可由

$$
[A\ B]=[E\ E]\begin{bmatrix}A&O\\O&B\end{bmatrix}
$$

得到。

对任意大小的矩阵 $A,B$，以及尺寸相容的 $C$，有

$$
\boxed{r\begin{bmatrix}A&O\\O&B\end{bmatrix}=r(A)+r(B)},
$$

$$
\boxed{r\begin{bmatrix}A&O\\C&B\end{bmatrix}\ge r(A)+r(B)}.
$$

证明思路是分别把 $A,B$ 化为相抵标准形：两块中的主元都可以保留；左下角 $C$ 可能增加秩，因此后一式一般只能写“$\ge$”。

### 14.3 伴随矩阵的秩只有三种情况

设 $A$ 为 $n\ge2$ 阶方阵，则

$$
\boxed{
r(A^*)=\begin{cases}
n,&r(A)=n,\\
1,&r(A)=n-1,\\
0,&r(A)\le n-2.
\end{cases}}
$$

证明分三类：

1. $r(A)=n$ 时，$A^*=(\det A)A^{-1}$ 可逆，故秩为 $n$。
2. $r(A)\le n-2$ 时，所有 $(n-1)$ 阶子式都为零，因此 $A^*=O$。
3. $r(A)=n-1$ 时，至少一个 $(n-1)$ 阶子式非零，所以 $A^*\ne O$；又 $AA^*=O$，故 $r(A)+r(A^*)\le n$，于是 $1\le r(A^*)\le1$。

### 14.4 幂等矩阵的秩

若 $A^2=A$，称 $A$ 为幂等矩阵。此时

$$
\boxed{r(A)+r(E-A)=n}.
$$

证明通常采用“夹住”的思路：

$$
n=r(E)=r(A+(E-A))\le r(A)+r(E-A),
$$

而 $A(E-A)=O$ 又给出

$$
r(A)+r(E-A)\le n.
$$

上下界相同，所以得到等式。

### 14.5 无解的证明可以写成另一个方程组

对 $A\in F^{m\times n},b\in F^m$，以下两个方程组**恰有一个有解**：

$$
Ax=b,
\qquad
\begin{bmatrix}A^T\\b^T\end{bmatrix}y
=\begin{bmatrix}0_n\\1\end{bmatrix}.
$$

后一个方程组即 $A^Ty=0$、$b^Ty=1$。

**直观解释（补充）：**若找到这样的 $y$，把原方程的第 $i$ 行乘 $y_i$ 后相加，左边变成 $0$，右边变成 $1$，就得到矛盾 $0=1$。因此 $y$ 给出了原方程组无解的一个可核验依据。

反过来为何无解时一定能找到 $y$？课件用秩说明：令

$$
r=r(A),\qquad s=r([A\mid b]),
$$

则 $r\le s\le r+1$，且

$$
r\begin{bmatrix}A^T\\b^T\end{bmatrix}=s,\qquad
r\begin{bmatrix}A^T&0_n\\b^T&1\end{bmatrix}=r+1.
$$

原方程有解当且仅当 $s=r$；后一个方程组有解当且仅当 $s=r+1$。两种情况恰好互补。

## 15. 做题方法与易错点

### 15.1 根据目标选择工具

| 目标 | 优先考虑的方法 |
| --- | --- |
| 计算矩阵乘积 | 检查维度，再逐行逐列计算；零块多时分块 |
| 计算 $A^k$ | 寻找 $A^2=cA$、$uv^T$、$\lambda E+N$、分块对角等结构 |
| 计算数字行列式 | 倍加消元化三角形，记录换行与倍乘因子 |
| 计算参数行列式 | 找相同行列和、公因子、Vandermonde 或分块结构 |
| 证明矩阵可逆 | 证明行列式非零、满秩，或构造单边逆 |
| 求二阶逆矩阵 | 直接用二阶公式 |
| 求一般逆矩阵 | 对 $[A\mid E]$ 作初等行变换 |
| 解 $AX=B$ | 直接对 $[A\mid B]$ 作行变换 |
| 求矩阵的秩 | 化为阶梯形，数非零行 |
| 讨论含参数方程组 | 化简到关键主元，再分类讨论是否为零 |
| 证明秩等式 | 用秩不等式分别给出上、下界 |
| 大矩阵含低秩乘积 | 尝试 $\det(E+AB)=\det(E+BA)$ 降阶 |

### 15.2 易错点速查

| 容易写错的说法 | 正确理解 |
| --- | --- |
| $AB=BA$ 总成立 | 乘法一般不交换，甚至两边类型可能不同 |
| $AB=O$ 意味着某因子为零 | 非零矩阵可能相乘为零 |
| $A\ne O$ 就可以约掉 $A$ | 需要可逆等足够条件，并从正确一侧乘逆 |
| $(A+B)^2=A^2+2AB+B^2$ 总成立 | 一般中间项为 $AB+BA$ |
| $(AB)^T=A^TB^T$ | 正确顺序是 $B^TA^T$ |
| $(AB)^{-1}=A^{-1}B^{-1}$ | 正确顺序是 $B^{-1}A^{-1}$ |
| $\det(cA)=c\det A$ | 对 $n$ 阶方阵是 $c^n\det A$ |
| $\det(A+B)=\det A+\det B$ | 行列式只对单独一行或一列有线性性 |
| 初等变换总保持行列式 | 倍加保持值；换行变号；倍乘改变倍数 |
| $\det A=0$ 就说明 $Ax=b$ 无解 | 也可能有无穷多解，需比较两个秩 |
| 伴随矩阵就是代数余子式原样排列 | 还要转置，即 $A_{ij}$ 放在 $(j,i)$ 位置 |
| 分块求逆就是每块求逆 | 只在分块对角等特定结构下成立 |
| 所有二阶子式为零，所以矩阵为零 | 只能推出秩至多为 $1$ |
| $r(A+B)=r(A)+r(B)$ | 一般只有上界，可能出现抵消 |
| $r(AB)=\min(r(A),r(B))$ | 一般只有上界，乘积的秩还可能更小 |
| 增广矩阵任意列变换都保持原解 | 列变换涉及未知量换元，需要跟踪其含义 |

## 16. 复习清单与自测

### 复习清单

- [ ] 能判断矩阵运算是否有意义，并写出结果的维度。
- [ ] 能用反例说明矩阵乘法不交换、有零因子、不能随意约去非零因子。
- [ ] 能解释 $Ax$ 是列向量的线性组合，以及 $AB$ 表示先 $B$ 后 $A$。
- [ ] 会用结构求矩阵幂，使用二项式展开前会检查可交换性。
- [ ] 掌握转置、对称与反对称分解、迹的常用性质。
- [ ] 能区分矩阵、行列式、子式、余子式与代数余子式。
- [ ] 能正确记录行列式在初等变换下的变化。
- [ ] 会计算 Vandermonde 行列式，理解其与插值唯一性的关系。
- [ ] 会使用 Cramer 法则，并说明适用条件。
- [ ] 能构造伴随矩阵，解释 $AA^*=(\det A)E$。
- [ ] 会用单边检验、行列式或秩证明可逆。
- [ ] 能正确进行分块乘法、转置与特殊结构的求逆。
- [ ] 理解初等矩阵左乘、右乘分别对应怎样的操作。
- [ ] 能用阶梯形求秩，并比较系数矩阵与增广矩阵的秩。
- [ ] 会用 $[A\mid E]$ 求逆、用 $[A\mid B]$ 解矩阵方程。
- [ ] 会用分块消元推导降阶公式，核对乘法顺序与可逆条件。
- [ ] 能用秩不等式证明伴随矩阵的秩公式与幂等矩阵的秩等式。

### 自测（补充）

1. 若 $A$ 为 $2\times3$、$B$ 为 $3\times4$，$AB$ 和 $BA$ 哪些有意义？结果是什么大小？
2. 令 $N=\begin{bmatrix}0&1\\0&0\end{bmatrix}$。计算 $N^2$ 和 $(E+N)^k$，$k$ 为非负整数。
3. 对 $A=\begin{bmatrix}1&2\\3&4\end{bmatrix}$，写出其对称部分与反对称部分。
4. $A$ 是四阶方阵，$\det A=3$。求 $\det(2A)$、$\det(A^{-1})$、$\det(A^*)$。
5. 计算 $\begin{vmatrix}1&1&1\\1&2&4\\1&4&16\end{vmatrix}$。
6. 若 $A^2-3A+2E=O$，证明 $A$ 可逆并写出 $A^{-1}$。
7. 令 $P=\begin{bmatrix}1&c\\0&1\end{bmatrix}$。说明 $PA$ 与 $AP$ 分别对二阶矩阵 $A$ 作什么操作。
8. 四个未知数的方程组，若 $r(A)=r([A\mid b])=2$，有几个自由变量？若增广矩阵的秩改为 $3$ 呢？
9. 求 $\begin{bmatrix}1&2&3\\2&4&6\end{bmatrix}$ 的秩，并给出一个 $\alpha\beta^T$ 分解。
10. 若 $B,C$ 可逆，写出 $\begin{bmatrix}B&D\\O&C\end{bmatrix}$ 的逆矩阵。
11. $A$ 为 $3\times5$、$B$ 为 $5\times4$，且 $r(A)=3,r(B)=4$。$r(AB)$ 可以是多少？
12. 若 $u=(1,2,3)^T$，求 $\det(E_3+uu^T)$。
13. 五阶方阵 $A$ 的秩分别为 $5,4,3$ 时，$A^*$ 的秩分别是多少？
14. 对方程组 $x_1+x_2=0$、$x_1+x_2=1$，找出列向量 $y$，使 $A^Ty=0$、$b^Ty=1$。

### 答案简要

1. $AB$ 有意义，是 $2\times4$；$BA$ 的内侧维度 $4\ne2$，没有意义。
2. $N^2=O$，所以 $(E+N)^k=E+kN=\begin{bmatrix}1&k\\0&1\end{bmatrix}$。
3. 对称部分为 $\begin{bmatrix}1&5/2\\5/2&4\end{bmatrix}$；反对称部分为 $\begin{bmatrix}0&-1/2\\1/2&0\end{bmatrix}$。
4. 分别为 $2^4\cdot3=48$、$1/3$、$3^3=27$。
5. 是 $x_1=1,x_2=2,x_3=4$ 的 Vandermonde 行列式，值为 $(2-1)(4-1)(4-2)=6$。
6. $A(3E-A)=2E$，故 $A^{-1}=(3E-A)/2$。
7. $PA$ 对应 $R_1\leftarrow R_1+cR_2$；$AP$ 对应 $C_2\leftarrow C_2+cC_1$。
8. 前者相容，有 $4-2=2$ 个自由变量；后者两秩不同，无解。
9. 秩为 $1$；可取 $\alpha=(1,2)^T$、$\beta=(1,2,3)^T$。
10. $\begin{bmatrix}B^{-1}&-B^{-1}DC^{-1}\\O&C^{-1}\end{bmatrix}$。
11. 由 $3+4-5\le r(AB)\le\min(3,4)$，可为 $2$ 或 $3$。
12. $1+u^Tu=1+1+4+9=15$。
13. 分别为 $5,1,0$。
14. $A=\begin{bmatrix}1&1\\1&1\end{bmatrix}$、$b=(0,1)^T$，取 $y=(-1,1)^T$。它对应“第二式减第一式”，得到矛盾 $0=1$。
