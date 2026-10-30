# Assignment 6 学习笔记：分块矩阵的乘法与求逆

课程：MATH1205 线性代数  
题目：[原题 PDF](assignment6.pdf)  
提交版：[中文解答 PDF](assignment6-sol.pdf)  
源文件：[LaTeX 主文件](solution/main.tex)  
配套讲义：[§2.4 分块矩阵](../../lectures/LA2.4.pdf)

原题包含 **40(2)、41(1)、41(2)、55(1)、55(2)**。第 40 题截图只有第 (2) 小问的乘法表达式，按矩阵乘积计算。以下 $E_k$ 是 $k$ 阶单位矩阵，$O$ 是相应尺寸的零矩阵。

## 1. 这次作业要学会什么

分块的目的，是把“大矩阵的运算”变成“少量小矩阵的运算”。先观察零元素的位置，再决定怎么划分。

| 题号 | 结构 | 方法 | 最容易出错的地方 |
| --- | --- | --- | --- |
| 40(2) | 两个矩阵都有大片零元素 | 选择相容的分块，只计算两个非零块乘积 | 误以为每个块必须是方阵 |
| 41(1) | 块对角矩阵 | 对角块分别求逆 | 忘记中间的标量块 $2$ 也要取逆 |
| 41(2) | 非零块位于右上、左下，阶数分别为 2、3 | 非零块求逆后交换位置 | 原矩阵与逆矩阵的行、列分块尺寸不同 |
| 55(1) | 带权的循环移位 | 分成一个对角块和一个标量块 | $a_1^{-1}$ 应在右上角 |
| 55(2) | 外层是交叉块，内层有反对角矩阵 | 先处理外层，再求各块的逆 | 反对角线系数取倒数后还要倒序 |

## 2. 分块运算的两条基本规则

### 2.1 相乘时检查“左边的列”与“右边的行”

若

$$
A=\begin{pmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{pmatrix},
\qquad
B=\begin{pmatrix}B_{11}&B_{12}\\B_{21}&B_{22}\end{pmatrix},
$$

且 $A$ 的各列块宽度分别等于 $B$ 的对应行块高度，则

$$
AB=
\begin{pmatrix}
A_{11}B_{11}+A_{12}B_{21}&A_{11}B_{12}+A_{12}B_{22}\\
A_{21}B_{11}+A_{22}B_{21}&A_{21}B_{12}+A_{22}B_{22}
\end{pmatrix}.
$$

乘法的外形与二阶矩阵相同，但每个“元素”现在是一个矩阵，因此：

- 要逐项检查尺寸是否允许相乘、相加。
- 乘法顺序不能改，通常 $A_{ij}B_{jk}\ne B_{jk}A_{ij}$。
- 子块可以是长方形，整体的行分法和列分法也可以不同。
- 只把 $AB$ 的结果写对还不够，分块过程必须与实际行列位置对应。

### 2.2 两类常用的分块求逆

**块对角矩阵。** 若 $P,Q$ 都可逆，则

$$
\begin{pmatrix}P&O\\O&Q\end{pmatrix}^{-1}
=\begin{pmatrix}P^{-1}&O\\O&Q^{-1}\end{pmatrix}.
$$

乘回去直接得到 $\operatorname{diag}(E,E)$。有更多对角块时同理。

**非零块在右上、左下的矩阵。** 若 $P$ 为 $p$ 阶可逆矩阵，$Q$ 为 $q$ 阶可逆矩阵，则

$$
M=
\begin{pmatrix}
O_{p\times q}&P\\
Q&O_{q\times p}
\end{pmatrix},
\qquad
\boxed{
M^{-1}=
\begin{pmatrix}
O_{q\times p}&Q^{-1}\\
P^{-1}&O_{p\times q}
\end{pmatrix}}.
$$

| 矩阵 | 行的分块尺寸 | 列的分块尺寸 |
| --- | --- | --- |
| $M$ | $p+q$ | $q+p$ |
| $M^{-1}$ | $q+p$ | $p+q$ |

为什么要交换两个逆块的位置？直接相乘：

$$
MM^{-1}
=\begin{pmatrix}PP^{-1}&O\\O&QQ^{-1}\end{pmatrix}
=\operatorname{diag}(E_p,E_q),
$$

$$
M^{-1}M
=\begin{pmatrix}Q^{-1}Q&O\\O&P^{-1}P\end{pmatrix}
=\operatorname{diag}(E_q,E_p).
$$

两者都是 $(p+q)$ 阶单位矩阵，只是分块边界不同。这个公式不要求 $p=q$，也不要求 $PQ=QP$；当 $p\ne q$ 时，$PQ$ 甚至没有定义。

也可以从方程理解。设 $x=(u^T,v^T)^T$，其中 $u$ 有 $q$ 个分量，$v$ 有 $p$ 个分量，则

$$
Mx=\begin{pmatrix}Pv\\Qu\end{pmatrix}.
$$

若输出是 $(y^T,z^T)^T$，恢复输入就需要 $u=Q^{-1}z$、$v=P^{-1}y$，所以逆块会交换位置。

## 3. 第 40(2) 题：让零块替你省计算

原式为 $AB$，其中

$$
A=\begin{pmatrix}
1&-1&0&0\\3&-1&0&0\\0&1&0&0\\0&0&2&-1
\end{pmatrix},
\quad
B=\begin{pmatrix}
1&0&0&0\\-1&0&0&0\\0&1&3&-1\\0&2&1&4
\end{pmatrix}.
$$

观察 $A$：前三行只用到前两列，最后一行只用到后两列。因此按行 $3+1$、列 $2+2$ 分块：

$$
A=\begin{pmatrix}C&O\\O&D\end{pmatrix},
\qquad
C=\begin{pmatrix}1&-1\\3&-1\\0&1\end{pmatrix},
\quad D=\begin{pmatrix}2&-1\end{pmatrix}.
$$

为了匹配 $A$ 的列分法，$B$ 的行必须也按 $2+2$ 分。再按列 $1+3$ 分块：

$$
B=\begin{pmatrix}F&O\\O&G\end{pmatrix},
\qquad
F=\begin{pmatrix}1\\-1\end{pmatrix},
\quad G=\begin{pmatrix}1&3&-1\\2&1&4\end{pmatrix}.
$$

于是只剩

$$
CF=
\begin{pmatrix}
1\cdot1+(-1)(-1)\\
3\cdot1+(-1)(-1)\\
0\cdot1+1(-1)
\end{pmatrix}
=\begin{pmatrix}2\\4\\-1\end{pmatrix},
$$

$$
DG=(2\cdot1-2,\ 2\cdot3-1,\ 2(-1)-4)=(0,5,-6).
$$

最终

$$
\boxed{AB=
\begin{pmatrix}
2&0&0&0\\
4&0&0&0\\
-1&0&0&0\\
0&0&5&-6
\end{pmatrix}}.
$$

**自检：** $AB$ 的第 2 列应为零，因为 $B$ 第 2 列为 $(0,0,1,2)^T$，故它等于 $A$ 的第 3 列加上第 4 列的两倍，恰好为零。

这里 $C$ 是 $3\times2$、$F$ 是 $2\times1$，所以 $CF$ 是 $3\times1$。这正说明“分块乘法”不要求块是方阵。

## 4. 第 41(1) 题：块对角矩阵逐块求逆

原矩阵是

$$
M_1=\operatorname{diag}(P,2,Q),
\quad
P=\begin{pmatrix}2&3\\-3&-5\end{pmatrix},
\quad
Q=\begin{pmatrix}8&5\\3&2\end{pmatrix}.
$$

二阶求逆公式为

$$
\begin{pmatrix}a&b\\c&d\end{pmatrix}^{-1}
=\frac1{ad-bc}\begin{pmatrix}d&-b\\-c&a\end{pmatrix},
\qquad ad-bc\ne0.
$$

因此

$$
|P|=-10+9=-1,\qquad
P^{-1}=\frac1{-1}\begin{pmatrix}-5&-3\\3&2\end{pmatrix}
=\begin{pmatrix}5&3\\-3&-2\end{pmatrix},
$$

$$
|Q|=16-15=1,\qquad
Q^{-1}=\begin{pmatrix}2&-5\\-3&8\end{pmatrix}.
$$

连同标量块 $2^{-1}=\frac12$，

$$
\boxed{
M_1^{-1}=
\begin{pmatrix}
5&3&0&0&0\\
-3&-2&0&0&0\\
0&0&\frac12&0&0\\
0&0&0&2&-5\\
0&0&0&-3&8
\end{pmatrix}}.
$$

**为什么不算五阶行列式和全部代数余子式？** 因为块对角结构已经把问题分成了两个二阶求逆和一个数的倒数。检验也只需检验 $PP^{-1}=E_2$、$QQ^{-1}=E_2$ 和 $2\cdot\frac12=1$。

## 5. 第 41(2) 题：逆矩阵会交换块的尺寸

原矩阵按行 $2+3$、列 $3+2$ 分块后是

$$
M_2=\begin{pmatrix}O_{2\times3}&R\\S&O_{3\times2}\end{pmatrix},
\quad
R=\begin{pmatrix}1&2\\2&3\end{pmatrix},
\quad
S=\begin{pmatrix}1&1&0\\0&1&1\\0&0&1\end{pmatrix}.
$$

先算小块：

$$
R^{-1}=\begin{pmatrix}-3&2\\2&-1\end{pmatrix}.
$$

对于 $S$，有两种理解方式。

**方法一：从方程回代。** 若 $Sx=y$，则

$$
x_3=y_3,\qquad x_2=y_2-y_3,\qquad
x_1=y_1-y_2+y_3.
$$

所以

$$
S^{-1}=\begin{pmatrix}1&-1&1\\0&1&-1\\0&0&1\end{pmatrix}.
$$

**方法二：联系上次作业的有限几何级数。** 写成

$$
S=E_3+N,\qquad
N=\begin{pmatrix}0&1&0\\0&0&1\\0&0&0\end{pmatrix},
\qquad N^3=O.
$$

因为

$$
(E_3+N)(E_3-N+N^2)=E_3+N^3=E_3,
$$

故 $S^{-1}=E_3-N+N^2$。特别注意右上角的 $1$ 来自 $N^2$，不能漏掉。

代入交叉块求逆公式：

$$
M_2^{-1}=
\begin{pmatrix}
O_{3\times2}&S^{-1}\\R^{-1}&O_{2\times3}
\end{pmatrix}
=
\boxed{\begin{pmatrix}
0&0&1&-1&1\\
0&0&0&1&-1\\
0&0&0&0&1\\
-3&2&0&0&0\\
2&-1&0&0&0
\end{pmatrix}}.
$$

**尺寸自检：** 逆矩阵上面有 3 行，下面有 2 行；左边有 2 列，右边有 3 列。若仍沿用原矩阵的分块尺寸，会连块都放不进去。

## 6. 第 55(1) 题：带权循环移位如何倒过来

设原矩阵为 $T$。它的非零元素为

$$
T_{i,i+1}=a_{i+1}\quad(1\le i<n),\qquad
T_{n,1}=a_1,
$$

其中每个 $a_i\ne0$。

### 6.1 按题意用分块求逆

令 $D=\operatorname{diag}(a_2,\ldots,a_n)$，则

$$
T=
\begin{pmatrix}
O_{(n-1)\times1}&D\\
a_1&O_{1\times(n-1)}
\end{pmatrix}.
$$

它与第 41(2) 题结构相同，只是两个非零块的阶数变成了 $n-1$ 与 $1$。因此

$$
\boxed{
T^{-1}=
\begin{pmatrix}
O_{1\times(n-1)}&a_1^{-1}\\
D^{-1}&O_{(n-1)\times1}
\end{pmatrix}}.
$$

写成元素形式最不容易误解省略号：

$$
(T^{-1})_{ij}=
\begin{cases}
1/a_1,&(i,j)=(1,n),\\
1/a_i,&2\le i\le n,\ j=i-1,\\
0,&\text{其余位置}.
\end{cases}
$$

也就是右上角放 $1/a_1$，次对角线从上到下放 $1/a_2,\ldots,1/a_n$。

### 6.2 从坐标变化理解位置

$$
Tx=
\begin{pmatrix}
a_2x_2\\a_3x_3\\\vdots\\a_nx_n\\a_1x_1
\end{pmatrix}.
$$

若 $Tx=y$，则

$$
x_1=\frac{y_n}{a_1},\qquad
x_i=\frac{y_{i-1}}{a_i}\quad(2\le i\le n).
$$

第一行必须从 $y$ 的最后一个分量恢复 $x_1$，所以 $1/a_1$ 位于右上角，不能留在左下角。

例如 $n=3$ 时：

$$
T=\begin{pmatrix}0&a_2&0\\0&0&a_3\\a_1&0&0\end{pmatrix},
\qquad
T^{-1}=\begin{pmatrix}0&0&1/a_1\\1/a_2&0&0\\0&1/a_3&0\end{pmatrix}.
$$

### 6.3 另一种自检：走完一圈

每做一次 $T$，坐标位置循环移动一次并乘上相应权重；做 $n$ 次回到原位置，且每个坐标都乘过全部权重，所以

$$
T^n=(a_1a_2\cdots a_n)E_n.
$$

这给出另一个等价表达式：

$$
T^{-1}=\frac{T^{n-1}}{a_1a_2\cdots a_n}.
$$

它适合验证结构，但本题使用分块法能更直接地写出每个元素。

## 7. 第 55(2) 题：先求外层，再处理反对角线

原矩阵为

$$
K=\begin{pmatrix}O&A_1\\A_2&O\end{pmatrix},
$$

其中 $A_1$ 的反对角线从右上到左下依次为 $a_1,\ldots,a_n$，其余元素为零，而

$$
A_2=\operatorname{diag}(a_{n+1},\ldots,a_{2n}).
$$

### 7.1 外层结构

先直接写出

$$
\boxed{K^{-1}=
\begin{pmatrix}O&A_2^{-1}\\A_1^{-1}&O\end{pmatrix}}.
$$

因为所有系数非零，$A_2$ 的逆为

$$
A_2^{-1}=\operatorname{diag}\left(
\frac1{a_{n+1}},\ldots,\frac1{a_{2n}}\right).
$$

剩下的问题是求 $A_1^{-1}$。

### 7.2 内层反对角矩阵

令 $J$ 为倒序置换矩阵，即 $J_{i,n+1-i}=1$，其余元素为零。它把一个向量的坐标顺序颠倒。颠倒两次恢复原状，所以

$$
J^2=E_n,\qquad J^{-1}=J.
$$

令 $D_1=\operatorname{diag}(a_1,\ldots,a_n)$，则

$$
A_1=D_1J,\qquad
\boxed{A_1^{-1}=JD_1^{-1}}.
$$

这里用到了 $(AB)^{-1}=B^{-1}A^{-1}$。乘积求逆要反序，所以不能写成 $D_1^{-1}J$。

展开后，$A_1^{-1}$ 仍只有反对角线非零，但从右上到左下依次为

$$
\boxed{\frac1{a_n},\frac1{a_{n-1}},\ldots,\frac1{a_1}}.
$$

精确的元素公式为

$$
(A_1^{-1})_{ij}=
\begin{cases}
1/a_{n+1-i},&j=n+1-i,\\
0,&j\ne n+1-i.
\end{cases}
$$

**用 $n=3$ 看清楚倒序：**

$$
A_1=\begin{pmatrix}0&0&a_1\\0&a_2&0\\a_3&0&0\end{pmatrix},
\qquad
A_1^{-1}=\begin{pmatrix}0&0&1/a_3\\0&1/a_2&0\\1/a_1&0&0\end{pmatrix}.
$$

例如乘积左上角是 $a_1\cdot(1/a_1)=1$。如果只在原位置取倒数而不倒序，左上角会变成 $a_1/a_3$，一般不等于 1。

所以最后的完整块形式为

$$
\boxed{
K^{-1}=
\begin{pmatrix}
O&\operatorname{diag}(a_{n+1}^{-1},\ldots,a_{2n}^{-1})\\
JD_1^{-1}&O
\end{pmatrix}}.
$$

这是一份完整答案：每个非零块的尺寸和所有元素都已经确定，无需展开成庞大的 $2n$ 阶矩阵。

## 8. 常见误区与检查办法

| 误区 | 为什么不对 | 正确做法 |
| --- | --- | --- |
| 给大矩阵随意画横线和竖线后直接相乘 | 左矩阵列块宽度可能不等于右矩阵行块高度 | 先写尺寸，再做分块乘法 |
| 认为每个块必须为方阵 | 分块只是一种记法，长方形块也可相乘 | 检查内侧维数即可 |
| 把一般矩阵的每个非零元素取倒数 | 元素倒数不是矩阵逆运算 | 检验乘积是否为单位矩阵 |
| 交叉块求逆后仍放在原位置 | 相乘时不能得到所需的对角单位块 | $P^{-1}$ 放左下，$Q^{-1}$ 放右上 |
| 求反对角矩阵的逆时只取倒数 | 逆运算要把输入、输出位置交换回来 | 反对角线上的系数倒序后取倒数 |
| 把 $S^{-1}$ 写成 $E_3-N$ | $(E_3+N)(E_3-N)=E_3-N^2$，而本题 $N^2\ne O$ | 补上 $N^2$ |
| 只凭矩阵稀疏就认定可逆 | 零系数可能导致整行或整列为零 | 检查所需非零块是否可逆 |

## 9. 自测题

先独立尝试，再看答案。

1. 若 $P$ 为 2 阶可逆矩阵，$Q$ 为 4 阶可逆矩阵，$M=\begin{pmatrix}O_{2\times4}&P\\Q&O_{4\times2}\end{pmatrix}$，写出 $M^{-1}$ 的分块尺寸。
2. 求 $\begin{pmatrix}0&2&0\\0&0&3\\5&0&0\end{pmatrix}$ 的逆。
3. 求 $\begin{pmatrix}0&0&2\\0&3&0\\5&0&0\end{pmatrix}$ 的逆，解释它与上一题的区别。
4. 若第 55(1) 题中某个 $a_i=0$，矩阵为什么不可逆？

<details>
<summary>参考答案</summary>

1. 逆矩阵按行 $4+2$、列 $2+4$ 分块：

   $$
   M^{-1}=\begin{pmatrix}O_{4\times2}&Q^{-1}\\P^{-1}&O_{2\times4}\end{pmatrix}.
   $$

2. 答案为

   $$
   \begin{pmatrix}0&0&1/5\\1/2&0&0\\0&1/3&0\end{pmatrix}.
   $$

   原矩阵做循环移位，逆矩阵做反向移位并除以对应权重。

3. 答案为

   $$
   \begin{pmatrix}0&0&1/5\\0&1/3&0\\1/2&0&0\end{pmatrix}.
   $$

   原矩阵做倒序排列并缩放；逆矩阵仍是反对角形，但权重倒序取倒数。

4. 原矩阵每行、每列恰好只有一个可能非零的元素。某个 $a_i=0$ 后，相应行和列都成为零，行列式为零，所以不可逆。

</details>

## 10. 如何编译提交版

在仓库根目录运行以下命令，使用 XeLaTeX 和 TeX Live 自带的 Fandol 中文字体：

~~~bash
mkdir -p tmp/pdfs/assignment6-build
xelatex -interaction=nonstopmode -halt-on-error \
  -output-directory=tmp/pdfs/assignment6-build \
  -jobname=assignment6-sol assignments/assignment6/solution/main.tex
xelatex -interaction=nonstopmode -halt-on-error \
  -output-directory=tmp/pdfs/assignment6-build \
  -jobname=assignment6-sol assignments/assignment6/solution/main.tex
cp tmp/pdfs/assignment6-build/assignment6-sol.pdf \
  assignments/assignment6/assignment6-sol.pdf
~~~

解答的数值矩阵已用精确分数运算核对；两个含参数的结论在提交版中均通过分块乘积证明，并额外检查了 $n=2,\ldots,8$ 的非零有理数样例。样例用于辅助检查下标，不能代替一般证明。
