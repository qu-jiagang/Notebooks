# 基本思想

有限元方法首先将偏微分方程的弱形式构造出来，然后通过对这个弱形式进行适当的近似和离散化，得到一个线性代数方程组。求解这个线性代数方程组得到的解通常被认为是原始偏微分方程的近似解，即它是原始方程的弱解。

## 基函数（shape function）和试探函数（trail/test function）

### 试探函数

试探函数在有限元的上下文中用于构建方程的弱形式。当我们将一个微分方程转化为其弱形式时，通常涉及到将微分方程乘以一个试探函数然后在整个域上积分。这样做的目的是为了确保方程的积分形式在更广泛的函数空间内成立，降低对解的光滑性要求。

### 基函数

基函数，有时也称为形状函数，是在每个单元内用来表示未知函数（如位移、温度等）的近似解的函数。一般情况下，基函数数量等于节点数

### 关系

在有限元分析的实际应用中，基函数通常同时用作试探函数，这种方法被称为 Galerkin 方法。这意味着，在构造弱形式和求解方程时，同一组基函数既用来近似未知函数，也用来作为试探函数。这种做法简化了问题的表达和求解过程，因为只需要考虑一组函数的性质。

## 示例：泊松方程

### 理论推导

$$-\frac{d^2u}{dx^2}=f(x)$$

给定$f(x)=-1$，且采用 Dirichlet 边界$u(0)=0,u(1)=0$
选取试探函数$\phi_i(x)$，可得方程弱解形式

$$\int_0^1-\frac{d^2u}{dx^2}\phi_i(x)dx = \int_0^1f(x)\phi_i(x)dx$$

对方程左边采用分部积分得到

$$
    \int_0^1\frac{d^2u}{dx^2}\phi_i(x)dx = \left[ \frac{du}{dx}\phi_i(x) \right]_0^1 - \int_0^1\frac{du}{dx}\frac{d\phi_i}{dx}dx
$$

由于边界条件为 Dirichlet 边界，因此边界项为 0，得到原方程的弱解形式

$$
    \int_0^1\frac{du}{dx}\frac{d\phi_i}{dx}dx = \int_0^1f(x)\phi_i(x)dx
$$

选取基函数于试探函数相同$\phi_i(x)$，可得

$$
    u(x)\approx \sum_{i=1}^n u_i\phi_i(x)
$$

其中$u_i$是未知系数
代入得

$$
    \sum_{i=1}^n u_i\int_0^1\frac{d\phi_i}{dx}\frac{d\phi_j}{dx}dx = \int_0^1f(x)\phi_j(x)dx
$$

令：

$$
K_{ij} = \int_0^1\frac{d\phi_i}{dx}\frac{d\phi_j}{dx}dx \\
F_j = \int_0^1f(x)\phi_j(x)dx
$$

则上式简化为

$$Ku=F$$

### 具体示例

在计算域[0,1]上均匀划分 3 个网格，网格大小$h=1/2=0.5$

基函数选取为：

$$
\phi_1(x) = 1 - \frac{x}{h} \\
\phi_2(x) = \frac{x}{h}
$$

将原方程的弱解形式代入到每个单元中

$$
\sum_{i=1}^n u_i\int_{kh}^{(k+1)h}\frac{d\phi_i}{dx}\frac{d\phi_j}{dx}dx = \int_{kh}^{(k+1)h}f(x)\phi_j(x)dx
$$

这里 i 是基函数，j 是试探函数

$$
    \begin{bmatrix}
    \int_{kh}^{(k+1)h}\frac{d\phi_1}{dx}\frac{d\phi_1}{dx}dx, \int_{kh}^{(k+1)h}\frac{d\phi_2}{dx}\frac{d\phi_1}{dx}dx \\
    \int_{kh}^{(k+1)h}\frac{d\phi_1}{dx}\frac{d\phi_2}{dx}dx, \int_{kh}^{(k+1)h}\frac{d\phi_2}{dx}\frac{d\phi_2}{dx}dx
    \end{bmatrix}
    \begin{bmatrix}
    u_1 \\
    u_2
    \end{bmatrix}
    =
    \begin{bmatrix}
    \int_{kh}^{(k+1)h}f(x)\phi_1(x)dx \\
    \int_{kh}^{(k+1)h}f(x)\phi_2(x)dx
    \end{bmatrix}
$$

这里

$$
\frac{d\phi_1}{dx} = -\frac{1}{h},\ \frac{d\phi_2}{dx} = \frac{1}{h}
$$

$$f(x)=-2$$

获得单元刚度矩阵$K_{local}$和装配向量$F_{local}$

$$
    K_{local} =
    \begin{bmatrix}
    \int_{kh}^{(k+1)h}\frac{d\phi_1}{dx}\frac{d\phi_1}{dx}dx, \int_{kh}^{(k+1)h}\frac{d\phi_2}{dx}\frac{d\phi_1}{dx}dx \\
    \int_{kh}^{(k+1)h}\frac{d\phi_1}{dx}\frac{d\phi_2}{dx}dx, \int_{kh}^{(k+1)h}\frac{d\phi_2}{dx}\frac{d\phi_2}{dx}dx
    \end{bmatrix}
    =
    \frac{1}{h}
    \begin{bmatrix}
    1, -1 \\
    -1, 1
    \end{bmatrix}
$$

$$
    F_{local}=
    \begin{bmatrix}
    \int_{kh}^{(k+1)h}f(x)\phi_1(x)dx \\
    \int_{kh}^{(k+1)h}f(x)\phi_2(x)dx
    \end{bmatrix}
    =
    -2\cdot \frac{h}{2}
    \begin{bmatrix}
    1 \\
    1
    \end{bmatrix}
    =
    -0.5
    \begin{bmatrix}
    1 \\
    1
\end{bmatrix}
$$

组装全局刚度矩阵和装备向量
将所有单元刚度矩阵和专配向量按节点所在位置叠加即可

$$
K= \frac{1}{h}
\begin{bmatrix}
1, -1, 0 \\
-1, 2, -1 \\
0, 1, -1
\end{bmatrix}
$$

$$
F=
\begin{bmatrix}
-0.5 \\
-1.0 \\
-0.5
\end{bmatrix}
$$

应用边界条件，即$u(0)=0,u(1)=0$，需要从全局矩阵中移除对应的固定节点的行和列

$$K=[4],F=[-1.0]$$

求解得到线性方程组的解

$$u=[-0.25]$$

填充边界条件后

$$
u=
\begin{bmatrix}
0 \\
-0.25 \\
0
\end{bmatrix}
$$

### 代码

```
import numpy as np
import matplotlib.pyplot as plt

# 定义问题参数
n_nodes = 11  # 节点数
element_indices = []
for i in range(n_nodes-1):
    element_indices.append((i, i))
length = 1.0  # 定义域长度
f_value = -2  # 源项常数

# 划分网格
h = length / (n_nodes - 1)
node_coords = np.linspace(0, 1, n_nodes)

# 初始化刚度矩阵和载荷向量
K = np.zeros((n_nodes, n_nodes))
F = np.zeros(n_nodes)

# 构建刚度矩阵和载荷向量
for (i, j) in element_indices:
    k_local = np.array([[1, -1], [-1, 1]]) * (1 / h)
    f_local = np.array([f_value * h / 2, f_value * h / 2])
    K[i:i+2, j:j+2] += k_local
    F[i:i+2] += f_local

# 应用边界条件
K = K[1:-1, 1:-1]
F = F[1:-1]

# 求解线性系统
u = np.linalg.solve(K, F)

# 填充边界条件
u = np.hstack([0, u, 0])

# 绘制结果
plt.plot(node_coords, u, '-o', label='Numerical Solution')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Finite Element Solution of 1D Poisson Equation')
plt.legend()
plt.grid(True)
plt.show()
```
