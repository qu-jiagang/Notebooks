# 特征线方法

## 基础
### 以连续性方程为例

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol u) = 0
$$

这里先只考虑一维情况，且速度为常数$a$， 也就是

$$
\frac{\partial \rho}{\partial t} + a\frac{\partial \rho}{\partial x} = 0
$$

对于Cauchy问题，初始条件为 $\rho(x, 0) = \rho_0(x)$，

1. 观察到，当 $\frac{dx}{dt}=a$ 以上PDE可以由转换成ODE
   
$$
\frac{d\rho}{dt} = \frac{\partial \rho}{\partial t} + \frac{\partial \rho}{\partial x}\frac{\partial x}{\partial t} = 0
$$
   
   我们把 $\frac{dx}{dt}=a$ 这个方程成为`特征线方程`，对其进行求解，它的解就叫特征线 $x(t, c) = at + c$。
2. 求解ODE，可以得到解为 $\rho = \rm Const$，也就是说 $\rho$ 是沿着特征线是不变的
3. 根据初始条件 $\rho(x, 0) = \rho_0(x)$，在沿着特征线方向上有：
   
$$
\rho(x(t,c), t) = \rho(x(0,c),0) = \rho(x(c), 0) = \rho_0(c)
$$

将 $x=at+c$ 代入可以得到

$$
\rho(x, t) = \rho_0(x-at)
$$

### 总结
我们整理下特征线法求解偏微分方程的步骤

1. 求特征线 $x = x(t, c)$

2. 沿特征线将原方程化为关于 $\rho = \rho(x(t,c), c t)$ 的ODE方程，并求出 $\rho = u(t,c)$

3. 从特征线方程解出 $c = \phi(x ,t)$，则所求的解为 $\rho = u(t, \phi(x, t))$

