# POD-Galerkin Model

> update 2024-10-26

Consider a dynamical system which evolves in a Hilbert space H, In particular, for $\boldsymbol u(t)\in \boldsymbol H$, $\boldsymbol u(t)$ satisfies

$$
    \dot{\boldsymbol u} = \boldsymbol X(\boldsymbol u)
$$

$\boldsymbol u(t)$ can be write with Proper Orthgonal Decomposition (POD),

$$
    \boldsymbol u(t) = \sum_k a_k(t) \boldsymbol\phi_k
$$

Inserting this into the dynamical system, we get

$$
    \dot{\boldsymbol u} = \sum_k \dot{a}_k(t) \boldsymbol\phi_k = \boldsymbol X (\boldsymbol u (t))
$$

In order to get the reduced order dynamical system of $a_k(t)$, by projecting the equation into the basis $\boldsymbol \phi_k$.

$$
    \langle\sum_k \dot{a}_k(t) \boldsymbol\phi_k, \boldsymbol\phi_l\rangle = \langle\boldsymbol X (\boldsymbol u (t)), \boldsymbol\phi_l\rangle
$$

Assuming that the basis $\boldsymbol \phi_k$ are orthgonal, 

$$
    \dot{a}_l(t) = \langle\boldsymbol X (\boldsymbol u (t)), \boldsymbol\phi_l\rangle
$$

## Governing Equations
3-dimensional Navier-Stokes equations is given by

$$
    \frac{\partial \rho}{\partial t} + (\boldsymbol u \cdot \nabla)\rho + \nabla \cdot \boldsymbol{u} = 0
$$

$$
    \rho \left(\frac{\partial \boldsymbol u}{\partial t} + (\boldsymbol u \cdot \nabla) \boldsymbol u \right ) = -\nabla p + \nu \nabla^2\boldsymbol u
$$

$$
    \rho\left( \frac{\partial T}{\partial t} + (\gamma - 1)T\nabla\cdot \boldsymbol u \right) = \nabla\cdot\boldsymbol (\rho \boldsymbol u T) + \frac{\gamma\nu}{Pr}\nabla^2T
$$

These equations are written in the form:

$$
    \frac{\partial \rho}{\partial t} = R(\boldsymbol q) , \quad
    \rho \frac{\partial u_i}{\partial t} = \boldsymbol U_1(\boldsymbol q) , \quad
    \rho \frac{\partial T}{\partial t} = \Theta(\boldsymbol q)
$$

where $\boldsymbol q = (\rho, u_1, u_2, u_3, T)$, $\boldsymbol R, \boldsymbol U_i, \boldsymbol \Theta$ are the nonlinear differential operators given by

$$
    \boldsymbol R(\boldsymbol q) = -(\boldsymbol u \cdot \nabla)\rho - \nabla \cdot \boldsymbol{u}
$$

$$
    \boldsymbol U_i(\boldsymbol q) = -\rho (\boldsymbol u \cdot \nabla) \boldsymbol u -\nabla p + \nu \nabla^2\boldsymbol u
$$

$$
    \boldsymbol \Theta(\boldsymbol q) = -(\gamma - 1)\rho T\nabla\cdot u + \nabla\cdot\boldsymbol (\rho \boldsymbol u T) + \frac{\gamma\nu}{Pr}\nabla^2T
$$

written more concisely as 

$$
    \boldsymbol A (\boldsymbol q)\dot{\boldsymbol q} = \boldsymbol f(\boldsymbol q)
$$

where $\boldsymbol A = diag(1, \rho, \rho, \rho, \rho)$, $\boldsymbol f (\boldsymbol q) = (\boldsymbol R(\boldsymbol q), \boldsymbol U_1(\boldsymbol q), \boldsymbol U_2(\boldsymbol q), \boldsymbol U_3(\boldsymbol q), \boldsymbol \Theta(\boldsymbol q))$, 

$\boldsymbol A$ is affine in $\boldsymbol q$, and may be written 

$$
    \boldsymbol A (\boldsymbol q) = \boldsymbol B + \boldsymbol L (\boldsymbol q) = diag(1, 0, 0, 0, 0) + diag(0, \rho, \rho, \rho, \rho)
$$

Furthermore, $\boldsymbol f$ is cubic and may be written

$$
    \boldsymbol f(\boldsymbol q) = \boldsymbol f_1(\boldsymbol q) + \boldsymbol f_2(\boldsymbol q, \boldsymbol q) + \boldsymbol f_3(\boldsymbol q, \boldsymbol q, \boldsymbol q) 
$$

where $\boldsymbol f_{1,2,3}$ are multilinear functions (linear in each argument).

## Galerkin projection

$$
    \boldsymbol f_1(\boldsymbol q) = 
    \begin{bmatrix}
        u_x + v_y + w_z \\
        - p_x + \nu(u_{xx} + u_{yy} + u_{zz}) \\
        - p_y + \nu(v_{xx} + v_{yy} + v_{zz}) \\
        - p_z + \nu(w_{xx} + w_{yy} + w_{zz}) \\
        \frac{\gamma\nu}{Pr}(T_{xx} + T_{yy} + T_{zz})
    \end{bmatrix}
$$

$$
    \boldsymbol f_2(\boldsymbol q^1, \boldsymbol q^2) = 
    \begin{bmatrix}
        u^1\rho^2_x + v^1\rho^2_y + w^1\rho^2_z \\
        0 \\
        0 \\
        0 \\
        0
    \end{bmatrix}
$$

$$
    \boldsymbol f_3(\boldsymbol q^1, \boldsymbol q^2, \boldsymbol q^3) = 
    \begin{bmatrix}
        0 \\
        -\rho^1 (u^2 u^3_x + v^2 u^3_y + w^2 u^3_z) \\
        -\rho^1 (u^2 v^3_x + v^2 v^3_y + w^2 v^3_z) \\
        -\rho^1 (u^2 w^3_x + v^2 w^3_y + w^2 w^3_z) \\
        -(\gamma - 1)\rho^1T^2(u^3_x + v^3_y + w^3_z) + \rho^1_x u^2T^3 + \rho^1 u^2_x T^3 + \rho^1 u^2 T^3_x + \rho^1_y v^2T^3 + \rho^1 v^2_y T^3 + \rho^1 v^2 T^3_y + \rho^1_z w^2T^3 + \rho^1 w^2_z T^3 + \rho^1 w^2 T^3_z
    \end{bmatrix}
$$

let $\boldsymbol \phi_k$ be a basis for the function space containing $\boldsymbol q$, 

$$
    \boldsymbol q(x, t) = \bar{\boldsymbol q}(x) + \sum_k a_k(t)\boldsymbol\phi_k
$$

Inserting this expression into the governing equation,

$$
    \left[\boldsymbol B + \boldsymbol L\left( \sum_l a_l\boldsymbol\phi_l\right)\right] \sum_k\dot{a}_k\boldsymbol\phi_k = \boldsymbol f (\boldsymbol q)
$$

Taking an inner product with $\boldsymbol\phi_j$ then gives

$$
    \sum_k \dot{a}_k\left( \langle\boldsymbol\phi_j, \boldsymbol B\phi_k\rangle + \sum_l a_l \langle\boldsymbol\phi_j, \boldsymbol L(\boldsymbol\phi_l)\phi_k\rangle \right) = \langle\boldsymbol\phi_j, \boldsymbol f(\boldsymbol q)\rangle
$$

The equation above can be written in matrix form as

$$
    \boldsymbol M (\boldsymbol a) \dot{\boldsymbol a} = \boldsymbol F (\boldsymbol a)
$$

where $\boldsymbol a = (a_1, a_2, ..., a_n)$, and

$$
    M_{kl} = \langle \boldsymbol\phi_k, \boldsymbol B\phi_l \rangle + 
    \sum_l a_l \langle\boldsymbol\phi_k, \boldsymbol L(\boldsymbol\phi_l)\phi_k\rangle
$$

$$
    F_k = \langle\boldsymbol\phi_k, \boldsymbol f(\boldsymbol q)\rangle = 
    \sum_l a_l \langle\boldsymbol\phi_k, \boldsymbol f_1(\boldsymbol q)\rangle + 
    \sum_{l, m}a_la_m\langle\boldsymbol\phi_k, \boldsymbol f_2(\boldsymbol\phi_l, \boldsymbol\phi_m)\rangle + 
    \sum_{l, m, n}a_la_ma_n\langle\boldsymbol\phi_k, \boldsymbol f_3(\boldsymbol\phi_l, \boldsymbol\phi_m, \boldsymbol\phi_n)\rangle
$$

where the details:

$$
    \langle\boldsymbol\phi_k, \boldsymbol f_1(\boldsymbol q)\rangle = \langle \boldsymbol\phi_k, \boldsymbol f_1(\bar{\boldsymbol q}) \rangle + \sum_l a_l \langle \boldsymbol\phi_k, \boldsymbol f_1(\boldsymbol\phi_l) \rangle 
$$

$$
    \langle\boldsymbol\phi_k, \boldsymbol f_2(\boldsymbol q, \boldsymbol q)\rangle = \langle \boldsymbol\phi_k, \boldsymbol f_2(\bar{\boldsymbol q}, \bar {\boldsymbol q}) \rangle + \sum_l a_l \langle \boldsymbol\phi_k, \boldsymbol f_2(\bar{\boldsymbol q}, \boldsymbol\phi_l) + \boldsymbol f_2(\boldsymbol\phi_l, \bar{\boldsymbol q}) \rangle + \sum_{l,m} a_la_m \langle \boldsymbol\phi_k, \boldsymbol f_2(\boldsymbol\phi_l, \boldsymbol\phi_m) \rangle 
$$

$$
    \langle\boldsymbol\phi_k, \boldsymbol f_3(\boldsymbol q, \boldsymbol q, \boldsymbol q)\rangle = \langle \boldsymbol\phi_k, \boldsymbol f_3(\bar{\boldsymbol q}, \bar {\boldsymbol q}, \bar{\boldsymbol q}) \rangle + \sum_l a_l \langle \boldsymbol\phi_k, \boldsymbol f_3(\bar{\boldsymbol q}, \bar{\boldsymbol q}, \boldsymbol\phi_l) + \boldsymbol f_3(\bar{\boldsymbol q}, \boldsymbol\phi_l, \bar{\boldsymbol q}) + \boldsymbol f_3(\boldsymbol\phi_l, \bar{\boldsymbol q}, \bar{\boldsymbol q}) \rangle + \sum_{l,m} a_la_m \langle \boldsymbol\phi_k, \boldsymbol f_3(\bar{\boldsymbol q}, \boldsymbol\phi_l, \boldsymbol\phi_m) + \boldsymbol f_3(\boldsymbol\phi_l, \bar{\boldsymbol q}, \boldsymbol\phi_m) + \boldsymbol f_3(\boldsymbol\phi_l, \boldsymbol\phi_m, \bar{\boldsymbol q}) \rangle + \sum_{l,m,n} a_la_ma_n \langle \boldsymbol\phi_k, \boldsymbol f_3(\boldsymbol\phi_l, \boldsymbol\phi_m, \boldsymbol\phi_n) \rangle
$$

The resulting Galerkin equations are given by

$$
    \dot{a}_k = b_k^1+ b_k^2 + b_k^3 + \sum_l a_l (L_{kl}^1 + L_{kl}^2 + L_{kl}^3) + \sum_{l,m} (a_la_m Q_{klm}^2 + a_la_m Q_{klm}^3) + \sum_{l,m,n} a_la_ma_n C_{klmn}^3 
$$

$$
    b_k^1 = \langle \boldsymbol\phi_k, \boldsymbol f_1(\bar{\boldsymbol q}) \rangle 
$$

$$
    b_k^2 = \langle \boldsymbol\phi_k, \boldsymbol f_2(\bar{\boldsymbol q}, \bar{\boldsymbol q}) \rangle 
$$

$$
    b_k^3 = \langle \boldsymbol\phi_k, \boldsymbol f_3(\bar{\boldsymbol q}, \bar{\boldsymbol q}, \bar{\boldsymbol q}) \rangle
$$

$$
    L_{kl}^1 = \langle \boldsymbol\phi_k, \boldsymbol f_1(\boldsymbol\phi_l) \rangle
$$

$$
    L_{kl}^2 = \langle \boldsymbol\phi_k, \boldsymbol f_2(\bar{\boldsymbol q}, \boldsymbol\phi_l) + \boldsymbol f_2(\boldsymbol\phi_l, \bar{\boldsymbol q}) \rangle
$$

$$
    L_{kl}^3 = \langle \boldsymbol\phi_k, \boldsymbol f_3(\bar{\boldsymbol q}, \bar{\boldsymbol q}, \boldsymbol\phi_l) + \boldsymbol f_3(\bar{\boldsymbol q}, \boldsymbol\phi_l, \bar{\boldsymbol q}) + \boldsymbol f_3(\boldsymbol\phi_l, \bar{\boldsymbol q}, \bar{\boldsymbol q}) \rangle
$$

$$
    Q_{klm}^2 = \langle \boldsymbol\phi_k, \boldsymbol f_2(\boldsymbol\phi_l, \boldsymbol\phi_m) \rangle
$$

$$
    Q_{klm}^3 = \langle \boldsymbol\phi_k, \boldsymbol f_3(\bar{\boldsymbol q}, \boldsymbol\phi_l, \boldsymbol\phi_m) + \boldsymbol f_3(\boldsymbol\phi_l, \bar{\boldsymbol q}, \boldsymbol\phi_m) + \boldsymbol f_3(\boldsymbol\phi_l, \boldsymbol\phi_m, \bar{\boldsymbol q}) \rangle
$$

$$
    C_{klmn}^3 = \langle \boldsymbol\phi_k, \boldsymbol f_3(\boldsymbol\phi_l, \boldsymbol\phi_m, \boldsymbol\phi_n) \rangle
$$
