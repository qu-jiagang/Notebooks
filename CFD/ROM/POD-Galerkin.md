# POD-Galerkin Model

> update 2024-10-25

Consider a dynamical system which evolves in a Hilbert space H, In particular, for $\bm u(t)\in \bm H$, $\bm u(t)$ satisfies
$$
    \dot{\bm u} = \bm X(\bm u)
$$

$\bm u(t)$ can be write with Proper Orthgonal Decomposition (POD),
$$
    \bm u(t) = \sum_k a_k(t) \bm\phi_k
$$

Inserting this into the dynamical system, we get
$$
    \dot{\bm u} = \sum_k \dot{a}_k(t) \bm\phi_k = \bm X (\bm u (t))
$$

In order to get the reduced order dynamical system of $a_k(t)$, by projecting the equation into the basis $\bm \phi_k$.
$$
    \langle\sum_k \dot{a}_k(t) \bm\phi_k, \bm\phi_l\rangle = \langle\bm X (\bm u (t)), \bm\phi_l\rangle
$$
Assuming that the basis $\bm \phi_k$ are orthgonal, 
$$
    \dot{a}_l(t) = \langle\bm X (\bm u (t)), \bm\phi_l\rangle
$$

## Governing Equations
3-dimensional Navier-Stokes equations is given by
$$
    \frac{\partial \rho}{\partial t} + (\bm u \cdot \nabla)\rho + \nabla \cdot \bm{u} = 0 \\
    \rho \left(\frac{\partial \bm u}{\partial t} + (\bm u \cdot \nabla) \bm u \right ) = -\nabla p + \nu \nabla^2\bm u \\
    \rho\left( \frac{\partial T}{\partial t} + (\gamma - 1)T\nabla\cdot \bm u \right) = \nabla\cdot\bm (\rho \bm u T) + \frac{\gamma\nu}{Pr}\nabla^2T
$$

where $\bm q = (\rho, u_1, u_2, u_3, T)$, $\bm R, \bm U_i, \bm \Theta$ are the nonlinear differential operators given by
$$
    \bm R(\bm q) = -(\bm u \cdot \nabla)\rho - \nabla \cdot \bm{u} \\
    \bm U_i(\bm q) = -\rho (\bm u \cdot \nabla) \bm u -\nabla p + \nu \nabla^2\bm u \\
    \bm \Theta(\bm q) = -(\gamma - 1)T\nabla\cdot u + \nabla\cdot\bm (\rho \bm u T) + \frac{\gamma\nu}{Pr}\nabla^2T
$$

written more concisely as 
$$
    \bm A (\bm q)\dot{\bm q} = \bm f(\bm q)
$$

where $\bm A = diag(1, \rho, \rho, \rho, \rho)$, $\bm f (\bm q) = (\bm R(\bm q), \bm U_1(\bm q), \bm U_2(\bm q), \bm U_3(\bm q), \bm \Theta(\bm q))$, 

$\bm A$ is affine in $\bm q$, and may be written 
$$
    \bm A (\bm q) = \bm B + \bm L (\bm q) = diag(1, 0, 0, 0, 0) + diag(0, \rho, \rho, \rho, \rho)
$$

Furthermore, $\bm f$ is cubic and may be written
$$
    \bm f(\bm q) = \bm f_1(\bm q) + \bm f_2(\bm q, \bm q) + \bm f_3(\bm q, \bm q, \bm q) 
$$
where $\bm f_{1,2,3}$ are multilinear functions (linear in each argument).

## Galerkin projection

$$
    \bm f_1(\bm q) = 
    \begin{bmatrix}
        u_x + v_y + w_z \\
        - p_x + \nu(u_{xx} + u_{yy} + u_{zz}) \\
        - p_y + \nu(v_{xx} + v_{yy} + v_{zz}) \\
        - p_z + \nu(w_{xx} + w_{yy} + w_{zz}) \\
        \frac{\gamma\nu}{Pr}(T_{xx} + T_{yy} + T_{zz})
    \end{bmatrix} \\ 
    \bm f_2(\bm q^1, \bm q^2) = 
    \begin{bmatrix}
        u^1\rho^2_x + v^1\rho^2_y + w^1\rho^2_z \\
        0 \\
        0 \\
        0 \\
        -(\gamma - 1)T^1(u^2_x + v^2_y + w^2_z)
    \end{bmatrix} \\
    \bm f_3(\bm q^1, \bm q^2, \bm q^3) = 
    \begin{bmatrix}
        0 \\
        -\rho^1 (u^2 u^3_x + v^2 u^3_y + w^2 u^3_z) \\
        -\rho^1 (u^2 v^3_x + v^2 v^3_y + w^2 v^3_z) \\
        -\rho^1 (u^2 w^3_x + v^2 w^3_y + w^2 w^3_z) \\
        \rho^1_x u^2T^3 + \rho^1 u^2_x T^3 + \rho^1 u^2 T^3_x + \rho^1_y v^2T^3 + \rho^1 v^2_y T^3 + \rho^1 v^2 T^3_y + \rho^1_z w^2T^3 + \rho^1 w^2_z T^3 + \rho^1 w^2 T^3_z
    \end{bmatrix} \\
$$

let $\bm \phi_k$ be a basis for the function space containing $\bm q$, 
$$
    \bm q(x, t) = \bar{\bm q}(x) + \sum_k a_k(t)\bm\phi_k
$$
Inserting this expression into the governing equation,
$$
    \left[\bm B + \bm L\left( \sum_l a_l\bm\phi_l\right)\right] \sum_k\dot{a}_k\bm\phi_k = \bm f (\bm q)
$$

Taking an inner product with $\bm\phi_j$ then gives
$$
    \sum_k \dot{a}_k\left( \langle\bm\phi_j, \bm B\phi_k\rangle + \sum_l a_l \langle\bm\phi_j, \bm L(\bm\phi_l)\phi_k\rangle \right) = \langle\bm\phi_j, \bm f(\bm q)\rangle
$$

The equation above can be written in matrix form as
$$
    \bm M (\bm a) \dot{\bm a} = \bm F (\bm a)
$$
where $\bm a = (a_1, a_2, ..., a_n)$, and
$$
    M_{jk} = \langle \bm\phi_j, \bm B\phi_k \rangle + \sum_l a_l \langle\bm\phi_j, \bm L(\bm\phi_l)\phi_k\rangle \\
    F_j = \langle\bm\phi_j, \bm f(\bm q)\rangle = \sum_l a_l \langle\bm\phi_j, \bm f_1(\bm q)\rangle + \sum_{l, m}a_la_m\langle\bm\phi_j, \bm f_2(\bm\phi_l, \bm\phi_m)\rangle + \sum_{l, m, n}a_la_ma_n\langle\bm\phi_j, \bm f_3(\bm\phi_l, \bm\phi_m, \bm\phi_n)\rangle
$$

where the details:
$$
    \langle\bm\phi_j, \bm f_1(\bm q)\rangle = \langle \bm\phi_j, \bm f_1(\bar{\bm q}) \rangle + \sum_l a_l \langle \bm\phi_j, \bm f_1(\bm\phi_l) \rangle \\
    \langle\bm\phi_j, \bm f_2(\bm q, \bm q)\rangle = \langle \bm\phi_j, \bm f_2(\bar{\bm q}, \bar {\bm q}) \rangle + \sum_l a_l \langle \bm\phi_j, \bm f_2(\bar{\bm q}, \bm\phi_l) + \bm f_2(\bm\phi_l, \bar{\bm q}) \rangle + \sum_{l,m} a_la_m \langle \bm\phi_j, \bm f_2(\bm\phi_l, \bm\phi_m) \rangle \\
    \langle\bm\phi_j, \bm f_3(\bm q, \bm q, \bm q)\rangle = \langle \bm\phi_j, \bm f_3(\bar{\bm q}, \bar {\bm q}, \bar{\bm q}) \rangle + \sum_l a_l \langle \bm\phi_j, \bm f_3(\bar{\bm q}, \bar{\bm q}, \bm\phi_l) + \bm f_3(\bar{\bm q}, \bm\phi_l, \bar{\bm q}) + \bm f_3(\bm\phi_l, \bar{\bm q}, \bar{\bm q}) \rangle + \sum_{l,m} a_la_m \langle \bm\phi_j, \bm f_3(\bar{\bm q}, \bm\phi_l, \bm\phi_m) + \bm f_3(\bm\phi_l, \bar{\bm q}, \bm\phi_m) + \bm f_3(\bm\phi_l, \bm\phi_m, \bar{\bm q}) \rangle + \sum_{l,m,n} a_la_ma_n \langle \bm\phi_j, \bm f_3(\bm\phi_l, \bm\phi_m, \bm\phi_n) \rangle
$$

The resulting Galerkin equations are given by
$$
    \dot{\bm a} = b_k^1+ b_k^2 + b_k^3 + \sum_l L_{kl}^1 a_l + \sum_{l,m} L_{kl}^2 a_la_m + \sum_{l,m,n} L_{kl}^3 a_la_ma_n \\
$$

$$
    b_k^1 = \langle \bm\phi_k, \bm f_1(\bar{\bm q}) \rangle \\
    b_k^2 = \langle \bm\phi_k, \bm f_2(\bar{\bm q}, \bar{\bm q}) \rangle \\
    b_k^3 = \langle \bm\phi_k, \bm f_3(\bar{\bm q}, \bar{\bm q}, \bar{\bm q}) \rangle \\
$$

$$
    L_{ik}^1 = \langle \bm\phi_i, \bm f_1(\bm\phi_k) \rangle \\
    L_{ik}^2 = \langle \bm\phi_i, \bm f_2(\bar{\bm q}, \bm\phi_k) + \bm f_2(\bm\phi_k, \bar{\bm q}) \rangle \\
    L_{ik}^3 = \langle \bm\phi_i, \bm f_3(\bar{\bm q}, \bar{\bm q}, \bm\phi_k) + \bm f_3(\bar{\bm q}, \bm\phi_k, \bar{\bm q}) + \bm f_3(\bm\phi_k, \bar{\bm q}, \bar{\bm q}) \rangle \\
$$

$$
    Q_{ijk}^2 = \langle \bm\phi_i, \bm f_2(\bm\phi_j, \bm\phi_k) \rangle \\
    Q_{ijk}^3 = \langle \bm\phi_i, \bm f_3(\bar{\bm q}, \bm\phi_j, \bm\phi_k) + \bm f_3(\bm\phi_j, \bar{\bm q}, \bm\phi_k) + \bm f_3(\bm\phi_j, \bm\phi_k, \bar{\bm q}) \rangle \\
$$

$$
    C_{ijkl}^3 = \langle \bm\phi_i, \bm f_3(\bm\phi_j, \bm\phi_k, \bm\phi_l) \rangle
$$
