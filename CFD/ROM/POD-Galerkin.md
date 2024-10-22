# Compressible 3D NS equations

## POD

## Governing Equations
3-dimensional Navier-Stokes equations is given by
$$
    \frac{\partial \rho}{\partial t} + (\bm u \cdot \nabla)\rho + \nabla \cdot \bm{u} = 0 \\
    \rho \left(\frac{\partial \bm u}{\partial t} + (\bm u \cdot \nabla) \bm u \right ) = -\nabla p + \nu \nabla^2\bm u \\
    \rho\left( \frac{\partial T}{\partial t} + (\gamma - 1)T\nabla\cdot u \right) = \gamma\nu\nabla^2\bm u + \frac{\gamma\nu}{Pr}\nabla^2T
$$

where $\bm q = (\rho, u_1, u_2, u_3, T)$, $\bm R, \bm U_i, \bm \Theta$ are the nonlinear differential operators given by
$$
    \bm R(\bm q) = -(\bm u \cdot \nabla)\rho - \nabla \cdot \bm{u} \\
    \bm U_i(\bm q) = -\rho (\bm u \cdot \nabla) \bm u -\nabla p + \nu \nabla^2\bm u \\
    \bm \Theta = -(\gamma - 1)T\nabla\cdot u + \gamma\nu\nabla^2\bm u + \frac{\gamma\nu}{Pr}\nabla^2T
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
$$
    \bm f_1 = 
    \begin{bmatrix}
        u_x + v_y + w_z \\
        - p_x + \nu(u_{xx} + u_{yy} + u_{zz}) \\
        - p_y + \nu(v_{xx} + v_{yy} + v_{zz}) \\
        - p_z + \nu(w_{xx} + w_{yy} + w_{zz}) \\
        \gamma\nu(T_{xx} + T_{yy} + T_{zz})
    \end{bmatrix}
$$


## Galerkin projection
let $\bm \phi_k$ be a basis for the function space containing $\bm q$, 
$$
    \bm q(t) = \sum_k a_k(t)\bm\phi_k
$$
Inserting this expression into the governing equation,
$$
    \left[\bm B + \bm L\left( \sum_l a_l\bm\phi_l\right)\right]\sum_k\dot{a}_k\bm\phi_k = \bm f (\bm q)
$$

