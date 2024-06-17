# RogueWaveInfiniteNLS.jl
Julia package for computing the rogue waves of infinite order for the focusing nonlinear Schrodinger equation.
This is a special solution of the focusing nonlinear equation denoted by $\Psi(X,T;\mathbf{G},\beta)$, which depends parametrically on a $2\times 2$ matrix and a scalar $\beta>0$. The matrix $\mathbf{G}$ is given by $$\mathbf{G}=\mathbf{G}(a,b)=\frac{1}{\sqrt{|a|^2+|b|^2}}\begin{bmatrix} a & b^* \\ -b & a^* \end{bmatrix},$$ where $c^*$ denotes the complex conjugate of $c\in\mathbb{C}$.

The main routine is `psi`.
