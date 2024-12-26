# 2024_master_slave
# Source code to reproduce the results in the paper:
Mí­guez, J., Molina-Bulla, H., & Mariño, I.P. (2024). Master-slave coupling scheme for synchronization and parameter estimation in the generalized Kuramoto-Sivashinsky equation. Physical Review E 110(5), 054206


## List of matlab scripts
- *kuramoto.m*: main script
- *graphs.m*: plots of signals & synchronisation errors using the outputs of *kuramoto.m*

## List of matlab functions
- *KSeuler.m*: Euler scheme for the Kuramoto-Sivashinsky equation (master equation)
- *KScoupled_euler.m*: coupled slave equation with static parameters
-  *KScoupled.m*: coupled slave equation with time-varying parameters
- *KSparupd.m*: parameter & variable update step
- *KSdudt2.m*: time derivatives of the Fourier coefficients

## Initialisation data

- *init_a0_N32.mat*: initial Fourier coefficient values with $N=32$.
- *init_a0_N64.mat*: initial Fourier coefficient values with $N=64$.
