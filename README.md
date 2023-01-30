## License

Copyright (c) 2022, Roland Griesmaier, Marvin Knöller, Rainer Mandel

This software is released under GNU General Public License, Version 3.
The full license text is provided in the file LICENSE included in this repository 
or can be obtained from http://www.gnu.org/licenses/

## Content of this project
This is a guide to generate the figures that have been used in the work

**„Inverse medium scattering for a nonlinear Helmholtz equation“** 

_By Roland Griesmaier, Marvin Knöller and Rainer Mandel._

You find all needed Matlab files to generate the figures. 

## An overview
- [ ] _evaluategfun_z_ generates the Herglotz density g, dependent on a possible shift z\in R^2
- [ ] _Finalplots_ plots the figures at the end of the computation
- [ ] _funhandle_zAbs_ evaluates the function handle corresponding to the factorization method
- [ ] _funhanlde_zReal_ evaluates the function handle corresponding to the monotonicity method
- [ ] _getc_ and _ToepPhi_ are used to evaluate the Toeplitz matrix in order to evaluate the 2d convolution from the nonlinear Lippmann Schwinger equation. Convolution is performed by using the 2d Fourier transform.
- [ ] _getUi_z_ generates incoming Herglotz fields, dependent on a possible shift z\in R^2
- [ ] _mycon_ is the constraint used in the optimization
- [ ] _NLHH evaluates_ the far field given an incoming field. This function uses a fixed point iteration arising from the nonlinear Lippmann Schwinger equation.
- [ ] _nonlinear_qh2_scaled_ gives the (scaled) function handle corresponding to a kite made of fused silica. 

The scripts _Numerical_Example_Fac.m_ and _Numerical_Example_Mon.m_ start the reconstruction of the kite using the factorization and the monotonicity method, respectively.

## Requirements and additional information
The computations have been carried out on a Cluster using 32 Cores.
Generating an example from scratch takes approximately 4 days.
Computations have been carried out using the Matlab 2018a version.

The code uses parallelization from the Matlab Parallelization Toolbox.
The code uses optimization from the Matlab Optimization Toolbox.
