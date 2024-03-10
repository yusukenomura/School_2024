# Copyright (c) 2024 Yusuke Nomura
# Usage of RBM solver

# Programs 

 ./src/exact_sampling/SD/  => RBM solver using using steepest descent (SD), exact sampling
 ./src/exact_sampling/SR/  => RBM solver using using stochastic reconfiguration (SR) method, exact sampling   
 ./src/MC_sampling/SR/     => RBM solver using using stochastic reconfiguration (SR) method, sampling using Monte Carlo method 


# Compilation 

 1. Go to ./src/exact_sampling/SD/ or ./src/exact_sampling/SR or ./src/MC_sampling/SR  

     cd ./src/exact_sampling/SD/
       or     
     cd ./src/exact_sampling/SR/    
       or     
     cd ./src/MC_sampling/SR/    

 2. Edit Makefile 

 3. Type "make"
   
     If the executable file RBM_solver_SD(R).x exists, the compilation is successful     


# Examples
 
 ./examples/1D_AF_chain/exact_sampling/08site  => 1D antiferromagnetic Heisenberg model (8site)  
 ./examples/1D_AF_chain/exact_sampling/16site  => 1D antiferromagnetic Heisenberg model (16site)
 ./examples/1D_AF_chain/MC_sampling/16site     => 1D antiferromagnetic Heisenberg model (16site)


# Execution

 1. copy executable files to working directory 

 2. ./RBM_solver_SD(R).x > log

    For exact sampling, RBM_solver_SD(R).x requires RBM.input and spin_configurations.txt as inputs
    For Monte Carlo sampling, RBM_solver_SR.x requires only RBM.input as input
    log is standard output    


