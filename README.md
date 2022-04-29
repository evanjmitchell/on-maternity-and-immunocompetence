This respository contains all of the Matlab code used in the paper "On maternity and the stronger immune system of women", by E. Mitchell, A. Graham, F. Ubeda, and G. Wild.

# System Requirements

All files require core Matlab, version R2019a or higher. This code makes use of the Parallel Computing Toolbox. System requirements are described in detail by the manufacturer at <https://www.mathworks.com/support/requirements/matlab-system-requirements.html>.

# Installation Guide - Section 3 of the Supplement

A Matlab license can be purchased from the manufacturer at <https://www.mathworks.com/store/>, but academic institutions often have site licenses available to faculty and students free of charge.

Core Matlab can be installed following the manufacturer's instructions at <https://www.mathworks.com/help/install/index.html>. Installation time will vary depending on hardware and operating system.

The user should download the following files from this repository and save them to their working directory:

(1) ResEquil.m
- Function that finds an asymptotically stable equilibrium solution to system (2) of the Supplement iteratively. Takes model parameters as inputs and returns an equilibrium.

(2) Wp.m
- Function that calculates invasion fitness of a mutant pathogen as the spectral radius of the matrix in equation (7) of the Supplement. Takes model parameters and an estimate of the equilibrium state of the resident population as inputs and returns invasion fitness of pathogen.

(3) Wh.m
- Function that calculates invasion fitness of a mutant pathogen as the spectral radius of the matrix in equation (16) of the Supplement. Takes model parameters and an estimate of the equilibrium state of the resident population as inputs and returns invasion fitness of host.

# Installation Guide - Section 4 of the Supplement

A Matlab license can be purchased from the manufacturer at <https://www.mathworks.com/store/>, but academic institutions often have site licenses available to faculty and students free of charge.

Core Matlab can be installed following the manufacturer's instructions at <https://www.mathworks.com/help/install/index.html>. Installation time will vary depending on hardware and operating system.

The user should download the following files from this repository and save them to their working directory:

(1) resident_VT.m 
- Function to define the system of differential equations governing the resident population.

(2) findEE_VT.m 
- Function to approximate the endemic equilibrium of the resident system.

(3) fitnessAlpha_VT.m 
- Function to compute the pathogen fitness.

(4) fitnessGamma_VT.m 
- Function to compute the host fitness.

(5) sgradAlphaF_VT.m 
- Function to define the selection gradient for alpha in female hosts.

(6) sgradAlphaM_VT.m 
- Function to define the selection gradient for alpha in male hosts.

(7) sgradGammaF_VT.m 
- Function to define the selection gradient for gamma in female hosts.

(8) sgradGammaM_VT.m 
- Function to define the selection gradient for gamma in male hosts.

(9) checkESSAlpha_VT.m 
- Function to check that the equilibrium values of alpha satisfy the evolutionary stability condition.

(10) checkESSGamma_VT.m 
- Function to check that the equilibrium values of gamma satisfy the evolutionary stability condition.

(11) findCSS_VT.m 
- Function to approximate the convergence stable evolutionary equilibrium. Serves as a wrapper function for scripts (1)-(10).

# Installation Guide - Section 5 of the Supplement

A Matlab license can be purchased from the manufacturer at <https://www.mathworks.com/store/>, but academic institutions often have site licenses available to faculty and students free of charge.

Core Matlab can be installed following the manufacturer's instructions at <https://www.mathworks.com/help/install/index.html>. Installation time will vary depending on hardware and operating system.

The user should download the following files from this repository and save them to their working directory:

(1) ResEquilAct.m
- Function that finds an asymptotically stable equilibrium solution to system (17) of the Supplement iteratively. Takes model parameters as inputs and returns an equilibrium.

(2) WpAct.m
- Function that calculates invasion fitness of a mutant pathogen as the spectral radius of the matrix in equation (19) of the Supplement. Takes model parameters and an estimate of the equilibrium state of the resident population as inputs and returns invasion fitness of pathogen.

(3) WhAct.m
- Function that calculates invasion fitness of a mutant pathogen as the spectral radius of the matrix in equation (20) of the Supplement. Takes model parameters and an estimate of the equilibrium state of the resident population as inputs and returns invasion fitness of host.

# Demo - Section 3 of the Supplement

The file maindemo.m is a Matlab script that demos how we estimate the joint ESS quadruple alpha_f*, alpha_m*, gamma_f*, and gamma_m* following Algorithm 1 of the Supplement. The script also assesses evolutionary stability by checking the second-order conditions described in the Supplement. The demo relies on functions ResEquil.m, Wp.m, and Wh.m and so these files should appear in the same working directory as maindemo.m. The output of the script is written to a comma-delimited file maindemo.csv and can be generated in about 0.5 sec on a standard laptop.

# Demo - Section 4 of the Supplement

The file demo.m is a demo of how to use the core Matlab files to generate data sets. Make sure that demo.m is saved in the same directory as the above Matlab files. To run the demo, open demo.m and click the green "Run" arrow in the top menu bar. This file generates a small test data set containing 20 observations and should match with the data_demo.xlsx file provided in the repository. It takes approximately 764.37 sec (12.74 min) to run on a laptop equipped with an Intel Core i5-11400H CPU and 8GB of RAM.

# Demo - Section 5 of the Supplement

The file maindemoAct.m is a Matlab script that demos how we estimate the joint ESS quadruple alpha_f*, alpha_m*, gamma_f*, and gamma_m* following Algorithm 1 of the Supplement. The script also assesses evolutionary stability by checking the second-order conditions described in the Supplement. The demo relies on functions ResEquilAct.m, WpAct.m, and WhAct.m and so these files should appear in the same working directory as maindemoAct.m. The output of the script is written to a comma-delimited file maindemoAct.csv and can be generated in about 0.5 sec on a standard laptop.

# Instructions for Use

Any of the demo files can be modified to generate larger data sets by adding more vertical transmission and cost difference values. Other parameters can be modified where they are introduced into the scripts.
