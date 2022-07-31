This respository contains all of the code used in the paper "On maternity and the co-evolution of a stronger immune system in women", by E. Mitchell, A. Graham, F. Ubeda, and G. Wild.

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
- Function that calculates invasion fitness of a mutant host as the spectral radius of the matrix in equation (14) of the Supplement. Takes model parameters and an estimate of the equilibrium state of the resident population as inputs and returns invasion fitness of host.

# Installation Guide - Section 4 of the Supplement

A Matlab license can be purchased from the manufacturer at <https://www.mathworks.com/store/>, but academic institutions often have site licenses available to faculty and students free of charge.

Core Matlab can be installed following the manufacturer's instructions at <https://www.mathworks.com/help/install/index.html>. Installation time will vary depending on hardware and operating system.

The user should download the following files from this repository and save them to their working directory:

(1) ResEquilAct.m
- Function that finds an asymptotically stable equilibrium solution to system (15) of the Supplement iteratively. Takes model parameters as inputs and returns an equilibrium.

(2) WpAct.m
- Function that calculates invasion fitness of a mutant pathogen as the spectral radius of the matrix in equation (17) of the Supplement. Takes model parameters and an estimate of the equilibrium state of the resident population as inputs and returns invasion fitness of pathogen.

(3) WhAct.m
- Function that calculates invasion fitness of a mutant host as the spectral radius of the matrix in equation (18) of the Supplement. Takes model parameters and an estimate of the equilibrium state of the resident population as inputs and returns invasion fitness of host.

# Demo - Section 3 of the Supplement

The file maindemo.m is a Matlab script that demos how we estimate the joint ESS quadruple alpha_f*, alpha_m*, gamma_f*, and gamma_m* following Algorithm 1 of the Supplement. The script also assesses evolutionary stability by checking the second-order conditions described in the Supplement. The demo relies on functions ResEquil.m, Wp.m, and Wh.m and so these files should appear in the same working directory as maindemo.m. The output of the script is written to a comma-delimited file maindemo.csv and can be generated in about 0.5 sec on a standard laptop.

# Demo - Section 4 of the Supplement

The file maindemoAct.m is a Matlab script that demos how we estimate the joint ESS quadruple alpha_f*, alpha_m*, gamma_f*, and gamma_m* following Algorithm 1 of the Supplement. The script also assesses evolutionary stability by checking the second-order conditions described in the Supplement. The demo relies on functions ResEquilAct.m, WpAct.m, and WhAct.m and so these files should appear in the same working directory as maindemoAct.m. The output of the script is written to a comma-delimited file maindemoAct.csv and can be generated in about 0.5 sec on a standard laptop.

# Instructions for Use

Any of the demo files can be modified to generate larger data sets by adding more vertical transmission and cost difference values. Other parameters can be modified where they are introduced into the scripts.

# Figure Source Code

The file SourceCode.zip contains the source code needed to create all of the figures in both the main text and supplement of the paper. To view this code, the user should download the zip file to their local machine. The code associated with each figure is contained in a separate subdirectory within the zip file, accompanied by a README document which outlines each of the files. The files to generate the figures are written in the Python language.
