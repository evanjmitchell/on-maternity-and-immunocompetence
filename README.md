This respository contains all of the Matlab code used in the paper "On maternity and the stronger immune system of women", by E. Mitchell, F. Ubeda, and G. Wild.

# System Requirements
All files require core Matlab, version R2019a or higher. No special-purpose toolboxes are required. System requirements are described in detail by the manufacturer at <https://www.mathworks.com/support/requirements/matlab-system-requirements.html>.

# Installation Guide
A Matlabe license can be purchased from the manufacturer at <https://www.mathworks.com/store/>, but academic institutions often have site licenses available to faculty and students free of charge.

Core Matlab can be installed following the manufacturer's instructions at <https://www.mathworks.com/help/install/index.html>. Installation time will vary depending on hardware and operating system.

The user should download the following files from this repository and save them to their working directory:

(1) resident_VT.m 
- function to define the system of differential equations governing the resident population

(2) findEE_VT.m 
- function to approximate the endemic equilibrium of the resident system

(3) fitnessAlpha_VT.m 
- function to compute the pathogen fitness

(4) fitnessGamma_VT.m 
- function to compute the host fitness

(5) sgradAlphaF_VT.m 
- function to define the selection gradient for alpha in female hosts

(6) sgradAlphaM_VT.m 
- function to define the selection gradient for alpha in male hosts

(7) sgradGammaF_VT.m 
- function to define the selection gradient for gamma in female hosts

(8) sgradGammaM_VT.m 
- function to define the selection gradient for gamma in male hosts

(9) checkESSAlpha_VT.m 
- function to check that the equilibrium values of alpha satisfy the evolutionary stability condition

(10) checkESSGamma_VT.m 
- function to check that the equilibrium values of gamma satisfy the evolutionary stability condition

(11) findCSS_VT.m 
- function to approximate the convergence stable evolutionary equilibrium

# Demo

Info on demo file for generating small data set.

# Instructions for Use

Info on generating figures from data.