# Simulating Shape-Assisted Self-Assembly 

## Problem Statement

Use Python to simulate the shape-assisted self-assembly of carpyridines into 1D fibers and 2D nanosheets in an aqueous environment. Then, use a machine-learning algorithm to analyse the relationship between the substituents in the carpyridine core and the quality and nature of stacking. 

## Background Information

Carpyridines form 1D nanofibers via pi-pi stacking and H-bonding interactions in an aqueous environment. These electronic and structural factors are the primary factor in the inducement of stacking/formation of supramolecular structures; however, the local topography and curvature of the molecule are shape-based features that assist in the formation of supramolecular structures by acting as saddle-shaped rotational locks and limiting freedom in the resultant structure. Hydrogen-bridging interactions also play a major role in this process; thus, the choice of solvent and ambient conditions are extremely important to the ooutcome of self-assembly.

## Specifics

__1) Collecting Structures__

The various car-pyridine structures were built by hand, except for 2H-CAR-Ph and 2H-CAR-C6 which were obtained from the CCSD Structural Database. The subsequent structures were then built by hand in Avogadro and saved as .pdb files.

__2) Simulation__

The code aims to simulate a solvated system of 50 self-assembling molecules with 20000 water molecules in a 10 nm box size. The topologies were created using OpenFF and the OpenFF-2.1.0 Force Field. HMR was also used to increase integrator timestep to 4 fs while maintaining stability and reducing NaN errors.

The simulation was conducted in three stages - a short NVT equilibriation of 0.2 ns, an NPT equilibriation of 0.4 ns and the final production run; which lasts between 5-100 ns depending on the final outcome expected.

If the data from stacking is to be used in ML models, the 5 ns time period is used to extract basic descriptors of stacking quality; while, for true visualisation of stacking, a long 100 ns run is required which i limited by my computational resources.

__3) Visualisation__

The final results are loaded into VMD to visualise the stacking process and ensure the programs works as expected



## Necessary Libraries

1) OpenMM

2) OpenFF

3) RDKit

4) NumPy

5) SciPy

6) MDTraj

## Basic Setup

cd projects

git clone 

__**-Note :Change carp.pdb to structures/carp.pdb**__

python modified.py

__Load assembly.dcd and topology.pdb into VMD__


## Problems Faced

1) The singular biggest problem I am currently facing is the computational resource bottleneck that makes it incredibly taxing to run a full 100 ns sim with 20000 waters; so the ML part of it has not been tested yet on a proper dataset as generating the data itself takes a huge time investment.

2) Some of the assembler structures are not verified via crystal structures, thus are not ideal for simulation.

## Improvements

1) Writing the program in Fortran/C should significantly improve speed.

2) Using the HIP/Metal platform in this code should lead to a 5-18x speedup however, I'm facing some installation issues at the current time.

3) GPU Acceleration can significantly boost processing times.

4) More efficient code

5) Better quality assembler structures/better forcefields

## File Structure

### Version 1

- simulation.py deals with creating the assembler+solvent topology and running the simulation. The visualisation of this is further done in VMD

### Version 2

-Involves running a 5 ns sim for 9 diff assemblers and extracting trainable features from all of them; will be added eventually


## Images/Videos 

The following are more-so toy models with around 50 assemblers and 400 waters in a 10 nm box over "CHECK EXACT" nm

