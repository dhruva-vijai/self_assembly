# Simulating Shape-Assisted Self-Assembly 

## Problem Statement

## Background Information

## Specifics

## Necessary Libraries

## Basic Setup

## Problems Faced

## Improvements

## File Structure

## Images/Videos 

-This project is a Python based MD simulation that describes shape-assisted self assembly of car-pyridine derivatives into 1D nano-fibers in an aqueous environment. 
-The project is written in Python using OpenMM and OpenFF to characterise a carpyridine and understand shape-assisted self-assembly.

- Note : carpyridine.pdb is the basic 2H-Car-Ph model, and to simulate other car-pyridine derivatives, the corresponding pdb and smiles files need to be downloaded first.
- The files : assembly.dcd and topology.pdb must be loaded into VMD to visualise the stacking process.

- The specific degree of stacking and 1D fiber formation can be measured by using the Tk Console in VMD to znslyse center-to-center distance and interplanar angle between the aromatic rings

- The program is run from the Command-Line/Terminal in a virtual environment created using MiniConda and Python3 with the following libraries installed : openmm/openff.toolkit/openff.units/rdkit/numpy

- Current simulation was capped at 10ns for me due to time constraints and speed limitations and took around 2 hrs
- Speed Improvements : 1) Writing the program in Fortran/C should significantly improve speed.
                       2) Using the HIP/Metal platform in this code should lead to a 5-18x speedup however, I'm facing some                                    installation issues at the current time.
                       3) GPU Acceleration can significantly boost processing times.
                       4) More efficient code

-Note : To analyse different car-pyridine derivatives, first; it has to be built in Avogadro and then optimised or could be downloaded as a pdb from any research sources
