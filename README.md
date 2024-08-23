# Zipfian Paper Simulations

This repository contains the code to run and then visualize the simulations.

The simulations are labeled sim1, sim2, and sim3, corresponding to the left, middle, and right panels of figure 2 in the paper.

sim1,sim2 and sim3 are implemented in simulationFunctions.py. simulationFunctions.py depends only on the numpy package.

To run each simulation, in the directory of simulationFunctions.py, run 

python3 simX.py 

in the terminal. This will create a csv file of the data output from the simulation and store it in the data/simX.py file

analysis.R creates the visualizations for the 3 simulations. After creating the 3 simulation data files, the simulation visuals are created by opening R (or RStudio) and running analysis.R

The simulations for the paper were created on an Apple M1 Macbook Pro with OS Ventura using Python 3.9.6 and R 4.2.2
