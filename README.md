# Lid-Driven Cavity Flow Simulation

This repository contains the source code for a 2D incompressible Navier-Stokes solver using the **SIMPLE algorithm** to simulate lid-driven cavity flow. The code implements the Finite Volume Method (FVM) on a staggered grid.

## Features

* Solves steady-state incompressible Navier-Stokes equations.
* Implements the SIMPLE algorithm for pressure-velocity coupling.
* Uses a staggered grid for velocity and pressure variables.
* Demonstrates grid independence by simulating with three different grid sizes.
* Visualizes results :
    * Contour plots of x-velocity, y-velocity, and pressure.
    * Line plots of velocity and pressure along centerlines.
    * Streamline plots for different grid sizes.
