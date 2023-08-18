# FireCat

[![DOI](https://zenodo.org/badge/679916980.svg)](https://zenodo.org/badge/latestdoi/679916980)

FireCat performs multiscale simulations by coupling ab initio microkinetic simulations from CatMAP with two-dimensional continuum transport simulations from EchemFEM (Firedrake). This coupling is performed in an iterative fashion. The examples provided here are specific to electrochemical CO2 reduction to CO on Au electrodes.

## Getting started

Please install the open-source finite element library [Firedrake](https://www.firedrakeproject.org/download.html).

Within the firedrake venv, install [EchemFEM](https://github.com/LLNL/echemfem) and this [fork of CatMAP](https://github.com/sringe/catmap-1).

In order to run the 2D examples, [GMSH](https://gmsh.info) is needed to the create the <tt>`.msh`</tt> mesh files from the given <tt>`.geo`</tt> files.
