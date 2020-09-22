# Order Parameter Distribution of 2D Ising Model

<br>This repository contains source code for Monte Carlo simulations of 2D Ising model.

## Introduction

The program calculates the statistical distribution of the order parameter of a [2D Ising model](https://en.wikipedia.org/wiki/Ising_model) at a given size L and a temperature (Ising coupling J) using Monte Carlo method. The difficulties will arise in large systems and especially at low temperature that the system gets stuck at local minima. To overcom, a Successive Umbrella Sampling is implemented that samples the entire range of order parameter in overlapping windows successively [Virnau and Müller, 2004](https://link.springer.com/chapter/10.1007/3-540-26565-1_18). To give a sence of what is the order parameter, one can think of magnetization in ferromagnetic matterials, the order parameter is +1 when all spins are up, -1 when all down, and 0 in case of 50-50 up/down spins.<br>

## Requirements

This program has no special requirements.<br>

## Compilation

Compile the program with **make**:

```bash
make clean
make
```

or run the command:

```bash
g++ -O2 -Wall main.cpp -o susIsing
```

To run the program for different sizes and temperature edit ``run.sh`` and run:

```bash
./run.sh
```
This creates a ``data`` directory and saves the output data therein.<br>

## Configuration

### Parameters:
The system size and the temperature (Ising coupling) are given by default and can be hardcoded as system setting in the **main.cpp**. They may optionally be set through input arguments by flags ``-L`` and ``-J``, e.g.:

```bash
./susIsing -L 100 -J 0.5
```

### Hyperparameters:
The number of steps of the Successive Umbrella Sampling in each window is given experimentaly based upon the desired accuracy and the available resources.<br>

## Output Data

* **prob\*.dat**: The probability distribution over the range of order paramter.

* **spin\*.dat**: A snapshot of the system at 50-50 composition.
