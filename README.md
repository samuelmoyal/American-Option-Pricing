# American-Option-Pricing


For a full description of this project, please refer to the detailed [project documentation](report.pdf).

# Summary


This repository contains a comprehensive analysis and implementation of various methods for pricing American options. The project is designed to provide a clear understanding of American option pricing mechanisms and to implement these using mathematical and computational techniques.

## Table of Contents

- [Introduction](#introduction)
- [Methods](#methods)
- [Usage](#usage)
- [Results](#results)



## Introduction

American options are financial derivatives that give the holder the right to exercise the option at any time before or at the expiration date. This flexibility makes American options more complex to price than European options, as they require consideration of optimal exercise strategies.

This project explores different methods to price American options, comparing their efficiency, accuracy, and computational demands.

## Methods

The repository implements and analyzes the following pricing methods:

1. **Binomial Trees**: A recursive method where each step divides the possible asset prices into an up or down move.
2. **Discretization Method**: A custom discretization approach tailored for this option pricing problem.
3. **Monte Carlo Simulation**: A stochastic method that simulates multiple paths of asset prices to approximate the option's value.

Each method includes explanations, code implementations, and analysis of its advantages and disadvantages.

## Usage

After installation, you can run the different pricing models using the following commands:

- **Binomial Trees**: Run the notebook `Binomial_tree.ipynb`
- **Discretization Method**: Run the notebook `DiscretisationVF.ipynb`
- **Monte Carlo Simulation**: Run the notebook `Monte_carlo.ipynb`

Each notebook includes configurable parameters for option characteristics and simulation details. Adjust these parameters as needed.

## Results

The project compares the results obtained from each method, highlighting their effectiveness in terms of accuracy and computational efficiency. Detailed results and analysis can be found in the `Report.pdf`.


