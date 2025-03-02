# Gravity Experiment Simulation

This project simulates a series of experiments to analyze gravitational forces using a mathematical model. The simulation generates data points for each experiment, applies transformations, and calculates differences in forces.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)

## Overview

The project simulates gravitational experiments by generating data points, applying transformations, and calculating differences in forces. It uses mathematical models to simulate the effects of gravity at different latitudes and heights.

## Features

- Simulates multiple experiments with configurable parameters.
- Generates data points with added Gaussian noise.
- Calculates and compares average forces from different experiments.
- Provides matrix inversion functions for 2x2 and 3x3 matrices.

## Installation

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. Compile the code using a C compiler, for example:
   ```bash
   gcc -o gravity_simulation main.c -lm
   ```

## Usage

Run the compiled program:
```bash
./gravity_simulation
```

The program will generate data files for each experiment and print the results to the console.

## Project Structure

- `main.c`: The main source file containing the implementation of the simulation.
- `experiment_*.txt`: Generated data files for each experiment.

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request for any improvements or bug fixes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.