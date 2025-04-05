# Plasma Surface Recombination

A Monte Carlo simulation for the recombination of
atoms in a plasma with a surface.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Extra Files](#usage)
- [GUI window](#gui-window)  
  - [Select Reactions](#select-reactions)  
  - [Simulation Parameters](#simulation-parameters)  
    - [General Parameters](#general-parameters)  
    - [Reaction Constants](#reaction-constants)  
    - [Concentrations](#concentrations) 
    - [Stop Time](#stop-time)
- [Output](#output)  
    - [Terminal](#terminal)  
    - [Plots](#plots)  
      - [Concentration Evolution](#concentration-evolution)
      - [Reaction Rate Evolution](#reaction-rate-evolution)
- [Authors](#authors)

## Installation

1. Clone the repository:
   ```sh
   git clone https://github.com/JoseMariano247/Plasma-Surface-Recombination
   cd Plasma-Surface-Recombination
   ```

2. For the code to run properly, the user needs to install
the following program from this website:

https://sourceforge.net/projects/vcxsrv/

3. Run the .exe file installed.

4. Run the XLaunch application.

5. On the application window, select the following:

multiple windows
display number: 0
----> click next
start no client
----> click next
disable access control
----> click next
save

6. Change the IP associated:

   ```sh
   export DISPLAY=localhost:0.0
   ```

7. To check if the installation was successful, run:

   ```sh
   xclock
   ```


## Usage

To run the GUI, enter in the terminal:

   ```sh
   python3 GUI.py
   ```

## Extra Files

Apart from the GUI, we have added files to verify results with a RK4 solver. Furthermore, you can compile these with make, but you will need to create a file inside the folder ´build´ with ´mkdir obj´.

## GUI Window

If the installation worked as expected, there should be a pop up

### Select Reactions

`Basic` is the example reaction provided. It is simply the
transformation of some atom A to B and vice-versa. Since this
is not a standard reaction, it cannot be used in conjunction
with other reactions.
`Physisorption` is the action of the Van der Waals to recombinate
`Chemisorption` is the formation of chemical bonds between the
surface and the atoms.
`Surface Diffusion` is the Eley-Rideal (E-R) mechanism
`Langmuir-Hinshelwood recombination` regenerates chemisorption and
physisorption sites, creating molecules in the process.

### Simulation Parameters

#### General Parameters

`Tw` is the temperature of the surface.
`Tg` is the temperature of the gas.
`M` is the molar mass of the atom.

#### Reaction Constants

`Sticking Probability for Physisorption` is a parameter for the probability of, if physisorption is to happen, it actually happens.
`Desorption Frequency` is a parameter to measure how often desorption happens.
`Desorption Barrier` is a parameter to measure how hard it is
for desorption to happen.
`Sticking Probability for Chemisorption` is a parameter for the probability of, if chemisorption is to happen, it actually happens.
`Pre-Exponential Factor for Recombination` is a factor for the
recombination to happen.
`Recombination Barrier` is a parameter that measures how hard it is for recombination to happen.
`Diffusion Frequency` is a parameter that measures how often diffusion happens.
`Diffusion Barrier` is a parameter that measures how hard it is for diffusion to happen.
`Recombination Barrier for Two Physisorbed Atoms` is a parameter to measure how hard it is for two physisorbed atoms to recombine.


#### Concentrations

`A` is the concentration of atoms in the gas.
`Af` is the number of physisorption sites already filled.
`As` is the number of chemisorption sites already filled.
`Fv` is the number of physisorption sites empty.
`Sv` is the number of chemisorption sites empty.
`A2` is the concentration of molecules in the gas.

#### Stop Time

`Stop Time` is a parameter that says when the simulation will end.

## Output

If all the selected parameters and given constants are valid, the simulation will be run.

### Terminal

On your terminal, there should be a progress bar telling how much the code has run. In the end, the data will be saved to an output.txt file and a window will pop up telling the user the simulation has been run successfully

### Plots

After closing the previously mentioned window, a few plots will appear in the order shown in this README.

#### Concentration Evolution

The first plot that will appear will be the one regarding the evolution of the concentrations of each species. This will be normalized to the maximum value of each concentration.

#### Reaction Rate Evolution

The second plot that will appear will be the one regarding the evolution of the reaction rates of each reaction. This will be normalized to the maximum value of each reaction rate.

## Authors

This projected was made by:

Daniela Estaço

José Mariano


Under the guidance of:

Professor Vasco Guerra