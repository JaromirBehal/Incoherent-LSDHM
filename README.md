# Incoherent-LSDHM - Data and Function Structure
Reproduction of main results present in paper: J. Běhal and M. Ježek, "Incoherent lateral shearing digital holographic microscopy," Optica 12, 620-625 (2025)
https://doi.org/10.1364/OPTICA.555786

The repository is organized as follows:

## Main MATLAB files
These are the main MATLAB scripts used for loading and processing the data for the final presentation.

- **`Fig2.m`** — processes camera images acquired during the experiments with USAF 1951 resolution target and cheek cells.
- **`Fig3abcd.m`** — processes camera images acquired during the experiments with diatoms.
- **`Fig3efg.m`** — processes camera images acquired during the experiments with yeasts.
- **`Fig3hij.m`** — processes camera images acquired during the experiments with yeasts.

## `data/`

Contains the raw camera images acquired during the experiments with USAF 1951 resolution target, cheek cells, diatoms, and yeast cells.

## `functions/`

Contains the MATLAB functions used by the main scripts:

- **`apodiz_SG.m`** — higher-order Gaussian function used for apodization of the input data.
- **`apodization_for_propag.m`** — function for applying apodization prior to the Fourier transform of an image.
- **`off_axis_reconstruction_general.m`** — function for holographic reconstruction.
