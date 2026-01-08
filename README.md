# Dark Matter Mapping from Weak Gravitational Lensing Data (KiDS)

## Overview  
This project reconstructs the **projected dark matter (convergence) field** from **weak gravitational lensing shear data** using observations from the **Kilo-Degree Survey (KiDS)**. By statistically inverting distortions in galaxy shapes, we infer the underlying mass distribution using the **Kaiser–Squires method**.  

The analysis uses real survey data that includes shear measurements for over **200,000 galaxies**. This demonstrates a complete weak-lensing pipeline from raw catalogs to convergence maps. 

---
## Results

![Weak-lensing convergence map](kappa_map.png)

**Figure 1:** Reconstructed convergence kappa map from KiDS weak-lensing data, tracing the projected dark matter distribution in the selected sky patch. The map is smoothed for visualization.

---
## Scientific Background  
Weak gravitational lensing comes from the bending of light by large groups of matter. This results in small but clear distortions in the shapes of distant galaxies. We can measure these distortions using the shear components **gamma_1** and **gamma_2**. By inverting these measurements, we can find the **convergence field** kappa, which shows the projected dark matter density.

The **Kaiser–Squires (KS) inversion** offers a method in Fourier space to rebuild **kappa** from shear measurements, assuming a flat sky. This makes it a common tool in observational cosmology.

---
## Data
- Survey: **Kilo-Degree Survey (KiDS)**
- Input: Galaxy shear catalog
- Fields used:
  - Right ascension (RA)
  - Declination (Dec)
  - Shear components **e_1, e_2**
  - Per-galaxy weights
- Dataset size: **2 * 10^5** galaxies (after quality cuts)

A rectangular sky patch is selected to ensure the validity of the flat-sky approximation.

---
## Methodology
The weak-lensing reconstruction pipeline consists of the following steps:

1. **Catalog filtering**
   - Selection of a contiguous sky patch
   - Removal of galaxies with invalid shear measurements or zero weight

2. **Projection to a Cartesian grid**
   - Tangent-plane approximation around the patch center
   - Binning galaxy shear measurements onto a regular grid using weighted averaging

3. **Shear preprocessing**
   - Mean subtraction to suppress DC modes
   - Handling of empty pixels and masking

4. **Kaiser–Squires inversion**
   - Fourier-space inversion of shear to convergence
   - Reconstruction of the projected dark matter field

5. **Apodization and smoothing**
   - Edge apodization to reduce FFT artifacts
   - Gaussian smoothing for visualization and noise suppression

6. **E/B-mode decomposition**
   - Separation of physical E-modes from spurious B-modes
   - Use of B-modes as a systematic diagnostic
---
