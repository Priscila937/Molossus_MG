# Geometric Morphometric (GM) Analysis of Bat Skulls in the Genus *Molossus*

This repository contains the script and associated information for the geometric morphometric (GM) analysis of bat skulls in the genus *Molossus*. The primary objective is to investigate the variation and evolution of skull shape and size, integrating genetic data for *Molossus* species.

---

## Methodology

### 1. Data
- **Images**: Two-dimensional images of skulls in ventral and dorsal views for each specimen and species.
- **Landmarks**: Type I landmarks were plotted: 17 on the dorsal view and 18 on the ventral view.
- **File Format**: Coordinates were extracted and saved in TPS files.

### 2. Procrustes Superimposition
- **Function**: `procSym` from the **Morpho v.2.12** package.
- **Purpose**: Removed size, position, and orientation effects for subsequent analyses.

### 3. Integration Analysis
- **Function**: `two.b.pls` from the **geomorph v.4.0.6** package.
- **Purpose**: Performed two-block partial least squares (PLS) analysis to assess integration between dorsal and ventral views.

### 4. Geometric Morphometric (GM) Analysis
- **Centroid Size (CS)**: Calculated as a proxy for skull size using superimposed coordinates.
- **Visualization**: Differences in CS across species visualized using violin box plots.
- **PCA**: Principal component analysis reduced data dimensionality.

### 5. Linear Regression and Allometry
- **Regression**: Tested the effect of CS on shape using linear regression.
- **Residuals**: Used to remove residual allometry in subsequent analyses.

### 6. Canonical Variate Analysis (CVA)
- **Function**: `CVA` from the **Morpho v.2.12** package.
- **Validation**: Cross-validated with 9,999 permutations.
- **Purpose**: Assessed cranial shape variation across species using raw and allometry-free data.

### 7. Phylomorphospace and Phylogenetic Signal
- **Functions**: `physignal` and `gm.prcomp` from the **geomorph** package.
- **Purpose**: Visualized morphospaces and tested for phylogenetic signal.
- **Tree Pruning**: Phylogenetic tree pruned to include only species with GM data.

### 8. Ancestral Character Reconstruction
- **Functions**: `contMap` from the **phytools** package and `fastAnc` from the **ape** package.
- **Purpose**: Visualized continuous evolution of cranial traits and reconstructed ancestral states using squared-change parsimony.

---

## Functions Used
- **procSym** (*Morpho*): Procrustes superimposition.
- **two.b.pls** (*geomorph*): Two-block PLS analysis for integration assessment.
- **CVA** (*Morpho*): Canonical variate analysis.
- **gm.prcomp** (*geomorph*): Principal component analysis.
- **physignal** (*geomorph*): Phylogenetic signal testing.
- **contMap** (*phytools*): Visualization of continuous traits on phylogenetic trees.
- **fastAnc** (*ape*): Ancestral character reconstruction.

---

## Expected Results
1. Identification of integration patterns between cranial views.
2. Visualization of shape and size differences among *Molossus* species.
3. Correlation between phylogeny and cranial shape in morphospace.
4. Reconstruction of cranial trait evolution in a phylogenetic context.

---

## References
- Adams, D. C., Collyer, M. L., & Kaliontzopoulou, A. (2022). *Geomorph: Software for Geometric Morphometric Analyses.* R Package version 4.0.4. [CRAN](https://cran.r-project.org/package=geomorph).
- Adams, D. C., Collyer, M. L. (2019). Comparing the strength of modular signal, and evaluating alternative modular hypotheses, using covariance ratio effect sizes with morphometric data. *Evolution, 73*(12), 2352–67. [DOI](https://doi.org/10.1111/evo.13867).
- Bookstein, F. L. (1991). *Morphometric Tools for Landmark Data.* Cambridge: Cambridge University Press.
- Blomberg, S. P., Garland, T., & Ives, A. R. (2003). Testing for phylogenetic signal in comparative data. *Evolution, 57*, 717–45. [DOI](https://doi.org/10.1111/j.0014-3820.2003.tb00285.x).
- Cardini, A., Jansson, A., & Elton, S. (2007). A geometric morphometric approach to study clinal variation in vervet monkeys. *Journal of Biogeography, 34*(10), 1663–78. [DOI](https://doi.org/10.1111/j.1365-2699.2007.01731.x).
- Collyer, M. L., Baken, E. K., & Adams, D. C. (2022). A standardized effect size for evaluating phylogenetic signal. *Methods in Ecology and Evolution, 13*(1), 367–82. [DOI](https://doi.org/10.1111/2041-210X.13749).
- Revell, L. J. (2012). Phytools: An R package for phylogenetic comparative biology. *Methods in Ecology and Evolution, 3*, 217–23. [DOI](https://doi.org/10.1111/j.2041-210X.2011.00169.x).
- Rohlf, F. J. (2015). The TPS series of software. *Hystrix, 26*(1), 1–4. [DOI](https://doi.org/10.4404/hystrix-26.1-11264).
- Schlager, S. (2017). Morpho and Rvcg - Shape Analysis in R. In Zheng G., Li S., Szekely G. (Eds.), *Statistical Shape and Deformation Analysis.* New York: Academic Press, 217–56.

---

## How to Use
1. Clone this repository:  
   ```bash
   git clone https://github.com/your-repo-name.git
