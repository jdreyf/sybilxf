# sybilxf
R package for metabolic modeling of Seahorse extracellular flux (XF) data. Used in the publication "Integrating Extracellular Flux Measurements and Genome-Scale Modeling Reveals Differences between Brown and White Adipocytes" (2017) by AK Ramirez, ..., JM Dreyfuss (http://www.cell.com/cell-reports/fulltext/S2211-1247(17)31717-5).

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

The package is on GitHub, so you can install it with `remotes`:
```
# install.packages("remotes")
library("remotes")
remotes::install_github("jdreyf/sybilxf")
```
See the vignette for a tutorial.

Notes:
- The LICENSE applies to the code, but not to the metabolic models in the data directory.
- This package is at major version zero (0.y.z), which is for initial development, as per https://semver.org. Anything may change at any time. The public API should not be considered stable.
