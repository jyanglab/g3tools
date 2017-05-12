[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# g3tools

This is an R packages for **Genomics**, **quantGen**, and **popGen** studies, especially for crop species. This package was intended for internal lab usage. It has not been extensively tested. Use at your own risk.

## Installation

Install [devtools](https://github.com/hadley/devtools) first, and then use `devtools` to install `g3tools` from github.

```R
#install.packages(devtools)
devtools::install_github("jyanglab/g3tools")
library(g3tools)
```

List all the functions in the package and find help.

```R
ls(getNamespace("g3tools"), all.names=FALSE)
```

## List of Functions

### Phenotype Simulations from genotypic data
1. `simcrop`: this function can simulate phenotypes with user specified distribution. (normal, gamma, or uniform).
2. `sim_qtl_pheno`: simulate phenotype for given number of QTLs.

### Quantative Genetics
1. `get_BLUP`: Estimate best linear unbiased predictor (BLUP) of the phenotypic values. 
2. `get_H2`: get broad sense heritability from BLUP model.



## Documentation

Stay tuned. It is under development!

## License
It is a free and open source software, licensed under [GPLv3](LICENSE).

