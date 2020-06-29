This repository is part of the supplementary material for the paper:
“Distinguishing Cause from Effect Using Quantiles: Bivariate Quantile Causal Discovery”

The file R/bqcd.R contains an implementation of our causal discovery method.

EXAMPLES
———————————————————
See R/benchmark-xyz for examples.
Benchmarks can be run from root directory.


CONTENTS
———————————————————
.
|-- R
|   |-- baselines
|   |   |-- GRAN_EMD
|   |   `-- Slope
|   |       |-- Slope.R
|   |       |-- resit
|   |-- baselines.R
|   |-- benchmark-ANLSMN.R
|   |-- benchmark-SIM.R
|   |-- benchmark-Tuebingen.R
|   |-- generate-ANLSMN-data.R
|   |-- bqcd.R
|   `-- utils
|       |-- data_generators.R
|       `-- read_tueb.R
|-- README
|-- code.Rproj
|-- data
|   |-- ANLSMN_pairs
|   |   |-- AN
|   |   |-- AN-s
|   |   |-- LS
|   |   |-- LS-s
|   |   `-- MN-U
|   |-- SIM_pairs
|   |   |-- SIM
|   |   |-- SIM-G
|   |   |-- SIM-c
|   |   `-- SIM-ln
|   `-- tuebingen_benchmark
`-- results

PACKAGES REQUIRED
———————————————————
install.packages(c("rvinecopulib", "statmod", "qrnn", "quantregForest")) # for bQCD
install.packages("MASS") # for the data generators
install.packages(c("CAM", "Hmisc", "np", "energy", "data.table", 
                   "boot", "energy", "gptk", "clue")) # for the baselines


