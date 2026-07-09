# R-codes-UPRA

This code package accompanies the manuscript

**A General U-Statistic Framework for High-Dimensional Multiple Change-Point Analysis**

## Supplementary R Code

This supplementary code package provides an implementation of the proposed moving-window U-statistic framework for high-dimensional multiple change-point analysis. The code includes testing, initial estimation, projection-based refinement, and confidence interval construction for multiple change-point locations.

## Overview

The code package contains:

-   Single-bandwidth implementation of the proposed procedure;
-   Multiscale bandwidth implementation for more stable practical performance;
-   Linear-kernel and sign-kernel versions of the moving-window U-statistic;
-   Initial change-point estimation based on bootstrap critical values;
-   Projection-based refinement using estimated active coordinates;
-   Confidence interval construction based on the argmax limiting distribution;
-   Data generation and auxiliary functions for simulation studies.

## File Structure

mosum.single.R \# Single-bandwidth implementation mosum.multiscale.R \# Multiscale bandwidth implementation

my.function.single.R \# Auxiliary functions for the single-bandwidth version my.function.multiscale.R \# Auxiliary functions for the multiscale version

argmax.R \# Simulation of the argmax limiting distribution covariance.est.R \# Estimation of covariance matrices for inference

single.bandwith.example.R \# Example script for the single-bandwidth procedure multiscale.bandwith.example.R \# Example script for the multiscale procedure

## Main Functions

The single-bandwidth procedure is implemented in:

``` r
mosum()
```

This function performs:

-   bootstrap-based change-point testing;
-   initial change-point estimation;
-   projection-based refinement;
-   confidence interval construction.

Example usage:

```         
source("my.function.single.R")
source("argmax.R")
source("covariance.est.R")
source("mosum.single.R")

dat <- Data_gen_multicpt(
  n = 600,
  p = 20,
  q0_list = c(0.3, 0.5, 0.7),
  cp_exit = TRUE,
  light_tail = TRUE,
  degree = 3,
  outlier = FALSE
)

res <- mosum(
  dat$X,
  bandwith = 80,
  B = 50,
  type = "sign",
  paral = "F"
)
```

The multiscale procedure is implemented in:

``` r
mosum.multiscale()
```

This function aggregates information from multiple bandwidths and then performs initial estimation, refinement, and confidence interval construction.

Example usage:

```         
source("my.function.multiscale.R")
source("argmax.R")
source("covariance.est.R")
source("mosum.multiscale.R")

dat <- Data_gen_multicpt(
  n = 600,
  p = 20,
  q0_list = c(0.3, 0.5, 0.7),
  cp_exit = TRUE,
  light_tail = TRUE,
  degree = 3,
  outlier = FALSE
)

res <- mosum.multiscale(
  dat$X,
  bandwidths = c(60, 80, 100),
  B = 50,
  type = "sign",
  paral = "F"
)
```

## Kernel Options

The argument type controls the kernel used in the moving-window U-statistic:

```         
type = "linear"    # linear kernel h(x,y) = y - x
type = "sign"      # sign kernel h(x,y) = sign(y - x)
```

The linear kernel is suitable for mean-type changes, while the sign kernel provides a rank-based robust alternative for heavy-tailed data or outlier contamination.

For variance-change simulations, one may first transform the data by taking componentwise squares and then apply the same mean-change procedure to the transformed data.

## Output

The main functions return a list containing:

```         
mosum.sta          # testing statistic
mosum.reject       # rejection indicator
cpt.est.initial    # initial estimated change-point locations
cpt.est.refine     # refined estimated change-point locations
cpt.conf           # confidence intervals for refined change-point locations
```

For the multiscale version, the output additionally includes:

```         
anchors            # anchor candidates across bandwidths
pooled.candidates  # pooled candidates from all bandwidths
merged.initial     # merged initial candidates used for refinement
```

## Reproducing Example Runs

To run the single-bandwidth example:

``` bash
Rscript single.bandwith.example.R
```

To run the multiscale bandwidth example:

``` bash
Rscript multiscale.bandwith.example.R
```

The example scripts generate simulated data, apply the proposed method, and return estimated change-point locations, refined estimates, and confidence intervals.

## Required Packages

The following R packages are required:

-   MASS
-   parallel

## Notes

-   The argument B denotes the number of bootstrap replications.
-   The argument bandwith in the single-bandwidth function specifies the moving-window bandwidth.
-   The argument bandwidths in the multiscale function specifies a collection of bandwidths, for example c(60, 80, 100).
-   The argument paral controls whether parallel computation is used.
-   Random seeds can be set before running the example scripts for reproducibility.

## Reproducibility

The provided scripts are designed so that:

-   each run corresponds to one simulation setting;
-   the single-bandwidth and multiscale procedures can be run separately;
-   repeating the scripts with the same random seed reproduces the corresponding numerical results.

## Correspondence

Yufeng Liu John D. MacArthur Professor of Statistics Department of Statistics University of Michigan 461 West Hall 1085 South University Ann Arbor, MI 48109-1107 USA

E-mail: [yufliu\@umich.edu](mailto:yufliu@umich.edu){.email}
