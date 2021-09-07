
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EC

## Overview

This package calculates the exact-corrected confidence interval [Hawila
& Berg, 2021](https://arxiv.org/pdf/2104.04660.pdf) and other confidence
intervals for noninferiorty binomial trials with various methods.
Additionally, as described in [Hawila & Berg
(2021)](https://arxiv.org/pdf/2104.04660.pdf), functions are provided to
calculate maximal sizes and p-values corresponding to different
confidence interval methods.

## Installation

``` r
library(devtools)
install_github("NourHawila/EC")
```

## Load the package

``` r
library(EC)
```

## Usage

We consider a noninferiority trial with treatment group (*T*) and
control/standard-of-care group (*C*) having a binary endpoint
representing whether or not an outcome is observed. Here we assume the
outcome is a positive event (e.g. resolution of a disease).

Let *P*<sub>*T*</sub> and *P*<sub>*C*</sub> be the probabilities the
outcome is observed, and let *δ* = *P*<sub>*T*</sub> − *P*<sub>*C*</sub>
represent the risk difference. We consider the following hypotheses:

*H*<sub>0</sub> : *δ* ≤  − *δ*<sub>0</sub>: “inferior trial”; *T* is
inferior to *C*

*H*<sub>1</sub> : *δ* &gt;  − *δ*<sub>0</sub>: “non-inferior trial”; *T*
is not inferior to *C*

Let `x.T`, `x.C`, `N.T`, `N.C`, `delta0`, and `alpha` be the
user-defined values for the noninferiority trial.

We use small sample sizes in the example below as the Chan & Zhang
confidence interval method can take a while to compute.

``` r
x.T=5
x.C=2
N.T=6
N.C=6  
delta0=.12
alpha=.05
```

### Confidence intervals

Note that only the exact corrected confidence interval, `ci_EC`, uses
the noninferiority margin *δ*<sub>0</sub> in constructing the confidence
interval. Note that *H*<sub>0</sub> is rejected if the lower bound of
the confidence interval is larger than  − *δ*<sub>0</sub>.

``` r
ci_EC(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,alpha=alpha,delta0=delta0)
#> $ci.lower
#> [1] -0.1090285
#> 
#> $ci.upper
#> [1] 0.7909168
#> 
#> $count
#> [1] 35
ci_CZ(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,alpha=alpha)
#> $ci.lower
#> [1] -0.14288
#> 
#> $ci.upper
#> [1] 0.8879
ci_MN(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,alpha=alpha)
#> $ci.lower
#> [1] -0.05788002
#> 
#> $ci.upper
#> [1] 0.8210233
#> 
#> $count
#> [1] 36
ci_Wald(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,alpha=alpha)
#> $ci.lower
#> [1] 0.01916852
#> 
#> $ci.upper
#> [1] 0.9808315
```

### p-values

As described in [Hawila & Berg,
2021](https://arxiv.org/pdf/2104.04660.pdf), we correspond a p-value
with each confidence interval method. Note that *H*<sub>0</sub> is
rejected if the p-value is smaller than *α*/2.

``` r
pval_Chan(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,delta0=delta0)
#> [1] 0.02271896
pval_CZ(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,delta0=delta0)
#> [1] 0.02994009
pval_MN(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,delta0=delta0)
#> [1] 0.01438484
pval_Wald(x.T=x.T,x.C=x.C,N.T=N.T,N.C=N.C,delta0=delta0)
#> [1] 0.005748107
```

### Maximal sizes

As described in [Hawila & Berg,
2021](https://arxiv.org/pdf/2104.04660.pdf), for a given value of
*N*<sub>*T*</sub>, *N*<sub>*C*</sub>, and *α*, we can also calculate the
maximal size for each confidence interval method. The size of the
exact-corrected method is bounded by *α*/2, but the sizes of the other
methods can be larger than *α*/2.

``` r
size_EC(alpha=alpha,N.T=N.T,N.C=N.C,delta0=delta0)
#> $level
#> [1] 0.02238953
#size_CZ(alpha=alpha,N.T=N.T,N.C=N.C,delta0=delta0)
size_MN(alpha=alpha,N.T=N.T,N.C=N.C,delta0=delta0)
#> $level
#> [1] 0.02986227
size_Wald(alpha=alpha,N.T=N.T,N.C=N.C,delta0=delta0)
#> $level
#> [1] 0.4644041
```

## Getting help

``` r
help(package="EC")
```
