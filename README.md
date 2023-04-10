
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hlatools

<!-- badges: start -->
<!-- badges: end -->

-   `hlatools` is a companion to the repository
    <https://github.com/BenSolomon/hla_benchmark>
-   The `hla_benchmark` repository contains all the code necessary to
    reproduce the analyses in the manuscript [“Prediction of HLA
    genotypes from single-cell transcription
    data”](https://www.biorxiv.org/content/10.1101/2022.06.09.495569v1)
-   In contrast, the `hlatools` package provides convenient
    implementation of the two models described in the manuscript
    -   Prediction of HLA-DRB345 copy numbers from read mapping
        statistics
        -   This corresponds to the `model_HLAD_DRB_kNN` notebook in the
            `hla_benchmark` repository
    -   Generation of composite “best” HLA genotype from the predictions
        of multiple genotyping algorithms
        -   This corresponds to the `model_genotype_composite` notebook
            in the `hla_benchmark` repository

## Installation

You can install the development version of hlatools from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BenSolomon/hlatools")
```

#### Issues with `kknn` implementation

`predict_drb345` uses a k-nearest neighbor model implemented using
`kknn`. However, namespace issues result in the following error whenever
the `CRAN` or `dev` versions of `kknn` are installed:

``` r
Error in get(ctr, mode = "function", envir = parent.frame()): object 'contr.dummy' of mode 'function' was not found
```

[A pull request has been
submitted](https://github.com/KlausVigo/kknn/pull/24) that fixes this
issue. The requirements for `kknn` in this package’s `DESCRIPTION` file
point to this remote so that it should install `kknn` with the proposed
fix. However, if you still receive this error, consider manually
installing the `kknn` pull request with:

``` r
devtools::install_github("https://github.com/KlausVigo/kknn/", ref = devtools::github_pull("24"))
```

## `predict_drb345`

This function predicts the copy number of HLA-DRB345 alleles based on
sequence fragment mapping statistics.

It takes input data in the following form:

``` r
library(hlatools)
data("drb345_data")
head(drb345_data)
#>        sample  DRB1  DRB3  DRB4 DRB5
#> 1 INCOV005-BL 34855  9486     1    0
#> 2 INCOV003-AC 66561  6892 12368    1
#> 3 INCOV003-BL 69175  7437 15819    2
#> 4 INCOV006-BL 48887     6    23    3
#> 5 INCOV005-AC 84603 23983     2    1
#> 6 INCOV007-AC 55972 12399     0    2
```

And is used to generating the following:

``` r
library(hlatools)
data("drb345_data")
drb345_pred <- predict_drb345(drb345_data)
head(drb345_pred)
#> # A tibble: 6 × 4
#>   sample       DRB3  DRB4  DRB5
#>   <chr>       <dbl> <dbl> <dbl>
#> 1 INCOV005-BL     2     0     0
#> 2 INCOV003-AC     1     1     0
#> 3 INCOV003-BL     1     1     0
#> 4 INCOV006-BL     1     0     0
#> 5 INCOV005-AC     2     0     0
#> 6 INCOV007-AC     2     0     0
```

These copy number predictions can then be used to filter genotype
predictions for HLA-DRB345 to keep only the relevant number of allele
copies.

## `predict_composite`

This function predicts the most accurate allele prediction for a given
HLA locus and field level when predictions from multiple genotyping
algorithms are available.

It takes input data in the following form:

``` r
library(hlatools)
data("composite_data")
head(composite_data)
#> # A tibble: 6 × 5
#>   sample      locus field   genotyper allele
#>   <chr>       <chr> <chr>   <chr>     <chr> 
#> 1 INCOV003-AC A     field_1 arcasHLA  26    
#> 2 INCOV003-AC A     field_1 arcasHLA  24    
#> 3 INCOV003-AC A     field_1 hlaminer  24    
#> 4 INCOV003-AC A     field_1 hlaminer  26    
#> 5 INCOV003-AC A     field_1 optitype  <NA>  
#> 6 INCOV003-AC A     field_1 phlat     24
```

And is used to generating the following:

``` r
library(hlatools)
data("composite_data")
composite_pred <- predict_composite(composite_data)
head(composite_pred)
#>        sample locus   field best_genotyper
#> 1 INCOV003-AC     A field_1          phlat
#> 2 INCOV003-AC     A field_2          phlat
#> 3 INCOV003-AC     A field_3          phlat
#> 4 INCOV003-AC     B field_1          phlat
#> 5 INCOV003-AC     B field_2          phlat
#> 6 INCOV003-AC     B field_3       arcasHLA
```
