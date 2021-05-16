scINSIGHT for interpreting single cell gene expression in biologically heterogeneous data
================
Kun Qian, Wei Vivian Li
2021-05-16

<!-- README.md is generated from README.Rmd. Please edit that file -->
<img src="https://github.com/Vivianstats/data-pkg/raw/main/img/scINSIGHT.png" height="200" align="right" />

## Latest News

> 2020/05/14:

-   Version 0.1.0 released!

## Introduction

scINSIGHT uses a novel matrix factorization model to jointly analyze multiple single-cell gene expression samples from biologically heterogeneous sources, such as different diseasephases, treatment groups, or developmental stages. It assumes that each gene pathway is asparse and non-negative linear combination of genes, and each cell is jointly defined by theexpression of common and condition-specific pathways. Given multiple gene expression sam-ples from different biological conditions, scINSIGHT aims to simultaneously identify common andcondition-specific gene pathways and quantify their expression levels in each sample in a lower-dimensional space .

Any suggestions on the package are welcome! For technical problems, please report to [Issues](https://github.com/Vivianstats/scINSIGHT/issues). For suggestions and comments on the method, please contact Kun (<kun_qian@foxmail.com>) or Vivian (<vivian.li@rutgers.edu>).

## Installation

You can install `scINSIGHT` from [CRAN](https://cran.r-project.org/web/packages/scINSIGHT/index.html) with:

``` r
install.packages("scINSIGHT")
```

## Usage

Please refer to the [package vignette](https://github.com/Vivianstats/scINSIGHT/wiki/scINSIGHT-vignette) for examples about how to use the package functions.
