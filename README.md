synchrony R package
===================

Download and Install
--------------------

To download the development version of the package, type the following at the R command line:

``` r
install.packages("devtools")
devtools::install_github("tgouhier/synchrony")
```

To download the release version of the package on CRAN, type the following at the R command line:

``` r
install.packages("synchrony")
```

About synchrony
---------------

The synchrony package contains methods for computing spatial, temporal, and spatiotemporal statistics
such as:

 - empirical univariate, bivariate and multivariate variograms
 - fitting variogram models
 - phase locking and synchrony analysis
 - generating autocorrelated and cross-correlated matrices

The package is fully described in Gouhier and Guichard (2014).

References
----------

Gouhier, T. C., and F. Guichard. 2014. **Synchrony: quantifying variability in space and time**. Methods in Ecology and Evolution 5:524–533.
