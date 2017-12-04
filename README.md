# densityEst

This package considers the problem of estimating a time series of density functions. A data set comprising many observations is recorded at each time period, and the associated probability density function is to be estimated for each time period. The principal motivation of this package arises due to the non-existence of previous studies on nonparametric estimation of a time series of densities, taking account of the time ordering of the densities. Previous studies ignored the time ordering of the densities and estimated the densities independently for every particular year, ignoring all other years. This package will address this issue by nonparametrically estimating a time series of density functions using logspline approach considering the time ordering of the densities.


# Installation

You can install the development version from Github

# install.packages("devtools")
devtools::install_github("ThilakshaSilva/densityEst")

