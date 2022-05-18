# R package Phub: A Regularization Approach to Component Selection in Hub Models for Network Inference

The package provides the algorithms of estimations for a series of hub models, and hub set selection approach for the hub model when hub set is unknown.
For details of hub models, see the paper ''Network Inference Using the Hub Models and Variants'' (available soon). An earlier version of the paper is available from [Arxiv](https://arxiv.org/pdf/2004.09709.pdf).

# Installation

``` r
# R CRAN version will be available soon

# Currently, install from Github
# install.packages("devtools")
devtools::install_github("zhibinghe/Phub")
```

# Usage

The main function is 'phub()', which is modified EM algorithm to perform component selection for hub set.
The function 'data.hub()' is a simple example showing the usage of 'phub()'. 
In addition, the package includes two real datasets: 'bakery5000' and 'passerines'.

``` r
#' Data Analysis using penalized hub model
#'
#' @param data  observed group data, the first \eqn{M} colunmns are corresponding to hub nodes
#' @param M size of potential hub set 
#' @param Lam a vector of tunning parameter \eqn{\lambda} 
#' @param bootstrap number of bootstrap times

#' @return a matrix of selected hub nodes for different \code{Lam}

#' @examples
#' data(bakery5000)
#' Lam = seq(0.015,0.035,0.005)
#' data.hub(bakery5000,40,Lam,0) # full data
#' data.hub(bakery5000,40,Lam,100) # bootstrapped data
data.hub(bakery5000,40,Lam,0)
```
