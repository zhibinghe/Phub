# R package Phub: A Regularization Approach to Component Selection in Hub Models for Latent Network Inference

The package provides the algorithms of estimations for  hub model and its variants, and hub set selection approach for the hub models when hub set is unknown. For details of hub models, see the papers ''[Network inference from grouped observations using hub models](https://www3.stat.sinica.edu.tw/statistica/oldpdf/A29n112.pdf)''  and "[Network Inference Using the Hub Model and Variants](https://www.tandfonline.com/doi/full/10.1080/01621459.2023.2183133)".

# Installation

``` r
# install 'Phub' from Github
# install.packages("devtools")
devtools::install_github("zhibinghe/Phub")
```

# Usage

The main function is 'phub()', which is modified EM algorithm to perform component selection for hub set.
The function 'data.hub()' is a simple example showing the usage of 'phub()'. 
In addition, the package includes two real datasets: 'bakery5000' and 'passerines'.

``` r
#' @examples
set.seed(2020)
n0 = 5; n=100; T=1000
A0 = GenA(n,n0,0.4,0.1,rep(0.05,n))
G0 = GenG(A0,T,c(0.2,rep(0.8/n0,n0)))
M = 10 
A = matrix(runif((M+1)*n),nrow=(M+1)); diag(A[-1,]) = 1
rho = runif(M+1); rho = rho/sum(rho)
phub(G0,A,rho,0.030,pen.type="log")$rho
phub(G0,A,rho,0.030,pen.type="plog")$rho
phub(G0,A,rho,0.8,pen.type="plasso")$rho

#' Data Analysis using penalized hub model
#' @param data  observed group data, the first \eqn{M} colunmns are corresponding to hub nodes
#' @param M size of potential hub set 
#' @param Lam a vector of tunning parameter \eqn{\lambda} 
#' @param bootstrap number of bootstrap times
#' @return a matrix of selected hub nodes for different lambda

#' @examples
data(bakery5000)
Lam = seq(0.015,0.035,0.005)
data.hub(bakery5000,40,Lam,0) # full data
data.hub(bakery5000,40,Lam,100) # bootstrapped data
data.hub(bakery5000,40,Lam,0)
```
