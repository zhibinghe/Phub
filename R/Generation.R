#' Simulate a simple pattern 'Adjacent matrix'
#'
#' @param n  number of nodes
#' @param n0  number of true hub nodes
#' @param pi a n-dimensional vector containg the probability of appearance in nonleader group
#' @param pf probability of following
#' @param pr probality of rejection

#' @return a (n0+1)*n matrix; the first row is \code{pi}
#' @export
#'
#' @examples
#' pi = rep(0.05,100)
#' A0 = GenA(100,5,0.4,0.1,pi)
#'
GenA = function(n,n0,pf,pr,pi){
  A = matrix(NA,n0,n)
  A_lab = vector("list",n0)
  # set of followers except hub itself
  for(i in 1:n0) A_lab[[i]] = seq(n0+i,n,n0)
  for(i in 1:n0){
    A[i,A_lab[[i]]] = pf
    A[i,-A_lab[[i]]] = pr
  }
  diag(A) = 1
  return(rbind(pi,A))
}
#' Simulate grouped data
#'
#' @param A (n0+1)*n-dimensional adjacent matrix
#' @param T number of observations
#' @param rho a (n0+1) dimensional vector containing component(hub) weight
#' @return \code{T}*n matrix
#' @export
#' @seealso \code{\link{GenA}}
#'
#' @examples
#' n0 = 5
#' A0 = GenA(100,n0,0.4,0.1,rep(0.05,100))
#' rho = c(0.2,rep(0.8/n0,n0))
#' G0 = GenG(A0,1000,rho)
#'
GenG = function(A,T,rho){
  n0=dim(A)[1]; n=dim(A)[2]
  G = matrix(0,T,n)
  for(t in 1:T){
    hub = sample(0:(n0-1),1,prob=rho)
    G[t,]=rbinom(n,1,prob=A[hub+1,])
  }
  return(G)
}
