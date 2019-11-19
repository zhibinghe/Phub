#' Multiple initial points (rho) version of function \code{\link{EM.hub}} (Deprecated now)
#'
#' @inheritParams  EM.hub
#' @param M  initial user-specifed number of components
#' @inheritParams  EM.hub
#' @inheritParams  EM.hub
#' @inheritParams  EM.hub
#'
#' @return a list of components
#' \item{A}{M*n matrix containg estimated correlation among nodes}
#' \item{rho}{M-dimensional vector containing estimated component weight}
#' \item{l}{loglikelihood}
#' \item{iteration}{number of iteration steps used to converge}
#' @export
#'
#' @examples
#' n0 = 5; n=100
#' A0 = GenA(100,n0,0.4,0.1,rep(0.05,n))
#' G0 = GenG(A0,1000,c(0.2,rep(0.8/n0,n0)))
#' M = 10; EMM.hub(G0,M,0.04)
#' M = 20; EMM.hub(G0,M,0.035)
#' \donttest{M = 50; EMM.hub(G0,M,0.019)}
#'
EMM.hub = function(G,M,lam,iter=1000,tol=0.01){
  T = dim(G)[1]; n = dim(G)[2]
  A = matrix(NA,M,n)
  A[1,] = colMeans(G) # nonleader group
  for(m in 2:M){
    if(sum(G[,(m-1)])==0) A[m,] = 0
    else A[m,] = colSums(G[,(m-1)]*G)/sum(G[,(m-1)])
  }
  diag(A[-1,]) = 1
  ## repeatation
  outp = list(10); cur = rep(NA,10)
  for(r in 1:10){
    rho = runif(M); rho = rho/sum(rho)
    outp[[r]] = EM.hub(G,A,rho,lam)
    cur[r] = outp[[r]]$l
  }
  maxid = which.max(cur)
  return(outp[[maxid]])
}
