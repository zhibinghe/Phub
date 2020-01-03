#' Multiple initial points (rho) version of function \code{\link{EM.hub}} (Deprecated now)
#'
#' @inheritParams  EM.hub
#' @param M  initial user-specifed number of components
#' @inheritParams  EM.hub
#' @param rep repeatation
#' @inheritparam EM.hub
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
#' M = 10; EMM.hub(G0,M,0.045)$rho
#' M = 20; EMM.hub(G0,M,0.035)$rho
#' \donttest{M = 50; EMM.hub(G0,M,0.019)}
#'
EMM.hub = function(G,M,lam,rep=10,nohub=TRUE,iter=1000,tol=10^(-6)){
  T = dim(G)[1]; n = dim(G)[2]
  # repeatation
  outp = list(rep); cur = rep(NA,rep)
  for(r in 1:rep){
    rho = runif(M+1); rho = rho/sum(rho)
    A = matrix(runif((M+1)*n),nrow=(M+1))
    diag(A[-1,])=1
    outp[[r]] = EM.hub(G,A,rho,lam)
    cur[r] = outp[[r]]$l
  }
  maxid = which.max(cur)
  return(outp[[maxid]])
}
