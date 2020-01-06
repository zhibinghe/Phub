#' \code{\link{phub}} with multiple initial points, both \code{A} and \code{rho} are randomly drawn from Unif(0,1)  
#'
#' @inheritParams  phub
#' @param M  initial user-specifed number of components
#' @inheritParams  phub
#' @param rep repeatation times
#' @inheritParams phub
#' @inheritParams  phub
#' @inheritParams  phub
#'
#' @return a list of components
#' \item{A}{M*n matrix containg estimated correlation among nodes}
#' \item{rho}{M-dimensional vector containing estimated component weight}
#' \item{l}{loglikelihood}
#' \item{iteration}{number of iteration steps used to converge}
#' @export
#' @examples
#' n0 = 5; n=100; T=1000
#' A0 = GenA(100,n0,0.4,0.1,rep(0.05,n))
#' G0 = GenG(A0,T,c(0.2,rep(0.8/n0,n0)))
#' M = 10; mphub(G0,M,0.045)$rho
#' \donttest{M = 20; mphub(G0,M,0.035)}
#' \donttest{M = 50; mphub(G0,M,0.019)}
#'
mphub = function(G,M,lam,rep=10,pen.type=c("plog","plasso","log"),iter.max=1000,tol=1e-4){
  T = dim(G)[1]; n = dim(G)[2]
  # repeatation
  outp = list(rep); cur = rep(NA,rep)
  for(r in 1:rep){
    rho = runif(M+1); rho = rho/sum(rho)
    A = matrix(runif((M+1)*n),nrow=(M+1))
    diag(A[-1,])=1
    outp[[r]] = phub(G,A,rho,lam)
    cur[r] = outp[[r]]$l
  }
  maxid = which.max(cur)
  return(outp[[maxid]])
}
