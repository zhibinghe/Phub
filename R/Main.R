#' Modified EM algorithm for hub model component selection
#'
#' @param G  observed group data, a T*n matrix
#' @param A  a matrix containg correlation among nodes; first row is the nonleader case
#' @param rho a vector composed of component weight; prior probability of components.
#' @param lam tuning parameter for component selection, lambda*Df; if lam=0, becomes a normal EM without penalty
#' @param nohub True/False indicates if rho0 should be penalized
#' @param iter maximum iteration steps
#' @param tol tolerated

#' @return a list of components
#' \item{A}{a matrix containg estimated correlation among nodes}
#' \item{rho}{a vector containing estimated component weight}
#' \item{l}{log-likelihood}
#' \item{iteration}{number of iterations used to converge}
#' @export
#'
#' @examples
#' n0 = 5; n=100
#' A0 = GenA(100,n0,0.4,0.1,rep(0.05,n))
#' G0 = GenG(A0,1000,c(0.2,rep(0.8/n0,n0)))
#' M = 10
#' A = matrix(runif((M+1)*n),nrow=(M+1)); diag(A[-1,])=1
#' rho = runif(M+1); rho = rho/sum(rho)
#' EM.hub(G0,A,rho,0.045)$rho
#'
EM.hub = function(G,A,rho,lam,nohub=TRUE,iter=1000,tol=10^(-6)){
  ## f(G(t)|Z(t)=k,A)
  cond.f = function(G,A){
    T = dim(G)[1]; q = dim(A)[1]
    Pr_cond = matrix(NA,T,q)
    for(t in 1:T) Pr_cond[t,] = apply(t(t(A)^G[t,])*t(t(1-A)^(1-G[t,])),1,prod)
    return(Pr_cond)
  }
  ## posterior prob. Htm and loglikelihood
  posterior = function(rho,Pr_cond){
    T = dim(Pr_cond)[1]; q = dim(Pr_cond)[2]
    H = matrix(NA,T,q)
    Htm = t(t(Pr_cond)*rho)
    H = Htm/rowSums(Htm)
    loglik = sum(log(rowSums(Htm)))
    return(list(H=H,l=loglik))
  }
  ## update rho: some rho will shrink to 0
  hatrho = function(H,Mi,lam){
    temp = sapply((colMeans(H) - lam)/(1-Mi*lam),function(x) max(0,x))
    # whether to penalized rho0
    if(nohub==TRUE) temp[1] = mean(H[,1])/(1-Mi*lam)
    return(temp)
  }
  ## update A: Am. will be 0 if rho_m=0
  hatA = function(G,H){
    q = dim(H)[2]; n = dim(G)[2]
    HatA = matrix(NA,q,n)
    for(m in 1:q){
      if (all(H[,m]==0)) HatA[m,] = 0
      else HatA[m,] = colSums(H[,m]*G)/sum(H[,m])
    }
    diag(HatA[-1,]) = 1 # No need to estimate A(m+1)m
    return(HatA)
  }
  ## EM Iteration
  count = 0;diff = 1
  L = vector();L[1] = 1
  while(count < iter & diff > 0.01){
    count = count + 1
    # E Step
    post = posterior(rho,cond.f(G,A))
    H = post$H
    L[count+1] = post$l
    if(is.infinite(L[count+1])) break
    # M Step
    Mi = sum(rho!=0)
    if(Mi*lam==1) stop("Please input a new lam avoiding 1-M*lam is 0")
    rho = hatrho(H,Mi,lam)
    A = hatA(G,H)
    diff = abs(L[count+1]-L[count])
  }
  ## delete rho when too small
  A[(rho<tol),] = 0
  rho[rho<tol] = 0
  return(list(A=A,rho=rho,l=L[count+1],iteration=count))
}




