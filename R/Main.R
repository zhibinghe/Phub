#' Modified EM algorithm for hub model component selection
#'
#' @param G  Observed group data, a T*n matrix
#' @param A  M*n matrix containg correlation among nodes; first row is the nonleader case
#' @param rho M-dimensional vector composed of component weight; prior probability of components.
#' @param lam tunning parameter for component selection, lambda*Df
#' @param iter maximum iteration steps

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
#' M = 11
#' A = matrix(runif(M*n),nrow=M); diag(A[-1,])=1
#' rho = runif(M); rho = rho/sum(rho)
#' pi = colMeans(G0)
#' EM.hub(G0,A,rho,0.03)
#'
EM.hub = function(G,A,rho,lam,iter=1000){
  ## f(G(t)|Z(t)=k,A)
  cond.f = function(G,A){
    T = dim(G)[1]; M = dim(A)[1]
    Pr_cond = matrix(NA,T,M)
    for(t in 1:T) Pr_cond[t,] = apply(t(t(A)^G[t,])*t(t(1-A)^(1-G[t,])),1,prod)
    return(Pr_cond)
  }
  ## posterior prob. Htm and loglikelihood
  posterior = function(rho,Pr_cond){
    T = dim(Pr_cond)[1];M = dim(Pr_cond)[2]
    H = matrix(NA,T,M)
    Htm = t(t(Pr_cond)*rho)
    # if Gm(t)=0 then m cannot be leader in t(th) group
    ind = which(rowSums(Htm)==0)
    if(length(ind)==0) H[,] = Htm/rowSums(Htm)
    else
      H[ind,] = 0
      H[-ind,] = Htm[-ind,]/rowSums(Htm[-ind,])
    loglik = sum(log(rowSums(Htm)))
    loglik[loglik==-Inf] = -9*10^5 # avoid likelihood being -inf
    return(list(H=H,l=loglik))
  }
  ## update rho: some rho will shrink to 0
  hatrho = function(H,Mi,lam) sapply((colMeans(H) - lam)/(1-Mi*lam),function(x) max(0,x))
  ## update A: Am. will be 0 if rho_m=0
  hatA = function(G,H){
    M = dim(H)[2]; n = dim(G)[2]
    HatA = matrix(NA,M,n)
    for(m in 1:M){
      if (all(H[,m]==0)) HatA[m,] = 0
      else HatA[m,] = colSums(H[,m]*G)/sum(H[,m])
    }
    diag(HatA[-1,]) = 1 # No need to estimate A(m+1)m
    return(HatA)
  }
  ## EM Iteration
  count = 0;diff = 1
  L = vector();L[1] = 1
  while((count < iter & diff > 10^(-6))){
    count = count + 1
    # E Step
    post = posterior(rho,cond.f(G,A))
    H = post$H
    L[count+1] = post$l
    # M Step
    Mi = sum(rho!=0)
    if(Mi*lam==1) stop("Please input a new lam avoiding 1-M*lam is 0")
    rho = hatrho(H,Mi,lam)
    if(all(rho==0)) break
    A = hatA(G,H)
    diff = abs((L[count+1]-L[count])/L[count])
  }
  return(list(A=A,rho=rho,l=L[count+1],iteration=count))
}





