#' Modified EM algorithm for hub model component selection
#'
#' @param G  observed group data, a T*n matrix
#' @param A  a matrix containg correlation among nodes; first row is the nonleader case
#' @param rho a vector composed of component weight; prior probability of components.
#' @param lam tuning parameter for component selection, degenrate to a standard EM without penalty if lam=0
#' @param pen.type type of penalty, including 'log': penalization for logarithm of all components;
#' 'plog': penalization for logarithm of partial components, rho0 (nonhub component) is always kept;
#' 'plasso': penalization for lasso form of partial compents
#' @param iter.max maximum iteration steps
#' @param tol threshold for shrinking rho to 0

#' @return a list of components
#' \item{A}{a matrix containg estimated correlation among nodes}
#' \item{rho}{a vector containing estimated component weight}
#' \item{l}{log-likelihood}
#' \item{iteration}{number of iterations used to converge}
#' @export
#' @import Rsolnp
#' @examples
#' n0 = 5; n=100; T=1000
#' A0 = GenA(100,n0,0.4,0.1,rep(0.05,n))
#' G0 = GenG(A0,T,c(0.2,rep(0.8/n0,n0)))
#' M = 8
#' A = matrix(runif((M+1)*n),nrow=(M+1)); diag(A[-1,])=1
#' rho = runif(M+1); rho = rho/sum(rho)
#' phub(G0,A,rho,0.035,pen.type="log")$rho
#' phub(G0,A,rho,0.035,pen.type="plog")$rho
#' phub(G0,A,rho,0.08,pen.type="plasso")$rho
phub = function(G,A,rho,lam,pen.type=c("plog","plasso","log"),iter.max=1000,tol=1e-4){
  ## f(G(t)|Z(t)=k,A)
  pen.type = match.arg(pen.type)
  posterior = function(G,A,rho){
    T = dim(G)[1]; q = dim(A)[1]
    Pr_cond = matrix(NA,T,q)
    for(t in 1:T) Pr_cond[t,] = apply(t(t(A)^G[t,])*t(t(1-A)^(1-G[t,])),1,prod)
  ## posterior probability H and loglikelihood
    id = which(rowSums(t(t(Pr_cond)*rho))==0)
    if(length(id)!=0){
      Pr_cond = Pr_cond[-id,]
      G = G[-id,] # delete bad smaple points
    } 
    Htm = t(t(Pr_cond)*rho)
    H = Htm/rowSums(Htm)
    loglike = sum(log(rowSums(Htm)))
    return(list(H=H,l=loglike,G=G))
  }
  ## update rho: some rho will shrink to 0
  if(pen.type=="log"){
    hatrho = function(rho,H){
      Mi = sum(rho!=0)
      if(Mi*lam==1) stop("\nPlease input a new lam avoiding 1-M*lam is 0\n",call.=FALSE)
      pars = sapply((colMeans(H) - lam)/(1-Mi*lam),function(x) max(0,x))
      pars[pars<tol] = 0
      return(pars)
    } 
  } 
  else{
    hatrho = function(rho,H){
      T = dim(G)[1]
      Hx = colSums(H)
      ind = which(rho!=0); len = length(ind)
      eqft = function(t) sum(t)
      if(pen.type=="plasso"){
        # convex function
        ft = function(t) -sum(Hx[ind]*log(t)) + T*lam*(1-t[1])
        pars = Rsolnp::solnp(pars=rho[ind],fun=ft,eqfun=eqft,eqB=1,LB=rep(0,len),
                       UB=rep(1,len),control=list(trace=0))$pars
      }
      else{
        epsi = 1e-8
        ans = list()
        ft = function(t) -sum(Hx[ind]*log(t)) + T*lam*sum(log(epsi+t[-1]))
        fun = function(i){
          if(i==1) par = rho[ind]
          else {par = runif(len); par = par/sum(par)} # normalized
          sol = Rsolnp::solnp(pars=par,fun=ft,eqfun=eqft,eqB=1,LB=rep(0,len),
                        UB=rep(1,len),control=list(trace=0))
          ans$pars = sol$pars
          ans$values = tail(sol$values,1)
          return(ans)
        }
        # adding one initial values rho from the last step.
        fval = lapply(as.list(1:3),fun) # repeatation
        best = sapply(fval,function(x) x$values)
        pars = fval[[which.min(best)]]$pars
      }
      pars[pars<tol] = 0
      pars[1] = 1 - sum(pars[-1]) # sum is 1
      rho[ind] = pars
      return(rho)
    }
  }
  ## update A: Am. will be 0 if rho_m=0
  hatA = function(G,H){
    q = dim(H)[2]; n = dim(G)[2]
    HatA = matrix(NA,q,n)
    for(m in 1:q){
      if (all(H[,m]==0)) HatA[m,] = 0
      else HatA[m,] = colSums(H[,m]*G)/sum(H[,m])
    }
    diag(HatA[-1,]) = 1 # Amm must be 1
    return(HatA)
  }
  ## EM Iteration
  count = 0;diff = 1
  L = vector();L[1] = 1
  while(count < iter.max & diff > 1e-6){
    count = count + 1
    # E Step
    post = posterior(G,A,rho)
    H = post$H
    G = post$G # some rows may be deleted
    L[count+1] = post$l
    if(is.infinite(L[count+1])) break
    # M Step
    rho = hatrho(rho,H)
    A = hatA(G,H)
    diff = abs((L[count+1]-L[count])/L[count+1])
  }
  return(list(A=A,rho=rho,l=L[count+1],iter=count))
}




