## scenario #1
library(Phub)
library(foreach)
n0 = 10; n=100; T=1000
pf = 0.4; pr = 0.1
M = 80 
Iter = 50
lam = seq(0.015,0.018,0.001)
# Bic selection with optimal lambda
bic.sele = function(lam){
  model = lapply(lam,phub,G=G0,A=A,rho=rho,pen.type="plog")
  sel = matrix(NA,length(lam),(M+1))
  best = rep(NA,length(lam))
  bic = function(x,hatM) x-0.5*log(T)*(n+1)*hatM # BIC criteria
  for(i in 1:length(lam)){
    sel[i,] = model[[i]]$rho
    tc = sum(model[[i]]$rho!=0)
    best[i] = bic(model[[i]]$l,tc)
  }
  out = cbind(lam,best,sel)
  colnames(out) = c("lambda","bic",paste("rho",0:M,sep=""))
  return(out)
} 
# parallel computing
c = 50 # number of cpu cores
cl = parallel::makeCluster(c)
doParallel::registerDoParallel(cl)
out = foreach::foreach(r=1:Iter,.packages="Phub") %dopar% {
  A0 = GenA(n,n0,pf,pr,rep(0.05,n))
  G0 = GenG(A0,T,c(0.2,rep(0.8/n0,n0)))
  # initial A
  A = matrix(NA,M,n)
  for(i in 1:M)A[i,] = colSums(G0[,i]*G0)/sum(G0[,i])
  diag(A) = 1; A = rbind(runif(n),A)
  # initial rho
  rho = c(0.2,rep(0.8/M,M))
  bic.sele(lam)
}
parallel::stopCluster(cl)
save.image("s1.RData")










