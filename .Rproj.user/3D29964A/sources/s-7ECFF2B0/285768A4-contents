#' Parallel computing fdr and power of change points estimation for different \code{gamma} and \code{nu}
#'
#' @param c number of cpu cores used for parallel computing
#' @inheritParams GenDY
#' @param Nu a vector of different \code{nu}s
#' @param Gamma a vector of different \code{gamma}s
#' @inheritParams Fdr
#' @inheritParams Fdr
#' @inheritParams ch.est
#' @inheritParams which.cp
#' @param iter iteration times for each combination of \code{gamma} and \code{nu}
#' @return a list of matrix with the same length as \code{Nu}, FDR and Power for different \code{Gamma} are displayed within each matrix
#' @export
#' @import foreach
#' @import doParallel
#' @import parallel
#' @examples \donttest{
#' size=12000
#' a = 1
#' A = a*(1:119)
#' H = seq(100,11900,100)
#' mu = GenMu(A,H,size=size)
#' z = GenZ(nu=2,size=size)
#' Gamma = seq(1,5,1)
#' Nu = seq(0,2,0.5)
#' model = fdr.gam(2,mu,Gamma,Nu,8,H,iter=100)
#'}
fdr.gam = function(c,mu,Gamma,Nu,b,th,B=100,level=0.1,iter=100){
  nu = NULL
  size = length(mu)
  cl = parallel::makeCluster(c)
  doParallel::registerDoParallel(cl)
  stem = foreach::foreach(nu = Nu) %:%
    foreach::foreach(gamma=Gamma,.combine='rbind',.packages="mSTEM") %dopar% {
      chest = ch.est(nu,gamma,size,B)
      temp = matrix(NA,iter,2);colnames(temp)=c("FDR","Power")
      for (i in 1:iter){
        y1 = GenDY(mu,GenZ(nu,size),gamma)
        cp = which.cp(y1,chest,level)
        temp[i,]=unlist(Fdr(uh=c(cp$peak,cp$vall),b=b,th=th))
      }
      colMeans(temp)
    }
  parallel::stopCluster(cl)
  return(stem)
}

