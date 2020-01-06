#' Generate a piecewise constant sequence starting from 0
#'
#' @param x a vector containing all values of change points
#' @param pos positions of change points, corresponding to \code{x}
#' @param size sample size
#' @return a piecewise constant sequence
#' @export
#' @seealso \code{\link{GenDY}}
#' @examples
#' GenMu(x=1:10,pos=seq(10,100,10),size=150)
GenMu = function(x,pos,size){
  obj = rep(0,size)
  len = length(x)
  for (i in 1:(len-1)) obj[pos[i]:(pos[i+1]-1)] = x[i]
  obj[pos[len]:size] = x[i]
  return(obj)
}
## Generate smooth random error
kern = function(x,v=1){
  temp = exp(-x^2/(2*v^2))
  return(temp/sum(temp))
}
## differential of kernal function
kern.d = function(gamma){
  x = (-4*gamma):(4*gamma)
  return(-x*kern(x,v=gamma)/gamma^2)
}
#' Generate Gaussian autocorrelated random error sequence based on White-noise and Gaussian kernal
#'
#' @param nu bandwidth of Gaussian kernal applied to White-noise, Whitenoise error if \code{nu} = 0
#' @inheritParams GenMu
#' @return a vector of random error
#' @export
#' @seealso \code{\link{GenDY}}
#' @examples
#' GenZ(nu=2,size=1000)
GenZ = function(nu,size){
  if(nu==0) return(rnorm(size))
  else{
    x = (-4*nu):(4*nu)
    return(conv(rnorm(size),kern(x,v=nu)))}
}
#' Generate first-order differential of a smoothed sequence Y
#'
#' @param mu a vector of piecewise constant
#' @param z a vector of stationary Gaussian random error
#' @param gamma bandwidth of nonparameter smoothing
#' @return a vector of the differential of Y
#' @export
#' @seealso \code{\link{GenMu}}/\code{\link{GenZ}}
#' @examples
#' mu = GenMu(x=1:10,pos=seq(10,100,10),size=150)
#' z = GenZ(nu=2,size=150)
#' GenDY(mu=mu,z=z,gamma=4)
GenDY = function(mu,z,gamma){
  w = kern.d(gamma)
  mu_gam = conv(mu,w,shape="same")
  z_gam = conv(z,w,shape="same")
  return(mu_gam+z_gam)
}
#' Estimate \strong{s2},\strong{lambda2},\strong{lambda4},\strong{Delta}
#'
#' @inheritParams GenZ
#' @inheritParams  GenDY
#' @inheritParams GenMu
#' @param B Montelarlo iteration times
#' @references Multiple Testing of Local Extrema for Detection of Change Points \url{https://arxiv.org/abs/1504.06384}
#' @return a list of s2,lambda2,lambda4,Delta
#' @export
#' @seealso \code{\link{which.cp}}
#' @examples ch.est(nu=2,gamma=4,size=1000,B=100)
ch.est = function(nu,gamma,size,B=100){
  z1 = matrix(NA,size,B)
  for (j in 1:B) z1[,j] = conv(GenZ(nu,size),kern.d(gamma))
  s2 = mean(apply(z1,2,var))
  lambda2 = mean(apply(z1,2,function(x) var(diff(x))))
  lambda4 = mean(apply(z1,2,function(x) var(diff(diff(x)))))
  Delta = s2*lambda4 - lambda2^2
  return(list(s2=s2,lambda2=lambda2,lambda4=lambda4,Delta=Delta))
}
#' Find locations of change points
#'
#' @param y1 a vector of the differential of sequence Y
#' @param chest output of function \code{\link{ch.est}}
#' @param level FDR control level
#'
#' @return a list of components
#' \item{peak}{a vector of peaks location}
#' \item{vall}{a vector of valleys location}
#' \item{pval}{a scalar of adjusted p-value based on FDR control}
#' \item{thresh}{a scalar of threshold for \code{y1}}
#' @export
#' @seealso \code{\link{ch.est}}/\code{\link{fdrBH}}
#'
#' @examples
#' mu = GenMu(x=1:10,pos=seq(10,100,10),size=150)
#' z = GenZ(nu=2,size=150)
#' y1 = GenDY(mu,z,gamma=4)
#' chest = ch.est(nu=2,gamma=8,size=150,B=100)
#' which.cp(y1,chest,level=0.1)
which.cp = function(y1,chest,level=0.1){
  size = length(y1)
  s2 = chest$s2; lambda2 = chest$lambda2
  lambda4 = chest$lambda4; Delta = chest$Delta
  Ly_max = which.peaks(y1)
  Ly_min = which.peaks(y1,decreasing=T)
  Ty = c(y1[Ly_max],-y1[Ly_min])
  pval = 1-pnorm(Ty*sqrt(lambda4/Delta))+sqrt((2*pi*lambda2^2)/(lambda4*s2))*dnorm(Ty/sqrt(s2))*pnorm(Ty*lambda2/sqrt(s2*Delta))
  pthresh = fdrBH(pval,level)
  if(pthresh==0){
    Tthresh = Inf
    index_max = NULL
    index_min = NULL
  }
  else{
    Tthresh = approx(pval,Ty,pthresh,method="linear")$y
    index_max = intersect(which(y1>=Tthresh),Ly_max)
    index_min = intersect(which(y1<=-Tthresh),Ly_min)
  }
  return(list(peak=index_max,vall=index_min,pval=pthresh,thresh=Tthresh))
}
#' Evaluate performance of estimated change points
#'
#' @param uh a vector of estimated change points locations
#' @param b a scalar of location tolerance, specified by user
#' @param th a vector of true change points locations
#'
#' @return a list of vector of \code{FDR} and \code{Power}
#' \item{FDR}{a scalar of fdr (false discovery rate)}
#' \item{Power}{a scalar of power (true positive rate)}
#' @export
#' @seealso \code{\link{which.cp}}
#' @examples
#' Fdr(uh=c(7,15,32,47),b=4,th=c(10,20,30,40,50))
Fdr = function(uh,b,th){
  confband = function(th,b) unique(unlist(lapply(th,function(x){seq(ceiling(x-b),floor(x+b))})))
  if (length(uh)==0) {
    FDR = 0
    Power = 0}
  else{
    n.tp = sum(uh %in% confband(th,b))
    FDR = 1 - n.tp/length(uh)
    Power = min(n.tp/length(th),1)
  }
  return(list(FDR=FDR,Power=Power))
}
