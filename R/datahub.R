#' Data Analysis using penalized hub model
#'
#' @param data  observed group data, the first \eqn{M} colunmns are corresponding to hub nodes
#' @param M size of potential hub set 
#' @param Lam a vector of tunning parameter \eqn{\lambda} 
#' @param bootstrap number of bootstrap times

#' @return a matrix of hub nodes weight for different \code{Lam}
#' @export
#' @import Rsolnp
#' @examples
#' data(bakery5000)
#' Lam = seq(0.015,0.035,0.005)
#' data.hub(bakery5000,40,Lam,0) # full data
#' data.hub(bakery5000,40,Lam,2) # bootstrapped data
#'   
data.hub = function(data,M,Lam,bootstrap=0){
  data = as.matrix(data)
  n = ncol(data)
  # 
  A = matrix(NA,M,n)
  for (i in 1:M) A[i,] = colSums(data[,i]*data)/sum(data[,i])
  diag(A) = 1; A = rbind(runif(n),A)
  rho = c(0.2,rep(0.8/M,M))
  #
  coef = matrix(NA,length(Lam),M+1)
  colnames(coef) = paste("H",0:M,sep="")
  rownames(coef) = paste("Lam",Lam,sep="")
  for (i in 1:length(Lam)){
    if (bootstrap == 0){
      coef[i,] = (phub(data,A,rho,Lam[i],pen.type="plog")$rho > 0) + 0
    }
    else{
      out = matrix(0,bootstrap,M+1)
      for(j in 1:bootstrap){
        sub = data[sample(1:nrow(data),nrow(data),replace=T),]
        out[j,] = (phub(sub,A,rho,Lam[i],pen.type="plog")$rho > 0) + 0 
      }
      coef[i,] = colMeans(out)
    }
  }
  return(coef)
}

