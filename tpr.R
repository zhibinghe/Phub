## Bic
load("scenario1.RData")
#
temp = lapply(out,function(x) x[which.max(x[,2]),-c(1,2)])
#temp = lapply(out,function(x) x[3,-c(1,2)])
fun = function(x){
  tp = sum(which(x!=0) %in% c(2:(n0+1)))
  fp = sum(x!=0) -tp -1
  tpr = tp/n0
  #fpr = fp/(n-n0)
  fpr = fp/(M-n0)
  return(c(tpr,fpr))
}

tt = lapply(temp,fun)
xx = do.call(rbind,tt)
