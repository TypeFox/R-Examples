Power <-
function(x){
   n.sim <- x$n.sim[1]
   n.ss <- x$n.ss   
   x$SeriesPerInd <- x$Series/x$Individuals
   SimVCVInd <- as.data.frame(matrix(unlist(x$SimVCVInd),n.ss,4, byrow=TRUE))
####
   PI <- unlist(lapply(x$pwI, sum))/x$n.sim[1]
   PS <- unlist(lapply(x$pwS, sum))/x$n.sim[1]
   Power1 <- data.frame(VarIDIntercept=SimVCVInd[,1], PowerVarIDIntercept=PI, VarIDSlope=SimVCVInd[,4], PowerVarIDSlope=PS, Individuals=x$Individuals, Series=x$Series, SeriesPerInd=x$SeriesPerInd, n.obs=unlist(x$n.obs))
   print(Power1)   
}
