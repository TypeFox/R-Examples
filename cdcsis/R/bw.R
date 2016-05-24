bw <-
function(x, y, z,index=1){ 
  cdcov_mean<-function(width) {mean(cdcor(x, y, z,width,index)$mcdcor)}
  width=optimize(cdcov_mean,interval=c(0,1))$minimum
  return(width)
}
