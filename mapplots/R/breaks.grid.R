breaks.grid <-
function(grd,quantile=0.975,ncol=12,zero=TRUE){
  
  # in order to deal with multigrid as well as grid, make it into a list
  if(is.list(grd)==F) grd <- list(grd)
  
  # work out the quantile value
  qua <- max(c(unlist(lapply(grd,quantile,probs=quantile,na.rm=T)),0),na.rm=T)
  
  # replace any values higher than quantile with quantile
  if(qua>0) grd <- lapply(grd,function(x) ifelse(x>qua,qua,x))
  
  # breakpoints
  hi <- max(unlist(lapply(grd,max,na.rm=T)))
  if(zero) {
    len <- ncol
    breaks <- c(0,seq(0,hi,length=len))
  } else {
    len <- ncol+1
    breaks <- seq(0,hi,length=len)
  }
  
  return(breaks)
}

