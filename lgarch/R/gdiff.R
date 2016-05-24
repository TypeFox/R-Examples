gdiff <-
function(x, lag=1, pad=TRUE, pad.value=NA)
{
  #check arguments:
  if(lag < 1) stop("lag argument cannot be less than 1")

  #zoo-related:
  iszoo.chk <- is.zoo(x)
  x <- as.zoo(x)
  x.index <- index(x)
  x <- coredata(x)
  isvector <- is.vector(x)
  x <- cbind(x)
  x.n <- NROW(x)
  x.ncol <- NCOL(x)

  #do the differencing:
  xdiff <- x - glag(x, k=lag, pad=TRUE, pad.value=NA)

  #pad operations:
  if(!pad){
    xdiff <- na.trim(as.zoo(xdiff))
    xdiff <- coredata(xdiff)
  }else{
    whereisna <- is.na(xdiff)
    xdiff[whereisna] <- pad.value
  }

  #transform to vector?:
  if(x.ncol==1 && isvector==TRUE){
    xdiff <- as.vector(xdiff)
  }

  #if(is.zoo(x)):
  if(iszoo.chk){
    if(pad){
      xdiff <- zoo(xdiff, order.by=x.index)
    }else{
      xdiff <- zoo(xdiff, order.by=x.index[c(lag+1):x.n])
    } #end if(pad)
  } #end if(iszoo.chk)

  #out:
  return(xdiff)
}
