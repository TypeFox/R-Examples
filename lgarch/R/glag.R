glag <-
function(x, k=1, pad=TRUE, pad.value=NA)
{
  #check arguments:
  if(k < 1) stop("Lag order k cannot be less than 1")

  #zoo-related:
  iszoo.chk <- is.zoo(x)
  x <- as.zoo(x)
  x.index <- index(x)
  x <- coredata(x)
  isvector <- is.vector(x)
  x <- cbind(x)
  x.n <- NROW(x)
  x.ncol <- NCOL(x)

  #do the lagging:
  x.nmink <- x.n - k
  xlagged <- matrix(x[1:x.nmink,], x.nmink, x.ncol)
  if(pad){
    xlagged <- rbind( matrix(pad.value,k,x.ncol) , xlagged)
  }

  #transform to vector?:
  if(x.ncol==1 && isvector==TRUE){
    xlagged <- as.vector(xlagged)
  }

  #if(is.zoo(x)):
  if(iszoo.chk){
    if(pad){
      xlagged <- zoo(xlagged, order.by=x.index)
    }else{
      xlagged <- zoo(xlagged, order.by=x.index[c(k+1):x.n])
    } #end if(pad)
  } #end if(iszoo.chk)

  #out:
  return(xlagged)
}
