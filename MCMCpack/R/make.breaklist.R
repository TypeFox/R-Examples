"make.breaklist" <- function(BF, threshold=3){
  eBF <- exp(BF)
  N <- nrow(BF)
  out <- rep(NA, N)
  for (i in 1:N){
    order.i <- order(eBF[i,], decreasing=TRUE)
    if(sum(is.na(eBF[i,]))>0){
      out[i] <- 1
    }
    else{
      if(eBF[i, order.i[1]] / eBF[i, order.i[2]] > threshold){
        out[i] <- order.i[1]
      }
      else if (eBF[i, order.i[1]] / eBF[i, order.i[2]] < threshold & order.i[1]<=order.i[2]){
        out[i] <- order.i[1]
      }
      else if (eBF[i, order.i[1]] / eBF[i, order.i[2]] < threshold & order.i[1]>order.i[2]){
        out[i] <- order.i[2]
      }
      else{
        cat("\n Error occurs at i = ", i)
      }
    }
  }
  return(out-1)
}
