rd <- function(x, digits=2, add=TRUE, max=(digits+3)){
  y <- round(x, digits=digits)
  yk <- format(y, nsmall=digits)
  nzero <- sum(unlist(y)==0)
  if(add==TRUE){
    while(nzero>0){
      zeros <- y==0
      digits <- digits+1
      y[zeros] <- round(x, digits=digits)[zeros]
      yk[zeros] <- format(y[zeros], nsmall=digits)
      nzero <- sum(y==0)
      if(digits>(max-1))
        nzero <- 0
    }
  }
  z <- sub("^([-]?)0[.]","\\1.", gsub(" +", "", yk))
  z
}
