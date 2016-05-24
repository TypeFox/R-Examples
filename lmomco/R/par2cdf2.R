"par2cdf2" <-
function(x, leftpara, rightpara, weight=NULL, ...) {
    end.min <-   .Machine$double.eps
    end.max <- 1-.Machine$double.eps
    qua.min <- par2qua2(end.min,leftpara,rightpara,...)
    qua.max <- par2qua2(end.max,leftpara,rightpara,...)
    f <- vector(mode="numeric", length=length(x))
    for(i in seq(1,length(x))) {
      QUAx <- x[i]
      if(QUAx <= qua.min) { f[i] <- end.min; next }
      if(QUAx >= qua.max) { f[i] <- end.max; next }

      "fn" <- function(F) {
        qua <- par2qua2(F, leftpara, rightpara, weight, ...)
        val <- QUAx - qua
        return(val)
      }

      root <- uniroot(fn,c(end.min,end.max))
      f[i] <- root$root
    }
    return(f)
}

