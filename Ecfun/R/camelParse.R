camelParse <- function(x, except=c('De', 'Mc', 'Mac')){
##
## 1.  strsplit
##
  x. <- strsplit(x, "")
##
## 2.  Find lower followed by upper
##
  nx <- length(x)
  out <- vector('list', length=nx)
  names(out) <- names(x)
  for(ix in 1:nx){
    xi <- x.[[ix]]
    lower <- (xi %in% letters)
    upper <- (xi %in% LETTERS)
    ni <- length(xi)
    camel <- which(lower[-ni] & upper[-1])
    begin <- c(1, camel+1)
    end <- c(camel, ni)
    X <- substring(x[ix], begin, end)
    for(ex in except){
        ei <- regexpr(ex, X)
        ej <- (ei+2-nchar(X))
        ej[ei<0] <- -1
        ek <- which(ej>0)
        for(ik in rev(ek)){
            X[ik] <- paste(X[ik], X[ik+1], sep='')
            X <- X[-(ik+1)]
        }
    }
    out[[ix]] <- X
  }
  out
}

