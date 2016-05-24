trMean <- function(b,blen) {
##  mytr <- function(x) x^2; dmytr <- function(x) 2*x
  mytr <- dmytr <- exp
  if (blen==0) return(b) 
  k <- length(b)
  Bidx <- seq_len(blen)+(k-blen)
  b[Bidx[1]] <- mytr(b[Bidx[1]])
  D <- diag(nrow=k)
  D[Bidx[1]:k,Bidx[1]] <- b[Bidx[1]]
  for (i in Bidx[-1]) {
    D[i:k,i] <- dmytr(b[i])
    b[i] <- b[i-1]+mytr(b[i])
  }
  attributes(b)$D <- D
  attributes(b)$idx <- Bidx
  return(b)
}


##' @export
plot.bptwin <- function(x,n=50,rg=range(x$B[,1]),xlab="Time",ylab="Concordance",...) {
  if (x$Blen>0) {
    ##    rg <- range(x$B[,1])
    t <- seq(rg[1],rg[2],length.out=n)
    B0 <- bs(t,degree=x$Blen)
    b0. <- coef(x)[x$midx0]
    b1. <- coef(x)[x$midx1]
    b0 <- trMean(b0.,x$Blen)
    b1 <- trMean(b1.,x$Blen)
    b00 <- tail(b0,x$Blen)
    b11 <- tail(b1,x$Blen)
    pr0 <- sapply(as.numeric(B0%*%b00+b0[1]), function(z)
                  pbvn(upper=rep(z,2),sigma=x$Sigma0))
    pr1 <- sapply(as.numeric(B0%*%b11+b1[1]), function(z)
                  pbvn(upper=rep(z,2),sigma=x$Sigma1))
    plot(pr0~t,type="l", xlab=xlab, ylab=ylab,...)
    lines(pr1~t,type="l",lty=2)
  }
  return(invisible(x))
}
