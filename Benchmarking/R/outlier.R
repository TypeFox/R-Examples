# $Id: outlier.R 140 2015-05-15 21:48:02Z B002961 $


outlier.ap <- function(X, Y, NDEL = 3, NLEN = 25, TRANSPOSE = FALSE)  
{
   xy <- cbind(X,Y)
   S <- det(t(xy)%*%xy)
   R <- NDEL
   K <- dim(X)[1]   # number of firms/observations
   Kset <- seq.int(K)
   n <- length(Kset)
   last = min(NLEN, K)   # gem 25 mindste vaerdier af RX
   ratio <- array(Inf, c(last,R) )
   imat <- matrix(NA, nrow=R, ncol=R)
   rmin <- array(Inf, R)
   for ( r in 1:R )  {
      # print( paste("Remove ",r," observations",sep=""), quote=FALSE )
      # remove r observationer
      RX <- rep(Inf,last)
      rrlast <- Inf
   
      if (n < r) 
          stop("n < r")
      e <- 0
      h <- r
      a <- 1L:r
      del <- Kset[a]
   
      Dxy = xy[-del,]
      rr <- det( t(Dxy) %*% Dxy )/S
      RX[1] <- rr
      rrlast <- RX[last]
   
      # count <- as.integer(round(choose(n, r)))
      # print( paste("Number of combinations", count), quote=FALSE )
      # i <- 2L
      nmmp1 <- n - r + 1L
      while (a[1L] != nmmp1) {
          if (e < n - h) {
              h <- 1L
              e <- a[r]
              j <- 1L
          }
          else {
              e <- a[r - h]
              h <- h + 1L
              j <- 1L:h
          }
          a[r - h + j] <- e + j
          del <- Kset[a]
   
          Dxy = xy[-del,]
          rr <- det( t(Dxy) %*% Dxy )/S
          if ( rr < rrlast ) {  # gem de |last| mindste
             if ( rr < min(RX) ) {
                imat[r,1:r] <- del
             }
             RX[last] <- rr
             RX <- sort(RX)
             rrlast <- RX[last]
          }
   
          # i <- i + 1L
      }
      # print(i)
   
      rmin[r] <- min(RX)
      Rratio <- log(RX/min(RX))
      ratio[,r] <- Rratio
      # if ( r==1 ) {
      #    # plot(rep(r,last),Rratio ,ylim=c(0,.6),xaxt="n")
      #    plot(1:R,rmin,type="n", ylim=c(0,.6),xaxt="n", ylab="Log ratio",xlab="r")
      #     axis(1, at=1:R, labels=c(1:R))
      # } 
      # points(rep(r,last),Rratio)
   } 
   # lines(1:R, ratio[2,],lty="dashed")
   return( list(ratio=ratio, imat=imat, r0=rmin) )
   
}  # outlier.ap



outlier.ap.plot <- function(ratio, NLEN = 25, xlab="r", ylab="Log ratio", 
                ..., ylim)  
{
   nlen <- min(NLEN, dim(ratio)[1])
   R <- dim(ratio)[2]
   if ( missing(ylim) ) ylim <- range(ratio[1:nlen,])
   ry <- matrix(1:R,nrow=nlen, ncol=R, byrow=TRUE)

   plot(1:R, rep(0,R), ylim=ylim, xaxt="n", ylab=ylab,xlab=xlab, ...)
   axis(1, at=1:R, labels=c(1:R))
   points(ry, ratio[1:nlen,])
   lines(1:R, ratio[2,],lty="dashed")
}  # outlier.ap.plot

