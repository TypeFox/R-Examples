# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Oct. 2012
# Version 1.0
# Licence GPL v3

.getFilterLag<-function(r,lag,n) {
  d1 <- (n -1) * lag
  d2 <- d1 + lag
  .filter(r=r,d1=d1,d2=d2)
}

if (!isGeneric("Variogram")) {
  setGeneric("Variogram", function(x,lag,cutoff,cells,size=100)
    standardGeneric("Variogram"))
}  

setMethod ('Variogram' ,signature(x='RasterLayer'),
           function (x,lag,cutoff,cells,size=100) {
             if (missing(cutoff)) cutoff<- sqrt((xmin(x)-xmax(x))^2+(ymin(x)-ymax(x))^2) / 3
             if (missing(lag)) lag <- res(x)[1]
             else if (lag < res(x)[1]) lag <- res(x)[1] 
             if (cutoff < lag) stop("cutoff should be greater than lag size")
             nlag <- ceiling(cutoff / lag)
             re <- res(x)[1]
             if (missing(cells)) {
               cells <-c(1:ncell(x))[which(!is.na(x[1:ncell(x)]))]
               if (length(cells) > size) cells <- cells[sample(1:length(cells),size)]
             }
             n <- length(cells)
             x <- as(x,"matrix")
             tbl <- matrix(nrow=n,ncol=nlag)
             for (nl in 1:nlag) {
               filter <- .getFilterLag(re,lag,nl)
               nf <- ncol(filter)
               out <- rep(NA,n)
               for (c in 1:n) {
                 xi <-  t(x)[cells[c]]
                 if (!is.na(xi)) {
                   xn <- .neighborRowCol(x,cells[c],nf)
                   xn[,1] <- xn[,1] * as(filter,"vector")
                   xn <- xn[!is.na(xn[,1]),]
                   xn <- unlist(lapply(1:nrow(xn),function(r) {x[xn[r,1],xn[r,2]]}))
                   xn <- xn[!is.na(xn)]
                   out[c] <- mean((xi - xn)^2)/2
                 } else out[c] <- NA
               }
               tbl[,nl] <- out
             }
             v <- new("RasterVariogram")
             v@lag <- lag
             v@nlags <- nlag
             v@variogramCloud <- tbl
             v@variogram <- data.frame(distance=seq(lag,lag*nlag,lag) - (lag/2),gamma=apply(tbl,2,mean,na.rm=TRUE))
             v
            }
           )
