Arrows <- function(x0, y0, x1 = x0, y1 = y0, length = 0.25, angle = 30,
       code = 2, col = par("fg"), lty = par("lty"),
       lwd = par("lwd"), warnZeroLength=FALSE, ...){
##
## 1. set up
##
    dat <- data.frame(x0, y0, x1, y1, length, angle,
                      code, col, lty, lwd, ...,
                      stringsAsFactors=FALSE)
##
## 2.  plot
##
    N <- nrow(dat)
    L2 <- with(dat, (x1-x0)^2 + (y1-y0)^2)
    L2[is.na(L2)] <- 0 
    wZL <- (warnZeroLength | (L2>0))
    for(i in 1:N){
#      dNA <- (dropNA & any(is.na(dat[i,])))
#      if(dNA){
#        wZL[i] <- FALSE
#      } 
      if(wZL[i]){
        do.call(arrows, dat[i,])
      }
    }
}
