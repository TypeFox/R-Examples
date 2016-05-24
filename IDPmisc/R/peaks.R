## peaks.R

peaks <- function(x, y = NULL, minPH, minPW, thr, stepF = 0.49){
  ## Author: Rene Locher
  ## Version: 2006-07-03

  xy <- getXY(x, y)
  xx <- xy$x
  yy <- xy$y
  
  min.yy <- min(yy, na.rm=TRUE)
  max.yy <- max(yy, na.rm=TRUE)

  ## make sure that curve is ascending at the beginning!
  yy <- c(min.yy,yy)
  xx <- c(xx[1],xx)

  if (missing(thr)) thr <- min(yy, na.rm=TRUE)
  if (missing(minPH)) minPH <- diff(range(yy), na.rm=TRUE)/10
  if (missing(minPW)) minPW <- 0

  ## The smaller stepF, the longer the calculation,
  ## the preciser the values of the peak widths
  if (stepF>=0.5) stop("'stepF' must be smaller than 0.5")

  peak.x <- peak.y <- peak.w <- NULL
  lev <- thr - stepF*minPH
  len.yy <- length(yy)
  
  repeat { ## find maximum above level lev
    lev <- lev + stepF*minPH
    if (lev >= max.yy) break
    hi <- ok(yy>lev)
    start <- which(diff(c(F,hi))>0)
    end <- which(diff(c(hi,F))<0)
    len <- length(start)
    if(len==0) next

    for (ii in 1:len){
      ## find maximum in each continuous part above lev
      x <- xx[start[ii]:end[ii]]
      y <- yy[start[ii]:end[ii]]
      i <- which.max(y)

      ## is local maximum already defined?
      if (is.element(x[i],peak.x)) next

      miny <- min(yy[max(1,start[ii]-1):min(end[ii]+1,len.yy)],
                  na.rm=TRUE)

      ## calculate height of Peak
      PH <- y[i]-miny

      ## calculate width at half peak height
      PW <- paste(as.numeric(y > (miny+PH/2)),collapse="")
      PW <- max(attr(gregexpr("1+", PW)[[1]],"match.length"))

      ## PW is calculated correctly only when the positions
      ## within the observed window is constant
      ## in rare cases PW might stay -1
      PW <-  if(PW>=0) abs(x[PW]-x[1])
      
      if (PH >= minPH && PW >= minPW) {
        peak.x <- c(peak.x,x[i])
        peak.y <- c(peak.y,y[i])
        peak.w <- c(peak.w,PW)
      }
    } ## for
  } ## repeat

  res <- data.frame(x=peak.x, y=peak.y, w=peak.w)
  res <- res[order(res$x),]
  return(res)
} ## peaks

