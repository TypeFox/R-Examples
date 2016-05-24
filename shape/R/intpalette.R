##==============================================================================
## intpalette    : interpolates a palette
##==============================================================================

intpalette <- function(inputcol, numcol=length(x.to), x.from= NULL,
   x.to  = NULL) {

  if (length(inputcol)<=1)
    return(rep(inputcol,numcol))
  ifelse (is.numeric(inputcol),
    rgb.col <- inputcol, rgb.col<-t(col2rgb(c(inputcol))) )
  if (is.null(x.from))
    x.from <- seq (0,1,length.out=nrow(rgb.col))
  if (is.null(x.to  ))
    x.to   <- seq (0,1,length.out=numcol)

  if (min(x.to) < min(x.from) | max(x.to) > max(x.from))
    stop ("intpalette: cannot interpolate; ranges x.to > ranges x.from")

  outcol   <- matrix(ncol=3,nrow=numcol)
  for (i in 1:3)
    outcol[,i] <- round(approx(x.from,rgb.col[,i],xout=x.to)$y)

  outcol[outcol<0]   <- 0
  outcol[outcol>255] <- 0
  color <- rgb(outcol[,1],outcol[,2],outcol[,3],maxColorValue = 255)

  return(color)

}
