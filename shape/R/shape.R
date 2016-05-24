##==============================================================================
##
## Spheres/Cylindres/Ellipses/Arrows
##    Karline Soetaert
##
##==============================================================================

##==============================================================================
## INTERNAL FUNCTION
##==============================================================================

val2col <- function (values,zlim,col)  {

  if (! is.null(values)) {  # a matrix of radius, z-values

    if (min(values[,1])<0)
      stop ("cannot draw shape: radiusses in first column of *values* are not positive")
    values <- values [sort(values[,1],index=TRUE)$ix,]   # sort on first column (radiusses)
     
    if (is.null(zlim)) {
      zlim<-range(values[,2])
    } else {
      values[,2]<-pmin(values[,2],zlim[2])
      values[,2]<-pmax(values[,2],zlim[1])
    }
     
    x.to   <- (values[,2]-zlim[1])/(diff(zlim))
    Col    <- intpalette (inputcol=col,x.to = x.to)
    nrad   <- nrow(values)
    values[,1] <- values[,1]/values[nrad,1]
    intrad <- c(0,values[,1])
  } else {
    Col <- col
    nrad <- length(Col)
    intrad <- c(0, 1:nrad)/nrad

    ncol <- length(col)
    if (ncol < nrad) 
      Col <- intpalette(col, nrad)
  }
  return(list(Col=Col,intrad=intrad,nrad=nrad))
}
