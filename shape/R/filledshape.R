##==============================================================================
## filledshape     : draws and colors a shape   color depend on radius
##==============================================================================

filledshape <- function(xyouter, xyinner=colMeans(xyouter),
  col=femmecol(100), values=NULL, zlim=NULL, lcol=NA, lwd=2, ...) {


  expand<-function(mat,npoints) {
    mat   <- matrix(ncol=2,mat)
    nfrom <- nrow(mat)
    nin   <- round(npoints/nfrom)
    nin   <- rep(nin,nfrom-1)
    nin   <- c(nin,npoints-sum(nin))
    mat   <-rbind(mat,mat[1,])
    out   <- NULL
    for (i in 1:nfrom) {
      x  <- approx(x=mat[c(i,i+1),1],n=nin[i])$y
      y  <- approx(x=mat[c(i,i+1),2],n=nin[i])$y
      out <- rbind(out,cbind(x,y) )
    }
    out
  }

  vv     <- val2col(values,zlim,col)
  intrad <- vv$intrad
  Col    <- vv$Col
  nrad   <- vv$nrad

  ## check inner and outer points
  npoint <- nrow(xyouter)
  nmid   <- length(xyinner)/2
  middle <- xyinner
  extern <- xyouter

  if ( nrad == 1 & nmid == 1) {
    polygon(xyouter[,1],xyouter[,2],col=Col,border=Col,...)
  } else if (nrad==1) {
    if (nmid < npoint)
      for (i in (nmid+1):npoint)
        middle <- rbind(middle ,xyinner)
    polygon(c(xyouter[,1], rev(middle[,1])), c(xyouter[,2], rev(middle[,2])),
            col=Col, border=Col, ...)

  } else {
    if (nmid < npoint)
      middle <- expand(middle,npoint)
    if (nmid > npoint)
      extern <- expand(extern,nmid  )
    ## start coloration
    inner  <- middle

    for (i in 1:nrad) {
      relrad    <- intrad[i+1]
      outer     <- inner
      inner[,1] <- middle [,1] + relrad * (extern[,1]-middle[,1])
      inner[,2] <- middle [,2] + relrad * (extern[,2]-middle[,2])
      polygon(c(outer[,1], rev(inner[,1])), c(outer[,2], rev(inner[,2])),
              col=Col[i], border=Col[i], ...)
    }
  }

  if (! is.na(lcol)) {
    lines(xyouter,lwd=lwd,col=lcol)
    if (length(xyinner)>2)
      lines(xyinner,lwd=lwd,col=lcol)
  }

} # of filledshape

