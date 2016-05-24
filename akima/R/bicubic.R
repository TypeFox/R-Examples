# R interface to Akimas regular grid spline interpolator
# (TOMS 464) improved version: TOMS 760
# http://www.netlib.org/toms/760

bicubic <- function(x,y,z,x0,y0){
  nx <- length(x)
  ny <- length(y)
  if(dim(z)[1]!=nx)
    stop("dim(z)[1] and length of x differs!")
  if(dim(z)[2]!=ny)
    stop("dim(z)[2] and length of y differs!")
  n0 <- length(x0)
  if(length(y0)!=n0)
    stop("length of y0 and x0 differs!")

  ret <- .Fortran("rgbi3p",
                  md=as.integer(1),
                  nxd=as.integer(nx),
                  nyd=as.integer(ny),
                  xd=as.double(x),
                  yd=as.double(y),
                  zd=as.double(z),                  
                  nip=as.integer(n0),
                  xi=as.double(x0),
                  yi=as.double(y0),
                  zi=double(n0),
                  ier=integer(1),
                  wk=double(3*nx*ny),
                  PACKAGE="akima")

    list(x=x0,y=y0,z=ret$zi) 
}

bicubic.grid <- function(x,y,z,xlim,ylim,dx,dy){
  nx <- length(x)
  ny <- length(y)
  if(dim(z)[1]!=nx)
    stop("dim(z)[1] and length of x differs!")
  if(dim(z)[2]!=ny)
    stop("dim(z)[2] and length of y differs!")

  xi <- seq(xlim[1],xlim[2],by=dx)
  nx <- length(xi)
  yi <- seq(ylim[1],ylim[2],by=dy)
  ny <- length(yi) 

  xmat <- matrix(rep(xi,ny),nrow=ny,ncol=nx,byrow=TRUE)
  ymat <- matrix(rep(yi,nx),nrow=ny,ncol=nx,byrow=FALSE)

  xy <- cbind(c(xmat),c(ymat))

  n0 <- nx*ny

  ret <- bicubic(x,y,z,xy[,1],xy[,2])
  
 # return cell boundaries
  list(x=xi,y=yi,z=t(matrix(ret$z,nrow=ny,ncol=nx,byrow=F)))
}
