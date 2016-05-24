## =============================================================================
## Maps a matrix u from (x,y) to (x2,y2) by 2-D linear interpolation
## =============================================================================
mapxy <- function(x, y, x2, y2, u) {
  Nx <- length(x)
  Ny <- length(y)

  if (min(x2) < min(x)) stop("'x' should embrace 'x2'")
  if (min(y2) < min(y)) stop("'y' should embrace 'y2'")

  dx   <- c(diff(x),1)  # 1= for last value
  dy <- c(diff(y),1)

  Du <- dim(u)
  if (Du[1] != Nx)
    stop("u and x not compatible")
  if (Du[2] != Ny)
    stop("u and y not compatible")

  Transf <- function (x2,y2,u) {

  # find embracing values : first interval
    ix <- findInterval(x2, x )
    iy <- findInterval(y2, y)

#    ix <- pmax(ix, 1)
#    iy <- pmax(iy, 1)

  # next interval
    ixp1 <- pmin(ix+1,Nx)
    iyp1 <- pmin(iy+1,Ny)

  # interpolation factor
    xfac <- (x2-x[ix])/dx[ix]
    yfac <- (y2-y[iy])/dy[iy]

  # interpolate
    (1-yfac)*((1-xfac)*u[cbind(ix,iy)]+xfac*u[cbind(ixp1,iy)]) +
    yfac*((1-xfac)*u[cbind(ix,iyp1)]+xfac*u[cbind(ixp1,iyp1)])

  } # end Transf

outer(x2, y2, FUN = Transf, u = u)
}

#mapxy(1:87, 1:61, seq(1,87,0.25), seq(1,61,0.25),volcano)->VOLCANO
#IM <- matrix(nrow=3,ncol=3,c(1,1,1,1,2,1,1,1,1))
#mapxy(1:3,1:3,seq(0,3,0.1),seq(0,3,0.1),IM)->im2

## =============================================================================
## From polar to cartesian coordinates
## =============================================================================

polar2cart <- function(out, r, theta, x = NULL, y = NULL)
{
  if (is.null(theta)) {   # one-D mapping to 2-D
    if (!is.vector(out))
      stop ("'out' should be a vector if 'theta' = NULL")
    afun <- approxfun(r, out, yleft = out[1]) 
    r2 <- matrix(nrow = length(x), ncol = length(y), 
            afun(outer (FUN = function(x, y) sqrt(x^2 + y^2), X = x, Y = y)))
    return(r2)
  
  }
  
  classout <- class(out)[1]
  if (!classout %in% c("steady2D", "deSolve", "matrix"))
    stop ("Class of 'out' should be one of steady2D, deSolve or matrix")

  ATTR <- attributes(out)
  if (length(ATTR$dimens)!= 2)
    stop ("input should be 2-dimensional")

  Nr   <- length(r)  -1
  Np   <- length(theta)-1
  maxr   <- max(r)
  maxtheta <- max(theta)
  minr   <- min(r)
  mintheta <- min(theta)

  # add leftmost boundary
  r      <- c(r[1],   0.5*(r[-1]+r[-(Nr+1)]), r[Nr+1])
  theta  <- c(theta[1], 0.5*(theta[-1]+theta[-(Np+1)]), theta[Np+1])
                                                # check dimensions...
  dr     <- c(diff(r),1)       # 1= for last value: no interpolation
  dtheta <- c(diff(theta),1)

  if (is.null (x))
    if (! is.null(y)) x<-y
  if (is.null (y)) {
    mr <- max(abs(r))
    x <- y <- seq(-mr, mr, len=Nr)
  }


  Transf <- function(x, y,  u) {

    Du <- dim(u)
    if (Du[1] != Nr)
      stop("u and r not compatible")
    if (Du[2] != Np)
      stop("u and theta not compatible")
  # augment u with boundary values
    u    <- rbind(c(u[1,1],u[1,], u[1,Np]) , 
                  cbind(u[,1],u, u[,Np]),
                  c(u[Nr,1],u[Nr,], u[Nr,Np]) )

    R   <- sqrt(x^2+y^2)
    Theta <- atan(y/x)
    Theta[x<0]           <- Theta[x < 0] + pi
    Theta[x>0 & y < 0]   <- Theta[x > 0 & y < 0] + 2*pi
    Theta[x==0 & y > 0]  <- pi/2
    Theta[x==0 & y < 0]  <- 3*pi/2
    Theta[x==0 & y == 0] <- 0

  # find embracing values : interval location
    iR <- findInterval(R, r )
    iP <- findInterval(Theta,theta)
    iR [iR == 0]   <- NA
    iR [iR > Nr+1] <- NA
    iP [iP == 0]   <- NA
    iP [iP > Np+1] <- NA

  # next interval
    iPp1 <- pmin(iP+1,Np+1)
    iRp1 <- pmin(iR+1,Nr+1)

  # interpolation factor
    iRfac <- (R-r[iR])/dr[iR]
    iPfac <- (Theta-theta[iP])/dtheta[iP]

  # interpolate inbetween 4 values
    (1-iPfac)*((1-iRfac)*u[cbind(iR,iP)]+iRfac*u[cbind(iRp1,iP)]) +
    iPfac*((1-iRfac)*u[cbind(iR,iPp1)]+iRfac*u[cbind(iRp1,iPp1)])

  } # end Transf


  Num   <- prod(ATTR$dimens)
  nspec <- ATTR$nspec
  if (classout[1] == "steady2D") {
    MAP   <- list()
    MAP$y <- NULL
    for (i in 1:nspec)  {
      ii <- (1+(i-1)*Num):(i*Num)
      U <- matrix(out$y[ii], nrow = ATTR$dimens[1], ncol= ATTR$dimens[2])
      MAP$y <- c(MAP$y, as.vector(outer(x, y, FUN = Transf, u = U)))
    }
  attributes(MAP) <- ATTR
  } else {               # Dynamic output
    MAP <- NULL
    nr = nrow(out)
    for ( j in 1: nr) {
      mm <- out[j,1]
      for (i in 1:nspec)  {
        ii <- ((1+(i-1)*Num):(i*Num))+1
        U <- matrix(out[j,ii], nrow = ATTR$dimens[1], ncol= ATTR$dimens[2])
        mm <- c(mm, as.vector(outer(x, y, FUN = Transf, u = U)))
      }
      MAP <- rbind(MAP, mm)
    }
  attributes(MAP)[c("istate","rstate","class","type","nspec","dimens")] <- 
    ATTR[c("istate","rstate","class","type","nspec","dimens")]
  }
  attributes(MAP)$dimens     <- c(length(x), length(y))
  attributes(MAP)$xgrid      <- x
  attributes(MAP)$ygrid      <- y
 
  MAP
}
