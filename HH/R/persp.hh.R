if.R(r={persp.setup <- function(...)
     .Defunct("persp.setup", package="HH", "'persp.setup' is not used in R.")},
     s={}   ## persp.setup is used in S-Plus, and is not needed or used in R
   )

persp.plane <- function(...)
  .Defunct("perspPlane", package="HH")

perspPlane <- function(x, y, z, persp.out, ...) {

  if.R(r=perspp <- trans3d, s={})

  lx <- length(x)
  ly <- length(y)

  pxx <- matrix(x,lx,2)
  pxy <- matrix(y[c(1,ly)],lx,2, byrow=TRUE)
  pxz <- z[1:lx, c(1,ly)]

  px1 <- perspp(pxx[,1],
                pxy[,1],
                pxz[,1],
                persp.out)
  px2 <- perspp(pxx[,2],
                pxy[,2],
                pxz[,2],
                persp.out)

  pyx <- matrix(x[c(1,lx)],ly,2, byrow=TRUE)
  pyy <- matrix(y,ly,2)
  pyz <- t(z[c(1,lx), 1:ly])

  py1 <- perspp(pyx[,1],
                pyy[,1],
                pyz[,1],
                persp.out)

  py2 <- perspp(pyx[,2],
                pyy[,2],
                pyz[,2],
                persp.out)

  segments(px1$x, px1$y, px2$x, px2$y, ...)
##  segments(px$x[,1], px$y[,1], px$x[,2], px$y[,2], ...)
  segments(py1$x, py1$y, py2$x, py2$y, ...)
##  segments(py$x[,1], py$y[,1], py$x[,2], py$y[,2], ...)
}

##floor
persp.floor <- function(...)
  .Defunct("perspFloor", package="HH")

perspFloor <- function(x, y, z, persp.out, ...) {
  perspPlane(x,
              y,
              matrix(min(z), length(x), length(y)),
              persp.out,
              ...)
}

##back wall along x
persp.back.wall.x <- function(...)
  .Defunct("perspBack.wall.x", package="HH")

perspBack.wall.x <- function(x, y, z, persp.out, ...) {
  perspPlane(x,
              rep(max(y), length(pretty(z))),
              matrix(pretty(z), length(x), length(pretty(z)), byrow=TRUE),
              persp.out,
              ...)
}

##back wall along y
persp.back.wall.y <- function(...)
  .Defunct("perspBack.wall.y", package="HH")

perspBack.wall.y <- function(x, y, z, persp.out, ...) {
  perspPlane(rep(max(x), length(pretty(z))),
              y,
              matrix(pretty(z), length(pretty(z)), length(y)),
              persp.out,
              ...)
}


## ##sample
## x <- 1:5
## y <- 1:10
## z <- matrix(1:50, 5, 10)
## persp.out <- persp(x, y, z)
## perspFloor(x, y, z, persp.out, col=2, lty=2)
## perspBack.wall.x(x, y, z, persp.out, col=3, lty=3)
## perspBack.wall.y(x, y, z, persp.out, col=4, lty=4)

## debugging tools
## perspSetup(col=c(0,0,0)) # set colors to background color: top, hidden, bottom
## perspSetup(restore=TRUE)    # restore default values
## 
## trace(perspPlane, exit=browser)
## 0
## untrace()
## 
## par()$usr
## 
## x
## y
## z
## lx
## ly
## px
## py
