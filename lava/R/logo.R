gfilter <- function(x,sigma=1) {
  gridfn <- function(fn,width,height,center=TRUE) {
    jx <- seq_len(height)
    jy <- seq_len(width)
    if (center) {
      jx <- jx/height-0.5
      jy <- jy/width-0.5
    }
    outer(jx, jy, FUN=fn)
  }
  width <- ncol(x); height <- nrow(x)
  oscunits <- gridfn(function(x,y) ((-1)^(x+y)),height=height,width=width,center=FALSE)
  x0 <- x*oscunits ## translate origo to center of image
  X <- fft(x0)
  d <- gridfn(function(x,y) (x^2+y^2),height=height,width=width,center=TRUE)
  Gn <- exp(-2*(base::pi*sigma)^2*d) # frequency response
  H <- X*Gn
  res <- Re(fft(H,inverse=TRUE))/(width*height)*oscunits
  return(res)
}

##' @export
lava <- function(seed,w=128,h=w,bw=4,sigma=5000,bg=20000,numcol=128,col=grDevices::heat.colors(numcol),...) {
  if (!missing(seed))
    set.seed(seed)
  x <- matrix(rnorm(w*h,bg,5000),nrow=h, ncol=w)
  x0 <- gfilter(x,sigma=4)
  y <- (x0-min(x0)+1)^1.2
  opt <- graphics::par(mai=c(0,0,0,0))
  graphics::image(y,axes=FALSE,col=col)
  graphics::par(opt)
  invisible(y)
}
