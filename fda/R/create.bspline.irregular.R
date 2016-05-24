create.bspline.irregular <- function (argvals,
      nbasis=max(norder, round(sqrt(length(argvals)))),
      norder=4,
      breaks=quantile(argvals, seq(0, 1, length.out=nbasis-norder+2)),
      dropind=NULL, quadvals=NULL, values=NULL,
      basisvalues=NULL, names="bspl", plot.=FALSE, ...){
##
## 1.  create.bspline.basis
##
#  k <- length(breaks)
#  bks <- breaks[-c(1, k)]
  bsp <- create.bspline.basis(range(argvals), nbasis=nbasis,
                              norder=norder, breaks=breaks,
                              dropind=dropind, quadvals=quadvals,
                              values=values, basisvalues=basisvalues,
                              names=names)
##
## 2.  plot
##
  if(plot.){
    dots <- list(...)
    if(!('x' %in% names(dots)))
        dots$x <- argvals
    if(!('type' %in% names(dots)))
        dots$type <- 'l'
    if(!('ylab' %in% names(dots)))
        dots$ylab <- 'argvals'
#  plot(argvals, type='l', ...)
    do.call(plot, dots)
#
    N <- length(argvals)
    Knots <- knots(bsp)
    Knots. <- c(bsp$rangeval[1], Knots, bsp$rangeval[2])
    qtiles <- sapply(Knots., function(x)mean(argvals<=x))

    points(N*qtiles, Knots., ...)
  }
##
## 3.  Done
##
  bsp
}
