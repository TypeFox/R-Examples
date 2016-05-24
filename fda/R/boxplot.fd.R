boxplot.fdSmooth <- function(x, z=NULL, ...){
  boxplot(x$fd, z, ...)
}

boxplot.fdPar <- function(x, z=NULL, ...){
  boxplot(x$fd, z, ...)
}

boxplot.fd <- function(x, z=NULL, ...){
  if(is.numeric(x)){
      fbplot(x, z, ...)
  } else {
      if(is.null(z)){
          rng <- getbasisrange(x$basis)
          z <- seq(rng[1], rng[2], length=101)
      }
      x. <- eval.fd(z, x)
#      x. <- predict(x, z)
      dots <- list(...)
      if(!('xlim' %in% names(dots)))xlim <- range(z)
      if(!('ylim' %in% names(dots)))
          ylim <- c(min(x.)-.5*diff(range(x.)),max(x.)+.5*diff(range(x.)))
      dots$fit <- x.
      dots$x <- z
      dots$xlim <- xlim
      dots$ylim <- ylim
      do.call(fbplot, dots)
#      fbplot(x., z, method=method, depth=depth, plot=plot, prob=prob,
#         color=color, outliercol=outliercol, barcol=barcol,
#         fullout=fullout, factor=factor, ...)
  }
}
