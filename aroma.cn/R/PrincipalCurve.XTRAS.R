setMethodS3("getBacktransforms", "PrincipalCurve", function(fit, dimensions=NULL, targetDimension=1L, range=NULL, length.out=100L, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dimensions':
  s <- fit$s;
  ndim <- ncol(s);
  if (is.null(dimensions)) {
    dimensions <- seq_len(ndim);
  }
  dimensions <- Arguments$getIndices(dimensions, max=ndim);

  # Argument 'range':
  if (is.null(range)) {
    range <- range(s, na.rm=TRUE);
  }
  range <- Arguments$getDoubles(range, length=c(2L,2L), disallow=c("NA", "NaN"));


  y <- seq(from=range[1L], to=range[2L], length.out=length.out);

  naValue <- as.double(NA);
  dim <- c(length(y), 2L, length(dimensions));
  XY <- array(naValue, dim=dim);

  for (kk in seq_along(dimensions)) {
    dim <- dimensions[kk];
    yN <- .backtransformPrincipalCurve(y, fit=fit, dimensions=dim,
                                      targetDimension=targetDimension);
    yN <- yN[,1L,drop=TRUE];
    xy <- cbind(y, yN);

    XY[,,kk] <- xy;
  } # for (kk ...)

  XY;
}) # getBacktransforms()


setMethodS3("plotBacktransforms", "PrincipalCurve", function(fit, ..., col="red", lwd=2, lty=1L, xlim=NULL, ylim=xlim, xlab="y", ylab="y*") {
  XY <- getBacktransforms(fit, ...);

  if (is.null(xlim)) {
    xlim <- range(XY, na.rm=FALSE);
  }
  if (is.null(ylim)) {
    ylim <- xlim;
  }

  ndim <- dim(XY)[3L];
  subplots(ndim);
  for (kk in seq_len(ndim)) {
    xy <- XY[,,kk];
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
    abline(a=0, b=1, lty=3L);
    lines(xy, col=col, lwd=lwd, lty=lty);
  }

  invisible(XY);
}) # plotBacktransforms()


############################################################################
# HISTORY:
# 2013-10-08
# o Added arguments 'col', 'lwd' and 'lty' to plotBacktransforms() and
#   argument 'xlim' now defaults to the maximum range of the data.
# 2010-01-14
# o Added getBacktransforms() and plotBacktransforms().
# o Created.
############################################################################
