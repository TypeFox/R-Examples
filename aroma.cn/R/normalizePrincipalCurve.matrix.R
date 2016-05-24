#########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizePrincipalCurve
#
# \encoding{latin1}
#
# @title "Normalizes data in K dimensions using principal curves"
#
# \description{
#   @get "title" such that afterward the data cluster (approximately
#   linearly) along the diagonal (in K dimensions).
# }
#
# @synopsis
#
# \arguments{
#  \item{x}{An NxK @matrix where the columns represent the
#      (K >= 2) dimensions.}
#  \item{...}{Additional arguments passed to
#      @see "aroma.light::fitPrincipalCurve" used for fitting the model.}
#  \item{center}{If @TRUE, normalized data is centered such that the median
#      signal in each dimension is at zero.}
#  \item{returnFit}{If @TRUE, the fitted principal curve parameters are
#      returned as an attribute.}
# }
#
# \value{
#   Returns an NxK @matrix.
# }
#
# \references{
#   [1] Hastie, T. and Stuetzle, W, \emph{Principal Curves}, JASA, 1989.
# }
#
# \seealso{
#   @see "aroma.light::fitPrincipalCurve" and
#   @see "aroma.light::backtransformPrincipalCurve".
# }
#
# @author "HB"
#
# @keyword internal
#*/#########################################################################
setMethodS3("normalizePrincipalCurve", "matrix", function(x, ..., center=TRUE, returnFit=FALSE) {
  # Fit principal curve
  fit <- .fitPrincipalCurve(x, ...);

  # Flip direction of 'lambda'?
  rho <- cor(fit$lambda, x[,1L], use="complete.obs");
  flip <- (rho < 0);

  # Sanity check
  stopifnot(identical(dim(fit$s), dim(x)));
  dx <- (fit$s - x);

  # Sanity check
  stopifnot(identical(dim(dx), dim(x)));
  stopifnot(identical(nrow(dx), length(fit$lambda)));
  xN <- fit$lambda + dx;
  stopifnot(identical(dim(xN), dim(x)));

  if (flip) {
    xN <- -xN;
  }

  if (center) {
    # Same center for each column
    for (cc in seq_len(ncol(x))) {
      mu <- median(x[,cc], na.rm=TRUE);
      muN <- median(xN[,cc], na.rm=TRUE);
      xN[,cc] <- xN[,cc] - (muN-mu);
    }
  }

  # Return fit?
  if (returnFit)
    attr(xN, "fit") <- fit;

  xN;
}) # normalizePrincipalCurve()



###########################################################################
# HISTORY:
# 2012-04-16
# o Now normalizePrincipalCurve() explicitly require aroma.light.
# 2012-03-30
# o Added Rdoc comments.
# 2008-10-08
# o Removed implementation for data.frame:s.
# 2008-05-27
# o Added normalizePrincipalCurve().
# o Created.  Will probably end up in aroma.light.
###########################################################################
