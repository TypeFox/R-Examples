###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitMultiDimensionalCone
#
# @title "Fits an affine transformation to multi-dimensional data"
#
# \description{
#  @get "title" using robust estimators.
# }
#
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric NxK @matrix with one column for each dimension and
#      where N is the number of data points.}
#   \item{alpha}{A @numeric @vector of decreasing values in (0,1).
#      This parameter "determines how far we are willing to press the
#      boundary of the [genotype cone]".  Lowering \code{alpha} expand
#      the cone.  When \code{alpha} goes to zero, all data points will
#      be on or inside the cone.}
#   \item{q,Q}{Percentiles in [0,100] for which data points that are
#      below (above) will be assigned zero weight in the fitting of
#      the parameters.}
#   \item{...}{Additional arguments passed to the \code{cfit()} function
#      of the \pkg{sfit} package.}
#   \item{flavor}{A @character string specifying what model/algorithm
#      should be used to fit the genotype cone.}
# }
#
# \value{
#   Returns the parameter estimates as a named @list with elements:
#    \item{M}{An estimate of the three vertices defining the genotype
#      triangle.  These three vertices are describes as an 2x3 @matrix
#      with column \code{origin}, \code{AA}, and \code{BB}.}
#    \item{Minv}{The inverse of \code{M}.}
#    \item{origin}{The estimate of the shift.}
#    \item{W}{The estimate of shear/rotation matrix with columns
#             \code{AA} and \code{BB}.}
#    \item{Winv}{The inverse of \code{W}.}
#    \item{params}{The parameters used for the fit, i.e.
#       \code{alpha}, \code{q}, \code{Q}, and  those passed in \code{...}.}
#    \item{dimData}{The dimension of the input data.}
# }
#
# \examples{
# if (require("sfit")) {
#  @include "../incl/fitMultiDimensionalCone.matrix.Rex"
# }
# }
#
# @author
#
# \seealso{
#  To backtransform data fitted using this method,
#  see @seemethod "backtransformMultiDimensionalCone".
#  Internally, the \code{cfit()} function the \pkg{sfit} package is used.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("fitMultiDimensionalCone", "matrix", function(y, alpha=c(0.10, 0.075, 0.05, 0.03, 0.01), q=2, Q=98, ..., flavor=c("sfit", "expectile")) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'alpha':
  if (flavor == "sfit") {
    requireNamespace("sfit") || throw("Package not loaded: sfit");
  } else if (flavor == "expectile") {
    # To please/fool R CMD check (in the case expectile is not installed)
    fitCone <- NULL; rm(list="fitCone");
    requireNamespace("expectile") || throw("Package not loaded: expectile");
    # Only final 'alpha' is needed by expectile::fitCone().
    alpha <- rev(alpha)[1];
  }


  # Fit simplex of (y_1, y_2, ..., y_K)
  if (flavor == "sfit") {
    fit <- sfit::cfit(y, alpha=alpha, q=q, Q=Q, ...);
  } else if (flavor == "expectile") {
    fit <- expectile::fitCone(y, alpha=alpha, ...);
    print(fit);
    throw("The rest is not implemented yet: ", flavor);
  }

  M <- fit$M;
  colnames(M) <- sprintf("dim%d", seq_len(ncol(M)));
  clazz <- class(M);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Re-arrange vertices in the order origin, 1st, 2nd, ..., Kth axis.
  # For bi-allele data, this means: (origin, A, B)
  # For resequence data, this means: (origin, A, C, G, T)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Identify the origin (0th vertex): It is the closest one to (0,0,...,0).
  dist <- apply(M, MARGIN=1, FUN=function(u) sum(u^2));
  idxOrigin <- which.min(dist);
  origin <- M[idxOrigin,];

  # Identify each of the dimensions: The vertex which has its greatest value
  # in the 1st position is the 1st dimension, and so on.
  M <- M[-idxOrigin,,drop=FALSE];
  dims <- apply(M, MARGIN=1, FUN=function(u) which.max(u));

  # Reorder the vertices accordingly
  o <- order(dims);
  M <- M[o,,drop=FALSE];
  rownames(M) <- sprintf("vertex%d", seq_len(nrow(M)));

  # Append the ordered vertices to the origin
  M <- rbind(origin, M);

  class(M) <- clazz;
  fit$M <- M;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate the backtransform matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift all vertices such that the 0th vertex is at (0,0,...,0).
  W <- M[-1,,drop=FALSE];
  W <- t(W) - origin;

  # Rescale such that the 1st dimension has value one in the first position.
  W <- W / W[1,1];

  # Find the inverse
  Winv <- solve(W);

  # Calculate the inverse
  Minv <- t(Winv %*% (t(M)-origin));
  class(Minv) <- clazz;
  fit$Minv <- Minv;

  W <- t(W);
  Winv <- t(Winv);

  fit$origin <- origin;
  fit$W <- W;
  fit$Winv <- Winv;

  fit$params <- list(
    alpha=alpha,
    q=q,
    Q=Q,
    ...
  );

  fit$dimData <- dim(y);

  fit;
}, private=TRUE) # fitMultiDimensionalCone()



############################################################################
# HISTORY:
# 2008-08-01
# o Generalized code to K dimensions.
# 2008-02-14
# o Added a self-contained example for fitGenotypeCone().
# 2007-09-08
# o Added 'dimData' to the return structure.
# 2007-06-12
# o Commented the code for re-arranging fit$X (only if retX=TRUE).
#   Code not really needed since backtransformGenotypeCone() is used.
# 2007-06-11
# o Calls sfit::cfit() explicitly.
# 2007-06-04
# o Added Rdoc comments.
# 2006-05-08
# o Created.
############################################################################
