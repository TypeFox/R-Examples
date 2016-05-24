###########################################################################/**
# @set "class=matrix"
# @RdocMethod matrixBlockPolish
#
# @title "Applies a polishing function to blocks of rows and columns repeatedly"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{z}{A @numeric KxN @matrix.}
#   \item{x}{A optional KxNx2 @array (or KxN @matrix).}
#   \item{blockSizes}{A positive @integer @vector of length two.}
#   \item{FUN}{A @function taking @numeric arguments \code{z} and
#      \code{x} and returns a @numeric object with either a scalar
#      or the same number of elements as in \code{z}.}
#   \item{...}{Additional arguments passed to the \code{FUN} @function.}
#   \item{tol}{A positive threshold specifying when the algorithm has
#      converged.}
#   \item{maxIter}{The maximum number of iterations.}
#   \item{returnEffects}{If @TRUE, the row and column effects are returned,
#      otherwise not.}
# }
#
# \value{
#   Returns a named @list.
# }
#
# @examples "../incl/matrixBlockPolish.matrix.Rex"
#
# @author
#
# \seealso{
#  @see "stats::medpolish".
#  \code{\link[aroma.light:medianPolish.matrix]{medianPolish}()}.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("matrixBlockPolish", "matrix", function(z, x=NULL, blockSizes=c(1,1), FUN, ..., tol=0.01, maxIter=10, returnEffects=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dim <- dim(z);

  # Argument 'x':
  if (is.null(x)) {
    x <- array(as.integer(NA), dim=c(dim, 2));
    for (dd in 1:2) {
      t <- matrix(seq_len(dim[dd]), nrow=dim[1], ncol=dim[2],
                                                             byrow=(dd == 2));
      x[,,dd] <- t;
    }
  } else if (is.matrix(x)) {
    if (any(dim(x) != dim)) {
      throw("Argument 'x' has a different dimension that 'z': ",
        paste(dim(x), collapse="x"), " != ", paste(dim, collapse="x"));
    }
    x <- array(x, dim=c(dim, 1));
  } else if (is.array(x)) {
    if (any(dim(x)[1:2] != dim)) {
      throw("The dimension of argument 'x' is incompatible with 'z': ",
        paste(dim(x), collapse="x"), " != ", paste(dim, collapse="x"));
    }
  } else {
    throw("Argument 'x' must be a matrix, array, or NULL: ", class(x)[1]);
  }

  # Argument 'blockSizes':
  blockSizes <- rep(as.integer(blockSizes), length.out=2);

  # Argument 'FUN':
  if (!is.function(FUN)) {
    throw("Argument 'FUN' is not a function: ", class(FUN)[1]);
  }



  ranges <- vector("list", length=2);
  for (dd in 1:2) {
    idxs1 <- seq(from=1, to=dim[dd], by=blockSizes[dd]);
    idxs2 <- c(idxs1[-1]-1, dim[dd]);
    ranges[[dd]] <- cbind(from=idxs1, to=idxs2);
  }
  # Not needed anymore
  idxs1 <- idxs2 <- NULL;

  if (returnEffects) {
    blockSizes <- sapply(ranges, FUN=function(r) r[,2]-r[,1]+1);
    maxBlockSizes <- sapply(blockSizes, FUN=max);
    effects <- vector("list", length=2);
    names(effects) <- c("rows", "columns");
    for (dd in 1:2) {
      nbrOfBlocks <- nrow(ranges[[dd]]);
      effects[[dd]] <- matrix(as.double(NA), nrow=nbrOfBlocks, ncol=maxBlockSizes[dd]*dim[-dd]);
    }
  }


  oldSum <- 0;
  for (ii in seq_len(maxIter)) {
    for (dd in 1:2) {
      range <- ranges[[dd]]
      froms <- range[,1];
      tos <- range[,2];
      nbrOfBlocks <- length(froms);

      for (kk in seq_len(nbrOfBlocks)) {
        idxs <- froms[kk]:tos[kk];

        # Get data
        if (dd == 1) {
          xB <- x[idxs,,-dd,drop=FALSE];
          zB <- z[idxs,,drop=FALSE];
        } else if (dd == 2) {
          xB <- x[,idxs,-dd,drop=FALSE];
          zB <- z[,idxs,drop=FALSE];
        }
        dim(xB) <- dim(xB)[1:2];

        # Polish data
        zB2 <- FUN(zB, xB, ...);
        # Not needed anymore
        xB <- NULL;
        if (returnEffects) {
          effects[[dd]][kk,] <- zB2;
        }

        zB <- zB - zB2;
        # Not needed anymore
        zB2 <- NULL;

        # Update data
        if (dd == 1) {
          z[idxs,] <- zB;
        } else if (dd == 2) {
          z[,idxs] <- zB;
        }

        # Not needed anymore
        zB <- NULL;
      } # for (kk ...)
    } # for (dd ...)

    newSum <- sum(abs(z), na.rm=TRUE);
    converged <- (newSum == 0 || abs(newSum - oldSum) < tol * newSum);
    if (identical(converged, TRUE))
      break;

    oldSum <- newSum;
  } # for (ii ...)

  res <- list(residuals=z, converged=converged, iter=ii);
  if (returnEffects) {
    effects <- lapply(effects, FUN=drop);
    res <- c(list(row=effects[[1]], col=effects[[2]]), res);
  }
  class(res) <- c("matrixBlockPolish");

  res;
})



############################################################################
# HISTORY:
# 2008-04-02
# o Verified against median polish.
# o Created.
############################################################################
