setMethodS3("fitLoessKD", "matrix", function(X, Y, ...) {
  # Argument 'X':
  X <- as.matrix(X);

  # Argument 'Y':
  Y <- as.matrix(Y);
  if (!all(dim(Y) == dim(X))) {
    dimXStr <- paste(dim(X), collapse="x");
    dimYStr <- paste(dim(Y), collapse="x");
    throw("The dimensions of argument 'Y' and 'X' do not match: ",
          dimYStr, " != ", dimXStr);
  }

  # Record input data
  data <- list(X=X, Y=Y);

  # Drop missing values
  ok <- rowAlls(is.finite(X) & is.finite(Y));
  X <- X[ok,,drop=FALSE];
  Y <- Y[ok,,drop=FALSE];

  # Estimate relationship for each dimension
  fitList <- list();
  for (cc in seq_len(ncol(X))) {
    fitList[[cc]] <- loess(Y[,cc] ~ X);
  }

  # Create predict function
  predictY <- function(X, ...) {
    # Argument 'X':
    if (ncol(X) != length(fitList)) {
      throw("The number of columns in argument 'X' does not match the number of fitted dimensions: ", ncol(X), " != ", length(fitList));
    }

    X <- as.matrix(X);
    naValue <- as.double(NA);
    Y <- array(naValue, dim=dim(X));
    ok <- rowAlls(is.finite(X));
    ok <- which(ok);
    X <- X[ok,,drop=FALSE];
    for (cc in seq_len(ncol(X))) {
      fit <- fitList[[cc]];
      yPred <- predict(fit, newdata=X);
      Y[ok,cc] <- yPred;
    }
    Y;
  } # predictY()

  fit <- list(data=data, fitList=fitList, predictY=predictY);

  class(fit) <- c("LoessKDFit", class(fit));

  fit;
});

 
setMethodS3("normalizeLoessKD", "matrix", function(X, fit, ...) {
  # Argument 'fit':
  fit <- Arguments$getInstanceOf(fit, "LoessKDFit");

  fit$predictY(X, ...);
})


setMethodS3("fitLoessKD", "data.frame", function(X, ...) {
  X <- as.matrix(X);
  fitLoessKD(X, ...);
})

setMethodS3("normalizeLoessKD", "data.frame", function(X, ...) {
  XN <- X;
  X <- as.matrix(X);
  XN <- normalizeLoessKD(X, ...);
  XN <- as.data.frame(XN);
  XN;
})


##############################################################################
# HISTORY
# 2010-09-19 [HB]
# o Created.
##############################################################################
