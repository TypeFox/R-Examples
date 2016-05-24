setMethodS3("fitC1C2Peaks", "PairedPSCBS", function(fit, ..., tol=0.05, onError=c("error", "warning", "skip"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'onError':
  onError <- match.arg(onError);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting the relationship between peaks in C1 and C2");

  d1d2 <- fitC1C2Densities(fit, orderBy="x");
  pList <- d1d2$pList;
  verbose && print(verbose, pList);

  D <- outer(pList$C1$x, pList$C2$x, FUN="-");
  idxs <- which(abs(D) < tol, arr.ind=TRUE);
  verbose && print(verbose, idxs);

  if (nrow(idxs) < 2) {
    msg <- sprintf("Cannot fit relationship between C1 and C2. Too few common peaks: %d", nrow(idxs));
    if (onError == "error") {
      plotC1C2Grid(fit);
      throw(msg);
    }
    if (onError == "warning") {
      warning(msg);
    }

    msg <- sprintf("Skipping fitting the relationship between peaks in C1 and C2. %s", msg);
    verbose && cat(verbose, msg);
    warning(msg);
    verbose && exit(verbose);
    return(NULL);
  }

  c1 <- pList$C1$x[idxs[,1]];
  c2 <- pList$C2$x[idxs[,2]];
  dd <- cbind(pList$C1$density[idxs[,1]], pList$C2$density[idxs[,2]]);
  dd <- rowSums(dd^2);
  w <- dd / sum(dd);
  verbose && print(verbose, cbind(c1=c1, c2=c2, weights=w));

  f <- lm(c2 ~ 1 + c1, weights=w);
  verbose && print(verbose, f);

  a <- coef(f)[1];
  b <- coef(f)[2];
  params <- list(a=a, b=b);
  verbose && print(verbose, params);
  res <- list(fit=f, params=params);

  verbose && str(verbose, res);

  verbose && exit(verbose);

  res;
}) # fitC1C2Peaks()


setMethodS3("fitC1C2Densities", "PairedPSCBS", function(fit, adjust=0.2, tol=0.05, orderBy=c("density", "x"), ...) {
  orderBy <- match.arg(orderBy);

  data <- extractC1C2(fit);
  n <- data[,4, drop=TRUE];
  n <- sqrt(n);
  w <- n/sum(n, na.rm=TRUE);
  adjust <- 0.2;

  dList <- list();
  for (cc in 1:2) {
    y <- data[,cc];
    ok <- is.finite(y) & is.finite(w);
    y <- y[ok];
    wt <- w[ok]/sum(w[ok]);
    d <- density(y, weights=wt, adjust=adjust);
    dList[[cc]] <- d;
  }
  names(dList) <- colnames(data)[1:2];

  type <- NULL; rm(list="type");  # To please R CMD check
  pList <- lapply(dList, FUN=function(d) {
    p <- .findPeaksAndValleys(d, tol=tol);
    p <- subset(p, type == "peak");
    p <- p[order(p[[orderBy]], decreasing=c("x"=FALSE, "density"=TRUE)[orderBy]),,drop=FALSE];
  });
  names(pList) <- names(dList);

  return(list(dList=dList, pList=pList));
}) # fitC1C2Densities()


##############################################################################
# HISTORY
# 2011-10-17 [HB]
# o ROBUSTNESS: Now calibrateC1C2() no longer gives an error if it cannot
#   fit diagonals due to lack of data points.
# 2011-10-16 [HB]
# o Now using getSegments(fit) instead of fit$output.
# 2011-07-10 [HB]
# o Updated code to work with the new column names in PSCBS v0.11.0.
# 2010-10-10 [HB]
# o Added fitC1C2Peaks().
# o Added calibrateC1C2() for PairedPSCBS.
##############################################################################
