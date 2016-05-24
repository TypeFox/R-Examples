testEqualityTcnByT <- function(dataL, dataR, alpha=0.02, ...) {
  nL <- length(dataL);
  nR <- length(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }
  print(c(nL=nL, nR=nR));

  # Test:
  #  H0: muL == muR
  #  H1: muL != muR
  fit <- t.test(dataL, dataR, paired=FALSE, var.equal=TRUE, 
                                            alternative="two.sided");
  t <- fit$statistic;
  p <- fit$p.value;
  fit$isSignificant <- (p < alpha);
  isEqual <- (!fit$isSignificant);
  attr(isEqual, "fit") <- fit;

  isEqual;
} # testEqualityTcnByT()


testEqualityTcnByMean <- function(dataL, dataR, alpha=0.02, ...) {
  nL <- length(dataL);
  nR <- length(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }

  muL <- mean(dataL, na.rm=TRUE);
  muR <- mean(dataR, na.rm=TRUE);
  delta <- abs(muL - muR);

  isEqual <- (delta <= alpha);

  isEqual;
} # testEqualityTcnByMean()


testEqualityC1C2ByMean <- function(dataL, dataR, alpha=0.02, ...) {
  # Sanity checks
  stopifnot(is.matrix(dataL));
  stopifnot(is.matrix(dataR));

  nL <- nrow(dataL);
  nR <- nrow(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }

  muL <- colMeans(dataL, na.rm=TRUE);
  muR <- colMeans(dataR, na.rm=TRUE);
  deltas <- abs(muL - muR);

  # Sanity checks
  stopifnot(length(deltas) == 2);

  # Eucledian distance
  delta <- sqrt(sum(deltas^2, na.rm=TRUE));

  isEqual <- (delta <= alpha);

  fit <- list(
    muL = muL,
    muR = muR,
    deltas = deltas,
    delta = delta,
    alpha = alpha
  );
  attr(isEqual, "fit") <- fit;

  isEqual;
} # testEqualityC1C2ByMean()


testEqualityC1orC2ByT <- function(dataL, dataR, alpha=0.02, ...) {
  # Sanity checks
  stopifnot(is.matrix(dataL));
  stopifnot(is.matrix(dataR));

  nL <- nrow(dataL);
  nR <- nrow(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }

  # Test C1 and C2 independently
  fitList <- lapply(1:2, FUN=function(cc) {
    t.test(dataL[,cc], dataR[,cc], paired=FALSE, var.equal=TRUE, 
                                              alternative="two.sided");
  });

  # Extract the significance values
  p <- sapply(fitList, FUN=function(fit) fit$p.value);

  fit <- list(
    fitList = fitList,
    p.value = p
  );
  fit$isSignificant <- (p < alpha);
  isEqual <- (!fit$isSignificant);
  isEqual <- all(isEqual);
  attr(isEqual, "fit") <- fit;

  isEqual;
} # testEqualityC1orC2ByT()


testEqualityC1C2ByChiSq <- function(dataL, dataR, alpha=0.02, ...) {
  # Sanity checks
  stopifnot(is.matrix(dataL));
  stopifnot(is.matrix(dataR));

  nL <- nrow(dataL);
  nR <- nrow(dataR);
  if (nL < 1 || nR < 1) {
    return(NA);    
  }

  muL <- colMeans(dataL, na.rm=TRUE);
  sdL <- colSds(dataL, na.rm=TRUE);
  muR <- colMeans(dataR, na.rm=TRUE);
  sdR <- colSds(dataR, na.rm=TRUE);
##  printf("nL=%d, muL=%.3f, sdL=%.3f\n", nL, muL, sdL);
##  printf("nR=%d, muR=%.3f, sdR=%.3f\n", nR, muR, sdR);

##   delta <- muL-muR;
##   s2 <- (nL-1)*sdL^2 + (nR-1)*sdR^2;
##   s2 <- s2 / (nL+nR);
##   print(s2);
##   print(colVars(rbind(dataL, dataR), na.rm=TRUE));
##   s <- sqrt(s2 / (1/nL+1/nR));
##   printf("delta=%.3f, s=%.3f\n", delta, s);

  # Test C1 and C2 independently
  fitList <- lapply(1:2, FUN=function(cc) {
    yL <- dataL[,cc,drop=TRUE];
    yR <- dataR[,cc,drop=TRUE];
##     str(yL); str(yR);
    t.test(yL, yR, paired=FALSE, var.equal=TRUE, alternative="two.sided");
  });

  # Extract the significance values
  t <- sapply(fitList, FUN=function(fit) fit$statistic);
  p <- sapply(fitList, FUN=function(fit) fit$p.value);
  printf("T[i] ~ t(df=.)-distribution: t=(%.2g,%.2g)\n", t[1], t[2]);
  printf("t-distribution: p=(%.2g,%.2g)\n", p[1], p[2]);

  # Back-transform the p-values into Gaussian N(0,1) quantiles
  x <- qnorm(1-p, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE);
  printf("X[1],X[2] ~ N(0,1)-distribution: p=(%.2g,%.2g)\n", p[1], p[2]);
  printf("X[1],X[2] ~ N(0,1)-distribution: x=(%.2g,%.2g)\n", x[1], x[2]);

  # Transform them to chi-square with 2 d.f.
  y <- x^2;
  printf("Y[i]=X[i]^2 ~ ChiSq(df=1)-distributions: y=(%.2g,%.2g)\n", y[1], y[2]);
  z <- sum(y);
  printf("Z = Y[1]+Y[2] ~ ChiSq(df=2)-distributions: z=%.2g\n", z);
  p <- 1-pchisq(z, df=2, ncp=0, lower.tail=TRUE, log.p=FALSE);
  printf("Z = Y[1]+Y[2] ~ ChiSq(df=2)-distributions: p=%.2g\n", p);

  fit <- list(
#    statistic = y,
    p.value = p
  );
  fit$isSignificant <- (p < alpha);
  isEqual <- (!fit$isSignificant);
  attr(isEqual, "fit") <- fit;

  isEqual;
} # testEqualityC1C2ByChiSq()



############################################################################
# HISTORY:
# 2011-01-23
# o Added testEqualityC1C2ByChiSq() stub.
# 2011-01-22
# o Added testEqualityC1orC2ByT().
# 2011-01-18
# o Added testEqualityC1C2ByMean().
# 2011-01-12
# o Added testEqualityTCN().
# o Created.
############################################################################
