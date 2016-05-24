setMethodS3("fit2d", "matrix", function(M, MARGIN=1, spar=0.7, h=20, ...) {
  ## aroma.light::robustSmoothSpline()
  use("aroma.light", how="load");

  nr <- nrow(M);
  nc <- ncol(M);

  if (MARGIN == 1) {
    # Normalize row by row
    n <- nr;
    byrow <- TRUE;
  } else if (MARGIN == 2) {
    # Normalize column by column
    n <- nc;
    byrow <- FALSE;
  }

  mu <- matrix(NA, nrow=nr, ncol=nc);
  x <- 1:n;
  nbands <- ceiling(n/h);
  for (kk in 1:nbands) {
    if (kk %% 10 == 0) print(kk);

    # Rows/columns in this band
    rr0 <- (kk-1)*h;
    rr <- seq(from=rr0+1, to=min(rr0+h,n));

    # Extract data
    if (MARGIN == 1) {
      Mb <- M[rr,,drop=FALSE];
      Mb <- colMedians(Mb, na.rm=TRUE);
    } else if (MARGIN == 2) {
      Mb <- M[,rr,drop=FALSE];
      Mb <- rowMedians(Mb, na.rm=TRUE);
    }

    # Fit smooth curve (1d)
    ok <- which(is.finite(Mb));
    fit <- aroma.light::robustSmoothSpline(x[ok], Mb[ok], spar=spar);
    # Not needed anymore
    ok <- NULL;
    Mp <- predict(fit, x=x)$y;
    Mp <- matrix(Mp, nrow=length(rr), ncol=length(Mb), byrow=byrow);

    if (MARGIN == 1) {
      mu[rr,] <- Mp;
    } else if (MARGIN == 2) {
      mu[,rr] <- Mp;
    }
  }

  mu;
}, private=TRUE)


setMethodS3("norm2d", "matrix", function(M, MARGIN=c(1,2), spar=0.7, h=20, ...) {
  n <- length(MARGIN);
  spar <- rep(spar, length.out=n);
  h <- rep(h, length.out=n);

  Mn <- M;
  for (kk in seq_len(n)) {
    mu <- fit2d(Mn, MARGIN=MARGIN[kk], spar=spar[kk], h=h[kk], ...);
    Mn <- Mn-mu;
  }

  Mn;
}, private=TRUE)


setMethodS3("calcMargins", "matrix", function(M, unshift=FALSE, ...) {
  if (unshift) {
    M <- M - median(M, na.rm=TRUE);
  }
  list(
    rows=rowMedians(M, na.rm=TRUE),
    cols=colMedians(M, na.rm=TRUE)
  );
}, private=TRUE)


############################################################################
# HISTORY:
# 2014-08-27
# o SPEEDUP: Now calcMargins() for matrix utilized colMedians().
# 2012-08-17
# o ROBUSTNESS: Now fit2d() for matrix utilizes use() for aroma.light.
# 2012-04-16
# o Now fit2d() explicitly require the 'aroma.light' package.
# o Dropped internal colMedians() from fit2d(); already in matrixStats.
# 2008-03-19
# o Created.
############################################################################
