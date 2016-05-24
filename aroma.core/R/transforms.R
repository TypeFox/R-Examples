log2center <- function(x) {
  # turn x into a *log2* ratio
  suppressWarnings({
    x <- log2(x);
  });
  # Remove the mean
  avg <- median(x, na.rm=TRUE);
  x <- x - avg;
  x;
}

log2neg <- function(x) {
  x <- -x;
  x[x <= 0] <- NA;
  log2(x);
}

log2pos <- function(x) {
  x[x <= 0] <- NA;
  log2(x);
}

log2abs <- function(x) {
  x <- abs(x);
  x[x <= 0] <- NA;
  log2(x);
}


sqrtcenter <- function(x) {
  # turn x into a *log2* ratio
  x <- sqrt(x);
  # Remove the mean
  avg <- median(x, na.rm=TRUE);
  x <- x - avg;
  x;
}

sqrtneg <- function(x) {
  x <- -x;
  x[x <= 0] <- NA;
  sqrt(x);
}

sqrtpos <- function(x) {
  x[x <= 0] <- NA;
  sqrt(x);
}

sqrtabs <- function(x) {
  x <- abs(x);
  x[x <= 0] <- NA;
  sqrt(x);
}

# "signed" sqrt()
sqrtsign <- function(x) {
  zero <- which(x == 0);
  neg <- which(x < 0);
  pos <- which(x > 0);
  x[zero] <- NA;
  x[pos] <-  sqrt( x[pos]);
  x[neg] <- -sqrt(-x[neg]);
  x;
}


############################################################################
# HISTORY:
# 2013-03-28
# o Added sqrtsign().
# 2008-03-17
# o Added log2center() and sqrtcenter().
# 2007-02-14
# o Created.
############################################################################
