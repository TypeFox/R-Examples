# The exact/permutation version of Jonckheere-Terpstra test
# asymptotic version is equivalent to cor.test for Kendall's tau

# this function computes the pdf of the statistic for small samples
# using Hardings generating function algorithm - blows up quickly
djonckheere <- function(gsize) {
  ng <- length(gsize)
  cgsize <- rev(cumsum(rev(gsize)))
  mxsum <- sum(gsize[-ng]*cgsize[-1]) + 1
  jrsum <- rep(0, mxsum)
  jrsum[1] <- 1
  zz <- .Fortran("djonck",
                 as.integer(mxsum),
                 jrsum=as.double(jrsum),
                 as.integer(ng),
                 as.integer(cgsize))
  zz$jrsum/sum(zz$jrsum)
}

# this function computes the pdf using the convolution by Mark van de Wiel
jtpdf <- function(gsize) {
  ng <- length(gsize)
  cgsize <- rev(cumsum(rev(gsize)))
  mxsum <- sum(gsize[-ng]*cgsize[-1]) + 1
  zz <- .Fortran("jtpdf",
                 as.integer(mxsum),
                 pdf=double(mxsum),
                 as.integer(ng),
                 as.integer(cgsize),
                 double(mxsum),
                 double(mxsum))
  zz$pdf
}

# now do the test as in wilcox.test
jonckheere.test <- function(x, g, alternative = c("two.sided", "increasing", "decreasing"), nperm=NULL) {
  if(!is.numeric(x)) stop("data values should be numeric")
  if(!is.numeric(g) & !is.ordered(g)) stop("group should be numeric or ordered factor")
  if (length(g) != length(x)) stop("lengths of data values and group don't match")
  alternative <- match.arg(alternative)
  METHOD <- "Jonckheere-Terpstra test"
  PERM <- !missing(nperm)
  # keep only the finite values
  ii <- is.finite(x) & is.finite(g)
  x <- x[ii]
  g <- g[ii]
  n <- length(x)
  if (n == 0) stop("either x (data) or g (group) missing for all observations")
  TIES <- length(unique(x)) != n
  gsize <- table(g)
  ng <- length(gsize)
  if (ng <= 1) stop("only one group has non-missing data")
  cgsize <- c(0,cumsum(gsize))
  x <- x[order(g)]
  jtrsum <- jtmean <- jtvar <- 0
  for(i in 1:(ng-1)) {
    na <- gsize[i]
    nb <- n-cgsize[i+1]
    jtrsum <- jtrsum + sum(rank(x[(cgsize[i]+1):n])[1:na]) - na*(na+1)/2
    jtmean <- jtmean + na*nb/2
    jtvar <- jtvar + na*nb*(na+nb+1)/12
  }
# this jtrsum will be small if data are increasing and large if decreasing
# to reverse this use 2*jtmean - jtrsum  
  jtrsum <- 2*jtmean - jtrsum
  STATISTIC <- jtrsum
  names(STATISTIC) <- "JT"
  if (PERM) {
    PVAL <- jtperm.p(x, ng, gsize, cgsize, alternative, nperm) 
  } else {
    if (n > 100 | TIES) {
      warning("Sample size > 100 or data with ties \n p-value based on normal approximation. Specify nperm for permutation p-value")
      zstat <- (STATISTIC-jtmean)/sqrt(jtvar)
      PVAL <- pnorm(zstat)
      PVAL <- switch(alternative,
                     "two.sided" = 2*min(PVAL, 1-PVAL, 0.5),
                     "increasing" = 1-PVAL,
                     "decreasing" = PVAL)
    } else {
      dPVAL <- sum(jtpdf(gsize)[1:(jtrsum+1)])
      iPVAL <- 1-sum(jtpdf(gsize)[1:(jtrsum)])
      PVAL <- switch(alternative,
                     "two.sided" = 2*min(iPVAL, dPVAL, 0.5),
                     "increasing" = iPVAL,
                     "decreasing" = dPVAL)
    }
  }
  RVAL <- list(statistic = STATISTIC,
               p.value = as.numeric(PVAL),
               alternative = alternative,
               method = METHOD)
  class(RVAL) <- "htest"
  RVAL
}

jtperm.p <- function(x, ng, gsize, cgsize, alternative, nperm) {
  n <- length(x)
  pjtrsum <- rep(0, nperm)
  for (np in 1:nperm){
    jtrsum <- 0
    for(i in 1:(ng-1)) {
      na <- gsize[i]
      nb <- n-cgsize[i+1]
# this jtrsum will be small if data are increasing and large if decreasing
      jtrsum <- jtrsum + sum(rank(x[(cgsize[i]+1):n])[1:na]) - na*(na+1)/2
    }
    pjtrsum[np] <- jtrsum
    # permute the data; this way the first value is the original statistic
    x <- sample(x)
  }
  # one-sided p-values
  # number of permuted values at least as small as original
  iPVAL <- sum(pjtrsum <= pjtrsum[1])/nperm
  # number of permuted values at least as large as original
  dPVAL <- sum(pjtrsum >= pjtrsum[1])/nperm
  # return p-value for the alternative of interest
  PVAL <- switch(alternative,
                 "two.sided" = 2*min(iPVAL, dPVAL, 0.5),
                 "increasing" = iPVAL,
                 "decreasing" = dPVAL)
  PVAL
}
