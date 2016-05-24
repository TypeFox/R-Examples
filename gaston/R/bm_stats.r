## Basic Statistics
set.stats <- function(x, set.ped_stats = TRUE, set.snps_stats = TRUE, set.p = TRUE, set.mu_sigma = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  if( is(x)!='bed.matrix' ) stop('x must be an object of class bed.matrix')
  if(!is.logical(set.p) | !is.logical(set.mu_sigma) | !is.logical(set.snps_stats) | !is.logical(set.ped_stats)) 
    stop('set.* arguments must be logical')

  if(set.ped_stats & set.snps_stats) { 
    st <- .Call('gg_geno_stats', PACKAGE = 'gaston', x@bed) 


    # completer snps
    st$snps$callrate <- 1-st$snps$NAs/nrow(x)
    n <- nrow(x) - st$snps$NAs;
    pp <- (2*st$snps$N2 + st$snps$N1)/(2*n);
    st$snps$maf <- pmin(pp,1-pp)
    st$snps$hz <- st$snps$N1/n

    # completer inds/ped
    st$inds$callrate <- 1-st$inds$NAs/ncol(x)
    n <- ncol(x) - st$inds$NAs
    st$inds$hz <- st$inds$N1/n

    # pour set_mu_sigma
    N <- nrow(x)
    s <- sqrt( (st$snps$N1 + 4*st$snps$N2 + (4*st$snps$NAs)*pp**2)/(N-1)  - (pp)**2*(4*N/(N-1)) );

    x@snps[, names(st$snps)] <- st$snps
    x@ped[ , names(st$inds)] <- st$inds
    if(verbose) cat("ped stats and snps stats have been set. \n") 
  } else if(set.snps_stats) {
    snps <- .Call('gg_geno_stats_snp', PACKAGE = 'gaston', x@bed)

    # completer snps
    snps$callrate <- 1-snps$NAs/nrow(x)
    n <- nrow(x) - snps$NAs;
    pp <- (2*snps$N2 + snps$N1)/(2*n);
    snps$maf <- pmin(pp,1-pp)
    snps$hz <- snps$N1/n

    # pour set_mu_sigma
    N <- nrow(x)
    s <- sqrt( (snps$N1 + 4*snps$N2 + (4*snps$NAs)*pp**2)/(N-1)  - (pp)**2*(4*N/(N-1)) );
    x@snps[, names(snps)] <- snps
    if(verbose) cat("snps stats have been set. \n") 
  } else if(set.ped_stats) { 
    inds <- .Call('gg_geno_stats_ind', PACKAGE = 'gaston', x@bed) 

    # completer inds/ped
    inds$callrate <- 1-inds$NAs/ncol(x)
    n <- ncol(x) - inds$NAs
    inds$hz <- inds$N1/n

    x@ped[ , names(inds)] <- inds
    if(verbose) cat("ped stats have been set. \n") 
  }

  if( all(c("N0", "N1", "N2", "NAs") %in% names(x@snps) )) {
    if(set.p | set.mu_sigma) {
      n <- nrow(x) - x@snps$NAs;
      pp <- (2*x@snps$N2 + x@snps$N1)/(2*n);
    }
    if(set.mu_sigma)
      s <- sqrt( (x@snps$N1 + 4*x@snps$N2)/(n-1) - n*(2*pp)**2/(n-1) );
    if(set.p) {
      x@p <- pp
      if(verbose) cat("'p' has been set. \n")
    }
    if(set.mu_sigma) {
      x@mu <- 2*pp; x@sigma <- s
      if(verbose) cat("'mu' and 'sigma' have been set.\n");
    }
  }
  x
}


set.hwe <- function(x, method = c("chisquare", "exact"), verbose = getOption("gaston.verbose",TRUE)) {
  if( !all(c("N0", "N1", "N2") %in% names(x@snps) )) {
    if(verbose) cat("Computing basic stats\n")
    x <- set.stats(x)
  }
  method <- match.arg(method)
  if(method == 'chisquare') {
    if(verbose) cat("Computing HW chi-square p-values\n")
    hwe_ <- .Call('gg_hwe_chi', PACKAGE = 'gaston', x@snps$N0, x@snps$N1, x@snps$N2)
  } else {
    if(verbose) cat("Computing HW exact test p-values\n")
    hwe_ <- .Call('gg_hwe', PACKAGE = 'gaston', x@snps$N0, x@snps$N1, x@snps$N2)
  }
  x@snps$hwe <- hwe_
  x
}

