LD.thin <- function(x, threshold, max.dist = 100e3, beg = 1, end = ncol(x), 
                    dist.unit = c("bases", "indices"), extract = TRUE, keep = c("left", "right", "random")) {
  if(is.null(x@mu) | is.null(x@sigma))
    stop("LD.thin needs mu and sigma to be set for LD computation (use set.stats)")

  dist.unit <- match.arg(dist.unit)
  if(dist.unit == "indices") x@snps$pos = seq_len(ncol(x))

  if( all(x@snps$pos == x@snps$pos[1]) )
    stop("Position of SNPs must be available")

  keep <- match.arg(keep)
  if(keep == "left") {
    w <- .Call("gg_ld_thin_left", x@bed, x@mu, x@sigma, threshold, as.integer(x@snps$pos), 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L)
  } else if (keep == "right"){
    w <- .Call("gg_ld_thin_right", x@bed, x@mu, x@sigma, threshold, as.integer(x@snps$pos), 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L)
  } else {
    w <- .Call("gg_ld_thin_random", x@bed, x@mu, x@sigma, threshold, as.integer(x@snps$pos), 
          as.integer(x@snps$chr), as.integer(max.dist), as.integer(beg)-1L, as.integer(end)-1L)
  }

  if(!extract) return(w)
  x[ , seq(beg,end)[w] ]
}

