setGeneric("cbind", signature="...")
setMethod("cbind", signature=c(...="bed.matrix"), definition= function(..., deparse.level=1)
    {
      L <- list(...)
      M <- lapply(L, function(x) x@bed)
      snps <- Reduce(rbind, lapply(L, function(x) x@snps))
      if(any(duplicated(snps$id))) 
        warning("Duplicated SNPs id's")
      if(any(duplicated(snps$pos)))
        warning("Duplicated SNPs pos's")
      bed <- .Call("gg_bind_snps",  PACKAGE="gaston", M)
      x <- new("bed.matrix", bed = bed, snps = snps, ped = L[[1]]@ped,
               p = NULL, mu = NULL, sigma = NULL, 
               standardize_p = FALSE, standardize_mu_sigma = FALSE )
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, set.snps_stats = FALSE, set.p = FALSE, set.mu_sigma = FALSE, verbose = FALSE)
      x
    }
)

setGeneric("rbind", signature="...")
setMethod("rbind", signature=c(...="bed.matrix"), definition= function(..., deparse.level=1)
    {
      L <- list(...)
      M <- lapply(L, function(x) x@bed)
      ped <- Reduce(rbind, lapply(L, function(x) x@ped))
      if(any( duplicated(ped$famid) & duplicated(ped$id) )) 
        warning("Duplicated individuals (same family and individual id)")
      bed <- .Call("gg_bind_inds",  PACKAGE="gaston", M)
      x <- new("bed.matrix", bed = bed, snps = L[[1]]@snps, ped = ped,
               p = NULL, mu = NULL, sigma = NULL, 
               standardize_p = FALSE, standardize_mu_sigma = FALSE )
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, set.ped_stats = FALSE, verbose = FALSE)
      x
    }
)
