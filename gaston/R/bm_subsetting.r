## Extraction
void_snps_stats <- function(x) {
  a <- names(x)[ !(names(x) %in% snpstatnames) ]
  x[,a]
} 

void_ped_stats <- function(x) {
  a <- names(x)[ !(names(x) %in% pedstatnames) ]
  x[,a]
} 

setMethod("[", signature(x="bed.matrix",i="numeric",j="missing", drop="missing"), 
    function( x, i, j) {
      if(any(i <= 0)) i <- (1:nrow(x))[i]
      x@bed <- .Call('gg_extract_inds_indices', x@bed, i) 
      x@ped <- x@ped[i,]
      x@snps <- void_snps_stats(x@snps)
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, set.ped_stats = FALSE, verbose = FALSE)
      x
    } );

setMethod("[", signature(x="bed.matrix",i="logical",j="missing", drop="missing"), 
    function( x, i, j) {
      if(any(is.na(i))) stop("NAs not allowed")   
      x@bed <- .Call('gg_extract_inds_bool', x@bed, i) 
      x@ped <- x@ped[i,]
      x@snps <- void_snps_stats(x@snps)
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, set.ped_stats = FALSE, verbose = FALSE)
      x
    } );
  
setMethod("[", signature(x="bed.matrix",i="missing",j="numeric", drop="missing"), 
    function( x, i, j) {    
      if(any(j <= 0)) j <- (1:ncol(x))[j]
      x@bed <- .Call('gg_extract_snps_indices', x@bed, j) 
      x@snps <- x@snps[j,]
      x@p <- x@p[j]
      x@mu <- x@mu[j]
      x@sigma <- x@sigma[j]
      x@ped <- void_ped_stats(x@ped)
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, set.snps_stats = FALSE, set.p = FALSE, set.mu_sigma = FALSE, verbose = FALSE)
      x
    } );
  
setMethod("[", signature(x="bed.matrix",i="missing",j="logical", drop="missing"), 
    function( x, i, j) {    
      if(any(is.na(j))) stop("NAs not allowed")   
      x@bed <- .Call('gg_extract_snps_bool', x@bed, j)
      x@snps <- x@snps[j,]
      x@p <- x@p[j]
      x@mu <- x@mu[j]
      x@sigma <- x@sigma[j]
      x@ped <- void_ped_stats(x@ped)
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, set.snps_stats = FALSE, set.p = FALSE, set.mu_sigma = FALSE, verbose = FALSE)
      x
    } );
  
setMethod("[", signature(x="bed.matrix",i="logical",j="logical", drop="missing"), 
    function( x, i, j) {    
      if(any(is.na(j))) stop("NAs not allowed")   
      x@bed <- .Call('gg_extract_snps_bool', x@bed, j)
      x@snps <- x@snps[j,]
      x@p <- x@p[j]
      x@mu <- x@mu[j]
      x@sigma <- x@sigma[j]

      if(any(is.na(i))) stop("NAs not allowed")   
      x@bed <- .Call('gg_extract_inds_bool', x@bed, i) 
      x@ped <- x@ped[i,]
      x@ped <- void_ped_stats(x@ped)
      x@snps <- void_snps_stats(x@snps)
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = FALSE)
      x
    } );
  

setMethod("[", signature(x="bed.matrix",i="logical",j="numeric", drop="missing"), 
    function( x, i, j) {    
      if(any(j <= 0)) j <- (1:ncol(x))[j]
      x@bed <- .Call('gg_extract_snps_indices', x@bed, j) 
      x@snps <- x@snps[j,]
      x@p <- x@p[j]
      x@mu <- x@mu[j]
      x@sigma <- x@sigma[j]

      if(any(is.na(i))) stop("NAs not allowed")   
      x@bed <- .Call('gg_extract_inds_bool', x@bed, i) 
      x@ped <- x@ped[i,]
      x@ped <- void_ped_stats(x@ped)
      x@snps <- void_snps_stats(x@snps)
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = FALSE)
      x
    } );

setMethod("[", signature(x="bed.matrix",i="numeric",j="logical", drop="missing"), 
    function( x, i, j) {    
      if(any(is.na(j))) stop("NAs not allowed")   
      x@bed <- .Call('gg_extract_snps_bool', x@bed, j)
      x@snps <- x@snps[j,]
      x@p <- x@p[j]
      x@mu <- x@mu[j]
      x@sigma <- x@sigma[j]

      if(any(i <= 0)) i <- (1:nrow(x))[i]
      x@bed <- .Call('gg_extract_inds_indices', x@bed, i) 
      x@ped <- x@ped[i,]
      x@ped <- void_ped_stats(x@ped)
      x@snps <- void_snps_stats(x@snps)
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = FALSE)
      x
    } );
  

setMethod("[", signature(x="bed.matrix",i="numeric",j="numeric", drop="missing"), 
    function( x, i, j) {    
      if(any(j <= 0)) j <- (1:ncol(x))[j]
      x@bed <- .Call('gg_extract_snps_indices', x@bed, j) 
      x@snps <- x@snps[j,]
      x@p <- x@p[j]
      x@mu <- x@mu[j]
      x@sigma <- x@sigma[j]

      if(any(i <= 0)) i <- (1:nrow(x))[i]
      x@bed <- .Call('gg_extract_inds_indices', x@bed, i) 
      x@ped <- x@ped[i,]
      x@ped <- void_ped_stats(x@ped)
      x@snps <- void_snps_stats(x@snps)
      if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = FALSE)
      x
    } );

