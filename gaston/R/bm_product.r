setMethod("%*%", signature(x="bed.matrix",y="matrix"), 
  function(x, y) {
    if(x@standardize_mu_sigma) 
      a <- .Call("gg_m4_loading_to_pc_ms", PACKAGE="gaston", x@bed, x@mu, x@sigma, y)
    else if(x@standardize_p)
      a <- .Call("gg_m4_loading_to_pc_p",  PACKAGE="gaston", x@bed, x@p, y)
    else
      stop("bed.matrix must be center/scaled for product")
    if(!is.null(x@ped$id)) {
      if(anyDuplicated(x@ped$id) == 0)
        rownames(a) <- x@ped$id
    }
    colnames(a) <- colnames(y)
    a
  }
);

setMethod("%*%", signature(x="matrix",y="bed.matrix"), 
  function(x, y) {
    if(y@standardize_mu_sigma)
      a <- t(.Call("gg_m4_pc_to_loading_ms", PACKAGE="gaston", y@bed, y@mu, y@sigma, t(x)))
    else if(y@standardize_p)
      a <- t(.Call("gg_m4_pc_to_loading_p",  PACKAGE="gaston", y@bed, y@p, t(x)))
    else
      stop("bed.matrix must be center/scaled for product")
    if(!is.null(y@snps$id)) {
      if(anyDuplicated(y@snps$id) == 0)
        colnames(a) <- y@snps$id
    }
    rownames(a) <- rownames(x)
    a
  }
);

setMethod("%*%", signature(x="bed.matrix",y="vector"), 
  function(x, y) {
    x %*% as.matrix(y)
  }
);

setMethod("%*%", signature(x="vector",y="bed.matrix"), 
  function(x, y) {
    matrix(x, nrow = 1) %*% y
  }
);


bed.loadings <- function(x, pc) {
  if(is.vector(pc)) dim(pc) <- c(length(pc),1)
  if(!is.matrix(pc)) stop("pc must be a matrix of PCs")
  if(!is(x, "bed.matrix")) stop("x must be a bed.matrix")
  if(!x@standardize_mu_sigma & !x@standardize_p) {
    if(!is.null(x@p)) 
      x@standardize_p <- TRUE
    else stop("Can't center/scale x for LD computation (use set.stat)\n")
 }

  if(x@standardize_mu_sigma)
    a <- .Call("gg_m4_pc_to_loading_ms", PACKAGE="gaston", x@bed, x@mu, x@sigma, pc)
  else if(x@standardize_p)
    a <- .Call("gg_m4_pc_to_loading_p",  PACKAGE="gaston", x@bed, x@p, pc)

  a <- a/sqrt(ncol(x))
  a <- sweep(a, 2L, sqrt(colSums(a**2)), "/")
  if(!is.null(x@snps$id)) {
    if(anyDuplicated(x@snps$id) == 0)
      rownames(a) <- x@snps$id
  }
  a
}


