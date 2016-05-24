GRM <- function(x, chunk = 1L) {
  if(!x@standardize_mu_sigma & !x@standardize_p) {
    if(!is.null(x@p)) x@standardize_p <- TRUE
    else stop("Can't center/scale x for LD computation (use set.stat)\n")
  }

  if(x@standardize_mu_sigma) {
    w <- ifelse(x@sigma == 0, 0, 1/x@sigma/sqrt(ncol(x)))   ### BEWARE ncol(x)-1 !!!
    K <- .Call('gg_Kinship', PACKAGE = 'gaston', x@bed, x@mu, w, chunk) 
  } else if(x@standardize_p) 
    K <- .Call('gg_Kinship_p', PACKAGE = 'gaston', x@bed, x@p, chunk)
  else 
    stop("bed.matrix must be centered/scaled for GRM computation")

  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0) 
      rownames(K) <- colnames(K) <- x@ped$id
  }

  K
}

reshape.GRM <- function(K, include = c(-Inf, +Inf), exclude) {
  diag(K) <- NA
  if(missing(exclude))
    w <- which(include[1] < K & K < include[2])
  else 
    w <- which(include[1] < K & K < include[2] & (K < exclude[1] | K > exclude[2]))
  I <- row(K)[w]
  J <- col(K)[w]
  R <- K[w]
  ww <- (I < J)
  i <- I[ww];
  j <- J[ww];
  data.frame(i = i, j = j, id_i = rownames(K)[i], id_j = colnames(K)[j], k = R[ww])
}
