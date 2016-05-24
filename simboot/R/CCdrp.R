CCdrp <-
function (x, cmat) 
{
  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("Argument 'x'must be a matrix or data.frame!")
  }
  ngroup <- ncol(x)
  Nsim <- nrow(x)
  chains <- x
  if (!is.matrix(cmat)) {
    stop("'cmat' must be a matrix, specifying the contrast coefficients")
  }
  if (ngroup != ncol(cmat)) {
    stop("ncol(cmat) must be the same as the number of means in muvec")
  }
  cs <- apply(cmat, 1, sum)
  if (any(cs != 0)) {
    warning("Rows of cmat do not sum up to zero. Are the contrasts appropriately defined?")
  }
  nchains <- apply(X = chains, MARGIN = 1, FUN = function(x) {
    cmat %*% x
  })
  if (nrow(cmat) == 1) {
    nchains <- matrix(nchains, nrow = 1)
  }
  rownames(nchains) <- rownames(cmat)
  out <- list(chains = t(nchains), x = x, cmat = cmat)
  return(out)
}

