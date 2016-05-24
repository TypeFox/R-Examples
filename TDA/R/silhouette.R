silhouette <-
function(Diag, p = 1, dimension = 1,
         tseq = seq(min(Diag[, 2:3]), max(Diag[, 2:3]), length = 500)) {

  if (((class(Diag) != "diagram" && class(Diag) != "matrix" &&
      !is.data.frame(Diag)) || NCOL(Diag) != 3) &&
      (!is.numeric(Diag) || length(Diag) != 3)) {
    stop("Diag should be a diagram, or a P by 3 matrix")
  }
  if (!is.numeric(dimension) || length(dimension) != 1 || dimension < 0) {
    stop("dimension should be an nonnegative integer")
  }
  if (!is.numeric(p) || any(p < 0)) {
    stop("p should be a vector of nonnegative number")
  }
  if (!is.numeric(tseq)) {
    stop("tseq should be a numeric vector")
  }

  if (is.numeric(Diag)) {
    Diag <- matrix(Diag, ncol = 3, dimnames = list(NULL, names(Diag)))
  }
  
  isNA <- length(which(Diag[, 1] == dimension))
  if (isNA == 0) {
    return(rep(0, length(tseq))) #in case there are no features with dimension "dimension"
  }
      
  Diag <- Diag[which(Diag[,1] == dimension), , drop = FALSE]

  left <- Diag[, 2]
  right <- Diag[, 3]
  Npoints <- length(left)
   
  ### Silhouette
  w <- outer(right - left, p, "^")
  w <- w %*% diag(1 / colSums(w), ncol = NCOL(w))
   
  ff <- matrix(0, nrow = length(tseq), ncol = length(p))

  for(i in seq_len(Npoints)) {
    tmp <- pmax(pmin(tseq - left[i], right[i] - tseq), 0)
    ff <- ff + tmp %*% w[i, , drop = FALSE]
  }
  return(ff)
}
