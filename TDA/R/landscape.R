landscape <-
function(Diag, dimension = 1, KK = 1,
         tseq = seq(min(Diag[, 2:3]), max(Diag[, 2:3]), length = 500)) {
    
  if (((class(Diag) != "diagram" && class(Diag) != "matrix" &&
      !is.data.frame(Diag)) || NCOL(Diag) != 3) &&
      (!is.numeric(Diag) || length(Diag) != 3)) {
    stop("Diag should be a diagram, or a P by 3 matrix")
  }
  if (!is.numeric(dimension) || length(dimension) != 1 || dimension < 0) {
    stop("dimension should be an nonnegative integer")
  }
  if (!is.numeric(KK) || any(KK <= 0)) {
    stop("KK should be a vector of positive integer")
  }
  if (!is.numeric(tseq)) {
    stop("tseq should be numeric")
  }

  if (is.numeric(Diag)) {
    Diag <- matrix(Diag, ncol = 3, dimnames = list(NULL, names(Diag)))
  }
  
  isNA <- length(which(Diag[, 1] == dimension))
  if (isNA == 0) {
    return(rep(0, length(tseq))) #in case there are no features with dimension "dimension"
  }
      
  Diag <- Diag[which(Diag[,1] == dimension), , drop = FALSE]
  
  Npoints <- nrow(Diag)

  fab <- matrix(NA, nrow = length(tseq), ncol = Npoints)
  lambda <- numeric()
  for (j in seq_len(Npoints)) {    
    fab[, j]  <- sapply(seq(along = tseq), FUN = function(i) {
        max(min(tseq[i] - Diag[j, 2], Diag[j, 3] - tseq[i]), 0)
      })
  }
  lambda <- sapply(seq(along = tseq), FUN = function(i) {
      sort(fab[i, ], decreasing = TRUE)[KK]
    })
  lambda[is.na(lambda)] <- 0
  if (length(KK) == 1) {
    lambda <- matrix(lambda)
  } else {
    lambda <- t(lambda)
  }
  return(lambda)
}