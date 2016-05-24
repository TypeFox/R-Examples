design.row <- function(d) {
  
  if (!is.matrix(d) || !is.numeric(d)) {
    stop("Please check your design matrix")
  }
  trt <- max(d)
  if (any(sort(unique(as.vector(d))) != 1:length(unique(as.vector(d))))) {
    stop("Please check your design matrix")
  }
  b <- nrow(d)
  k <- ncol(d)
  if (any(c(trt, b, k) == 1)) {
    stop("Please check your design matrix")
  }
  a1 <- floor(k/trt) + 1
  
  td <- Td(d)  # Treatment design matrix
  bd <- kronecker(diag(b), rep(1, k))  # Block design matrix
  occ <- diag(t(td) %*% td)  # Number of occurrences of treatments in design d
  rinc <- t(td) %*% bd  # Incidence matrix
  concur <- rinc %*% t(rinc)  # Concurrence matrix
  
  # Check type of design
  type <- rep(FALSE, 6)
  
  type[1] <- all(occ == (b * k/trt))  # TRUE if all treatments occur equally often in d
  type[2] <- (min(rinc) > a1 - 1.1) && (max(rinc) < a1 + 0.1)
  # TRUE if the design is binary w.r.t rows in the sense that each treatment
  # occurs a1 or a1-1 times in each row
  type[3] <- (length(unique(concur[upper.tri(concur)])) == 1)
  # TRUE if d is pairwise balanced
  type[4] <- (k < trt)  # TRUE, if d is incomplete w.r.t. rows
  type[5] <- (k == trt)  # TRUE, if d is complete w.r.t. rows
  type[6] <- !any(as.logical(as.vector(rinc) - rinc[1, 1]))  # TRUE, if d is uniform on the rows
  
  # Output
  
  names(occ) <- 1:trt
  rownames(rinc) <- 1:trt
  colnames(rinc) <- 1:b
  rownames(concur) <- 1:trt
  colnames(concur) <- 1:trt
  
  list(occ, rinc, concur, type)
} 
