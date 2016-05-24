isCbalanced <- function(d, preperiod = FALSE) {
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
  
  balance <- FALSE
  
  V <- matrix(0, nrow = trt, ncol = trt)
  
  if (preperiod) {
    V <- diag(k)[c(k, 1:(k - 1)), ]
  } else {
    V <- rbind(numeric(k), diag(k)[1:(k - 1), ])
  }
  
  M <- t(Td(d)) %*% kronecker(diag(b), t(V)) %*% Td(d)
  ## M=[m_ij] where i is left neighbour to j
  
  cat("\n")
  
  if (!any(as.logical(diag(M))) && length(unique(as.vector(M - 2 * b * k * diag(trt)))) ==  2) {
    balance <- TRUE
    cat("The design is (first order) carry-over balanced.", "\n", "\n")
  } else {
    cat("The design is not (first order) carry-over balanced.", "\n", "\n")
  }
  ## If M is completely symmetric and the diagonal of M is zero then d is
  ## neighbourbalanced. Hence, if the diagonal is zero and
  ## length(unique(as.vector(M-2*b*k*diag(trt)))) [negative diagonal!] is two then
  ## d is carry-over balanced.
  
  cat("Left neighbour incidence matrix M_ij (i is left neighbour of j)", "\n", "\n")
  print(M)
  invisible(list(balance, M))
} 
