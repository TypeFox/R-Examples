#Check to make sure the magnitude of the estimator must be no longer than |r|.
#Described in Appendeix:NORTA RVG Methods step 3.3 according to the reference paper.
check_y_estimator<- function (y_estimator, r_matrix, msave, mk,
                              invcdfnames, paramslists) {
  r_matrix[lower.tri(r_matrix)] <- 0
  upper_elements <- r_matrix[upper.tri(r_matrix)]
  upper_elements[upper_elements==1] <- 0.999
  upper_elements[upper_elements==-1] <- -0.999
  r_matrix[upper.tri(r_matrix)] <- upper_elements
  ndim <- ncol(r_matrix)
  y_estimator[lower.tri(y_estimator)] <- 0
  diag(y_estimator) <- 0
  indexmat <- which(abs(y_estimator) > abs(r_matrix), arr.ind = TRUE)
  if (length(indexmat) == 0) {
      y_estimator <- y_estimator + t(y_estimator)
      diag(y_estimator) <- 1
      return(list(y_estimator = y_estimator, mk = mk))
  }
  t <- 1
  while (any(abs(y_estimator) > abs(r_matrix))) {
  t <- t + 1
  w_bar <- matrix(rnorm(msave*ndim, mean = 0, sd = 1), nrow = msave)
  y_msave_estimator <- Perform_chol(r_matrix , w_bar, invcdfnames, paramslists)
  y_new_estimator <- ((t - 1) / t) * y_estimator + 1 / t * y_msave_estimator
  y_estimator <- y_new_estimator
  y_estimator[lower.tri(y_estimator)] <- 0
  diag(y_estimator) <- 0
  mk <- mk + msave
  }
  y_estimator <-  y_estimator + t( y_estimator)
  diag( y_estimator) <- 1
  return(list(y_estimator = y_estimator, mk = mk))
}
