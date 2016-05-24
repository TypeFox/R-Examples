clusterprob <- function(fit, newdata, ...) {
  if (!is(fit, "disclapmixfit")) stop("object must be a disclapmixfit")
  
  if (!is.matrix(newdata) || !is.integer(newdata)) {
    stop("newdata is not an integer matrix")
  }
  
  if (ncol(newdata) != ncol(fit$y)) {
    stop("newdata must have the same number of columns as the original model was fitted for")
  }
  
  wic <- rcpp_calculate_wic(newdata, fit$y, fit$disclap_parameters, fit$tau)      
  vic_matrix <- rcpp_calculate_vic(wic)

  return(vic_matrix)
}

