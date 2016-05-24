# @title Impute missing multinomial values via MLE
# @description Given an observation of multinomial values, \eqn{X= (X_1, ..., X_p)}, 
# impute the missing values via the MLE estimates of \eqn{\theta= (\theta_1, ..., \theta_p)}
# where \eqn{X_j \sim M(1, \theta_j)}.
# @param miss_val A vector / row of multinomial variables containing missing values.
# @param row_ind A vector of integers corresponding to the row indices of \code{MLEx_y}
# which match the missing values for \code{miss_val}.
# @param MLEx_y A \code{data.frame} returned from \code{\link{multinomial_em}}
# @return A complete observation \eqn{X} (ie-without missing values).
impute_multinomial <- function(miss_val, row_ind, MLEx_y) {
  ml_vals <- MLEx_y[row_ind,]
  ml_vals <- ml_vals[ml_vals$theta_y == max(ml_vals$theta_y),]
  
  if (nrow(ml_vals) == 1) {
    miss_val[which(is.na(miss_val))] <- ml_vals[which(is.na(miss_val))]
  } else {
    imp <- ml_vals[sample(1:nrow(ml_vals), size= 1),]
    miss_val[which(is.na(miss_val))] <- imp[which(is.na(miss_val))]
  }
  return(miss_val)
}

# this provides a wrapper to \code{\link{impute_multinomial}} above such 
# that all missing values are imputed.
impute_multinomial_all <- function(dat_miss, MLEx_y, p) {
  gc()
  marg_ind <- marg_complete_compare(dat_miss, MLEx_y[, 1:p], marg_to_complete= TRUE)
  
  for (i in 1:nrow(dat_miss)) {
    dat_miss[i,] <- impute_multinomial(dat_miss[i,], row_ind= marg_ind[[i]], MLEx_y= MLEx_y)
  }
  
  return(dat_miss)
}

