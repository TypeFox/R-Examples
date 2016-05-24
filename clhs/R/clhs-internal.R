## Objective function for cLHS
##
## (c) Pierre Roudier, Landcare Research, 2011.
##

.lhs_obj <- function(
  size,
  data_continuous_sampled,
  data_factor_sampled,
  continuous_strata,
  weights = list(numeric = 1, factor = 1, corelation = 1),
  cor_mat,
  factor_obj
  ) {

  # Continuous variables
  #
  n_cont_variables <- ncol(data_continuous_sampled)

  cont_data_strata <- lapply(1:n_cont_variables, function(i) list(data_continuous_sampled[, i], continuous_strata[, i]) )
  cont_obj_sampled <- lapply(cont_data_strata, function(x) hist(x[[1]], breaks = x[[2]], plot = FALSE)$counts)
  cont_obj_sampled <- matrix(unlist(cont_obj_sampled), ncol = n_cont_variables, byrow = FALSE)

  delta_obj_continuous <- rowSums(abs(cont_obj_sampled - 1))

  # Factor variables
  #
  n_factor_variables <- ncol(data_factor_sampled)
  if (n_factor_variables > 0) {
    factor_obj_sampled <- lapply(1:n_factor_variables, function(x) table(data_factor_sampled[, x])/nrow(data_factor_sampled))
    delta_obj_factor <- lapply(1:n_factor_variables, function(x) sum(abs(factor_obj_sampled[[x]] - factor_obj[[x]])))

    delta_obj_factor <- unlist(delta_obj_factor)#/length(delta_obj_factor) # do we need to ponder w/ the number of factors?
  }
  else
    delta_obj_factor <- 0

  # Correlation of continuous data
  #
  cor_sampled <- cor(data_continuous_sampled)
  cor_sampled[is.na(cor_sampled)] <- 1 # when there's only one observation, cor() throws NAs - we set these to 1

  delta_obj_cor <- sum(abs(cor_mat - cor_sampled))

  # Objective function
  #
  obj <- weights[[1]]*sum(delta_obj_continuous) + weights[[2]]*sum(delta_obj_factor) + weights[[3]]*delta_obj_cor

  # Returning results
  #
  list(obj = obj, delta_obj_continuous = delta_obj_continuous, delta_obj_factor = delta_obj_factor, delta_obj_cor = delta_obj_cor)
}
