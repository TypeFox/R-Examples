# Compute control-variate estimate ,that's, equation (7) in the reference paper.
# the initial sample size m1 is adjusted from 40 to 60.
y_control_estimator <- function (mk_correlated_standard_norm, r_mat,
                                invcdfnames, paramslists) {
  ndim <- ncol(mk_correlated_standard_norm)
  mk <- nrow(mk_correlated_standard_norm )
  if (class(mk_correlated_standard_norm) != "matrix") {
    stop("mk_correlated_standard_norm must be matrix !")
  }
  if (length(invcdfnames) != length(paramslists)) {
    stop("inversecdfs should have the same length paramslists !")
  }
  ndim <- ncol(mk_correlated_standard_norm)
  if (mk < 60) {
    stop("the number of observations must be not less than 60!")
  }
  transform_mat <- NULL
  for (i in 1:ndim) {
    funcall <- as.call(c(as.name(invcdfnames[i]),
                         list(pnorm(mk_correlated_standard_norm)[ ,i]), paramslists[[i]]))
    transform_mat <- cbind(transform_mat, eval(funcall))
  }
  res <- matrix(rep(0, ndim * ndim), nrow = ndim)
  diag(res) <- 1
  for (i in 1:(ndim-1))
    for (j in (i+1):ndim) {
      if (length(which(!duplicated(transform_mat[ ,i])[-1]))==0 ||
           length(which(!duplicated(transform_mat[ ,j])[-1]))==0)
      tmp_crude <- 0
      else
      tmp_crude <- cor(transform_mat[ ,i], transform_mat[ ,j])
      theta <- cor(mk_correlated_standard_norm[ ,i],
                                   mk_correlated_standard_norm[ ,j])
     r <- r_mat[i,j]
     # micro/macro replications, 20 macro replications
     micro_size <- floor(mk / 20)
     lowbound <- 0:19*micro_size+1
     upbound <- 1:20*micro_size
     tmp <- Map(function(x, y) {
          theta_micro <- cor(mk_correlated_standard_norm[ ,i][x:y],
                             mk_correlated_standard_norm[ ,j][x:y])
           if(length(which(!duplicated(transform_mat[ ,i][x:y])[-1]))==0 ||
              length(which(!duplicated(transform_mat[ ,j][x:y])[-1]))==0  )
           y_micro_crude_estimator <- 0
           else
           y_micro_crude_estimator <- cor(transform_mat[ ,i][x:y],
                                         transform_mat[ ,j][x:y])
           return(list(theta_micro =  theta_micro,
                  y_micro_crude_estimator = y_micro_crude_estimator))

           }, lowbound, upbound)
     tmp <- unlist(tmp)
     theta_macro <- tmp[names(tmp)=="theta_micro"]
     y_macro_crude_estimator <- tmp[names(tmp)=="y_micro_crude_estimator"]
     cov_estimator <- cov(y_macro_crude_estimator, theta_macro)
     var_estimator <- var(theta_macro)
     beta <- cov_estimator / var_estimator
     res[j,i] <- res[i,j] <- tmp_crude - beta * (theta - r)
    }
   res
}
