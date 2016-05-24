##' Significance test for bootstrapped coefficients
##'
##' For internal use only.
##' @param bc the parameter matrix as returned from correlation of the
##' bootstapped samples (n * 1000 matrix, with n = number of
##' parameters)
##' @param ci confidence level, one of c(0.01, 0.05, 0.1)
##' @param c0 for weibull: correlation vector for original model
##' @param method method for significance test
##' @return a logical vector indicating significance
##' @keywords internal
ptest <- function(bc, ci, c0 = NULL, method) {
  limits <- switch(as.character(ci),
                   "0.01" = c(5, 995),
                   "0.05" = c(25, 975),
                   "0.1" = c(50, 950))
  m <- dim(as.matrix(bc))[1]
  if (method == "weibull") {
    ## estimating p
    ## sorting values
    s <- t(apply(bc, 1, sort))
    ## matrix of ranks
    j <- matrix(rep(1:1000, m), nrow = 1000)

    ps <- numeric(m)

    for (i in 1:m) {
      if (c0[i] <= min(bc[i,])) {
        j1 <- 1
      } else {
        if (c0[i] >= max(bc[i,])) {
          j1 <- 1000
        } else {
          ## interpolate rank if within range of monte-carlo simulations
          j1 <- approx(s[i,], j[,i], c0[i])$y
        }
      }
      ## Weibull formula for n = 1000
      ps[i] <- j1 / 1001
    }
    is_sig <- ifelse(ps < limits[1]/1000 | ps > limits[2]/1000, TRUE,
                     FALSE)
    ci_lower <- NA
    ci_upper <- NA
    coefs <- c0
  } 
  if (method == "range") {
    coefs <- apply(bc, 1, median)
    ci_lower <- apply(bc, 1, function(x) sort(x)[limits[1]])
    ci_upper <- apply(bc, 1, function(x) sort(x)[limits[2]])
    is_sig <- logical(m)
    for (i in 1:m) {
      if (sign(ci_upper[i]) != sign(ci_lower[i])) {
        is_sig[i] <- FALSE
      } else {
        if (abs(coefs[i]) > abs((abs(ci_upper[i]) - abs(ci_lower[i]))/2)) {
          is_sig[i] <- TRUE
        } else {
          is_sig[i] <- FALSE
        }
      }
    }
  }
  if (method == "none") {
    coefs <- bc
    ci_lower <- ci_upper <- is_sig <- rep(NA, m)
  }
  data.frame(coef = coefs,
             significant = is_sig,
             ci_lower = ci_lower,
             ci_upper = ci_upper)
}
