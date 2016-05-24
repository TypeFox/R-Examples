aout.chisq <-
function(data, param, alpha = 0.1, hide.outliers = FALSE, ncp = 0, 
                       lower = auto.l, upper = auto.u, method.in = "Newton", 
                       global.in = "gline", 
                       control.in = list(sigma = 0.1, maxit = 1000,
                                         xtol = 1e-12, ftol = 1e-12, btol = 1e-4)){
  # check arguments
  if (!is.numeric(param) | !is.vector(param) | !identical(all.equal(length(param), 1), TRUE)) 
    stop("param must be a numeric vector of length 1.")
  if (!is.numeric(data) | !is.vector(data)) 
    stop("data must be a numeric vector.")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  dfr <- param
  a <- alpha
  auto.l <- qchisq(a/2, dfr, ncp)
  auto.u <- qchisq(1 - a/2, dfr, ncp)
  
  # fn2 defines the set of equations for which the roots are searched
  fn2 <- function(x, freedom, noncen, alpha){
    y <- numeric(2)
    y[1] <- 1 - alpha - pchisq(x[2], freedom, noncen) + pchisq(x[1], freedom, noncen)
    y[2] <- dchisq(x[1], freedom, noncen) - dchisq(x[2], freedom, noncen)
    y
  }
  
  if (dfr > 2 | (dfr == 2 & ncp > 2)){
    # if the following condition holds, the density is increasing-decreasing
    temp1 <- nleqslv(c(lower, upper), fn = fn2, alpha = a, freedom = dfr, noncen = ncp, 
                     method = method.in, global = global.in, control = control.in)
    if (temp1$termcd != 6) temp.region <- temp1$x
    # if the algorithm produces a regular jacobian, run the algorithm
    else { 
      # if the algorithm produces a singular jacobian...
      if (identical(all.equal(dchisq(0, dfr, ncp), 0), TRUE))
        # if this cond. holds, we get the singular jacobian because of bad initial values
        temp.region <- nleqslv(c(0.1, upper), fn = fn2, alpha = a, freedom = dfr, 
                               noncen = ncp, method = method.in, global = global.in, 
                               control = control.in)$x
      else temp.region <- c(0, qchisq(1 - a, dfr, ncp))
      # if lim_{x->0} f(x) != 0, take the upper (1-a) tail region as outlier region 
    }
  }
  else temp.region <- c(0, qchisq(1 - a, dfr, ncp))
  # this is the inlier region for strictly decreasing densities
  
  # give the results of the analysis
  temp <- data.frame(data = data, is.outlier = (data < temp.region[1] | 
                                                  data > temp.region[2]))
  if (hide.outliers == FALSE) temp
  else temp[temp[,2] == FALSE, 1]
}
