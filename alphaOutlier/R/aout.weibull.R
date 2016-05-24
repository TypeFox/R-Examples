aout.weibull <-
function(data, param, alpha = 0.1, hide.outliers = FALSE,
                         lower = auto.l, upper = auto.u, method.in = "Broyden",
                         global.in = "qline", 
                         control.in = list(sigma = 0.1, maxit = 1000, xtol = 1e-12,
                                           ftol = 1e-12, btol = 1e-4)){
  # check arguments
  if (!is.numeric(param) | !is.vector(param) | !identical(all.equal(length(param), 2), TRUE)) 
    stop("param must be a numeric vector of length 2.")
  if (!is.numeric(data) | !is.vector(data)) 
    stop("data must be a numeric vector.")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  bet <- param[1]
  lambda <- param[2]
  
  
  auto.l <- qweibull(alpha/2, bet, lambda)
  auto.u <- qweibull(1 - alpha/2, bet, lambda)
  
  fn3 <- function(x, betta, lammda, a) {
    # defines the set of equations for which the roots shall be found
    y <- numeric()
    y[1] <- 1 - a - pweibull(x[2], betta, lammda) + pweibull(x[1], betta, lammda)
    y[2] <- dweibull(x[1], betta, lammda) - dweibull(x[2], betta, lammda)
    y
  }
  
  if (bet > 1) # then f is increasing-decreasing density
  {
    temp.sol <- nleqslv(c(auto.l, auto.u), fn = fn3, a = alpha, betta = bet, 
                        lammda = lambda, method = method.in, global = global.in, 
                        control = control.in)
    # if the jacobian of temp.sol is non-singular, use this solution
    if (!identical(all.equal(temp.sol$termcd, 6), TRUE)) temp.region <- temp.sol$x
    else { # if the jacobian is singular
      warning("Lower bound of inlier region is only approximately zero")
      temp.region <- c(0, qweibull(1 - alpha, bet, lambda))
    }
  }
  else temp.region <- c(0, qweibull(1 - alpha, bet, lambda))
  # give the results of the analysis
  temp <- data.frame(data = data, is.outlier = (data < temp.region[1] | 
                                                  data > temp.region[2]))
  if (identical(all.equal(hide.outliers, FALSE), TRUE)) temp
  else temp[temp[,2] == FALSE, 1]
}
