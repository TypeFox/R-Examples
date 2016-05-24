aout.gandh <-
function(data, param, alpha = 0.1, hide.outliers = FALSE){
  # check arguments
  if (!is.numeric(param) | !is.vector(param) | !identical(all.equal(length(param), 4), TRUE)) 
    stop("param must be a numeric vector of length 4.")
  if (!is.numeric(data) | !is.vector(data)) 
    stop("data must be a numeric vector.")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  median2 <- param[1]
  scale2 <- param[2]
  g2 <- param[3]
  h2 <- param[4]
  # quantile function of g and h distribution
  qgandh <- function(p, median2, scale2, g2, h2){
    temp <- qnorm(p)
    if (g2 != 0) median2 + scale2 * (exp(g2 * temp) - 1) * exp(temp^2 * h2 / 2) / g2
    else median2 + scale2 * temp * exp(temp^2 * h2 / 2)
  }
  # objective function in minimization 
  qgandhmin <- function(vec, median2, scale2, g2, h2) 
    qgandh(vec[2], median2, scale2, g2, h2) - qgandh(vec[1], median2, scale2, g2, h2)
  # inequality constraint in minimization
  Fgreater1ma <- function(vec, median2, scale2, g2, h2) vec[2] - vec[1]
  
  res.opt <- solnp(pars = c(alpha/2, 1-alpha/2), fun = qgandhmin, ineqfun = Fgreater1ma, 
                   ineqLB = 1-alpha, ineqUB = 1, UB = c(alpha, 1), LB = c(0, 1-alpha), 
                   median2 = median2, scale2 = scale2, g2 = g2, h2 = h2)
  temp.region <- qgandh(res.opt$pars, median2, scale2, g2, h2)
  # give the results of the analysis
  temp <- data.frame(data = data, is.outlier = (data < temp.region[1] | 
                                                  data > temp.region[2]))
  if (identical(all.equal(hide.outliers, FALSE), TRUE)) temp
  else temp[temp[,2] == FALSE, 1]
}
