aout.binom <-
function(data, param, alpha = 0.1, hide.outliers = FALSE){
  # check arguments
  if (!is.numeric(param) | !is.vector(param) | !identical(all.equal(length(param), 2), TRUE)) 
    stop("param must be a numeric vector of length 2.")
  if (!identical(all.equal(param[1], round(param[1])), TRUE))
    stop("First element of param must be the sample size (an integer).")
  if (param[2] <= 0 | param[2] >= 1)
    stop("Second element of param must be the probability of success (between 0 and 1).")
  if (!is.numeric(data) | !is.vector(data)) 
    stop("data must be a numeric vector.")
  if (any(data > param[1]))
    stop("No element of data may be larger than ", param[1], ".")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  size <- param[1]
  prob <- param[2]
  x <- 0:size
  prob.vector <- dbinom(x, size, prob)
  temp.region <- order(prob.vector)[which(cumsum(sort(prob.vector)) < alpha)] - 1
  # give the results of the analysis
  temp <- data.frame(data = data, is.outlier = (data %in% temp.region))
  if (hide.outliers == FALSE) temp
  else temp[temp[,2] == FALSE, 1]
}
