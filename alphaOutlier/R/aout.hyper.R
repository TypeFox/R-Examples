aout.hyper <-
function(data, param, alpha = 0.1, hide.outliers = FALSE){
  # check arguments
  if (!is.numeric(param) | !is.vector(param) | !identical(all.equal(length(param), 3), TRUE)) 
    stop("param must be a numeric vector of length 3.")
  if (!is.numeric(data) | !is.vector(data)) 
    stop("data must be a numeric vector.")
  if (max(data) > param[3]) 
    stop("Impossible value for k, must be larger than max(data).")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  m <- param[1]
  n <- param[2]
  k <- param[3]
  x <- 0:k
  prob.vector <- dhyper(x, m, n, k)
  temp.region <- order(prob.vector)[which(cumsum(sort(prob.vector)) < alpha)] - 1
  # give the results of the analysis
  temp <-  data.frame(data = data, is.outlier = (data %in% temp.region))
  if (identical(all.equal(hide.outliers, FALSE), TRUE)) temp
  else temp[temp[,2] == FALSE, 1]
}
