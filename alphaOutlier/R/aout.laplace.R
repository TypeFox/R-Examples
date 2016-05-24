aout.laplace <-
function(data, param, alpha = 0.1, hide.outliers = FALSE){
  # check arguments
  if (!is.numeric(param) | !is.vector(param) | !identical(all.equal(length(param), 2), TRUE)) 
    stop("param must be a numeric vector of length 2.")
  if (param[2] <= 0)
    stop("sigma must be positive.")
  if (!is.numeric(data) | !is.vector(data)) 
    stop("data must be a numeric vector.")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  mu <- param[1]
  sigma <- param[2]
  temp.region <- c(mu + sigma * log(alpha), mu - sigma*log(alpha))
  # give the results of the analysis
  temp <- data.frame(data = data, is.outlier = (data < temp.region[1] | 
                                                  data > temp.region[2]))
  if (identical(all.equal(hide.outliers, FALSE), TRUE)) temp
  else temp[temp[,2] == FALSE, 1]
}
