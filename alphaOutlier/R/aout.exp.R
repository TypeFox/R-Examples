aout.exp <-
function(data, param, alpha = 0.1, hide.outliers = FALSE, theta = 0){
  # check arguments
  if (!is.numeric(param) | !is.vector(param) | !identical(all.equal(length(param), 1), TRUE)) 
    stop("param must be a numeric vector of length 1.")
  if (!is.numeric(theta) | !is.vector(theta) | !identical(all.equal(length(theta), 1), TRUE)) 
    stop("theta must be a numeric vector of length 1 (0 for central exponential distribution).")
  if (!is.numeric(data) | !is.vector(data)) 
    stop("data must be a numeric vector.")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  temp.region <- c(theta, theta - log(alpha)/param)
  # give the results of the analysis
  temp <- data.frame(data = data, is.outlier = (data < temp.region[1] | 
                                                  data > temp.region[2]))
  if (identical(all.equal(hide.outliers, FALSE), TRUE)) temp
  else temp[temp[,2] == FALSE, 1]
}
