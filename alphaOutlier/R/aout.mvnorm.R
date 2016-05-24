aout.mvnorm <-
function(data, param, alpha = 0.1, hide.outliers = FALSE){
  # check arguments
  if (!identical(all.equal(length(param), 2), TRUE)) 
    stop("param must be a list with two elements.")
  if (!is.vector(param[[1]]))
    stop("First entry of param must be the mean vector.")
  if (!is.matrix(param[[2]]))
    stop("Second entry of param must be the covariance matrix.")
  if (is.data.frame(data)) as.matrix(data)
  if (!is.matrix(data)) 
    stop("data must be a data.frame or matrix.")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  n <- nrow(data)
  d <- ncol(data)
  te1 <- numeric(n)
  for (i in 1:n){
    te1[i] <- (t(data[i,] - param[[1]])) %*% solve(param[[2]]) %*% (data[i,] - param[[1]])
  }
  temp.region <- qchisq(df = d, 1 - alpha)
  # give the results of the analysis
  temp <- data.frame(data = data, is.outlier = (te1 > temp.region))
  if (identical(all.equal(hide.outliers, FALSE), TRUE)) temp
  else temp[temp[,d+1] == FALSE, -(d+1)]
}
