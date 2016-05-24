aout.cg <- function(data, param, alpha = 0.1, hide.outliers = FALSE){
  # check arguments
  if (!is.list(param) | !identical(all.equal(length(param), 3), TRUE)) 
    stop("param must be a list of length 3.")
  if (!is.numeric(data) | !is.matrix(data) | !identical(all.equal(ncol(data), 2), TRUE)) 
    stop("data must be a numeric matrix with two columns.")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  if (!identical(all.equal(length(param[[1]]), length(param[[2]])), TRUE) | 
        !identical(all.equal(length(param[[2]]), length(param[[3]])), TRUE)) 
    stop("The vectors contained by param must have the same length.")
  if (!identical(all.equal(sum(param[[1]]), 1), TRUE)) 
    stop("The sum of the probabilities in the first element of param must be 1.")
  
  # end check arguments
  # determine the outlier region
  # helper functions:
  uppbound <- function(K.alpha){
    x <- -2 * sigma^2 * log(K.alpha * sqrt(2 * pi) * sigma / p)
    sqrt(ifelse(x >= 0, x, 0)) + mu
  }
  lowbound <- function(K.alpha){
    x <- -2 * sigma^2 * log(K.alpha * sqrt(2 * pi) * sigma / p)
    -sqrt(ifelse(x >= 0, x, 0)) + mu
  }
  findroot <- function(K.alpha){
    1 - t(pnorm(uppbound(K.alpha), mu, sigma) 
          - pnorm(lowbound(K.alpha), mu, sigma)) %*% p - alpha
  }
  
  p <- param[[1]]
  mu <- param[[2]]
  sigma <- param[[3]]
  K <- uniroot(findroot, interval=c(0, 0.3989423/min(sigma)), tol = 1e-10)$root 
  # K is between 0 and the maximum of the Gaussian curves by definition
  
  temp.region <- cbind(lowbound(K), uppbound(K))
  is.out <- logical(nrow(data))
  for (i in 1:nrow(data)){
    is.out[i] <- (data[i,2] < temp.region[data[i,1], 1] | 
                    data[i,2] > temp.region[data[i,1], 2])
  }
  # give the results of the analysis
  temp <- data.frame(data = data, is.outlier = is.out)
  if (identical(all.equal(hide.outliers, FALSE), TRUE)) temp
  else temp[temp[,3] == FALSE, -3]
}