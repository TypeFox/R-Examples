aout.nbinom <-
function(data, param, alpha = 0.1, hide.outliers = FALSE){
  # check arguments
  if (!is.numeric(param) | !is.vector(param) | !identical(all.equal(length(param), 2), TRUE)) 
    stop("param must be a numeric vector of length 2.")
  if (!is.numeric(data) | !is.vector(data)) 
    stop("data must be a numeric vector.")
  if (!identical(all.equal(length(alpha), 1), TRUE) | alpha <= 0 | alpha >= 1) 
    stop("alpha must be a real number between 0 and 1, but it is ", alpha, ".")
  # end check arguments
  # determine the outlier region
  size <- param[1]
  prob <- param[2]
  
  tu <- qnbinom(alpha, size, prob)
  to <- qnbinom(1-alpha, size, prob)
  
  P.now <- sum(dnbinom((tu+1):(to-1), size, prob))
  
  while(P.now <= 1-alpha){ # uses the ascending-descending character of the p.d.f.
    if (dnbinom(tu, size, prob) > dnbinom(to, size, prob)) {
      next.P <- dnbinom(tu, size, prob)
      tu <- tu - 1
    }
    else {
      next.P <- dnbinom(to, size, prob)
      to <- to + 1
    }     
    P.now <- P.now + next.P
  }
  temp.region <- (tu+1):(to-1) # is the inlier region
  # give the results of the analysis
  temp <-  data.frame(data = data, is.outlier = !(data %in% temp.region))
  if (identical(all.equal(hide.outliers, FALSE), TRUE)) temp
  else temp[temp[,2] == FALSE, 1]
}
