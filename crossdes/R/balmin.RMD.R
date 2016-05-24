balmin.RMD <- function(trt, n, p) {
  
  # check if input parameters make sense
  ifelse(any(as.logical(c(trt, n, p)%%1)), stop("Please check your design parameters."), 
    ifelse(p >= trt, stop("Please check your design parameters."), lambda <- n/trt))
  if ((lambda%%1)) {
    stop("There is no balanced minimal RMD for these parameters.")
  }
  if (lambda * (p - 1) != (trt - 1)) {
    stop("There is no balanced minimal RMD for these parameters.")
  }
  
  # construct c-vector that gives rise to a difference set
  vec <- numeric(trt)
  if (!(trt%%2)) {
    for (i in 1:(trt/2)) {
      vec[(2 * i - 1):(2 * i)] <- c(i, trt + 1 - i)
    }
  } else {
    for (i in 1:(trt + 1)%/%4) {
      vec[(2 * i - 1):(2 * i)] <- c(2 * i - 1, trt + 2 - 2 * i)
    }
  }
  
  dummy <- 2 * ((trt + 1)%/%4)
  if (!((trt - 1)%%4)) {
    vec <- c(vec[1:dummy], (trt + 1)/2, rev(vec[1:dummy]))
  }
  if (!((trt + 1)%%4)) {
    vec <- c(vec[1:dummy], rev(vec[1:(dummy - 1)]))
  }
  
  # construct design matrix (rows represent periods, columns represent subjects)
  des <- matrix(0, p, n)
  for (i in 1:lambda) {
    for (j in 1:trt) {
      des[, trt * (i - 1) + j] <- (vec[((p - 1) * (i - 1) + 1):((p - 1) * i + 1)] + j - 1)%%trt
    }
  }
  des[!des] <- trt  # rename treatment number 0 to treatment number trt.
  t(des)  # final design with subjects as rows
} 
