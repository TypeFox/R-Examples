
# Author: P. Poncet

tsybakov <-
function(x,                        # sample (the data)
         bw = NULL,                # bandwidth
         a,                        # parameter 'a'
         alpha = 0.9,              # parameter 'alpha' (if 'a' is missing, a <- (1:length('x'))^(-alpha))
         kernel = "triangular",    # kernel used
         dmp = TRUE,               # Should Djeddour et al. method be used?
         par = shorth(x))          # initial value of the algorithm
{
#####################################
# Tsybakov's recursive mode estimator
#####################################

  if (pmatch(tolower(kernel), "normal", nomatch = 0)) {
    kernel <- "gaussian"
  } else kernel <- match.arg(tolower(kernel), .kernelList) # '.kernelList' is defined in 'K.R'
  
  K <- paste(".kernel.d", kernel, sep = "") # The derivative of the kernel 'kernel' is used
  
  nx <- length(x)
    
  if (missing(a)) {
    a <- (1:nx)^(-alpha)
  }
  
  if (is.null(bw)) bw <- (1:nx)^(-1/7)
    
  ## Initialization
  M <- par
  M.dmp <- par

  b <- a/(bw^2)
  p <- bw^3
  p <- p/cumsum(p)
  
  for (n in 1:nx) {
    M <- M + b[n]*do.call(K, list((M-x[n])/bw[n]))$k #eval(call(K, (M-x[n])/bw[n]))$k
    M.dmp <- M.dmp + p[n]*(M - M.dmp)
  }
      
  ## Output
  return(ifelse(dmp, M.dmp, M))
}

#Tsybakov <- tsybakov
