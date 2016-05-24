Smean <-
function(y, N=Inf, level=0.95){
  
  # y = data (vector)
  # N = population size (default N=Inf without finite population correction)
  # level = coverage probability for confidence intervals N(default =0.95)
  
    n <- length(y)
    # input treatment
      if(level<0 | level>1) stop("Wrong input: ", sQuote("level")," has to be probability between 0 and 1.")
      if(n<=0) stop("Wrong input: ", sQuote("length(y)")," has to be positive integer.")
      if(!is.numeric(N)) stop("Wrong input: ", sQuote("N")," is not a number or ", sQuote("Inf"))
      if(n > N) stop("Wrong input: ", sQuote("length(y)"), " has to be smaller than population ", sQuote("N"))
      if(any(is.na(y))) y <- na.omit(y)
    q <- qnorm((1+level)/2)
    
    # for N large enough (without finite population correction)
    if(N==Inf) { 
      mean.est = mean(y)                                     
      var.est = (1/(n*(n-1))) * sum((y - mean.est)^2)
  #    test = var(y)
      lowerLimit =  mean.est -  q*sqrt(var.est)
      upperLimit =  mean.est +  q*sqrt(var.est)
    }  
    # for specified population size N (with finite population correction)
    else  {
      mean.est = mean(y)
      var.est = (N-n)/N * (1/(n*(n-1))) * sum((y - mean.est)^2)
      lowerLimit <- mean.est -  q*sqrt(var.est)
      upperLimit <-  mean.est +  q*sqrt(var.est)
    }
    # calculation of se
    if(var.est>0) se <- sqrt(var.est)
    else{
      se <- NA
      warning("Standard error is ", sQuote("NA"),", because calculations for variance of mean has been not positive. Confidence intervals may not be valid.")
    }  
    # return argument
    ret <- list()
    ret$call <- list(y=y,n=n,N=N,level=level)
    ret$mean <- mean.est
    ret$se <- se
    ret$ci <- c(lowerLimit,upperLimit)
    structure(ret, class= "Smean")
  }
