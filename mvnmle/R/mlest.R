mlest<-function(data,...)
  {
    # Takes MVN data with missing values and calculates the MLE of the mean vector and the var-cov matrix

    data<-as.matrix(data)
    sortlist<-mysort(data) # put data with identical patterns of missingness together
    
    nvars<-ncol(data)
    nobs<-nrow(data)
    if(nvars>50)
      stop("mlest cannot handle more than 50 variables.")

    startvals<-getstartvals(data) # find starting values

    lf<-getclf(data=sortlist$sorted.data, freq=sortlist$freq)
    mle<-nlm(lf,startvals,...)

    muhat<-mle$estimate[1:nvars] # extract estimates of mean
    del<-make.del(mle$estimate[-(1:nvars)]) # extract estimates of sigmahat
    factor<-solve(del,diag(nvars))
    sigmahat<-t(factor) %*% factor
    list(muhat=muhat, sigmahat=sigmahat, value=mle$minimum, gradient=mle$gradient,
         stop.code=mle$code, iterations=mle$iterations)
  }

