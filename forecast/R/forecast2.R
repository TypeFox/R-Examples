# Mean forecast
meanf <- function(x,h=10,level=c(80,95),fan=FALSE, lambda=NULL, biasadj=FALSE)
{
  xname <- deparse(substitute(x))
  n <- length(x)
  #if(!is.ts(x))
	if(!is.null(lambda))
	{
		origx <- x
		x <- BoxCox(x,lambda)
	}
  meanx <- mean(x, na.rm=TRUE)
  fits <- rep(meanx,length(x))
  res <- x-fits
  f <- rep(meanx,h)
  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  nconf <- length(level)
  lower <- upper <- matrix(NA,nrow=h,ncol=nconf)
  s <- sd(x,na.rm=TRUE)
  for(i in 1:nconf)
  {
    if(n > 1)
      tfrac <- qt( 0.5 - level[i]/200, n-1)
    else
      tfrac <- -Inf
    w <- -tfrac * s*sqrt(1+1/n)
    lower[,i] <- f-w
    upper[,i] <- f+w
  }
  colnames(lower) <- colnames(upper) <- paste(level,"%",sep="")
  if(is.ts(x))
  {
    fits <- ts(fits)
    res <- ts(res)
    tsp(fits) <- tsp(res) <- tsp(x)
    freq <- frequency(x)
    f <- ts(f,start=tsp(x)[2]+1/freq,frequency=freq)
    lower <- ts(lower,start=tsp(x)[2]+1/freq,frequency=freq)
    upper <- ts(upper,start=tsp(x)[2]+1/freq,frequency=freq)
    #fits <- ts(rep(NA,n))
    #if(n > 1)
    #{
    #  for(i in 2:n)
    #    fits[i] <- mean(x[1:(i-1)],na.rm=TRUE)
    #}
    #res <- x - fits	
  }
  
	if(!is.null(lambda))
	{
		fits <- InvBoxCox(fits,lambda)
		x <- origx
		f <- InvBoxCox(f,lambda)
		if(biasadj){
		  f <- InvBoxCoxf(x = list(level = level, mean = f, upper = upper, lower = lower), lambda = lambda)
		}
		lower <- InvBoxCox(lower,lambda)
		upper <- InvBoxCox(upper,lambda)
	}	

  junk <- list(method="Mean",level=level,x=x,xname=xname,mean=f,lower=lower,upper=upper,
        model=list(mu=f[1],mu.se=s/sqrt(length(x)),sd=s), lambda=lambda, fitted=fits, residuals=res)
  junk$model$call <- match.call()

  return(structure(junk,class="forecast"))
}

thetaf <- function(x,h=10,level=c(80,95),fan=FALSE)
{
    if(fan)
        level <- seq(51,99,by=3)
    else
    {
        if(min(level) > 0 & max(level) < 1)
            level <- 100*level
        else if(min(level) < 0 | max(level) > 99.99)
            stop("Confidence limit out of range")
    }
    fcast <- ses(x,h=h)
    tmp2 <- lsfit(0:(length(x)-1),x)$coef[2]/2
    alpha <- fcast$model$par["alpha"]
    n <- length(x)
    fcast$mean <- fcast$mean + tmp2*(0:(h-1) + (1-(1-alpha)^n)/alpha)
    fcast.se <- sqrt(fcast$model$sigma) * sqrt((0:(h-1))*alpha^2+1)
    nconf <- length(level)
    fcast$lower <- fcast$upper <- matrix(NA,nrow=h,ncol=nconf)
    for(i in 1:nconf)
    {
        zt <- -qnorm( 0.5 - level[i]/200)
        fcast$lower[,i] <- fcast$mean - zt*fcast.se
        fcast$upper[,i] <- fcast$mean + zt*fcast.se
    }
    fcast$level <- level
    fcast$method <- "Theta"
    fcast$model <- list(alpha=alpha,drift=tmp2,sigma=fcast$model$sigma)
    fcast$model$call <- match.call()
    return(fcast)
}

# Random walk
rwf <- function(x,h=10,drift=FALSE,level=c(80,95),fan=FALSE,lambda=NULL,biasadj=FALSE)
{
  xname <- deparse(substitute(x))
  n <- length(x)
  freq <- frequency(x)
  nn <- 1:h
  if(!is.ts(x))
    x <- ts(x)
  if(!is.null(lambda))
  {
    origx <- x
    x <- BoxCox(x,lambda)
  }
  if(drift)
  {
    fit <- summary(lm(diff(x) ~ 1,na.action=na.exclude))
    b <- fit$coefficients[1,1]
    b.se <- fit$coefficients[1,2]
    s <- fit$sigma
#    res <- ts(c(NA,residuals(fit)))
    method <- "Random walk with drift"
  }
  else
  {
    b <- b.se <- 0
    s <- sd(diff(x),na.rm=TRUE)
#    fits <- ts(x[-n],start=tsp(x)[1]+1/freq,frequency=freq)
#    res <- ts(c(NA,diff(x)))
    method <- "Random walk"
  }
  fits <- ts(c(NA,x[-n]) + b, start=tsp(x)[1], frequency=freq)
  res <- x - fits
  #tsp(res) <- tsp(fits) <- tsp(x)
  f <- x[n] + nn*b
  se <- sqrt((nn*s^2) + (nn*b.se)^2)

  if(fan)
    level <- seq(51,99,by=3)
  else
  {
    if(min(level) > 0 & max(level) < 1)
      level <- 100*level
    else if(min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  nconf <- length(level)
  z <- qnorm(.5 + level/200)
  lower <- upper <- matrix(NA,nrow=h,ncol=nconf)
  for(i in 1:nconf)
  {
    lower[,i] <- f - z[i]*se
    upper[,i] <- f + z[i]*se
  }
  lower <- ts(lower,start=tsp(x)[2]+1/freq,frequency=freq)
  upper <- ts(upper,start=tsp(x)[2]+1/freq,frequency=freq)
  colnames(lower) <- colnames(upper) <- paste(level,"%",sep="")
  fcast <- ts(f,start=tsp(x)[2]+1/freq,frequency=freq)
  #fits <- x - res
  if(!is.null(lambda))
  {
    x <- origx
    fcast <- InvBoxCox(fcast,lambda)
    if(biasadj){
      fcast <- InvBoxCoxf(x = list(level = level, mean = fcast, upper = upper, lower = lower), lambda = lambda)
    }
    fcast <- InvBoxCox(fcast,lambda)
    fits <- InvBoxCox(fits,lambda)
    upper <- InvBoxCox(upper,lambda)
    lower <- InvBoxCox(lower,lambda)
  }

  junk <- list(method=method,level=level,x=x,xname=xname,mean=fcast,lower=lower,upper=upper,
      model=list(drift=b,drift.se=b.se,sd=s), fitted = fits, residuals = res, lambda=lambda)
  junk$model$call <- match.call()

  return(structure(junk,class="forecast"))
}

BoxCox <- function(x,lambda)
{
  if(lambda < 0)
	x[x < 0] <- NA
  if(lambda==0)
    out <- log(x)
  else
    out <- (sign(x)*abs(x)^lambda - 1)/lambda
  if(!is.null(colnames(x)))
	colnames(out) <- colnames(x)
  return(out)	
}

InvBoxCox <- function(x,lambda)
{
	if(lambda < 0)
		x[x > -1/lambda] <- NA
    if(lambda==0)
        out <- exp(x)
    else
    {
        xx <- x*lambda + 1
        out <- sign(xx)*abs(xx)^(1/lambda)
    }
    if(!is.null(colnames(x)))
	  colnames(out) <- colnames(x)
    return(out)	
}

InvBoxCoxf <- function(x=NULL, fvar=NULL, lambda=NULL){
  if(is.null(lambda)){
    stop("Must specify lambda using lambda=numeric(1)")
  }
  if(is.null(fvar)){
    level <- max(x$level)
    if(NCOL(x$upper)>1 & NCOL(x$lower)){
      i <- match(level,x$level)
      x$upper <- x$upper[,i]
      x$lower <- x$lower[,i]
    }
    if(level>1){
      level <- level/100      
    }
    level <- mean(c(level,1))
    #Note: Use BoxCox transformed upper and lower values
    fvar <- ((x$upper-x$lower)/stats::qnorm(level)/2)^2
  }
  else{
    x <- list(mean=x)
  }
  if("matrix"%in%class(fvar)){
    fvar <- diag(fvar)
  }
  return(x$mean * (1 + 0.5*fvar*(1-lambda)/(x$mean)^(2*lambda)))
}

forecast.StructTS <- function(object,h=ifelse(object$coef["epsilon"]>1e-10, 2*object$xtsp[3], 10),level=c(80,95),fan=FALSE,lambda=NULL,biasadj=FALSE,...)
{
    xname <- deparse(substitute(x))
    x <- object$data
    pred <- predict(object,n.ahead=h)
    if(fan)
        level <- seq(51,99,by=3)
    else
    {
        if(min(level) > 0 & max(level) < 1)
            level <- 100*level
        else if(min(level) < 0 | max(level) > 99.99)
            stop("Confidence limit out of range")
    }
    nint <- length(level)
    lower <- matrix(NA,ncol=nint,nrow=length(pred$pred))
    upper <- lower
    for(i in 1:nint)
    {
        qq <- qnorm(0.5*(1+level[i]/100))
        lower[,i] <- pred$pred - qq*pred$se
        upper[,i] <- pred$pred + qq*pred$se
    }
    colnames(lower) = colnames(upper) = paste(level,"%",sep="")
    if(is.element("seas",names(object$coef)))
        method <- "Basic structural model"
    else if(is.element("slope",names(object$coef)))
        method <- "Local linear structural model"
    else 
        method <- "Local level structural model"
	
	fits <- x - residuals(object)
	if(!is.null(lambda))
	{
		fits <- InvBoxCox(fits,lambda)
		x <- InvBoxCox(x,lambda)
		pred$pred <- InvBoxCox(pred$pred,lambda)
		if(biasadj){
		  pred$pred <- InvBoxCoxf(x = list(level = level, mean = pred$pred, upper = upper, lower = lower), lambda = lambda)
		}
		lower <- InvBoxCox(lower,lambda)
		upper <- InvBoxCox(upper,lambda)
	}
		
	
    return(structure(list(method=method,model=object,level=level,mean=pred$pred,lower=lower,upper=upper,
        x=x,xname=xname,fitted=fits,residuals=residuals(object)),
        class="forecast"))
}

forecast.HoltWinters <- function(object, h=ifelse(frequency(object$x)>1,2*frequency(object$x),10),
  level=c(80,95), fan=FALSE, lambda=NULL, biasadj=FALSE,...)
{
    xname <- deparse(substitute(x))
    x <- object$x
    if(!is.null(object$exponential))
      if(object$exponential)
        stop("Forecasting for exponential trend not yet implemented.")

    if(fan)
        level <- seq(51,99,by=3)
    else
    {
        if(min(level) > 0 & max(level) < 1)
            level <- 100*level
        else if(min(level) < 0 | max(level) > 99.99)
            stop("Confidence limit out of range")
    }
    nint <- length(level)

    pred <- predict(object,n.ahead=h,prediction.interval=TRUE,level=level[1]/100)
    pmean <- pred[,1]
    upper <- lower <- matrix(NA,ncol=nint,nrow=length(pred[,1]))
    se <- (pred[,2]-pred[,3])/(2*qnorm(0.5*(1+level[1]/100)))
    for(i in 1:nint)
    {
        qq <- qnorm(0.5*(1+level[i]/100))
        lower[,i] <- pmean - qq*se
        upper[,i] <- pmean + qq*se
    }
    colnames(lower) = colnames(upper) = paste(level,"%",sep="")
	
	
	if(!is.null(lambda))
	{
		object$fitted[,1] <- InvBoxCox(object$fitted[,1],lambda)
		x <- InvBoxCox(x,lambda)
		pmean <- InvBoxCox(pmean,lambda)
		if(biasadj){
		  pmean <- InvBoxCoxf(x = list(level = level, mean = pmean, upper = upper, lower = lower), lambda = lambda)
		}
		lower <- InvBoxCox(lower,lambda)
		upper <- InvBoxCox(upper,lambda)
	}	
	
    return(structure(list(method="HoltWinters", model=object, level=level,
        mean=pmean, lower=lower, upper=upper,
        x=x, xname=xname, fitted=object$fitted[,1], residuals=residuals(object)),
        class="forecast"))
}


## CROSTON

croston <- function(x,h=10,alpha=0.1)
{
    if(sum(x<0) > 0)
        stop("Series should not contain negative values")
    out <- croston2(x,h,alpha)
    out$x <- x
    out$xname <- deparse(substitute(x))
    if(!is.null(out$fitted))
      out$residuals <- x-out$fitted
    out$method <- "Croston's method"
    return(structure(out,class="forecast"))
}

croston2 <- function(x,h=10,alpha=0.1,nofits=FALSE)
{
    x <- as.ts(x)
    y <- x[x>0]
    tsp.x <- tsp(x)
    freq.x <- tsp.x[3]
    start.f <- tsp.x[2] + 1/freq.x
    if(length(y)==0) # All historical values are equal to zero
    {
      fc <- ts(rep(0,h), start=start.f, frequency=freq.x)
      if(nofits)
        return(fc)
      else
        return(list(mean=fc, fitted=ts(x*0, start=tsp.x[1], frequency=freq.x)))
    }
    tt <- diff(c(0,(1:length(x))[x>0])) # Times between non-zero observations
    if(length(y)==1 & length(tt)==1) # Only one non-zero observation
    {
      y.f <- list(mean=ts(rep(y,h), start=start.f, frequency=freq.x))
      p.f <- list(mean=ts(rep(tt,h), start=start.f, frequency=freq.x))
    }
    else if(length(y)<=1 | length(tt)<=1) # length(tt)==0 but length(y)>0. How does that happen?
        return(list(mean=ts(rep(NA,h), start=start.f, frequency=freq.x)))
    else
    {
      y.f <- ses(y,alpha=alpha,initial="simple",h=h)
      p.f <- ses(tt,alpha=alpha,initial="simple",h=h)
    }
    ratio <- ts(y.f$mean/p.f$mean,start=start.f, frequency = freq.x)
    if(nofits)
        return(ratio)
    else
    {
        n <- length(x)
        junk <- x*NA
        if(n > 1)
        {
          for(i in 1:(n-1))
              junk[i+1] <- croston2(x[1:i],h=1,alpha=alpha,nofits=TRUE)
        }
        junk <- ts(junk)
        tsp(junk) <- tsp.x
        return(list(mean = ratio, fitted = junk, model=list(demand=y.f,period=p.f)))
    }
}


snaive <- function(x,h=2*frequency(x),level=c(80,95),fan=FALSE, lambda=NULL)
{
    fc <- forecast(Arima(x,seasonal=list(order=c(0,1,0),period=frequency(x)), lambda=lambda),h=h,level=level,fan=fan)
    # Remove initial fitted values and error
    m <- frequency(x)
    fc$fitted[1:m] <- NA
    fc$residuals[1:m] <- NA
    fc$method <- "Seasonal naive method"
    return(fc)
}

naive <- function(x,h=10,level=c(80,95),fan=FALSE, lambda=NULL)
{
    fc <- forecast(Arima(x,order=c(0,1,0),lambda=lambda),h,level=level,fan=fan)
    # Remove initial fitted values and error
    fc$fitted[1] <- NA
    fc$residuals[1] <- NA
    fc$method <- "Naive method"
    return(fc)
}
