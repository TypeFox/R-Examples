##
## Phillips-Perron Test
##
ur.pp <- function(x, type=c("Z-alpha", "Z-tau"), model=c("constant", "trend"), lags=c("short", "long"), use.lag=NULL){
  x <- na.omit(as.vector(x))
  n <- length(x)
  y <- x[-1]
  y.l1 <- x[-n]
  n <- n-1
  lags <- match.arg(lags)
  model <- match.arg(model)
  type <- match.arg(type)
  if(!(is.null(use.lag))){
    lmax <- as.integer(use.lag)
    if(lmax < 0){
      warning("\nuse.lag has to be positive and integer; lags='short' used.")
      lmax <- trunc(4*(n/100)^0.25)}
  }else if(lags == "short"){
    lmax <- trunc(4*(n/100)^0.25)
  }else if(lags == "long"){
    lmax <- trunc(12*(n/100)^0.25)}
  if(model=="trend"){
    cval <- as.matrix(t(c(-3.9638-8.353/n-47.44/(n^2), -3.4126-4.039/n-17.83/(n^2), -3.1279-2.418/n-7.58/(n^2))))
    colnames(cval) <- c("1pct", "5pct", "10pct")
    rownames(cval) <- "critical values"
    model <- "with intercept and trend"
    trend <- (1:n) - n/2
    test.reg <- summary(lm(y ~ y.l1 + trend))
    res <- residuals(test.reg)
    my.tstat <- coef(test.reg)[1, 3]
    beta.tstat <- coef(test.reg)[3, 3]
    res <- residuals(test.reg)
    s <- 1/n*(sum(res^2))
    myybar <- (1/n^2)*sum((y-mean(y))^2)
    myy <- (1/n^2)*sum(y^2)
    mty <- (n^(-5/2))*(t(1:n)%*%y)
    my <- (n^(-3/2))*sum(y)
    idx <- 1:lmax
    coprods <- sapply(idx, function(l) t(res[-c(1:l)])%*%(res[-c((n-l+1):n)]))
    weights <- 1 - idx/(lmax+1)
    sig <- s + (2/n)*(t(weights)%*%coprods)
    lambda <- 0.5*(sig-s)
    lambda.prime <- lambda/sig
    M <- (1-n^(-2))*myy - 12*mty^2 + 12*(1 + 1/n)*mty*my - (4 + 6/n + 2/n^2)*my^2
    my.stat <- sqrt(s/sig)*my.tstat - lambda.prime*sqrt(sig)*my/(sqrt(M)*sqrt((M+my^2)))
    beta.stat <- sqrt(s/sig)*beta.tstat - lambda.prime*sqrt(sig)*(0.5*my - mty)/(sqrt(M/12)*sqrt(myybar))
    aux.stat <- as.matrix(c(round(my.stat, 4), round(beta.stat, 4)))
    rownames(aux.stat) <- c("Z-tau-mu", "Z-tau-beta")
    colnames(aux.stat) <- "aux. Z statistics"
    if(type=="Z-tau"){
      tstat <- (coef(test.reg)[2, 1] - 1)/coef(test.reg)[2, 2]
      teststat <- sqrt(s/sig)*tstat-lambda.prime*sqrt(sig)/sqrt(M)
    }else if(type=="Z-alpha"){
      alpha <- coef(test.reg)[2, 1]
      teststat <- n*(alpha-1)-lambda/M
      cval <- as.matrix(t(c(NA, NA, NA)))
    }
  }else if(model=="constant"){
    cval <- as.matrix(t(c(-3.4335-5.999/n-29.25/(n^2), -2.8621-2.738/n-8.36/(n^2), -2.5671-1.438/n-4.48/(n^2))))
    colnames(cval) <- c("1pct", "5pct", "10pct")
    rownames(cval) <- "critical values"
    model <- "with intercept"
    test.reg <- summary(lm(y ~ y.l1))
    my.tstat <- coef(test.reg)[1, 3]
    res <- residuals(test.reg)
    s <- 1/n*(sum(res^2))
    myybar <- (1/n^2)*sum((y-mean(y))^2)
    myy <- (1/n^2)*sum(y^2)
    my <- (n^(-3/2))*sum(y)
    idx <- 1:lmax
    coprods <- sapply(idx, function(l) t(res[-c(1:l)])%*%(res[-c((n-l+1):n)]))
    weights <- 1 - idx/(lmax+1)
    sig <- s + (2/n)*(t(weights)%*%coprods)
    lambda <- 0.5*(sig-s)
    lambda.prime <- lambda/sig
    my.stat <- sqrt(s/sig)*my.tstat + lambda.prime*sqrt(sig)*my/(sqrt(myy)*sqrt(myybar))
    aux.stat <- as.matrix(round(my.stat, 4))
    rownames(aux.stat) <- "Z-tau-mu"
    colnames(aux.stat) <- "aux. Z statistics"
    if(type=="Z-tau"){
      tstat <- (coef(test.reg)[2, 1] - 1)/coef(test.reg)[2, 2]
      teststat <- sqrt(s/sig)*tstat-lambda.prime*sqrt(sig)/sqrt(myybar)
    }else if(type=="Z-alpha"){
      alpha <- coef(test.reg)[2, 1]
      teststat <- n*(alpha-1)-lambda/myybar
      cval <- as.matrix(t(c(NA, NA, NA)))
    }
  }
  new("ur.pp", y=y, type=type, model=model, lag=as.integer(lmax), cval=cval, teststat=as.numeric(teststat), testreg=test.reg, auxstat=aux.stat, res=res, test.name="Phillips-Perron")
}
