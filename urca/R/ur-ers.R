##
## Elliott, Rothenberg and Stock-Test
##
ur.ers <- function(y, type=c("DF-GLS", "P-test"), model=c("constant", "trend"), lag.max=4){
  type <- match.arg(type)
  model <- match.arg(model)
  lag.max <- as.integer(lag.max)
  if(lag.max < 0){
    warning("\nlag.max bust be greater or equal to one and integer; setting lag.max=4")
  lag.max <- 4}
  lag.max <- lag.max+1
  idx <- 2:lag.max
  y <- na.omit(as.vector(y))
  nobs <- length(y)
  if(nobs < 50){
    rowsel <- 1
  }else if(nobs < 100){
    rowsel <- 2
  }else if(nobs <= 200){
    rowsel <- 3
  }else if(nobs > 200){
    rowsel <- 4}
  if(model=="constant"){
    ahat <- 1 - 7.0/nobs
    ya <- c(y[1], y[2:nobs]-ahat*y[1:(nobs-1)])
    za1 <- c(1, rep(1-ahat, nobs-1))
    yd.reg <- summary(lm(ya ~ -1 + za1))
    yd <- y - coef(yd.reg)[1]
  }else if(model=="trend"){
    ahat <- 1 - 13.5/nobs
    ya <- c(y[1], y[2:nobs]-ahat*y[1:(nobs-1)])
    za1 <- c(1, rep(1-ahat, nobs-1))
    trd <- 1:nobs
    za2 <- c(1, trd[2:nobs]-ahat*trd[1:(nobs-1)])
    yd.reg <- summary(lm(ya ~ -1 + za1 + za2))
    yd <- y - coef(yd.reg)[1] - coef(yd.reg)[2]*trd 
  }
  what <- function(x, z=y){
    z.l <-  z[1:(nobs-1)]
    z.diff <- diff(z)
    z.dlags <- embed(diff(z), x)[, -1]
    data.what <- data.frame(cbind(z.diff[-(1:(x-1))], z.l[-(1:(x-1))], z.dlags))
    bic <- BIC(lm(data.what))
    return(bic)
  }
  if(type=="P-test"){
    cvals.ptest <- array(c(1.87, 1.95, 1.91, 1.99, 2.97, 3.11, 3.17, 3.26, 3.91, 4.17, 4.33, 4.48, 4.22, 4.26, 4.05, 3.96, 5.72, 5.64, 5.66, 5.62, 6.77, 6.79, 6.86, 6.89), c(4, 3, 2))
    res <- residuals(yd.reg)
    if(model=="constant"){
      null.res <- c(0, diff(y))
      cvals <- as.matrix(t(cvals.ptest[rowsel, , 1]))
      model <- "with intercept"
    }else if(model=="trend"){
      null.res <- c(0, diff(y))
      null.res <- null.res - mean(null.res)
      cvals <- as.matrix(t(cvals.ptest[rowsel, , 2]))
      model <- "with intercept and trend"
    }
  sig.null <- sum(null.res^2)
  sig.res <- sum(res^2)
  if(lag.max > 1){
    bic <- sapply(idx, what, z=y)
    BIC.opt <- which.min(bic)+1
    y.l <-  y[1:(nobs-1)]
    y.diff <- diff(y)
    y.dlags <- embed(diff(y), BIC.opt)[, -1]
    data.what <- data.frame(cbind(y.diff[-(1:(BIC.opt-1))], y.l[-(1:(BIC.opt-1))], y.dlags))
    what.reg <- summary(lm(data.what))
    npar <- nrow(what.reg$coef)
    sumlc <- sum(what.reg$coef[3:npar,1])
    lag.max <- BIC.opt-1
  }else if(lag.max <= 1){
    y.diff <- diff(y)
    y.l <- y[1:(nobs-1)]
    what.reg <- summary(lm(y.diff ~ y.l))
    sumlc <- 0
    lag.max <- lag.max-1
  }
  what.sq <- what.reg$sigma^2/(1-sumlc)^2
  teststat <- (sig.res - ahat*sig.null)/what.sq
  test.reg <- NULL
  }else if(type=="DF-GLS"){
    if(model=="constant"){
      cvals <- as.matrix(t(c(-2.5658-1.960/nobs-10.04/(nobs**2),-1.9393-0.398/nobs,-1.6156-0.181/nobs)))
      model <- "with intercept"
    }else if(model=="trend"){
      cvals.dfgls.tau <- matrix(-1*c(3.77, 3.58, 3.46, 3.48, 3.19, 3.03, 2.93, 2.89, 2.89, 2.74, 2.64, 2.57), nrow=4, ncol=3)
      cvals <- as.matrix(t(cvals.dfgls.tau[rowsel,]))
      model <- "with intercept and trend"
    }
    yd.l <-  yd[1:(nobs-1)]
    yd.diff <- diff(yd)
    if(lag.max > 1){
      yd.dlags <- embed(diff(yd), lag.max)[, -1]
      data.dfgls <- data.frame(cbind(yd.diff[-(1:(lag.max-1))], yd.l[-(1:(lag.max-1))], yd.dlags))
      colnames(data.dfgls) <- c("yd.diff", "yd.lag", paste("yd.diff.lag", 1:(lag.max-1), sep=""))
      dfgls.form <- formula(paste("yd.diff ~ -1 + ", paste(colnames(data.dfgls)[-1], collapse=" + ")))
    }else if(lag.max <=1){
      data.dfgls <- data.frame(cbind(yd.diff, yd.l))
      colnames(data.dfgls) <- c("yd.diff", "yd.lag")
      dfgls.form <- formula("yd.diff ~ -1 + yd.lag")
    }
    dfgls.reg <- summary(lm(dfgls.form, data=data.dfgls))
    teststat <- coef(dfgls.reg)[1,3]
    test.reg <- dfgls.reg
    lag.max <- lag.max-1
  }
  colnames(cvals) <- c("1pct", "5pct", "10pct")
  rownames(cvals) <- c("critical values")
  new("ur.ers", y=y, yd=yd, type=type, model=model, lag=as.integer(lag.max), cval=round(cvals, 2), teststat=teststat, testreg=test.reg, test.name="Elliot, Rothenberg and Stock")
}
