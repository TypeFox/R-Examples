##
## Phillips-Ouliaris Test
##
ca.po <- function(z, demean=c("none", "constant", "trend"), lag=c("short", "long"), type=c("Pu", "Pz"), tol=NULL){
  z <- na.omit(as.matrix(z))
  if(ncol(z)<2 || ncol(z)>6){
    stop("Please provide a matrix with at least two and maximal six columns")}
  demean <- match.arg(demean)
  lag <- match.arg(lag)
  type <- match.arg(type)
  nobs <- nrow(z)
  m <- ncol(z)
  zl <- z[2:nobs,]
  zr <- z[1:(nobs-1),]
  nobs <- nobs-1
  if(lag == "short"){
    lmax <- trunc(4*(nobs/100)^0.25)
  }else if(lag == "long"){
    lmax <- trunc(12*(nobs/100)^0.25)
  }
  if(demean=="none"){
    ari3 <- 1
    model <- "none"
    res <- residuals(lm(zl ~ zr - 1))
    if(type=="Pu"){
      resu <- residuals(lm(z[,1] ~ z[,-1] -1))
      test.reg <- summary(lm(z[,1] ~ z[,-1] -1))
    }else if(type=="Pz"){
      test.reg <- summary(lm(zl ~ zr - 1))}
  }else if(demean=="constant"){
    ari3 <- 2
    model <- "with constant only"
    res <- residuals(lm(zl ~ zr))
    if(type=="Pu"){
      resu <- residuals(lm(z[,1] ~ z[,-1]))
      test.reg <- summary(lm(z[,1] ~ z[,-1]))
    }else if(type=="Pz"){
      test.reg <- summary(lm(zl ~ zr))}
  }else if(demean=="trend"){
    ari3 <- 3
    model <- "with constant and linear trend"
    trd <- 1:nobs
    res <- residuals(lm(zl ~ zr + trd))
    if(type=="Pu"){
      trd <- 1:(nobs+1)
      resu <- residuals(lm(z[,1] ~ z[,-1] + trd))
      test.reg <- summary(lm(z[,1] ~ z[,-1] + trd))
    }else if(type=="Pz"){
      test.reg <- summary(lm(zl ~ zr + trd))}
  }
  xi2 <- 1/nobs*t(res)%*%res
  index <- 1:lmax
  wsl <- 1-index/(lmax+1)
  xi.mat <- sapply(index, function(x){
    wsl[x]*(t(res[-c(1:x),])%*%res[-c((nobs-x+1):nobs),] +  t(res[-c((nobs-x+1):nobs),])%*%res[-c(1:x),])
  })
  xi.mat <- array(xi.mat, c(m, m, lmax))
  smat <- matrix(0, m, m)
  for(i in 1:lmax)
    smat <- smat + xi.mat[,,i]
  omega <- xi2 + 1/nobs*smat
  if(type=="Pz"){
    Mzz <- (1/nobs*t(zl)%*%zl)
    Mzzinv <- solve(Mzz, tol=tol)
    teststat <- nobs*sum(diag(omega%*%Mzzinv))
    cvals <- array(c(33.9267, 62.1436, 99.2664, 143.0775, 195.6202, 40.8217, 71.2751, 109.7426, 155.8019, 210.2910, 55.1911, 89.6679, 131.5716, 180.4845, 237.7723, 47.5877, 80.2034, 120.3035, 168.8572, 225.2303, 55.2202, 89.7619, 132.2207, 182.0749, 241.3316, 71.9273, 109.4525, 153.4504, 209.8054, 270.5018, 71.9586, 113.4929, 163.1050, 219.5098, 284.0100, 81.3812, 124.3933, 175.9902, 234.2865, 301.0949, 102.0167, 145.8644, 201.0905, 264.4988, 335.9054), c(5, 3, 3))
    cval <- as.matrix(t(cvals[m-1, ,ari3]))
  }else if(type=="Pu"){
    w11 <- omega[1,1]
    w21 <- omega[2:m,1]
    O22 <- omega[-1,-1]
    O22inv <- solve(O22, tol=tol)
    w112 <- w11 - t(w21)%*%O22inv%*%w21
    ssqr <- 1/nobs*sum(resu^2)
    teststat <- as.numeric(nobs*w112/ssqr)
    cvals <- array(c(20.3933, 26.7022, 33.5359, 39.2826, 44.3725, 25.9711, 32.9392, 40.1220, 46.2691, 51.8614, 38.3413, 46.4097, 55.7341, 63.2149, 69.4939, 27.8536, 33.6955, 39.6949, 45.3308, 50.3537, 33.713, 40.5252, 46.7281, 53.2502, 57.7855, 48.0021, 53.8731, 63.4128, 71.5214, 76.7705, 41.2488, 46.1061, 52.0015, 57.3667, 61.6155, 48.8439, 53.8300, 60.2384, 65.8706, 70.7416, 65.1714, 69.2629, 78.3470, 84.5480, 91.0392)  , c(5, 3, 3))
    cval <- as.matrix(t(cvals[m-1, ,ari3]))
    res <- as.matrix(resu)
  }
  colnames(cval) <- c("10pct", "5pct", "1pct")
  rownames(cval) <- "critical values"
  new("ca.po", z=z, type=type, model=model, lag=as.integer(lmax), cval=cval, res=res, teststat=teststat, testreg=test.reg, test.name="Phillips and Ouliaris")
}
