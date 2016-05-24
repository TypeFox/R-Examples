##
## KPSS-Test
##
ur.kpss <- function(y, type=c("mu", "tau"), lags=c("short", "long", "nil"), use.lag=NULL){
  y <- na.omit(as.vector(y))
  n <- length(y)
  type <- match.arg(type)
  lags <- match.arg(lags)
  if(!(is.null(use.lag))){
    lmax <- as.integer(use.lag)
    if(lmax < 0){
      warning("\nuse.lag has to be positive and integer; lags='short' used.")
    lmax <- trunc(4*(n/100)^0.25)}
  }else if(lags == "short"){
    lmax <- trunc(4*(n/100)^0.25)
  }else if(lags == "long"){
    lmax <- trunc(12*(n/100)^0.25)
  }else if(lags == "nil"){
    lmax <- 0
  }
  if(type=="mu"){
    cval <- as.matrix(t(c(0.347, 0.463, 0.574, 0.739)))
    colnames(cval) <- c("10pct", "5pct", "2.5pct", "1pct")
    rownames(cval) <- "critical values"
    res <- y - mean(y)
  }else if(type=="tau"){
    cval <- as.matrix(t(c(0.119, 0.146, 0.176, 0.216)))
    colnames(cval) <- c("10pct", "5pct", "2.5pct", "1pct")
    rownames(cval) <- "critical values"
    trend <- 1:n
    res <- residuals(lm(y ~ trend))
  }
  S <- cumsum(res)
  nominator <- sum(S^2)/n^2
  s2 <- sum(res^2)/n
  if(lmax == 0){
    denominator <- s2
  }else{
    index <- 1:lmax
    x.cov <- sapply(index, function(x) t(res[-c(1:x)])%*%res[-c((n-x+1):n)])
    bartlett <- 1-index/(lmax+1)
    denominator <- s2 + 2/n*t(bartlett)%*%x.cov
  }
  teststat <- nominator/denominator
  new("ur.kpss", y=y, type=type, lag=as.integer(lmax), teststat=as.numeric(teststat), cval=cval, res=res , test.name="KPSS") 
}
