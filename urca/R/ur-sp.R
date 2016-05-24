##
## Schmidt-Phillips Test
##
ur.sp <- function(y, type=c("tau", "rho"), pol.deg=c(1, 2, 3, 4), signif=c(0.01, 0.05, 0.1)){
  y <- na.omit(as.vector(y))
  type <- match.arg(type)
  signif <- signif[1]
  signif.val <- c(0.01, 0.05, 0.1)
  if(!(signif %in% signif.val)){
    warning("\nPlease, provide as signif one of c(0.01, 0.05, 0.1); signif=0.01 used")
  signif <- 0.01}
  pol.deg <- pol.deg[1]
  if(!(pol.deg %in% c(1:4))){
    warning("\nPlease, provide as polynomial degree one of c(1, 2, 3, 4); po.deg=1 used")
  pol.deg <- 1}
  n <- length(y)
  lag <- trunc(12*(n/100)^0.25)
  idx <- 1:lag
  trend1 <- 1:n
  y.diff <- diff(y)
  if(pol.deg==1){
    yi.hat <- (y[n]-y[1])/(n-1)
    phi.y <- y[1]-yi.hat
    S.hat <- y - phi.y - yi.hat*trend1
    S.hat.l1 <- S.hat[-n]
    test.reg <- summary(lm(y.diff ~ 1 + S.hat.l1))
    sp.data <- data.frame(cbind(y[2:n],  y[1:(n-1)], trend1[2:n]))
    colnames(sp.data) <- c("y", "y.lagged", "trend.exp1")
    corr.reg <- summary(lm(sp.data))
    res <- residuals(corr.reg)
    sig.eps <- (1/n)*sum(res^2)
    coprods <- sapply(idx, function(x) t(res[-c(1:x)])%*%(res[-c((n-x):(n-1))]))
    weights <- (2*(lag-idx)/lag)^2
    sig <- sig.eps + (2/n)*(t(weights)%*%coprods)
    omega2.hat <- sig.eps/sig
  }else if(pol.deg==2){
    trend2 <- trend1^2
    S.hat <- c(0, cumsum(residuals(summary(lm(y.diff ~ trend1[2:n])))))
    test.reg <- summary(lm(y.diff ~ S.hat[-n] + trend1[2:n]))
    sp.data <- data.frame(cbind(y[2:n],  y[1:(n-1)], trend1[2:n], trend2[2:n]))
    colnames(sp.data) <- c("y", "y.lagged", "trend.exp1", "trend.exp2")
    corr.reg <- summary(lm(sp.data))
    res <- residuals(corr.reg)
    sig.eps <- (1/n)*sum(res^2)
    coprods <- sapply(idx, function(x) t(res[-c(1:x)])%*%(res[-c((n-x):(n-1))]))
    weights <- (2*(lag-idx)/lag)^2
    sig <- sig.eps + (2/n)*(t(weights)%*%coprods)
    omega2.hat <- sig.eps/sig
  }else if(pol.deg==3){
    trend2 <- trend1^2
    trend3 <- trend1^3    
    S.hat <- c(0, cumsum(residuals(summary(lm(y.diff ~ trend1[2:n] + trend2[2:n])))))
    test.reg <- summary(lm(y.diff ~ S.hat[-n] + trend1[2:n] + trend2[2:n]))
    sp.data <- data.frame(cbind(y[2:n],  y[1:(n-1)], trend1[2:n], trend2[2:n], trend3[2:n]))
    colnames(sp.data) <- c("y", "y.lagged", "trend.exp1", "trend.exp2", "trend.exp3")
    corr.reg <- summary(lm(sp.data))
    res <- residuals(corr.reg)
    sig.eps <- (1/n)*sum(res^2)
    coprods <- sapply(idx, function(x) t(res[-c(1:x)])%*%(res[-c((n-x):(n-1))]))
    weights <- (2*(lag-idx)/lag)^2
    sig <- sig.eps + (2/n)*(t(weights)%*%coprods)
    omega2.hat <- sig.eps/sig
  }else if(pol.deg==4){
    trend2 <- trend1^2
    trend3 <- trend1^3
    trend4 <- trend1^4
    S.hat <- c(0, cumsum(residuals(summary(lm(y.diff ~ trend1[2:n] + trend2[2:n] + trend3[2:n])))))
    test.reg <- summary(lm(y.diff ~ S.hat[-n] + trend1[2:n] + trend2[2:n] + trend3[2:n]))
    sp.data <- data.frame(cbind(y[2:n],  y[1:(n-1)], trend1[2:n], trend2[2:n], trend3[2:n], trend4[2:n]))
    colnames(sp.data) <- c("y", "y.lagged", "trend.exp1", "trend.exp2", "trend.exp3", "trend.exp4")
    corr.reg <- summary(lm(sp.data))
    res <- residuals(corr.reg)
    sig.eps <- (1/n)*sum(res^2)
    coprods <- sapply(idx, function(x) t(res[-c(1:x)])%*%(res[-c((n-x):(n-1))]))
    weights <- (2*(lag-idx)/lag)^2
    sig <- sig.eps + (2/n)*(t(weights)%*%coprods)
    omega2.hat <- sig.eps/sig
  }
  if(type=="rho"){
    rho <- n*coef(test.reg)[2,1]
    teststat <- rho/omega2.hat
    cval <- .spcv(obs=n, type="rho", pol.deg=pol.deg, signif=signif)
  }else if(type=="tau"){
    tau <- coef(test.reg)[2,3]
    teststat <- tau/sqrt(omega2.hat)
    cval <- .spcv(obs=n, type="tau", pol.deg=pol.deg, signif=signif)
  }
  new("ur.sp", y=y, type=type, polynomial=as.integer(pol.deg), teststat=as.numeric(teststat), cval=cval, signif=signif, res=res, testreg=corr.reg, test.name="Schmidt-Phillips")
}
