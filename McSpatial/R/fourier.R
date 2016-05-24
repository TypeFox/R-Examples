fourier <- function(form,q=1,minq=0,maxq=0,crit="gcv",data=NULL) {  
cat("Reminder:  first explanatory variable is used for fourier expansion","\n")

  mat <- model.frame(form,data=data)
  y <- mat[,1]
  z <- mat[,2]
  n = length(y)
  nx = ncol(mat)-2
  if (nx>0) {xmat <- as.matrix(mat[,3:ncol(mat)])}
  xnames <- colnames(mat)[3:ncol(mat)]
  minz = min(z)
  maxz = max(z)
  z <- 2*pi*(z-minz)/(maxz-minz)
  square <- z^2

  searchq = maxq>minq
  if (searchq==FALSE) {maxq = q}
  sinvar <- array(0,dim=c(n,maxq))
  cosvar <- array(0,dim=c(n,maxq))
  for (j in seq(1,maxq)) {
    sinvar[,j] <- sin(j*z)
    cosvar[,j] <- cos(j*z)
  }
  qstar = maxq

  newform <- y~z+square
  if (nx>0) {newform <- update(newform, ~. +xmat)}

  if (searchq==TRUE) {
    fit <- lm(newform)
    k = length(fit$coef) 
    sig2 <- mean(residuals(fit)^2)
    if (crit=="gcv") {critq = n*(n*sig2)/((n-k)^2) }
    if (crit=="sc")  {critq = log(sig2) + log(n)*k/n }
    if (crit=="aic") {critq = log(sig2) + 2*k/n }
    qstar = 0
    mincrit = critq
    cat("Information Criterion, Linear:    ", mincrit,"\n")

    newform <- update(newform, ~. +sinvar[,1:j]+cosvar[,1:j])
    for (j in seq(minq,maxq)) {
      fit <- lm(newform) 
      k = length(fit$coef) 
      sig2 <- mean(residuals(fit)^2)
      if (crit=="gcv") {critq = n*(n*sig2)/((n-k)^2) }
      if (crit=="sc")  {critq = log(sig2) + log(n)*k/n }
      if (crit=="aic") {critq = log(sig2) + 2*k/n }
      cat("Information Criterion, q =",j,":", critq,"\n")
      if (critq < mincrit) {
        qstar = j
        mincrit = critq 
      }
    }

    cat("Information Criterion minimizing q = ",qstar,"\n")
    if (qstar==maxq) {cat("Warning:  best q = maximum allowed; may want to try higher value for maxq","\n") }
   }
 

  maxq = ifelse(searchq==TRUE,qstar,q)
  newform <- y~z+square
  if (qstar>0) {
    sinvar <- as.matrix(sinvar[,1:maxq])
    cosvar <- as.matrix(cosvar[,1:maxq])
    newform <- update(newform, ~. +sinvar[,1:maxq]+cosvar[,1:maxq])
  }
  if (nx>0) {newform <- update(newform, ~. +xmat)}
  fit <- lm(newform)

  k = length(fit$coef) 
  sig2 <- mean(residuals(fit)^2)
  yhat <- fitted(fit)
  rss = sum(residuals(fit)^2)
  sig2 = rss/n
  aic = log(sig2) + 2*k/n
  sc = log(sig2) + log(n)*k/n
  gcv = n*rss/((n-k)^2)

  fourierhat <- yhat
  nx1 = 3+2*maxq+1
  nx2 = nx1 + nx -1
  if (nx>0) {
    bhat <- fit$coefficients[nx1:nx2]
    xbhat <- as.array(as.matrix(xmat)%*%as.matrix(bhat))
    fourierhat <- yhat - xbhat + mean(xbhat)
    names(fit$coefficients)[nx1:nx2] <- xnames
  }
  names(fit$coefficients)[1:3] <- c("Intercept", "z", "square")
  sinnames <- paste("sin(",c(1:maxq),sep="")
  sinnames <- paste(sinnames,"z)",sep="")
  cosnames <- paste("cos(",c(1:maxq),sep="")
  cosnames <- paste(cosnames,"z)",sep="")
  names(fit$coefficients)[4:(4+maxq-1)] <- sinnames
  names(fit$coefficients)[(4+maxq):(4+2*maxq-1)] <- cosnames
  fourierhat <- fourierhat - mean(fourierhat)  + mean(y)
  
  print(summary(fit))
  out <- list(yhat,rss,sig2,aic,sc,gcv,fit$coef,fourierhat,qstar)
  names(out) <- c("yhat","rss","sig2","aic","sc","gcv","coef","fourierhat","q")
  return(out)
}

