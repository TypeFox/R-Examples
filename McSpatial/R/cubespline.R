cubespline <- function(form,knots=1,mink=1,maxk=1,crit="gcv",data=NULL) { 
cat("Reminder:  first explanatory variable is used for spline","\n")

  mat <- model.frame(form,data=data)
  y <- mat[,1]
  z <- mat[,2]
  n = length(y)
  nx = ncol(mat)-2
  if (nx>0) {xmat <- as.matrix(mat[,3:ncol(mat)])}
  zname <- colnames(mat)[2]
  xnames <- colnames(mat)[3:ncol(mat)]
  square <- z^2
  cube <- z^3
  minz = min(z)
  maxz = max(z)
  searchk <- maxk>mink
  maxk = ifelse(searchk==TRUE,maxk,knots)
  kstar = knots

  newform <- y~z+square+cube
  if (nx>0) {newform <- update(newform, ~. +xmat)}

  if (searchk==TRUE) {
    fit <- lm(newform)
    k = length(fit$coef) 
    sig2 <- mean(residuals(fit)^2)
    if (crit=="gcv") {critk = n*(n*sig2)/((n-k)^2) }
    if (crit=="sc")  {critk = log(sig2) + log(n)*k/n }
    if (crit=="aic") {critk = log(sig2) + 2*k/n }
    kstar = 0
    mincrit = critk
    cat("Information Criterion, Linear:    ", mincrit,"\n")

    newform <- update(newform, ~. +kvar)
    for (i in seq(mink,maxk)) {
      knots <- seq(minz,maxz,length=(i+2))
      kvar <- array(0,dim=c(n,i))
      for (j in seq(1:i)) {
        kvar[,j] <- ifelse(z>=knots[j+1], (z-knots[j+1])^3,0)   
      }
      fit <- lm(newform)
      k = length(fit$coef) 
      sig2 <- mean(residuals(fit)^2)
      if (crit=="gcv") {critk = n*(n*sig2)/((n-k)^2) }
      if (crit=="sc")  {critk = log(sig2) + log(n)*k/n }
      if (crit=="aic") {critk = log(sig2) + 2*k/n }
      cat("Information Criterion, k =",i,":", critk,"\n")
      if (critk < mincrit) {
        kstar = i
        mincrit = critk 
      }
    }
    cat("Information Criterion minimizing k = ",kstar,"\n")
    if (kstar==maxk) {cat("Warning:  best k = maximum allowed; may want to try higher value for maxk","\n") } 
    maxk = kstar
  }  

  if (kstar>0) {
    knots <- seq(minz,maxz,length=(maxk+2))
    kvar <- array(0,dim=c(n,maxk))
    for (j in seq(1:maxk)) {
      kvar[,j] <- ifelse(z>=knots[j+1], (z-knots[j+1])^3,0)   
    }
  }
 
  newform <- y~z+square+cube
  if (kstar>0) {newform <- update(newform, ~. +kvar)}
  if (nx>0) {newform <- update(newform, ~. +xmat)}
  fit <- lm(newform)

  k = length(fit$coef)
  yhat <- fitted(fit)
  rss = sum(residuals(fit)^2)
  sig2 = rss/n
  aic = log(sig2) + 2*k/n
  sc = log(sig2) + log(n)*k/n
  gcv = n*rss/((n-k)^2)

  splinehat <- yhat
  if (nx>0) {
    bhat <- fit$coefficients[(k-nx+1):k]
    xbhat <- as.array(as.matrix(xmat)%*%as.matrix(bhat))
    splinehat <- yhat - xbhat + mean(xbhat)
    names(fit$coefficients)[(k-nx+1):k] <- xnames
  }
  names(fit$coefficients)[1:4] <- c("Intercept", zname, "square", "cube")
  

  print(summary(fit))
  out <- list(yhat,rss,sig2,aic,sc,gcv,fit$coef,splinehat,knots)
  names(out) <- c("yhat","rss","sig2","aic","sc","gcv","coef","splinehat","knots")
  return(out)
}

