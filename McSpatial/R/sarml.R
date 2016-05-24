sarml <- function(form,wmat=NULL,shpfile=NULL,wy=NULL,eigvar=NULL,startrho=NULL,print=TRUE,data=NULL) {
  library(spdep)
  newform <- as.formula(form,env=data)

  mat <- model.frame(newform,data=data)
  y <- mat[,1]
  xmat <- model.matrix(newform,data=data)
  namey = colnames(mat)[1]
  namex = colnames(xmat)
  
# don't need w matrix if wy and eigvar are both provided
# need w to calculate wy if eigvar provided but not wy
# need w to calculate eigvar if wy provided but not eigvar

  needw <- ((identical(wy,NULL))|(identical(eigvar,NULL)))&identical(wmat,NULL)
  if (identical(eigvar,NULL)&identical(shpfile,NULL)) {stop("Shape file needed for eigenvalue calculation")}
  if (needw==TRUE&identical(shpfile,NULL)) {stop("Shape file needed to calculate W")}
  if (needw==TRUE|identical(eigvar,NULL)) {neighbors <- poly2nb(shpfile,queen=TRUE)}

  if (needw==TRUE) {
    wmat <- nb2mat(neighbors,zero.policy=TRUE)
    nlist <- nb2listw(neighbors,zero.policy=TRUE) 
  }
  if (identical(wy,NULL)) {wy <- as.vector(wmat%*%y)}
  if (identical(eigvar,NULL)) {eigvar <- eigenw(nlist,quiet=TRUE)}
 
  xmat <- cbind(xmat,wy)
  nk = ncol(xmat)
  n = nrow(xmat)

  xdata <- data.frame(mat,wy)
  newform <- as.formula(newform,env=xdata)

# iterate concentrated likelihood function
  if (identical(startrho,NULL)) { 
    fit <- lm(newform, data=xdata)
    err0 <- residuals(fit)
    b0 <- coef(fit)
    formwy <- update(newform, wy~.)
    fit <- lm(formwy, data=xdata)
    err1 <- residuals(fit)
    b1 <- coef(fit)
    lc <- function(rho) {
      err2 <- (err0 - rho*err1)^2
      lc = n*log(mean(err2))/2 - sum(log(1-rho*eigvar))
      return(lc)
    }
    startrho = optimize(lc,lower=-.99,upper=.99)$minimum
  }

  e <- y-startrho*wy
  xdata <- data.frame(xdata,e)
  fit <- lm(update(form, e~.),data=xdata)
  startb <- fit$coef
  bmat <- c(startb,startrho)
  startsig2 <- summary(fit)$sigma^2

 # nlm minimizes function logl
  logl <- function(b) {
    s2 = b[nk+1]
    rho = b[nk]
    bmat <- b[1:nk]
    xb <- xmat%*%bmat
    e <- y - xb
    out <-  n*log(pi)/2 + n*log(s2)/2 + sum(e^2)/(2*s2) - sum(log(1 - rho*eigvar))
    return(out)
  }
  options(warn=-1)
  nlfit <- nlm(logl, c(startb,startrho,startsig2), hessian=TRUE)
  options(warn=0)

  b <- nlfit$estimate
  rho = b[nk]
  sig2 = b[nk+1]
  xx <- crossprod(xmat)
  xx[nk,nk] = xx[nk,nk] + sig2*sum((eigvar/(1-rho*eigvar))^2)
  vmat <- sig2*solve(xx)
  semat <- sqrt(diag(vmat))

  outmat <- cbind(bmat, semat, bmat/semat, 2*(1-pnorm(abs(bmat)/semat)) )
  colnames(outmat) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  rownames(outmat)[nk] <- "rho"
  outmat
  if (print==TRUE) {
    print(outmat)
    cat("sig2 = ", sig2, "\n")
  }

  outmat <- list(bmat[1:(nk-1)], rho, sig2, vmat, eigvar)
  names(outmat) <- c("beta", "rho", "sig2", "vmat", "eigvar")
  return(outmat)
}

