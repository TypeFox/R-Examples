mshapiro.test <-
function (x) {
  if (!is.matrix(x)) {x <- as.matrix(x)}
  x <- x[complete.cases(x),]
  x <- t(x)
  n <- ncol(x)
  if (n<3 || n>5000) {stop("sample size must be between 3 and 5000")}
  rng <- range(x)
  rng <- rng[2]-rng[1]
  if (rng==0) {stop("all `x[]' are identical")}
  Us <- apply(x,1,mean)
  R <- x-Us
  M.1 <- solve(R%*%t(R),tol=1e-18)
  Rmax <- diag(t(R)%*%M.1%*%R)
  C <- M.1%*%R[,which.max(Rmax)]
  Z <- t(C)%*%x
  result <- shapiro.test(Z)
  result$method <- "Multivariate Shapiro-Wilk normality test"
  result$data.name <- paste("(",paste(rownames(x),collapse=","),")",sep="")
  return(result)
}
