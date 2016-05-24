linear.fit<-function(x,y,weights,tol=1e-8,max.iter=25,verbose=FALSE){

  if (!is.matrix(x))
    x <- as.matrix(x)

  nvars <- NCOL(x)
  nobs <- NROW(x)
  if (missing(weights))
    weights=rep(1,nobs)

  fit <- .Fortran("dqrls", qr = x * weights, n = nobs,
          p = nvars, y = weights * y, ny = as.integer(1), tol = min(1e-07,
          glm.control()$epsilon/1000), coefficients = double(nvars),
          residuals = double(nobs), effects = double(nobs),
          rank = integer(1), pivot = 1:nvars, qraux = double(nvars),
          work = double(2 * nvars), PACKAGE = "CNVassoc")

  coeffic <- fit$coefficients
  mu <- drop(x%*%coeffic)
  sigma<-sqrt(sum((y-mu)^2*weights)/sum(weights))

  coeffic<-as.vector(coeffic)
  names(coeffic)<-colnames(x)

  return(list(coeffic=coeffic,sigma=sigma))

}
