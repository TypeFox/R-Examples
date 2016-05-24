TYLERshape <- function (X, location = TRUE, eps = 1e-06, maxiter = 100) 
{
  Xnames <- colnames(X)
  X <- as.matrix(X)
  p <- ncol(X)
  
  if (location==TRUE) {
        mu <- colMeans(X)
        X<-t(t(X)-mu)
        }
  if (location==FALSE) {
        mu <- rep(0,p)
        names(mu) <- Xnames
        }
  if (is.numeric(location)) {
        mu <- location
        X<-t(t(X)-mu)
        names(mu) <- Xnames
        }
    
res <- .Call("Tyler0", X, FALSE, eps, maxiter, PACKAGE = "fastM") 
Sigma <- res$S
Sigma <- Sigma / det(Sigma)^(1/p)
rownames(Sigma) <- colnames(Sigma) <- Xnames
  
if (res$nG > eps) {warning("TYLERshape did not converge", call. = FALSE)}
return(list(mu = mu, Sigma=Sigma, iter=res$iter))

}
