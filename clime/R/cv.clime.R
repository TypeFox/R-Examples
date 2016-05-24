cv.clime <- function(clime.obj, loss=c("likelihood", "tracel2"), fold=5) {
  x <- clime.obj$x
  if (is.null(x)) stop("No x in clime object.  Use x instead of sigma for computing clime!")
  n <- nrow(x)
  p <- ncol(x)
  
  part.list <- cv.part(n, fold)
  
  lambda <- clime.obj$lambda
  lpfun <- clime.obj$lpfun
  
  nlambda <- length(lambda)

  lossname <- match.arg(loss, c("likelihood", "tracel2"))
  lossfun <- match.fun(lossname )
  
  loss.re <- matrix(0, nrow = fold, ncol = nlambda)
  for (j in 1:fold) {
      x.train <- x[part.list$trainMat[,j],]
      clime.cv <- clime(x.train, lambda, standardize = FALSE, perturb = clime.obj$perturb, linsolver=lpfun)
      x.test <- x[part.list$testMat[,j],]
      ntest <- nrow(x.test) 
      for (jl in 1:nlambda) {
        loss.re[j,jl] <- loss.re[j,jl]  + lossfun( (cov(x.test)*(1-1/ntest)),  clime.cv$Omegalist[[jl]])
      }
  }

  loss.mean <- apply(loss.re, 2, mean)
  loss.sd <- apply(loss.re, 2, sd)
  
  lambdaopt <- lambda[which.min(loss.mean)]
  

  outlist <- list(lambdaopt=lambdaopt, loss=lossname, lambda=lambda, loss.mean=loss.mean, loss.sd = loss.sd, lpfun=lpfun)
  class(outlist) <- c("cv.clime")
  return(outlist)
}
