imputePCA <- function (X, ncp = 2, scale=TRUE, method=c("Regularized","EM"),row.w=NULL,coeff.ridge=1,threshold = 1e-6,seed = NULL,nb.init=1,maxiter=1000,...){

impute <- function (X, ncp = 4, scale=TRUE, method=NULL,threshold = 1e-6,seed = NULL,init=1,maxiter=1000,row.w=NULL,coeff.ridge=1,...){
    moy.p <- function(V, poids) {
        res <- sum(V * poids,na.rm=TRUE)/sum(poids[!is.na(V)])
    }
    ec <- function(V, poids) {
        res <- sqrt(sum(V^2 * poids,na.rm=TRUE)/sum(poids[!is.na(V)]))
    }

   nb.iter <- 1
   old <- Inf
   objective <- 0
   if (!is.null(seed)){set.seed(seed)}
   X <- as.matrix(X)
   ncp <- min(ncp,ncol(X),nrow(X)-1)
   missing <- which(is.na(X))
   mean.p <- apply(X, 2, moy.p,row.w)
   Xhat <- t(t(X)-mean.p)
   et <- apply(Xhat, 2, ec,row.w)
   if (scale) Xhat <- t(t(Xhat)/et)
   if (any(is.na(X))) Xhat[missing] <- 0
   if (init>1) Xhat[missing] <- rnorm(length(missing)) ## random initialization
   fittedX <- Xhat
   if (ncp==0) nb.iter=0
 
   while (nb.iter > 0) {
       Xhat[missing] <- fittedX[missing]
       if (scale) Xhat=t(t(Xhat)*et)
       Xhat <- t(t(Xhat)+mean.p)
       mean.p <- apply(Xhat, 2, moy.p,row.w)
       Xhat <- t(t(Xhat)-mean.p)
       et <- apply(Xhat, 2, ec,row.w)
       if (scale) Xhat <- t(t(Xhat)/et)

       svd.res <- svd.triplet(Xhat,row.w=row.w,ncp=ncp)
#       sigma2 <- mean(svd.res$vs[-(1:ncp)]^2)
	   sigma2  <- nrow(X)*ncol(X)/min(ncol(X),nrow(X)-1)* sum((svd.res$vs[-c(1:ncp)]^2)/((nrow(X)-1) * ncol(X) - (nrow(X)-1) * ncp - ncol(X) * ncp + ncp^2))
       sigma2 <- min(sigma2*coeff.ridge,svd.res$vs[ncp+1]^2)
       if (method=="em") sigma2 <-0
       lambda.shrinked=(svd.res$vs[1:ncp]^2-sigma2)/svd.res$vs[1:ncp]
       fittedX = tcrossprod(t(t(svd.res$U[,1:ncp,drop=FALSE]*row.w)*lambda.shrinked),svd.res$V[,1:ncp,drop=FALSE])
       fittedX <- fittedX/row.w
	   diff <- Xhat-fittedX
	   diff[missing] <- 0
       objective <- sum(diff^2*row.w)
#       objective <- mean((Xhat[-missing]-fittedX[-missing])^2)
       criterion <- abs(1 - objective/old)
       old <- objective
       nb.iter <- nb.iter + 1
       if (!is.nan(criterion)) {
         if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
         if ((objective < threshold) && (nb.iter > 5))  nb.iter <- 0
       }
       if (nb.iter > maxiter) {
         nb.iter <- 0
         warning(paste("Stopped after ",maxiter," iterations"))
       }
   }
   if (scale) Xhat <- t(t(Xhat)*et)
   Xhat <- t(t(Xhat)+mean.p)
   completeObs <- X
   completeObs[missing] <- Xhat[missing]
   if (scale) fittedX <- t(t(fittedX)*et)
   fittedX <- t(t(fittedX)+mean.p)

   result <- list()
   result$completeObs <- completeObs
   result$fittedX <- fittedX
   return(result) 
}

#### Main program
 method <- match.arg(method,c("Regularized","regularized","EM","em"),several.ok=T)[1]
 obj=Inf
 method <- tolower(method)
 if (ncp>min(nrow(X)-2,ncol(X)-1)) stop("ncp is too large")
 if (is.null(row.w)) row.w = rep(1,nrow(X))/nrow(X)
 for (i in 1:nb.init){
  if (!any(is.na(X))) return(X)
  res.impute=impute(X, ncp=ncp, scale=scale, method=method, threshold = threshold,seed=if(!is.null(seed)){(seed*(i-1))}else{NULL},init=i,maxiter=maxiter,row.w=row.w,coeff.ridge=coeff.ridge)
  if (mean((res.impute$fittedX[!is.na(X)]-X[!is.na(X)])^2) < obj){
    res <- res.impute
    obj <- mean((res.impute$fittedX[!is.na(X)]-X[!is.na(X)])^2)
  }
 }
return(res)
}

