gcrq.rq.cv <-
function(y, B, X, taus, interc=FALSE, monotone, ndx, lambda, deg, dif, var.pen=NULL, cv=TRUE, nfolds=10, 
foldid=NULL, eps=.0001){
#perform 'nfolds' cross-validation to select lambda.
#B e X devono essere (eventualmente) matrici!!!
    Rho <- function(u, tau) u * (tau - (u < 0))
    if(length(lambda)<=1) stop("lambda should be a vector")
    n<-length(y)
    if (is.null(foldid))
        foldid = sample(rep(seq(nfolds), length = n)) #sample(seq(nfolds), size=n, replace=TRUE)
    else nfolds = max(foldid)
    if (nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")
    CV<-matrix(NA, length(lambda), nfolds)
    for(j in 1:length(lambda)){
      lambda1<-lambda[j]
      for(i in seq(nfolds)) {
          which <- foldid == i
          y_sub = y[!which]
          B_sub<-list(B[!which, , drop = FALSE])
          if(ncol(X)<=0){
          fit<-ncross.rq.fitXB(y=y_sub, B=B_sub, X=NULL, taus=taus, monotone=monotone, ndx=ndx, lambda=lambda1, deg=deg,
              dif=dif, var.pen=var.pen)
          xreg<-B[which,,drop=FALSE]
          } else {
          fit<-ncross.rq.fitXB(y=y_sub, B=B_sub, X=X[!which, , drop = FALSE], taus=taus, interc=interc, monotone=monotone, ndx=ndx, 
              lambda=lambda1, deg=deg, dif=dif, var.pen=var.pen)
          xreg<-cbind(X,B)[which,,drop=FALSE]
          }
          fit.values<-predictQR(fit, xreg=xreg)
          fit.values<-if(is.matrix(fit.values)) rbind(as.numeric(colnames(fit.values)),fit.values) else matrix(c(taus, fit.values), ncol=1)
          rho.values<-apply(fit.values, 2, function(z) sum(Rho(y[which]-z[-1],z[1])))        
          CV[j,i]<-sum(rho.values)
          }
      }
    attr(CV, "foldid")<-foldid
    lambda.ok <- lambda[which.min(apply(CV, 1, mean))]
    if(cv) lambda.ok<-list(lambda.ok=lambda.ok, CV=CV)
    return(lambda.ok)
    }
