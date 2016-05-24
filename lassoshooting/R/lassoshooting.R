lassoshooting <- function(X=NULL, y=NULL, lambda, XtX = NULL, Xty = NULL,
	   thr = 1e-06, maxit = 10000, nopenalize = NULL,
	   penaltyweight = NULL, trace = 0, ...) {
  .External("ccd", trace=trace,X=X,y=y,XtX=XtX,Xty=Xty,lambda=lambda,thr=thr,maxit=maxit,nopenalize=nopenalize,penaltyweight=penaltyweight,..., PACKAGE="lassoshooting")
}
softthresh <- function(x,t) {
  .External("R_softthresh",x,t)
}
