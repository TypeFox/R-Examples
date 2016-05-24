binomControl<-function(relErr=1+1e-07,tol=.00001,pRange=c(1e-10,1-1e-10),minn=1,maxn=10^4,PRINT.STEPS=FALSE){
    if (relErr<=1) stop("relErr should be greater than 1")
    if (tol<=0) stop("tol should be greater than 0")
    if (pRange[1]>pRange[2] | pRange[1]<0 | pRange[1]>1 | pRange[2]<0 | pRange[2]>1) stop("pRange should be a range of possible values for the binomial parameter p")
if (minn>maxn) stop("minn should be less than maxn")
minn<-as.integer(minn)
maxn<-as.integer(maxn)
if (!is.logical(PRINT.STEPS)) stop("PRINT.STEPS must be logical")
    list(relErr=relErr,tol=tol,pRange=pRange,
         minn=minn,maxn=maxn,PRINT.STEPS=PRINT.STEPS)
}