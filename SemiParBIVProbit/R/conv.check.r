conv.check <- function(x){

e.v <- eigen(x$fit$hessian, symmetric=TRUE, only.values = TRUE)$values

cat("\nLargest absolute gradient value:",max(abs(x$fit$gradient)))

if(x$hess==TRUE) mv <- "Observed" else mv <- "Expected" 

if (min(e.v) > 0) cat("\n",mv," information matrix is positive definite\n",sep="") else cat("\n",mv," information matrix is NOT positive definite\n",sep="")

cat("Eigenvalue range: [",min(e.v),",",max(e.v),"]\n", sep = "")

if( (x$l.sp1!=0 || x$l.sp2!=0 || x$l.sp3!=0 || x$l.sp4!=0 || x$l.sp5!=0 || x$l.sp6!=0 || x$l.sp7!=0) && x$fp==FALSE ){

cat("\nTrust region iterations before smoothing parameter estimation:",x$iter.if)
cat("\nLoops for smoothing parameter estimation:",x$iter.sp) 
cat("\nTrust region iterations within smoothing loops:",x$iter.inner,"\n\n")

}else{cat("\nTrust region iterations:",x$iter.if,"\n\n")}

if(!is.null(x$conv.sp)){
if(x$iter.sp >= x$iterlimsp) cat("Smoothing algorithm reached the max. number of loops allowed.\n\n")
}


}
