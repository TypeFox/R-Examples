`Beta2parcor` <-
function(Beta,verbose=FALSE){
    Dummy=Beta*t(Beta)
    if (verbose==TRUE){
    cat("\nNumber of pairwise regression coefficients with conflicting signs:", (sum((Dummy) < 0))/2, "\n")
  cat("Number of partial correlation coefficients greater than 1 in absolute value:",(sum((Dummy) >1))/2 ,
  "\n\n")
  }
    Dummy[Dummy<0]=0 # if a product is <0 the partial correlation coefficient is set to 0
    Dummy[Dummy>1]=1 # partial correlation coefficients should be in the range of [-1,1]
    P=sign(Beta)*sqrt(Dummy)
    diag(P)=rep(1,ncol(P))
    return(P)
}
