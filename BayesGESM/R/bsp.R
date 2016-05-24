bsp <-
function(x, kn){
     if(!is.numeric(x)) stop("Variable in B-spline must be numeric!!",call.=FALSE)
as.matrix(x)
  if(missingArg(kn)){
 kn <-  floor(length(x)^(1/5))
}
 else{
if(kn<1 | kn != floor(kn))
stop("the number of knots must be a positive integer!!", call.=FALSE)
  }
  B <- bs(x,knots=quantile(x,prob=seq(1,kn,by=1)/(kn+1)),intercept=FALSE)
  attr(x,"B") <- B
  attr(x,"kn") <- kn
 x
}
