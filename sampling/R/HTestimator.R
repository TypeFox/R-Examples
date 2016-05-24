HTestimator<-function(y,pik) 
{
if(any(is.na(pik))) stop("there are missing values in pik")
if(any(is.na(y))) stop("there are missing values in y")
if(length(y)!=length(pik)) stop("y and pik have different sizes")
crossprod(y,1/pik)
}
