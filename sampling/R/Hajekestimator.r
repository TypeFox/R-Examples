Hajekestimator<-function(y,pik,N=NULL,type=c("total","mean")) 
{
if(any(is.na(pik))) stop("there are missing values in pik")
if(any(is.na(y))) stop("there are missing values in y")
if(length(y)!=length(pik)) stop("y and pik have different sizes")
if(missing(type) | is.null(N)) {
                   if(missing(type)) warning("the type estimator is missing")
                   warning("by default the mean estimator is computed")
                   est<-crossprod(y,1/pik)/sum(1/pik)}
else if(type=="total") est<-N*crossprod(y,1/pik)/sum(1/pik) 
est
}





