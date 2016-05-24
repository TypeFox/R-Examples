Hajekstrata<-function(y,pik,strata,N=NULL,type=c("total","mean"),description=FALSE) 
{
    str <- function(st, h, n) .C("str", as.double(st), as.integer(h), 
        as.integer(n), s = double(n), PACKAGE = "sampling")$s
    if(any(is.na(pik))) stop("there are missing values in pik")
    if(any(is.na(y))) stop("there are missing values in y")
    if(length(y)!=length(pik)) stop("y and pik have different sizes")
    if (is.matrix(y)) 
        sample.size <- nrow(y)
    else sample.size <- length(y)
    if(!is.vector(N)) N <- as.vector(N)
    h <- unique(strata)
    if(length(N)!=length(h)) stop("N should be a vector with the length equal to the number of strata")
    options(warn=-1)
    s1 <- 0
    for (i in 1:length(h)) {
        s <- str(strata, h[i], sample.size)
        est <- Hajekestimator(y[s == 1], pik[s == 1], type="mean")
        s1 <- s1 + est*N[i]
        if(description)
          if(type=="mean")
            cat("For stratum ",i,", the Hajek estimator is:",est,"\n")
          else cat("For stratum ",i,", the Hajek estimator is:",est*N[i],"\n")
             }
if(description)
  cat("The Hajek estimator is:\n")
if(type=="mean")
 return(s1/sum(N))
else return(s1)
}



