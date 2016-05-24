HTstrata<-function (y, pik, strata, description=FALSE) 
{
    str <- function(st, h, n) .C("str", as.double(st), as.integer(h), 
        as.integer(n), s = double(n), PACKAGE = "sampling")$s
    if(any(is.na(pik))) stop("there are missing values in pik")
    if(any(is.na(y))) stop("there are missing values in y")
    if(length(y)!=length(pik)) stop("y and pik have different sizes")
    if (is.matrix(y)) 
        sample.size <- nrow(y) else sample.size <- length(y)
    h <- unique(strata)
    s1 <- 0
    for (i in 1:length(h)) {
        s <- str(strata, h[i], sample.size)
        est<-HTestimator(y[s == 1], pik[s == 1])
        s1 <- s1 + est
        if(description)
         cat("For stratum",i,",the Horvitz-Thompson estimator is:",est,"\n")
    }
if(description)
  cat("The Horvitz-Thompson estimator is:\n")
s1
}
