ksdrift <- function(x,bw,n=512){
 len <- length(x)
 xval <- seq(min(x), max(x), length=n)
 if(missing(bw))
  bw <- len^(-1/5)*sd(x)
  y <- sapply(xval, function(xval) { 
  tmp <- dnorm(xval, x[1:(len-1)], bw)
  sum(tmp * diff(x)) / (deltat(x) * sum(tmp))})
  invisible(list(x=xval, y=y))
}
 
ksdiff <- function(x,bw,n=512){
 len <- length(x)
 xval <- seq(min(x), max(x), length=n)
 if(missing(bw))
  bw <- len^(-1/5)*sd(x)
  y <- sapply(xval, function(xval) { 
  tmp <- dnorm(xval, x[1:(len-1)], bw)
  sum(tmp * as.numeric(diff(x))^2) / (deltat(x) * sum(tmp))})
  invisible(list(x=xval, y=sqrt(y)))
}


ksdens <- function(x,bw,n=512){
   len <- length(x)
   if(missing(bw))
    bw <- len^(-1/5)*sd(x)
   invisible(density(x,bw=bw,n=n))
}
