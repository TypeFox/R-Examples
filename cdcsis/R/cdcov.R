#' Perform the measurement of the conditional independence between x and y given z.
#'
#' @title cdcov
#' @name cdcov
#' @param x a numeric vector (the name of a variable).
#' @param y a numeric vector (the name of a variable).
#' @param z a numeric vector (the names of the conditioning variables).
#' @param width bandwidth about z.
#' @param index a parametric of distance.
#' @usage cdcov(x,y,z,width,index=1)
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' z <- x + y
#' result=cdcov(x,y,z,0.25,1)
#' 

cdcov <-
function(x,y,z,width,index=1) {
  x<-as.matrix(x);y<-as.matrix(y);z<-as.matrix(z)
  dim_x<-dim(x); n<-dim_x[1]; p<-dim_x[2]; q<-dim(y)[2];d<-dim(z)[2]
  k<-numeric(n*n)
  CDCOV<-numeric(n)
  re <-.C("cdCOV", as.double(t(x)), as.double(t(y)), as.double(t(z)), as.integer(n), 
        as.integer(p), as.integer(q), as.integer(d), as.double(index),as.double(width), k=as.double(k),
        cd=as.double(CDCOV)) 
  cdc <- list(cdcov=re$cd, width=width, index=index)
  return(cdc)
}
