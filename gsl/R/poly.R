"gsl_poly" <- function(c_gsl, x){
  x.vec <- as.vector(x)
  attr <- attributes(x)
  jj <- .C("gsl_poly",
           as.double(c_gsl),
           as.integer(length(c_gsl)),
           as.double(x.vec),
           as.integer(length(x.vec)),
           ans=as.double(x.vec),
           PACKAGE="gsl"
           )
  ans <- jj$ans
  attributes(ans) <- attr
  return(ans)
}  
