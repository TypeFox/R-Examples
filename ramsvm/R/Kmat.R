Kmat <- function(x, 
                 y,  
                 kernel = "gaussian",  
                 kparam = 1.0) {

  if( kernel == "polynomial" ) {
    obj <- (x %*% t(y) + 1.0)^kparam
  } else if( kernel == "gaussian" ) {
    normx <- drop((x^2) %*% rep(1.0, ncol(x)))
    normy <- drop((y^2) %*% rep(1.0, ncol(y)))
    temp <- x %*% t(y)
    temp <- (-2.0 * temp + normx) + outer(rep(1.0, nrow(x)), normy, "*")
    obj <- exp(-temp / (kparam^2))
  } else {
    obj <- NULL
  }

  return(obj)

}
