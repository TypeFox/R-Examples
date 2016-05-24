replication <-
function(object, repperdr = 1, ...){ 
     subset = object$subset
     innov = object$innov
     m <- object$m
     theta <- as.matrix(object[[1]][,1:m])
     n <- nrow(theta)
     m <- ncol(theta)
     p <- object$p
     if (repperdr != as.integer(repperdr) || repperdr <= 0)
          stop("replication per draw must be positive integers")
     theta <- array(apply(theta, 2, rep, repperdr), c(m, n * repperdr))
     eps <- array(rnorm(m * n * repperdr), c(m, n * repperdr))
     if (innov == "normal"){
          Z <- object$Z
          theta <- theta + Z * eps  
     }
     else{
          sig2 <- as.matrix(object[[1]][,(m+p+2):(2*m+p+1)])
          sig2 <- array(apply(sig2, 2, rep, repperdr), c(m, n * repperdr))
          theta <- theta + sig2 * eps
     }
     theta
}
