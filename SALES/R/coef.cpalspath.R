coef.cpalspath <- function(object, s = NULL, type = 
    c("coefficients", "nonzero"), ...) {
    type <- match.arg(type)
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- rbind2(b0, object$beta)
    t0 <- t(as.matrix(object$t0))
    rownames(t0) <- "(Intercept)"
    ntheta <- rbind2(t0, object$theta)
    if (!is.null(s)) {
      lambda <- object$lambda
      lamlist <- lambda.interp(lambda, s)
      vnames <- dimnames(nbeta)[[1]]
      dimnames(nbeta) <- list(NULL, NULL)
      nbeta <- nbeta[,lamlist$left,drop=FALSE]%*%Diagonal(x=lamlist$frac) +
               nbeta[,lamlist$right,drop=FALSE]%*%Diagonal(x=1-lamlist$frac)
      dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
      vnames <- dimnames(ntheta)[[1]]
      dimnames(ntheta) <- list(NULL, NULL)
      ntheta <- ntheta[,lamlist$left,drop=FALSE]%*%Diagonal(x=lamlist$frac) +
               ntheta[,lamlist$right,drop=FALSE]%*%Diagonal(x=1-lamlist$frac)
      dimnames(ntheta) <- list(vnames, paste(seq(along = s)))
    }
    if (type == "coefficients") return(list(beta = nbeta, theta = ntheta))
    if (type == "nonzero") # return locations of nonzero coefficients for each lambda
      return(list(nzbeta = nonzero(nbeta[-1, , drop = FALSE], bystep = TRUE),
        nztheta = nonzero(ntheta[-1, , drop = FALSE], bystep = TRUE)))
} 
