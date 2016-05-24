predict.alspath <- function(object, newx, s = NULL, 
    type = c("response"), ...) {
    type <- match.arg(type)
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- rbind2(b0, object$beta)
    if (!is.null(s)) {
      vnames <- dimnames(nbeta)[[1]]
      dimnames(nbeta) <- list(NULL, NULL)
      lambda <- object$lambda
      lamlist <- lambda.interp(lambda, s)
      nbeta <- nbeta[,lamlist$left,drop=FALSE]%*%Diagonal(x=lamlist$frac) +
               nbeta[,lamlist$right,drop=FALSE]%*%Diagonal(x=1-lamlist$frac)
      dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    nfit <- as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
    nfit
} 
