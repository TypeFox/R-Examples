# arbitrary differentiable function

fSolver<-function(z, a, extra) {
    x <- z
    if ((is.null(extra$fobj)) || (is.null(extra$gobj))) stop("fSolver needs the additional arguments fobj and gobj!")
    fobj <- extra$fobj
    gobj <- extra$gobj
    n <- length(x)
    if (length(a)==0) { 
      indi <- diag(n)
    } else {
      if (is.vector(a)) a <- matrix(a,1,2)
      indi <- mkIndi(a,n)
    }
    z <- drop(crossprod(indi,x))
    p <- optim(z,
         fn=function(u) fobj(drop(indi%*%u)),
         gr=function(u) drop(crossprod(indi,gobj(drop(indi%*%u)))),
         method="BFGS")
    y <- drop(indi%*%(p$par))
    f <- p$value
    gy <- gobj(y)
    if (length(a)==0) lbd <- 0
        else lbd <- mkLagrange(a,gy)
    return(list(x = y, lbd = lbd, f = f, gx = gy))
}