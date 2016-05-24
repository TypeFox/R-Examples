.complexLog <-
function (x) {
  if ( is.double(x) & x>0 ) {
      y <- log(x)
  } else {
      if ( is.double(x) & x<0 )
          warning("Log of negative real number, using complex log!")
      y <- log(x+0i)
  }
  return ( y )
}
.dist2 <-
function (x, x2) {
  xdim <- dim(as.matrix(x))
  x2dim <- dim(as.matrix(x2))

  xMat <- array(apply(as.matrix(x*x),1,sum), c(xdim[1], x2dim[1]))
  x2Mat <- t(array(apply(as.matrix(x2*x2),1,sum), c(x2dim[1], xdim[1])))

  if ( xdim[2] != x2dim[2] )
    stop("Data dimensions are not matched.")

  n2 <-   xMat+x2Mat-2*tcrossprod(x, x2)

  return (n2)
}
.distfit <-
function(data, dist = "normal") {
  if (dist == "gamma") {
    cdf <- qgamma 
  }

  else if (dist == "normal") {
    cdf <- qnorm
  }

  else {
    stop("Unknown distribution.")
  }

  t <- optim(c(1, 1), fn=.distfit_obj, gr=NULL, data, cdf)

  return (t)
}
.distfit_obj <-
function(theta, y, cdf) {
  p <- c(.05, .25, .50, .75, .95)
  x <- cdf(p, theta[1], theta[2])
  r <- .5 * sum((x - y)^2)

  return (r)
}
.fn_line <-
function (linemin, fun, para0, direction, ...) {
  ## y = fn (x)
  func <- function(x, ...) fun(x, ...)
  
  para <- para0 + linemin * direction
  
  ans <- func(para, ...)
  
  return (ans)
}
.gradLnDiffErfs <-
function(x1, x2, fact1, fact2) {
  m <- pmin(as.matrix(x1)^2, as.matrix(x2)^2)
  dlnPart <- 2/sqrt(pi) * (exp(-x1^2 + m) * fact1 - exp(-x2^2 + m) * fact2)

  g <- list(dlnPart=dlnPart, m=m)
  return (g)

}
.jitChol <-
function ( M, Num=10, silent=FALSE ) {
  jitter <- 0
  jitter1 <- abs(mean(diag(M)))*1e-6
  eyeM <- diag( 1, nrow=length(M[,1]), ncol=length(M[1,]) )

  for ( i in 1:Num ) {
    ## clear the last error message
    try(stop(""),TRUE)

    Ch <- try( chol( M + jitter*eyeM ), silent=TRUE )

    nPos <- grep("not positive definite",  geterrmessage())

    if ( length(nPos) != 0 ) {
      jitter1 <- jitter1*10
      jitter <- jitter1

      if (! silent) {
        warnmsg <- paste("Matrix is not positive definite, adding",
                         signif(jitter,digits=4), "jitter!")
        warning(warnmsg)
      }
    }
    else break
  }

  return (list(chol=Ch, jitter=jitter))
}
.jitCholInv <-
function ( M, Num=10, silent=FALSE ) {
  jitter <- 0
  jitter1 <- abs(mean(diag(M)))*1e-6
  eyeM <- diag( 1, nrow=length(M[,1]), ncol=length(M[1,]) )

  for ( i in 1:Num ) {

    ## clear the last error message
    try(stop(""),TRUE)

    Ch <- try( chol( M + jitter*eyeM ), silent=TRUE )

    nPos <- grep("not positive definite",  geterrmessage())

    if ( length(nPos) != 0 ) {
      jitter1 <- jitter1*10
      jitter <- jitter1

      if (! silent) {
        warnmsg <- paste("Matrix is not positive definite, adding",
                         signif(jitter,digits=4), "jitter!")
        warning(warnmsg)
      }
    }
    else break
  }

  invCh <- try (solve( Ch, eyeM ), silent=TRUE)

  if ( class(invCh) == "try-error" ) {
    return (NaN)
  }
  else {
    invM <- invCh %*% t(invCh)

    if ( jitter == 0 ) {
      ans <- list(invM=invM, jitter=jitter, chol=Ch)
    }
    else ans <- list(invM=invM, jitM=M+jitter*eyeM , jitter=jitter, chol=Ch)

    return (ans)
  }
}
.kernFactors <-
function (kern, factorType) {
  factors <- list()

  if ( length(kern$transforms) > 0 ) {
    funcName <- paste(kern$type, "KernExtractParam", sep="")
    func <- get(funcName, mode="function")
    params <- func(kern)

    for (i in seq(along=kern$transforms)) {
      factors[[i]] <- list()
      factors[[i]]$index <- kern$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(kern$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        factors[[i]]$val <- func(params[factors[[i]]$index], factorType, kern$transformArgs[[i]])
      else
        factors[[i]]$val <- func(params[factors[[i]]$index], factorType)
    }
  }
  return (factors)
}
.kernTestCombinationFunction <-
function (kern1, kern2) {
  if (kern1$type == "selproj" && kern2$type == "selproj")
    funcName <- paste(kern1$comp[[1]]$type, "X", kern2$comp[[1]]$type, "KernCompute", sep="")
  else
    funcName <- paste(kern1$type, "X", kern2$type, "KernCompute", sep="")

  if ( !exists(funcName, mode="function") ) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}
.multiKernCacheBlock <-
function(kern, fhandle, i, j, x1, x2=NULL, arg1, arg2=NULL) {

  global_cache <- get("cache", envir = kern$cache)
  if ((length(global_cache) >= i) && (length(global_cache[[i]]) >= j))
    cache <- global_cache[[i]][[j]]
  else
    cache <- list()
  key <- c(x1, x2)

  for (k in seq(along=cache)) {
    if (length(key) == length(cache[[k]]$key) && all(key == cache[[k]]$key)) {
      #cat("multiKernCacheBlock: cache hit ", i, j, x1, x2, "\n")
      return (cache[[k]]$value)
    }
  }

  #cat("multiKernCacheBlock: cache miss ", i, j, x1, x2, "\n")
  # No match if we get all the way here
  if (is.null(arg2)) {
    if (is.null(x2))
      K <- fhandle(arg1, x1)
    else
      K <- fhandle(arg1, x1, x2)
  }
  else {
    if (is.null(x2))
      K <- fhandle(arg1, arg2, x1)
    else
      K <- fhandle(arg1, arg2, x1, x2)
  }
  cache <- append(cache, list(list(key=key, value=K)))
  if (length(global_cache) < i)
    global_cache[[i]] <- list()
  global_cache[[i]][[j]] <- cache
  assign("cache", global_cache, envir=kern$cache)
  return(K)
}
.multiKernComputeBlock <-
function (kern, i, j, x1, x2=NULL) {
  if ( i==j ) {
    funcName <- paste(kern$comp[[i]]$type, "KernCompute", sep="")
    transpose <- 0
    arg1 <- kern$comp[[i]]

    func <- get(funcName, mode="function")
    if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
      K <- .multiKernCacheBlock(kern, func, i, j, arg1=arg1, x1=x1, x2=x2)
    }
    else {
      if (is.null(x2))
        K <- func(arg1, x1)
      else
        K <- func(arg1, x1, x2)
    }
  } else {

    if ( j<i ) {
      funcName <- paste(kern$block[[i]]$cross[j], "KernCompute", sep="")
      transpose <- kern$block[[i]]$transpose[j]
    } else {
      funcName <- paste(kern$block[[j]]$cross[i], "KernCompute", sep="")
      transpose <- !kern$block[[j]]$transpose[i]
    }

    if ( transpose ) {
      arg1 <- kern$comp[[j]]
      arg2 <- kern$comp[[i]]
    } else {
      arg1 <- kern$comp[[i]]
      arg2 <- kern$comp[[j]]      
    }
    
    func <- get(funcName, mode="function")
    if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
      K <- .multiKernCacheBlock(kern, func, i, j, arg1=arg1, arg2=arg2, x1=x1, x2=x2)
    }
    else {
      if (is.null(x2))
        K <- func(arg1, arg2, x1)
      else
        K <- func(arg1, arg2, x1, x2)
    }
  }
  return (K)
}
.multiKernGradientBlock <-
function (kern, x, x2, covGrad, i, j) {
  if ( nargs()<6 ) {
    j <- i
    i <- covGrad
    covGrad <- x2
    x2 <- array()
  }

  if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
    if (i==j)
      return (0)
    else
      return (list(g1=0, g2=0))
  }

  if ( i==j ) {
    funcName <- paste(kern$comp[[i]]$type, "KernGradient", sep="")
    transpose <- 0
    arg1 <- kern$comp[[i]]
    factors <- .kernFactors(kern$comp[[i]], "gradfact")

    func <- get(funcName, mode="function")

    if ( is.na(x2) ) {
      g <- func(arg1, x, covGrad)
    } else {
      g <- func(arg1, x, x2, covGrad)
    }
    for (i in seq(along=factors))
      g[factors[[i]]$index] <- g[factors[[i]]$index]*factors[[i]]$val
    
  } else {
    if ( j<i ) {
      funcName <- paste(kern$block[[i]]$cross[j], "KernGradient", sep="")
      transpose <- kern$block[[i]]$transpose[j]
    } else {
      funcName <- paste(kern$block[[j]]$cross[i], "KernGradient", sep="")
      transpose <- kern$block[[j]]$transpose[i]
    }

    if ( transpose ) {
      arg1 <- kern$comp[[j]]
      factors1 <- .kernFactors(kern$comp[[j]], "gradfact")
      arg2 <- kern$comp[[i]]
      factors2 <- .kernFactors(kern$comp[[i]], "gradfact")
    } else {
      arg1 <- kern$comp[[i]]
      factors1 <- .kernFactors(kern$comp[[i]], "gradfact")      
      arg2 <- kern$comp[[j]]
      factors2 <- .kernFactors(kern$comp[[j]], "gradfact")
    }

    func <- get(funcName, mode="function")
    if ( any(is.na(x2)) ) {
      gList <- func(arg1, arg2, x, covGrad)
    } else {
      gList <- func(arg1, arg2, x, x2, covGrad)
    }

    g1 <- gList$g1
    g2 <- gList$g2
    
    for (i in seq(along=factors1))
      g1[factors1[[i]]$index] <- g1[factors1[[i]]$index]*factors1[[i]]$val

    for (i in seq(along=factors2))
      g2[factors2[[i]]$index] <- g2[factors2[[i]]$index]*factors2[[i]]$val

    if ( transpose ) {
      g <- g2
      g2 <- g1
      g1 <- g
    }
    g <- list(g1=g1, g2=g2)   
    
  }
  return (g)
}
