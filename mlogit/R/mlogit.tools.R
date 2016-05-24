num.gradient <- function(f, param, ...){
  m <- match.call(expand.dots = TRUE)
  m[[1]] <- as.name(m[[2]])
  m[[2]] <- NULL
  eps <- 1E-10
  nc <- c()
  for (i in 1:length(param)){
    params <- param
    parami <- param
    params[i] <- params[i] + eps
    parami[i] <- parami[i] - eps
    m$param <- params
    lnls <- eval(m, parent.frame())
    m$param <- parami
    lnli <- eval(m, parent.frame())
    nc <- c(nc, (lnls-lnli)/(2*eps))
  }  
  nc 
}

num.hessian <- function(f, param, ...){
  m <- match.call(expand.dots = TRUE)
  m[[1]] <- as.name(m[[2]])
  m[[2]] <- NULL
  m$gradient <- TRUE
  eps <- 1E-8
  K <- length(param)
  nc <- c()
  for (i in 1:K){
    params <- param
    parami <- param
    params[i] <- params[i] + eps
    parami[i] <- parami[i] - eps
    m$param <- params
    gs <- attr(eval(m, parent.frame()), "gradient")
    m$param <- parami
    gi <- attr(eval(m, parent.frame()), "gradient")
    nc <- c(nc, (gs-gi)/(2*eps))
  }  
  matrix(nc, K, K)
}

print.est.stat <- function(x, ...){
  et <- x$elaps.time[3]
  i <- x$nb.iter[1]
  halton <- x$halton
  method <- x$method
  if (!is.null(x$type) && x$type != "simple"){
    R <- x$nb.draws
    cat(paste("Simulated maximum likelihood with", R, "draws\n"))
  }
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  cat(paste(method, "method\n"))
  tstr <- paste(h, "h:", m, "m:", s, "s", sep="")
  cat(paste(i,"iterations,",tstr,"\n"))
  if (!is.null(halton)) cat("Halton's sequences used\n")
  if (!is.null(x$eps)) cat(paste("g'(-H)^-1g =", sprintf("%5.3G", as.numeric(x$eps)),"\n"))
  if (is.numeric(x$code)){
    msg <- switch(x$code,
                  "1" = "gradient close to zero",
                  "2" = "successive function values within tolerance limits",
                  "3" = "last step couldn't find higher value",
                  "4" = "iteration limit exceeded"
                  )
    cat(paste(msg, "\n"))
  }
  else cat(paste(x$code, "\n"))
}

suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))){
    d <- dim(x[[1]])
    s <- matrix(0,d[1],d[2])
    for (i in 1:n){
      x[[i]][is.na(x[[i]])] <- 0
      s <- s+x[[i]]
    }
  }
  else{
    s <- rep(0,length(x[[n]]))
    for (i in 1:n){
      x[[i]][is.na(x[[i]])] <- 0
      s <- s+x[[i]]
    }
  }
  s
}

numderiv <- function(f, param, ...){
  m <- match.call()
  m[[1]] <- as.name(m[[2]])
  m[[2]] <- NULL
  eps <- 1E-4
  nc <- c()
  for (i in 1:length(param)){
    params <- param
    parami <- param
    params[i] <- params[i]+eps
    parami[i] <- parami[i]-eps
    m$param <- params
    lnls <- eval(m, parent.frame())
    m$param <- parami
    lnli <- eval(m, parent.frame())
    nc <- c(nc, (lnls-lnli)/(2*eps))
  }
  nc
}

###################################
##  Test the gradient
##   nd <- f
##   nd[["f"]] <- nd[[1]]
##   nd[[1]] <- as.name("numderiv")
##   x <- eval(nd, parent.frame())
##   print(x)
###################################


mlogit.optim <- function(logLik, start,
                         method = c('bfgs', 'nr', 'bhhh'),
                         iterlim = 2000,
                         tol = 1E-06,
                         ftol = 1E-08,
                         steptol = 1E-10,
                         print.level = 0,
                         constPar = NULL,
                         ...){
  method <- match.arg(method)
  param <- start 
  callT <- match.call(expand.dots = TRUE)
  optimoptions <- c('iterlim', 'tol', 'method', 'print.level', 'constPar', 'ftol', 'steptol')
  chi2 <- 1E+10
  i <- 0
  K <- length(param)
  d <- rep(0, K)

  # construct a vector of fixed parameters
  fixed <- rep(FALSE, K)
  names(fixed) <- names(start)
  if (!is.null(constPar)) fixed[constPar] <- TRUE

  # construct a call for the function
  f <- callT
  # if the model is updated based on a multinomial logit model, change the method to bfgs
  # if (method == 'nr' && as.character(f[[1]]) != "lnl.mlogits") method <- 'bfgs'
  m <- match(optimoptions, names(callT), 0L)
  if (sum(m)) f <- f[-m]
  f[[1]] <- as.name(f[[2]])
  f$gradient <- TRUE
#  f$steptol <- steptol
  f$stptol <- steptol
  if (method == 'nr') f$hessian <- TRUE else f$hessian <- FALSE
  f[[2]] <- NULL
  names(f)[2] <- 'param'
  # eval a first time the function, the gradient and the hessian
  x <- eval(f, parent.frame())

  if (FALSE){
    print(attr(x, "gradient"))
    nd <- f
    nd[["f"]] <- nd[[1]]
    nd[[1]] <- as.name("numderiv")
    nd$gradient <- FALSE
    print(nd)
    bla <- eval(nd, parent.frame())
    print(bla)
    cat("____________\n");
#    stop()
    }

  if (print.level > 0)
    cat(paste("Initial value of the function :", as.numeric(x), "\n"))
  g <- attr(x, "gradient")
  if (method == 'nr')   H <- attr(x, "hessian")[!fixed, !fixed]
  if (method == 'bhhh') H <- crossprod(attr(x, "gradi")[, !fixed])
  if (method == 'bfgs') Hm1 <- solve(crossprod(attr(x, "gradi")[, !fixed]))
#  if (method == 'bfgs') Hm1 <- diag(length(start[!fixed]))

  repeat{
    # save the previous values of the function and of the gradient
    oldx <- x
    oldg <- g

    # Compute the direction, ie d = H^-1 g
    
    # For the predict method, I don't want the solve
    if (iterlim > 0){
      if (method == "bfgs") d[!fixed] <- - as.vector(Hm1 %*% g[!fixed])
      else d[!fixed] <- - as.vector(solve(H, g[!fixed]))
    }
    i <- i + 1
    
    if (i > iterlim){
      # exit if the iteration limit is reached
      code <- 4
      break
    }

    # indicate in the call the previous parameters vector, the
    # direction and the value of the function
    f$param <- param
    f$direction <- d
    f$initial.value <- oldx

    # eval the function and compute the gradient and the hessian
    x <- eval(f, parent.frame())
    if (is.null(x)){
      # x is null if steptol is reached
      code = 3
      break
    }
    if (abs(x - oldx) < ftol){
      code = 2
      break
    }
    g <- attr(x, "gradient")
    step <- attr(x, "step")
    param[!fixed] <- param[!fixed] + step * d[!fixed]
    if (method == 'nr')   H <- attr(x, "hessian")[!fixed, !fixed]
    if (method == 'bhhh') H <-  crossprod(attr(x, "gradi")[, !fixed])
    if (method == 'bfgs'){
      incr <- step * d
      y <- g - oldg
      Hm1 <- Hm1 +
        outer( incr[!fixed], incr[!fixed]) *
          (sum(y[!fixed] * incr[!fixed]) +
           as.vector( t(y[!fixed]) %*% Hm1 %*% y[!fixed])) /
             sum(incr[!fixed] * y[!fixed])^2 -
               (Hm1 %*% outer(y[!fixed], incr[!fixed])
                + outer(incr[!fixed], y[!fixed]) %*% Hm1)/
                  sum(incr[!fixed] * y[!fixed])
    }

    # compute the quadratic form of the gradient
    chi2 <- -  crossprod(d[!fixed], oldg[!fixed])
    
    # print some informations about the iteration
    if (print.level > 0){
      chaine <- paste("iteration ",i,", step = ",step,
                      ", lnL = ",round(x,8),", chi2 = ",
                      round(chi2,8),"\n",sep="")
    cat(chaine)
    }
    if (print.level > 1){
      resdet <- rbind(param = param, gradient = g)
      print(round(resdet,3))
      cat("--------------------------------------------\n")
    }

    if (abs(chi2) < tol){
      # exit if the quadratic form of the gradient is small enough
      code = 1
      break
    }

  }
  if (code == 3) x <- oldx
  names(attr(x, 'gradient')) <- colnames(attr(x, 'gradi')) <- names(param)
  attr(x, "fixed") <- fixed
  est.stat = structure(list(elaps.time = NULL, nb.iter = i, eps = chi2,
    method = method, code = code), class = 'est.stat')
  result <- list(optimum = x,
                 coefficients = param,
                 est.stat = est.stat
                 )
  result
}
