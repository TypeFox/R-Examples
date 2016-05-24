##
## Methods for estimating GO-GARCH models
## ==================================================
##
##
## Method definition for objects of class "Goestica"
## "Goestica" extends directly "GoGARCH"
##
setMethod(f = "goest", signature(object = "Goestica"), definition = function(object, initial, garchlist, ...){
  X <- object@X
  m <- ncol(X)
  n <- nrow(X)
  P <- object@P
  Id <- diag(m)
  Dsqr <- object@Dsqr
  ica <- fastICA(X, n.comp = m, ...)
  W <- ica$W
  Z <- P %*% Dsqr %*% t(P) %*% W
  Zinv <- solve(Z)
  Y <- X %*% Zinv
  fitted <- apply(Y, 2, function(x) do.call("garchFit", c(list(formula = object@garchf, data = quote(x)), garchlist)))
  H <- matrix(unlist(lapply(fitted, function(x) x@h.t)), ncol = m, nrow = n)
  Hdf <- data.frame(t(H))
  Ht <- lapply(Hdf, function(x) Z %*% diag(x) %*% t(Z))
  names(Ht) <- rownames(object@X)            
  result <- new("Goestica", ica = ica, estby = "fast ICA", U = W, Z = Z, Y = Y, H = Ht, models = fitted, X = object@X, P = object@P, Dsqr = object@Dsqr, V = object@V, garchf = object@garchf, name = object@name)
  return(result)  
})
##
## Method definition for objects of class "Goestmm"
## "Goestmm" extends directly "GoGARCH"
##
setMethod(f = "goest", signature(object = "Goestmm"), definition = function(object, lag.max, garchlist, ...){
  lag.max <- abs(as.integer(lag.max))
  X <- object@X
  m <- ncol(X)
  n <- nrow(X)
  P <- object@P
  Id <- diag(m)
  Dsqr <- object@Dsqr
  S <- P %*% Dsqr %*% t(P)
  Sinv <- solve(S)
  S <- X %*% Sinv
  if(lag.max < 1){
    U <- Id
    Umatched <- list(U)
    weights <- 1
  } else {
    SSI <- array(dim = c(m, m, n))
    for(i in 1:n){
      SSI[, , i] <- S[i, ] %*% t(S[i, ]) - diag(m)
    }
    Phil <- lapply(1:lag.max, function(x) cora(SSI, lag = x))
    evs <- lapply(Phil, function(x) eigen(x, symmetric = TRUE))
    evmin <- unlist(lapply(evs, function(x){
      sel <- combn(1:m, 2)
      diffs2 <- (x$values[sel[1, ]] - x$values[sel[2, ]])^2
      min(diffs2)
    }))
    denom <- sum(evmin)
    weights <- evmin / denom
    Ul <- lapply(evs, function(x) x$vectors)
    Ul[[1]] <- Umatch(Id, Ul[[1]])
    Sm <- matrix(0, nrow = m, ncol = m)
    for(i in 1:lag.max){
      Ul[[i]] <- Umatch(Ul[[1]], Ul[[i]])
      mmprod <- weights[i] * (Id - Ul[[i]]) %*% solve(Id + Ul[[i]])
      Sm <- Sm + mmprod
    }
    Umatched <- Ul
    U <- (Id - Sm) %*% solve(Id + Sm)
  }
  Y <- S %*% U
  Z <- P %*% Dsqr %*% t(P) %*% t(U)
  fitted <- apply(Y, 2, function(x) do.call("garchFit", c(list(formula = object@garchf, data = quote(x)), garchlist)))
  H <- matrix(unlist(lapply(fitted, function(x) x@h.t)), ncol = m, nrow = n)
  Hdf <- data.frame(t(H))
  Ht <- lapply(Hdf, function(x) Z %*% diag(x) %*% t(Z))
  names(Ht) <- rownames(object@X)          
  result <- new("Goestmm", weights = weights, Umatched = Umatched, estby = "Methods of Moments", U = U, Z = Z, Y = Y, H = Ht, models = fitted, X = object@X, P = object@P, Dsqr = object@Dsqr, V = object@V, garchf = object@garchf, name = object@name) 
  return(result)  
})
##
## Method definition for objects of class "Goestnls"
## "Goestnls" extends directly "GoGARCH"
##
setMethod(f = "goest", signature(object = "Goestnls"), definition = function(object, initial, garchlist, ...){
  d <- ncol(object@X)
  if(is.null(initial)){
    l <- d * (d + 1)/2
    initial <- rep(0.1, l)
  } else {
    l <- length(initial)
    if (l != d * (d + 1)/2) {
      stop(paste("\nLength of initial vector does not match length of vech(B).\n", "It should have length: ", d * (d + 1)/2, sep = ""))
    }
  }
  X <- object@X
  m <- ncol(X)
  n <- nrow(X)
  Dsqr <- object@Dsqr
  Dsqri <- diag(1 / diag(Dsqr))
  P <- object@P
  S <- X %*% P %*% Dsqri 
  SSI <- list()
  length(SSI) <- n
  for(i in 1:n){
    SSI[[i]] <- S[i, ] %*% t(S[i, ]) - diag(m)
  }
  SSI0 <- SSI[-1]
  SSI1 <- SSI[-n]
  SSI <- list(SSI0 = SSI0, SSI1 = SSI1)  
  nlsobj <- optim(par = initial, fn = gonls, SSI = SSI, ...)
  B <- unvech(nlsobj$par)
  U <- eigen(B)$vectors
  Z <- P %*% Dsqr %*% t(U)
  Y <- S %*% U
  fitted <- apply(Y, 2, function(x) do.call("garchFit", c(list(formula = object@garchf, data = quote(x)), garchlist)))
  H <- matrix(unlist(lapply(fitted, function(x) x@h.t)), ncol = m, nrow = n)
  Hdf <- data.frame(t(H))
  Ht <- lapply(Hdf, function(x) Z %*% diag(x) %*% t(Z))
  names(Ht) <- rownames(object@X)          
  result <- new("Goestnls", nls = nlsobj, estby = "non-linear Least-Squares", U = U, Z = Z, Y = Y, H = Ht, models = fitted, X = object@X, P = object@P, Dsqr = object@Dsqr, V = object@V, garchf = object@garchf, name = object@name) 
  return(result)  
})
##
## Method definition for objects of class "Goestml"
## "Goestml" extends directly "GoGARCH"
##
setMethod(f = "goest", signature(object = "Goestml"), definition = function(object, initial, garchlist, ...){
  d <- ncol(object@X)
  if(is.null(initial)){
    l <- d * (d - 1)/2
    initial <- seq(3.0, 0.1, length.out = l)
  } else {
    l <- length(initial)
    if (l != d * (d - 1)/2) {
      stop(paste("\nLength of initial vector does not match implied dimension of orthogonal matrix.\n", "It should have length: ", d * (d - 1)/2, sep = ""))
    }
  }
  llobj <- nlminb(start = initial, objective = gollh, object = object, garchlist = garchlist, lower = 1.5e-8, upper = pi/2, ...)
  gotheta <- gotheta(llobj$par, object, garchlist)
  result <- new("Goestml", opt = llobj, estby = "maximum likelihood", gotheta)
  return(result)  
})
