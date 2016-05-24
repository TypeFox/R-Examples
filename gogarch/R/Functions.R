##
## This file includes the following functions:
## ===========================================
##
## gogarch
## Umatch
## UprodR
## Rd2
## cora
## goinit
## gollh
## gonls
## gotheta
## unvech
##
## ============================================
##
##
## gogarch: main function for estimating GO-GARCH models
##
gogarch <-
function(data, formula, scale = FALSE, estby = c("ica", "mm", "ml", "nls"), lag.max = 1, initial = NULL, garchlist = list(init.rec = "mci", delta = 2, skew = 1, shape = 4, cond.dist = "norm", include.mean = FALSE, include.delta = NULL, include.skew = NULL, include.shape = NULL, leverage = NULL, trace = FALSE, algorithm = "nlminb", hessian = "ropt", control = list(), title = NULL, description = NULL), ...){
  estby <- match.arg(estby)
  Call <- match.call()
  gini <- goinit(X = data, garchf = formula, scale = scale)
  gomod <- new("GoGARCH", gini)
  if(estby == "ml"){
    goestml <- new("Goestml", gomod)
    gogarch <- goest(object = goestml, initial = initial, garchlist = garchlist, ...)
  }
  if(estby == "nls"){
    goestnls <- new("Goestnls", gomod)
    gogarch <- goest(object = goestnls, initial = initial, garchlist = garchlist, ...)
  }
  if(estby == "mm"){
    goestmm <- new("Goestmm", gomod)
    gogarch <- goest(object = goestmm, lag.max = lag.max, garchlist = garchlist, ...)
  }  
  if(estby == "ica"){
    goestica <- new("Goestica", gomod)
    gogarch <- goest(object = goestica, initial = initial, garchlist = garchlist, ...)
  }
  gogarch@CALL <- Call
  gogarch@name <- deparse(substitute(data))
  return(gogarch)
}
##
## Umatch: Matching of orthogonal matrices. This function is employed  
## whence GO-GARCH models are estimated by methods of moments
##
Umatch <-
function(from, to){
  cols <- ncol(from)
  mat <- matrix(0, nrow = cols, ncol = cols)
  for(i in 1:cols){
    inner <- abs(colSums(from[, i] * to))
    maxcol <- which.max(inner)
    mat[, i] <- to[, maxcol]
    to <- as.matrix(to[, -c(maxcol)])
  }
  signs <- matrix(sign(diag(mat)), nrow = cols, ncol = cols, byrow = TRUE)
  mat <-  signs * mat
  if(det(mat) < 0.0){
    colminus <- which.min(abs(colSums(from * mat)))
    mat[, colminus] <- -1.0 * mat[, colminus]
  }
  return(mat)
}
##
## UprodR: This function computes an orthogonal matrix as the product of two-dimensional rotation matrices.
##
UprodR <-
function(theta){
  theta <- as.vector(theta)
  l <- length(theta)
  d <- as.integer(0.5 + sqrt(0.5^2 + 2*l))  
  if(l != d * (d - 1) / 2){
    stop("\nLength of theta does not match implied dimension of U.\n")
  }
  Id <- diag(d)
  U <- Id
  rc <- combn(x = d, m = 2)
  idx <- seq(along.with = theta)
  Rs <- lapply(idx, function(x){
    tmp <- Id
    tmp[rc[, x], rc[, x]] <- Rd2(theta = theta[x])
    return(tmp)
  })
  for(i in 1:l) U <- U %*% Rs[[i]]
  result <- new("Orthom", M = U)
  return(result)
}
##
## Rd2: This function returns a two-dimensional rotation matrix for a given Euler angle.
##
Rd2 <-
function(theta){
  theta <- as.vector(theta)
  if(length(theta) > 1){
    stop("\nLength of argument 'theta' should be one.\n")
    }
  if((theta <= 0) | (theta > pi/2)){
    stop("\nTheta should be in the interval [0, pi/2).\n")
  }
  R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol = 2, nrow = 2)
  return(R)              
}
##
## cora: Computation of autocorrelations/autocovariances of a matrix process.
## This function is utilized whence a GO-GARCH model is estimated by the
## methods of moments
##
cora <-
function(SSI, lag = 1, standardize = TRUE){
  lag <- abs(as.integer(lag))
  dims <- dim(SSI)
  Gamma <- matrix(0, nrow = dims[1], ncol = dims[2])
  SSIp <- array(dim = dims)
  for(i in 1:dims[3]){
    SSIp[, ,i] <- SSI[, ,i] %*% SSI[, ,i]
    Gamma <- Gamma + SSIp[, , i]
  }
  Gamma <- Gamma / dims[3]
  Gsvd <- svd(Gamma)
  Gsqrtinv <- Gsvd$u %*% diag(1/sqrt(Gsvd$d)) %*% t(Gsvd$u)
  idx <- 1:dims[3]
  if(identical(lag, as.integer(0))){
    idx1 <- idx
    idx2 <- idx
  } else {
    idx1 <- idx[-c(1:lag)]
    idx2 <- rev(rev(idx)[-c(1:lag)])
  }
  nl <- length(idx1)
  Gamma <- matrix(0, nrow = dims[1], ncol = dims[2])
  SSIc <- array(dim = c(dims[1], dims[2], nl))
  for(i in 1:nl){
    SSIc[, , i] <- SSI[, , idx1[i]] %*% SSI[, , idx2[i]]
    Gamma <- Gamma + SSIc[, , i]
  }
  Gamma <- Gamma / nl
  if(standardize){
    cora <- Gsqrtinv %*% Gamma %*% Gsqrtinv
  } else {
    cora <- Gamma
  }
  cora <- (cora + t(cora)) / 2
  return(cora)
}
##
## goinit: Function for creating an object of class "Goinit"
##
goinit <-
function(X, garchf = ~ garch(1, 1), scale = FALSE){
  dname <- deparse(substitute(X))
  X <- as.matrix(X)
  if(ncol(X) > nrow(X)){
    stop("\nMatrix has more columns than rows.\n")
  }
  garchf <- as.formula(garchf)
  if(scale){
    X <- scale(X)
  }
  V <- t(X) %*% X / nrow(X)
  svd <- svd(V)
  P <- svd$u
  Dsqr <- diag(sqrt(svd$d))
  result <- new("Goinit", X = X, V = V, P = P, Dsqr = Dsqr, garchf = garchf, name = dname)
  return(result)
}
##
## gollh: The log-likelihood function of GO-GARCH models.
## This function is employed whence GO-GARCH models are estimated by Maximum Likelihood
##
gollh <-
function(params, object, garchlist){
  gotheta <- gotheta(theta = params, object = object, garchlist = garchlist)
  m <- ncol(object@X)
  n <- nrow(object@X)
  H <- matrix(unlist(lapply(gotheta@models, function(x) x@h.t)), ncol = m, nrow = n)
  Hinv <- 1.0 / H
  arg1 <- n * m * log(2 * pi)
  arg2 <- log(det(gotheta@Z %*% t(gotheta@Z))) * n
  arg3 <- sum(log(apply(H, 1, prod)))
  arg4 <- sum(rowSums(gotheta@Y * Hinv * gotheta@Y))
  ll <- -0.5 * (arg1 + arg2 + arg3 + arg4)
  negll <- -1.0 * ll
  return(negll)
}
##
## gonls: The target function to be minimized whence GO-GARCH models are estimated
## by non-linear least-squares.
##
gonls <-
function(params, SSI){
  B <- unvech(params)
  n <- length(SSI[[1]])
  fl <- list()
  length(fl) <- n
  for(i in 1:n){
    M <- (SSI[[1]][[i]] - B %*% SSI[[2]][[i]] %*% B)
    fl[[i]] <- M %*% M
  }
  f <- sum(unlist(lapply(fl, function(x) sum(diag(x))))) / n
  return(f)   
}
##
## gotheta: For a given vector of Euler angles, this function computes a GO-GARCH model.
## The function is called during estimation of GO-GARCH models by maximum likelihood.
##
gotheta <-
function(theta, object, garchlist = list(init.rec = "mci", delta = 2, skew = 1, shape = 4, cond.dist = "norm", include.mean = FALSE, include.delta = NULL, include.skew = NULL, include.shape = NULL, leverage = NULL, trace = FALSE, algorithm = "nlminb", hessian = "ropt", control = list(), title = NULL, description = NULL)){
  if(!any(inherits(object, what = c("Goinit", "GoGARCH", "Goestml")))) {
    stop("\nObject is neither of class 'Goinit', 'GoGARCH' or 'Goestml'.\n")
  }
  l <- length(theta)
  d <- as.integer(0.5 + sqrt(0.5^2 + 2 * l))
  if (l != d * (d - 1)/2) {
    stop(paste("\nLength of theta does not match implied dimension of orthogonal matrix.\n", "It should have length: ", d, sep = ""))
  }
  m <- ncol(object@X)
  n <- nrow(object@X)
  U <- UprodR(theta)@M
  Z <- object@P %*% object@Dsqr %*% t(U)
  Zinv <- solve(Z)
  Y <- object@X %*% Zinv
  fitted <- apply(Y, 2, function(x) do.call("garchFit", c(list(formula = object@garchf, data = quote(x)), garchlist)))
  H <- matrix(unlist(lapply(fitted, function(x) x@h.t)), ncol = m, nrow = n)
  Hdf <- data.frame(t(H))
  Ht <- lapply(Hdf, function(x) Z %*% diag(x) %*% t(Z))
  names(Ht) <- rownames(object@X)
  result <- new("GoGARCH", U = U, Z = Z, Y = Y, H = Ht, models = fitted, X = object@X, P = object@P, Dsqr = object@Dsqr, V = object@V, garchf = object@garchf, name = object@name, CALL = match.call())
  return(result)
}
##
## unvech: Reverts the vech-operator and returns a symmetric matrix
##
unvech <-
function(v){
  v <- as.vector(v)
  l <- length(v)
  n <- -(1 - sqrt(1 + 8 * l)) / 2
  if(n %% 1 != 0.0){
    stop("\nCannot produce symmetric matrix, check length of v.\n")
  }
  X <- matrix(0, ncol = n, nrow = n)
  X[lower.tri(X, diag = TRUE)] <- v
  Z <- X + t(X)
  diag(Z) <- diag(X)
 return(Z)
}
