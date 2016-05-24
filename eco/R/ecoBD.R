ecoBD <- function(formula, data = parent.frame(), N=NULL){
  mf <- match.call()
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  vnames <- attr(tt, "variables")
  vnamesR <- vnames[[2]]
  
  if (is.matrix(eval.parent(mf$data)))
    data <- as.data.frame(data)
  X <- model.matrix(tt, data)
  Y <- as.matrix(model.response(model.frame(tt, data = data)))
  N <- eval(mf$N, data)
  n.obs <- nrow(X)

  ## counts
  if (all(X>1) & all(Y>1)) {
    if (!is.null(N)) {
      if (!all(apply(X, 1, sum) == N))
        X <- cbind(X, N-apply(X, 1, sum))
      if (!all(apply(Y, 1, sum) == N))
        Y <- cbind(Y, N-apply(Y, 1, sum))
      if(any(X<0) || any(Y<0))
        stop("Invalid inputs for X, Y, or/and N")
    }
    else {
      if (!all(apply(X, 1, sum) == apply(Y, 1, sum)))
        stop("X and Y do not sum to the same number. Input N.")
      N <- apply(X, 1, sum)
    }
    C <- ncol(X)
    R <- ncol(Y)
    Wmin <- Wmax <- Nmin <- Nmax <- array(NA, c(n.obs, R, C))
    clab <- rlab <- NULL
    if (length(vnames) == 3)
      clab <- c(vnames[[3]], paste("not",vnames[[3]]))
    else {
      for (j in 1:C) {
        if ((j == C) & (length(vnames) < j+2))
          clab <- c(clab, "other")
        else
          clab <- c(clab, vnames[[j+2]])
      }
    }
    if (length(vnamesR) == 1)
      rlab <- c(vnamesR, paste("not",vnamesR))
    else {
      for (i in 1:R) {
        if ((i == R) & (length(vnamesR) < i+1))
          rlab <- c(rlab, "other")
        else
          rlab <- c(rlab, vnamesR[[i]])
      }
    }
    for (i in 1:R) {
      for (j in 1:C) {
        Nmin[,i,j] <- apply(cbind(0, X[,j]+Y[,i]-N), 1, max)
        Nmax[,i,j] <- apply(cbind(Y[,i], X[,j]), 1, min)
        Wmin[,i,j] <- Nmin[,i,j]/X[,j]
        Wmax[,i,j] <- Nmax[,i,j]/X[,j]
      }
    }
    dimnames(Wmin) <- dimnames(Wmax) <- dimnames(Nmin) <-
      dimnames(Nmax) <-
        list(if (is.null(rownames(X))) 1:n.obs else rownames(X),
             rlab, clab)
  }
  else { ## proportions
    if (any(apply(X, 1, sum) > 1.000000001))
      stop("invalid input for X")
    if (any(apply(X, 1, sum) < 0.9999999999))
      X <- cbind(X, 1-X)
    if (any(apply(Y, 1, sum) > 1.0000000001))
      stop("invalid input for Y")
    if (any(apply(Y, 1, sum) < 0.9999999999))
      Y <- cbind(Y, 1-Y)
    C <- ncol(X)
    R <- ncol(Y)
    Wmin <- Wmax <- array(NA, c(n.obs, R, C))
    clab <- rlab <- NULL
    if (length(vnames) == 3)
      clab <- c(vnames[[3]], paste("not",vnames[[3]]))
    else {
      for (j in 1:C) {
        if ((j == C) & (length(vnames) < j+2))
          clab <- c(clab, "other")
        else
          clab <- c(clab, vnames[[j+2]])
      }
    }
    if (length(vnamesR) == 1)
      rlab <- c(vnamesR, paste("not",vnamesR))
    else {
      for (i in 1:R) {
        if ((i == R) & (length(vnamesR) < i+1))
          rlab <- c(rlab, "other")
        else
          rlab <- c(rlab, vnamesR[[i]])
      }
    }
    for (i in 1:R) {
      for (j in 1:C) {
        Wmin[,i,j] <- apply(cbind(0, (X[,j]+Y[,i]-1)/X[,j]), 1, max)
        Wmax[,i,j] <- apply(cbind(1, Y[,i]/X[,j]), 1, min)
      }
    }
    dimnames(Wmin) <- dimnames(Wmax) <-
      list(if (is.null(rownames(X))) 1:n.obs else rownames(X),
           rlab, clab)
    colnames(X) <- clab
    colnames(Y) <- rlab
    if (!is.null(N)) {
      Nmin <- Nmax <- array(NA, c(n.obs, R, C), dimnames =
                            dimnames(Wmin))
      for (i in 1:R) 
        for (j in 1:C) {
          Nmin[,i,j] <- Wmin[,i,j]*X[,j]*N
          Nmax[,i,j] <- Wmax[,i,j]*X[,j]*N
        }
    }
    else
      Nmin <- Nmax <- NULL
  }

  ## aggregate bounds
  aggWmin <- aggWmax <- matrix(NA, R, C, dimnames =
                               list(dimnames(Wmin)[[2]], dimnames(Wmin)[[3]]))
  if (is.null(N))
    for (j in 1:C) {
      aggWmin[,j] <- apply(Wmin[,,j], 2, weighted.mean, X[,j])
      aggWmax[,j] <- apply(Wmax[,,j], 2, weighted.mean, X[,j])
    }
  else
    for (j in 1:C) {
      aggWmin[,j] <- apply(Wmin[,,j], 2, weighted.mean, X[,j]*N)
      aggWmax[,j] <- apply(Wmax[,,j], 2, weighted.mean, X[,j]*N)
    }

  if (!is.null(Nmin) & !is.null(Nmax)) {
    aggNmin <- aggNmax <- matrix(NA, R, C, dimnames =
                                 list(dimnames(Nmin)[[2]], dimnames(Nmin)[[3]]))
    for (j in 1:C) {
      aggNmin[,j] <- apply(Nmin[,,j], 2, sum)
      aggNmax[,j] <- apply(Nmax[,,j], 2, sum)
    }
  }
  else
    aggNmin <- aggNmax <- NULL
    
  ## output
  res <- list(call = mf, X = X, Y = Y, N = N, aggWmin = aggWmin,
              aggWmax = aggWmax, aggNmin = aggNmin, aggNmax = aggNmax,
              Wmin = Wmin, Wmax = Wmax, Nmin = Nmin, Nmax = Nmax)
  class(res) <- c("ecoBD", "eco")
  return(res)
}
