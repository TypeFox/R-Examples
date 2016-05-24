##############################
# Make Draws and Coefficients
##############################

## Halton function: copied from mlogit package
halton <- function(prime = 3, length = 100, drop = 10){
  halt <- 0
  t <- 0
  while(length(halt) < length + drop){
    t <- t + 1
    halt <- c(halt, rep(halt, prime - 1) + rep(seq(1, prime - 1, 1) / prime ^ t, each = length(halt)))
  }
  halt[(drop + 1):(length + drop)]
}

## Make draw function. Modified from mlogit package
make.draws <- function(R, Ka, haltons){
  # Create the matrix of random numbers
  if (!is.null(haltons)){
    length.haltons <- rep(R,Ka)
    prime <- c(3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
               47, 53, 59, 61, 71, 73, 79, 83, 89, 97, 101, 103,
               107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
               173, 179, 181, 191, 193, 197, 199)
    drop.haltons <- rep(100, Ka)
    if (!is.na(haltons) && !is.null(haltons$prime)){
      if (length(haltons$prime) != Ka){
        stop("wrong number of prime numbers indicated")
      }
      else{
        prime <- haltons$prime
      }
      if (!is.na(haltons) && !is.null(haltons$drop)){
        if (!length(haltons$drop) %in% c(1,Ka)) stop("wrong number of drop indicated")
        if (length(haltons$drop) == 1){
          drop.haltons <- rep(haltons$drop, Ka)
        }
        else{
          drop.haltons <- haltons$drop
        }
      }
    }
    random.nb <- numeric(0)
    i <- 0
    for (i in 1:Ka){
      random.nb <- cbind(random.nb,qnorm(halton(prime[i],R,drop.haltons[i])))
    }
  }
  else{
    random.nb <- matrix(rnorm(R*Ka), ncol = Ka, nrow = R)
  }
  t(random.nb)  
}

## Make Lower Triangular
makeL <- function(x){
  K <- (-1+sqrt(1 + 8 * length(x)))/2
  mat <- matrix(0, K, K)
  mat[lower.tri(mat, diag = TRUE)] <- x
  mat
}

## Make Random Coefficients Version 0.2: include mvar list, include Johnson Sb
## distribution
Make.rcoef <- function(beta, sigma, ranp, omega, correlation, Pi = NULL, Slist = NULL, mvar = NULL){
  names.r    <- names(beta)  
  Ka    <- nrow(omega) 
  R     <- ncol(omega)  
  br    <- matrix(NA, Ka, R) 
  d.mu  <- d.sigma <- br   
  
  beta  <- drop(beta)
  sigma <- drop(sigma)
  names(beta) <- names(sigma) <- names.r
  rownames(br) <- rownames(d.mu) <- rownames(d.sigma) <- rownames(omega) <- names.r
  
  if (!is.null(Pi)){
    d.pis           <- matrix(NA, length(Pi), R)
    names.pi        <- rep(names(mvar), sapply(mvar, length))
    names(Pi) <- rownames(d.pis) <- names.pi
  }
  
  if (correlation){
    rownames(omega) <- NULL
    L <- makeL(sigma)
    L.omega <- tcrossprod(L, t(omega))
    rownames(L.omega) <- names.r
    br <- beta + L.omega
    d.mu <- matrix(1, Ka, R)
    rownames(d.mu) <- names.r
    if (!is.null(Pi)) d.pis[,] <- 1
    d.sigma <- omega[rep(1:Ka, Ka:1), ]
    for (i in 1:Ka){
      var  <- names.r[i]
      distr <- ranp[var]
      sigi <- i + cumsum(c(0, (Ka-1):1))[1:i]
      if (any(names(Slist) == var)) {
        pic <- as.vector(Pi[names(Pi) == var])
        br[var, ] <- br[var, ] + drop(tcrossprod(Slist[[var]], t(pic)))
      }
      if (distr == "cn"){
        br[var, ]  <- pmax(br[var, ], 0)
        d.mu[var,] <- as.numeric(br[var, ] > 0)
        d.sigma[sigi,] <- repRows(as.numeric(br[var, ] > 0), length(sigi))  * d.sigma[sigi, ]
        if (any(names(Slist) == var)){
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <- repRows(as.numeric(br[var, ] > 0), hmt)
        } 
      }
      if (distr == "ln"){
        br[var, ] <- exp(br[var, ])
        d.mu[var, ] <- br[var, ]
        d.sigma[sigi, ] <- repRows(br[var, ], length(sigi)) * d.sigma[sigi, ]
        if (any(names(Slist) == var)){
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <-  repRows(br[var, ], hmt)
        } 
      }
    }
    ## no correlation ##
  } else {
    for (j in 1:Ka) {
      var  <- names.r[j]
      distr <- ranp[var]
      if (distr == "n"){
        if (any(names(Slist) == var)){
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * omega[var,, drop = F]
        } else {
          br[var, ] <- beta[var] + sigma[var] * omega[var,, drop = F]
        }
        d.mu[var, ]    <- 1
        d.sigma[var, ] <- omega[var, ]
        if (any(names(Slist) == var)) d.pis[rownames(d.pis) == var, ] <- 1   
      }
      if (distr == "ln"){
        if (any(names(Slist) == var)){
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- exp(beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * omega[var, , drop = F])
        } else {
          br[var, ] <- exp(beta[var] + sigma[var] * omega[var, , drop = F])
        } 
        d.mu[var, ]    <- br[var, , drop = F]
        d.sigma[var,]  <- d.mu[var, ] * omega[var, ]
        if (any(names(Slist) == var)){
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <- repRows(d.mu[var, ], hmt)
        } 
      }
      if (distr == "cn"){
        if (any(names(Slist) == var)){
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- pmax(beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * omega[var, , drop = F], 0)
        }else{
          br[var, ] <- pmax(beta[var] + sigma[var] * omega[var, , drop = F], 0)
        }
        d.mu[var, ] <- as.numeric(br[var, ] > 0)
        d.sigma[var, ] <- d.mu[var, ] * omega[var, ]
        if (any(names(Slist) == var)){
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <-  repRows(as.numeric(br[var, ] > 0), hmt)
        } 
      }
      if (distr == "u"){
        etauni <- 2 * pnorm(omega[var,, drop = F]) - 1
        if(any(names(Slist) == var)){
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- beta[var] + drop(tcrossprod(Slist[[var]], t(pic)))  + etauni * sigma[var]
        } else {
          br[var, ] <- beta[var] + etauni * sigma[var]
        }
        d.mu[var, ] <- 1
        d.sigma[var, ] <- etauni
        if (any(names(Slist) == var)) d.pis[rownames(d.pis) == var, ] <- 1
      }  
      if (distr == "t"){
        etauni <- pnorm(omega[var,, drop = F])
        eta05  <- etauni < 0.5
        d.mu[var, ] <- 1
        d.sigma[var, ] <- eta05 * (sqrt(2 * etauni) - 1) + 
          !eta05 * (1 - sqrt(2 * (1 - etauni)))
        if (any(names(Slist) == var)){
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * d.sigma[var, ] 
        } else {
          br[var, ] <- beta[var] + sigma[var] * d.sigma[var, ]
        }
        if (any(names(Slist) == var)) d.pis[rownames(d.pis) == var, ] <- 1 
      }
      if (distr == "sb"){
        if (any(names(Slist) == var)){
          pic <- as.vector(Pi[names(Pi) == var])
          btemp <- beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * omega[var, , drop = F]
        } else {
          btemp <- beta[var] + sigma[var] * omega[var, , drop = F]
        } 
        br[var, ] <- exp(btemp) / (1 + exp(btemp))
        d.mu[var, ]    <- br[var, ] - br[var, ]^2
        d.sigma[var,]  <- d.mu[var, ] * omega[var, ]
        if (any(names(Slist) == var)){
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <- repRows(d.mu[var, ], hmt)
        } 
      }
    }
  }
  if(!is.null(Pi)){
    list(br = br, d.mu = d.mu, d.sigma = d.sigma, d.pis = d.pis)
  }else{
    list(br = br, d.mu = d.mu, d.sigma = d.sigma)
  }
}






