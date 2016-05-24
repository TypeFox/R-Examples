## Log-Likelihood functions for version 0.2
## version 0.2:
## - Add panel data
## - more flexible formulation of the hierarchical variables

## Standard Poisson model
lnpoisson <- function(theta, y, X,
                      weights = NULL, ...){
  if (is.null(weights)) weights <- 1
  index <- tcrossprod(X, t(theta))
  index <- pmin(index, 700)
  mu    <- exp(index)
  pi    <- dpois(y, mu)
  pi    <- pmax(pi, .Machine$double.eps)
  ll    <- sum(weights * log(pi))
  
  ## gradient
  G      <- as.vector(y - mu) * X 
  colnames(G) <- names(theta)
  attr(ll, 'gradient') <- weights * G 
  
  ## hessian
  H       <- crossprod(drop(- mu) * X, X)
  colnames(H) <- rownames(H) <- names(theta)
  attr(ll,'hessian') <- H
  attr(ll,'probabilities') <- pi
  ll
}

## Poisson with random parameters
lnlpoisson.ran <- function(theta, y, X, S = NULL, ranp, R, correlation, link,
                           weights = NULL, haltons = NULL, seed = 123, 
                           make.estb = FALSE, id = NULL, gradient = TRUE, 
                           mvar, ...){
  ## get globals
  N    <- nrow(X)
  K    <- ncol(X)
  panel <- !is.null(id)
  if (panel){
    n <- length(unique(id))
    if (length(weights) == 1) weights <- rep(weights, N)
  }
  hier <- !is.null(S)
  
  ## get variables
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  fixed <- !(Kc == 0L)
  Xa   <- X[ , Vara, drop = F]                        
  Xc   <- X[ , Varc, drop = F]                        
  
  # get parameters
  if (fixed) gamma <- theta[colnames(Xc)]
  beta.bar <- theta[paste('mean', colnames(Xa), sep = '.')]
  names(beta.bar) <- colnames(Xa)
  if (hier){
    hname  <- unlist(lapply(names(mvar), function(x) outer(x, mvar[[x]], FUN = paste, sep = ".")))
    phi <- theta[hname]
    sigma <- theta[-c(1:(K + length(hname)))]
    thecolumns <- lapply(mvar, function (x) {sort(match(x, colnames(S)))}) 
    Slist  <- lapply(thecolumns, function(x) S[, x, drop = FALSE])
  } else {
    sigma <- theta[-c(1:K)] 
    phi <- NULL
  }               
  
  ## get random draws
  set.seed(seed)
  Omega    <- make.draws(R * ifelse(panel, n, N), Ka, haltons) 
  
  ## fixed part
  if (fixed) ZB <- as.vector(crossprod(t(Xc), gamma)) 
  ## Random Part
  nind <- ifelse(panel, n, N)
  if (panel) theIds <- unique(id)
  XB    <- matrix(NA, N, R)
  if (hier) XaS <- matrix(NA, N, length(phi))
  if (make.estb){
    id.name <- if (panel) theIds else rownames(X)
    Br <- array(NA, dim = c(nind, R, Ka), 
                dimnames = list(id.name, NULL, names(beta.bar))) 
  }
  for (i in 1:nind){
    if (panel){
      anid <- theIds[i]
      theRows <- which(id == anid)
    }
    else theRows <- i
    if (hier) selecS <- lapply(Slist, function(x) x[i, , drop = FALSE]) else selecS <- NULL
    beta.r   <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                           omega = Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE], 
                           correlation = correlation, Pi = phi, Slist = selecS, mvar = mvar)
    XB[theRows, ]  <- crossprod(t(Xa[theRows, , drop = FALSE]), beta.r$br)
    if (hier && gradient) XaS[theRows, ] <- Reduce(cbind, lapply(names(selecS), function(x) kronecker(selecS[[x]], 
                                                   rep(1, length(theRows))) * Xa[theRows, x]))
    if (make.estb) Br[i, , ] <- t(beta.r$br)
  }
  
  index <- if (fixed) ZB + XB else XB
  index <- pmin(index, 700)
  # note that exp(709) = Inf ==> gradient = Inf
  mu    <- exp(index)
  Pitr   <- dpois(y, lambda = mu)
  if (panel) Pir <- apply(Pitr, 2, tapply, id, prod) else Pir <- Pitr
  Pir   <- pmax(Pir, .Machine$double.eps)
  Pi    <- rowSums(Pir) / R
  if (make.estb) Qir <- Pir / (Pi * R)
  lls      <- if (panel) sum(log(Pi) * weights[!duplicated(id)]) else sum(log(Pi) * weights)
  
  ## gradient
  if (gradient){
    lambda <- y - mu
    Qir    <- Pir / (Pi * R)
    if (panel) Qir <- Qir[as.character(id), ]
    eta    <- Qir * lambda            
    dUdb <- matrix(NA, N, Ka)
    if (hier) dUdphi <- matrix(NA, N, length(phi))
    dUds <- if (correlation) matrix(NA, N, (0.5 * Ka * (Ka + 1)))  else matrix(NA, N, Ka)  
    for (i in 1:nind){
      if (panel){
        anid <- theIds[i]
        theRows <- which(id == anid)
      }
      else theRows <- i
      if (hier) selecS <- lapply(Slist, function(x) x[i, , drop = FALSE]) else selecS <- NULL
      beta.r   <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                             omega = Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE], 
                             correlation = correlation, Pi = phi, Slist = selecS, mvar = mvar)
      dUdb[theRows, ] <- tcrossprod(eta[theRows, , drop = FALSE], beta.r$d.mu)  
      dUds[theRows, ] <- tcrossprod(eta[theRows, , drop = FALSE], beta.r$d.sigma)
      if (hier) dUdphi[theRows, ] <- tcrossprod(eta[theRows, , drop = FALSE], beta.r$d.pis)
    }
    if (correlation){
      vecX <- c()
      for (i in 1:Ka){
        vecX <- c(vecX, i:Ka)
      }
      Xac <- Xa[,vecX]
    } else {
      Xac <- Xa  
    }
    gbarfi <- if (fixed) Xc  * rowSums(eta) else numeric()
    gbarmi <- Xa  * dUdb
    gbarvi <- Xac * dUds
    gbarphi <- if (hier)  XaS * dUdphi else numeric()
    gbari  <- cbind(gbarfi, gbarmi, gbarphi, gbarvi)
    colnames(gbari) <- names(theta)
    attr(lls, 'gradient') <- weights * gbari
  }
  if (make.estb){
    attr(lls,'bi') <- Br
    attr(lls,'Qir') <- Qir
    attr(lls,'probabilities') <- rowSums(Pitr) / R
  }    
  lls
}

## Standard binary
lnbinary <- function(theta, y, X, link, 
                   weights = NULL, ...){
  if (is.null(weights)) weights <- 1
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  ddfun <- switch(link,
                  "logit" = function(x) (1 - 2 * pfun(x)) * pfun(x) * (1 - pfun(x)),
                  "probit"= function(x) -x * dnorm(x))  
  mill  <- function(x) dfun(x) / pmax(pfun(x), .Machine$double.eps)
  millh <- function(x) ddfun(x) / pmax(pfun(x), .Machine$double.eps) - 
                       (dfun(x) / pmax(pfun(x), .Machine$double.eps))^2
  index  <- tcrossprod(X, t(theta))
  q      <- 2 * y - 1 
  pi     <- pfun(q * index)
  pi     <- pmax(pi, .Machine$double.eps)
  ll     <- sum(weights * log(pi))
  
  ## Gradient
  G <- as.vector(q * mill(q * index)) * X
  colnames(G) <- names(theta)
  attr(ll,'gradient') <- weights * G
  
  ## Hessian
  H <- crossprod(as.vector(millh(q * index)) * X, X)
  colnames(H) <- rownames(H) <- names(theta)
  attr(ll,'hessian') <- H
  attr(ll,'probabilities') <- pi
  ll
}

## Binary with Random parameters
lnlbinary.ran <- function(theta, y, X, S = NULL, ranp, R, correlation, link,
                          weights = NULL, haltons = NULL, seed = 123, make.estb = FALSE,
                          id = NULL, gradient = TRUE, mvar, ... ){
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  mill <- function(x) dfun(x) / pmax(pfun(x), .Machine$double.eps)
  
  ## get globals
  N    <- nrow(X)
  K    <- ncol(X)
  panel <- !is.null(id)
  if (panel){
    n <- length(unique(id))
    if (length(weights) == 1) weights <- rep(weights, N)
  }
  hier <- !is.null(S)
  ## get variables
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  fixed <- !(Kc == 0L)
  Xa   <- X[, Vara, drop = F]                        
  Xc   <- X[, Varc, drop = F]                        
  
  ## get vector of parameters
  if (fixed) gamma <- theta[colnames(Xc)]
  beta.bar <- theta[paste('mean', colnames(Xa), sep = '.')]
  names(beta.bar) <- colnames(Xa)
  if (hier){
    hname  <- unlist(lapply(names(mvar), function(x) outer(x, mvar[[x]], FUN = paste, sep = ".")))
    phi <- theta[hname]
    sigma <- theta[-c(1:(K + length(hname)))]
    thecolumns <- lapply(mvar, function (x) {sort(match(x, colnames(S)))}) 
    Slist  <- lapply(thecolumns, function(x) S[, x, drop = FALSE]) # S variables in a list
  } else {
    sigma <- theta[-c(1:K)] 
    phi <- NULL
  }                 
  
  ## make random draws
  set.seed(seed)
  Omega    <- make.draws(R * ifelse(panel, n, N), Ka, haltons) 
  
  ## fixed part of index
  if (fixed) ZB <- as.vector(crossprod(t(Xc), gamma)) 
  ## random part of index
  nind <- ifelse(panel, n, N)
  if (panel) theIds <- unique(id)
  XB <- matrix(NA, N, R)
  if (hier) XaS <- matrix(NA, N, length(phi))
  if (make.estb){
    id.name <- if (panel) theIds else rownames(X)
    Br <- array(NA, dim = c(nind, R, Ka), 
                dimnames = list(id.name, NULL, names(beta.bar))) 
  }
  for (i in 1:nind){
    if (panel){
      anid <- theIds[i]
      theRows <- which(id == anid)
    }
    else theRows <- i
    if (hier) selecS <- lapply(Slist, function(x) x[i, , drop = FALSE]) else selecS <- NULL
    beta.r   <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                           omega = Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE], 
                           correlation = correlation, Pi = phi, Slist = selecS, mvar = mvar)
    XB[theRows, ]  <- crossprod(t(Xa[theRows, , drop = FALSE]), beta.r$br)
    if (hier && gradient) XaS[theRows, ] <- Reduce(cbind, lapply(names(selecS), function(x) kronecker(selecS[[x]], 
                                                   rep(1, length(theRows))) * Xa[theRows, x]))
    if (make.estb) Br[i, , ] <- t(beta.r$br)
  }
  ## get probability and log-likelihhod
  q <- 2 * y - 1
  index <- if (fixed) ZB + XB else XB
  Pitr   <- pfun(q * index)
  if (panel) Pir <- apply(Pitr, 2, tapply, id, prod) else Pir <- Pitr
  Pir   <- pmax(Pir, .Machine$double.eps)
  Pi    <- rowSums(Pir) / R
  if (make.estb) Qir <- Pir / (Pi * R)
  lls      <- if (panel) sum(log(Pi) * weights[!duplicated(id)]) else sum(log(Pi) * weights)
  
  ## make gradient
  if (gradient){
    lambda <- q * mill(q * index)
    Qir    <- Pir / (Pi * R)
    if (panel) Qir <- Qir[as.character(id), ]
    eta    <- Qir * lambda            
    dUdb <- matrix(NA, N, Ka)
    if (hier) dUdphi <- matrix(NA, N, length(phi))
    dUds <- if (correlation) matrix(NA, N, (0.5 * Ka * (Ka + 1))) else  matrix(NA, N, Ka) 
    for (i in 1:nind){
      if (panel){
        anid <- theIds[i]
        theRows <- which(id == anid)
      }
      else theRows <- i
      if (hier) selecS <- lapply(Slist, function(x) x[i, , drop = FALSE]) else selecS <- NULL
      beta.r   <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                             omega = Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE], 
                             correlation = correlation, Pi = phi, Slist = selecS, mvar = mvar)
      dUdb[theRows, ] <- tcrossprod(eta[theRows, ], beta.r$d.mu)  
      dUds[theRows, ] <- tcrossprod(eta[theRows, ], beta.r$d.sigma)
      if (hier) dUdphi[theRows, ] <- tcrossprod(eta[theRows, , drop = FALSE], beta.r$d.pis)
    }    
    if (correlation){
      vecX <- c()
      for (i in 1:Ka){
        vecX <- c(vecX, i:Ka)
      }
      Xac <- Xa[, vecX]
    } else{
      Xac <- Xa  
    }
    gbarfi <- if (fixed) Xc  * rowSums(eta) else numeric()
    gbarmi <- Xa  * dUdb
    gbarphi <- if (hier)  XaS * dUdphi else numeric()
    gbarvi <- Xac * dUds
    gbari  <- cbind(gbarfi, gbarmi, gbarphi, gbarvi)
    colnames(gbari) <- names(theta)
    attr(lls, 'gradient') <- weights * gbari
  }
  if (make.estb){
    attr(lls,'bi') <- Br
    attr(lls,'Qir') <- Qir
    attr(lls,'probabilities') <- rowSums(Pitr) / R
  }    
  lls
}

## Standard Ordered Model
lnordered <- function(theta, y, X, link,
                    weights = NULL, ... ){
  if (is.null(weights)) weights <- 1
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis
  )
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis
  )
  ## it is assumed that y is a "factor" class
  m       <- sort(unique(y))
  m       <- unclass(m)
  J       <- length(levels(y))
  ## Just J-2 alpha are estimated, so the gradient is just for 2:(J-1) cuts
  m       <- as.matrix(m[2:(J - 1)], nrow = (J - 2))               
  y       <- unclass(y)
  delta   <- t(kronecker(m, t(rep(1, nrow(X)))))    
  ## indicator for alternative and previous alternative
  deltaj  <- delta == y
  deltak  <- delta == y - 1
  
  alpha    <- theta[1:(J - 2)]
  kappa    <- c(-Inf, cumsum(c(0, exp(alpha))) , +Inf)
  beta     <- theta[-c(1:(J - 2))]
  index    <- tcrossprod(X, t(beta))
  eta1     <- kappa[y + 1] - index
  eta2     <- kappa[y] - index
  pi       <- pfun(eta1) - pfun(eta2)
  pi       <- ifelse(pi <= 0, .Machine$double.eps, pi) 
  ll       <- sum(weights * log(pi))
  
  ## Gradient of parameters
  phi1     <- dfun(eta1) ; phi2 <- dfun(eta2)
  lambda   <- drop((phi2 - phi1) / pi)
  gbeta    <- as.vector(lambda) * X
  
  ## Gradient of kappa
  lamkap   <- (deltaj * drop(phi1) - deltak * drop(phi2)) / as.vector(pi)
  gkappa   <- lamkap %*% jacobian(alpha)
  G        <- cbind(gkappa, gbeta)
  colnames(G) <- names(theta)
  attr(ll,'gradient') <- weights * G  
  attr(ll,'probabilities') <- pi
  ll
}

## Ordered Model with Random Parameters
lnordered.ran <- function(theta, y, X, S = NULL, ranp, R, correlation, link,
                          weights = NULL, haltons = NULL, seed = 123, make.estb = FALSE,
                          id = NULL, gradient = TRUE, mvar, ...){
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  N    <- nrow(X)
  K    <- ncol(X)
  panel <- !is.null(id)
  if (panel){
    n <- length(unique(id))
    if (length(weights) == 1) weights <- rep(weights, N)
  }
  hier <- !is.null(S)
  
  ## Get variables
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  fixed <- !(Kc == 0L)
  Xa   <- X[, Vara, drop = F]                        
  Xc   <- X[, Varc, drop = F] 
  
  m       <- unclass(sort(unique(y)))
  J       <- length(levels(y))
  m       <- as.matrix(m[2:(J - 1)], nrow = (J - 2))               
  y       <- unclass(y)
  delta   <- t(kronecker(m, t(rep(1, nrow(X)))))    
  deltaj  <- delta == y
  deltak  <- delta == y - 1
  
  ## Parameters
  alpha    <- theta[1:(J - 2)]
  kappa    <- c(-Inf, cumsum(c(0, exp(alpha))), +Inf)
  if (fixed) gamma <- theta[colnames(Xc)]
  beta.bar <- theta[paste('mean', colnames(Xa), sep = '.')]
  names(beta.bar) <- colnames(Xa)
  if (hier){
    hname  <- unlist(lapply(names(mvar), function(x) outer(x, mvar[[x]], FUN = paste, sep = ".")))
    phi <- theta[hname]
    sigma <- theta[-c(1:(K + (J - 2) + length(hname)))]
    thecolumns <- lapply(mvar, function (x) {sort(match(x, colnames(S)))}) 
    Slist  <- lapply(thecolumns, function(x) S[, x, drop = FALSE]) # S variables in a list
  } else {
    sigma <- theta[-c(1:(K + (J - 2)))] 
    phi <- NULL
  }  
  
  ## Random Draws
  set.seed(seed)
  Omega <- make.draws(R * ifelse(panel, n, N), Ka, haltons) 
  
  ## Make Fixed Part
  if (fixed) ZB <- as.vector(crossprod(t(Xc), gamma)) 
  ## Make Random Part
  nind <- ifelse(panel, n, N)
  if (panel) theIds <- unique(id)
  XB    <- matrix(NA, N, R)
  if (hier) XaS <- matrix(NA, N, length(phi))
  if (make.estb){
    id.name <- if (panel) theIds else rownames(X)
    Br <- array(NA, dim = c(nind, R, Ka), 
                dimnames = list(id.name, NULL, names(beta.bar))) 
  }
  for (i in 1:nind){
    if (panel){
      anid <- theIds[i]
      theRows <- which(id == anid)
    }
    else theRows <- i
    if (hier) selecS <- lapply(Slist, function(x) x[i, , drop = FALSE]) else selecS <- NULL
    beta.r   <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                           omega = Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE], 
                           correlation = correlation, Pi = phi, Slist = selecS, mvar = mvar)
    XB[theRows, ]  <- crossprod(t(Xa[theRows, , drop = FALSE]), beta.r$br)
    if (hier && gradient) XaS[theRows, ] <- Reduce(cbind, lapply(names(selecS), function(x) kronecker(selecS[[x]], 
                                                   rep(1, length(theRows))) * Xa[theRows, x]))
    if (make.estb) Br[i, , ] <- t(beta.r$br)
  }
  
  index    <- if (fixed) ZB + XB else XB
  eta1     <- kappa[y + 1] - index
  eta2     <- kappa[y]   - index
  Pitr      <- pfun(eta1) - pfun(eta2)
  if (panel) Pir <- apply(Pitr, 2, tapply, id, prod) else Pir <- Pitr
  Pir      <- pmax(Pir, .Machine$double.eps)
  Pi       <- rowSums(Pir) / R
  if (make.estb) Qir <- Pir / (Pi * R)
  lls      <- if (panel) sum(log(Pi) * weights[!duplicated(id)]) else sum(log(Pi) * weights)
  
  ## Gradients
  if (gradient){
    phi1     <- pmax(dfun(eta1), .Machine$double.eps)
    phi2     <- pmax(dfun(eta2), .Machine$double.eps) 
    lambda   <- (phi2 - phi1) / pmax(Pitr, .Machine$double.eps)
    Qir      <- Pir / (Pi * R)
    if (panel) Qir <- Qir[as.character(id), ]
    eta      <- Qir * lambda
    
    gkappa   <- vector(mode = "list", length = (J - 2))
    for (j in 1:(J - 2)){
      gkappa[[j]] <- matrix(NA, N, R)
      gkappa[[j]] <- drop(deltaj[ , j]) * (phi1 / pmax(Pitr, .Machine$double.eps)) - 
        drop(deltak[ ,j]) * (phi2 / pmax(Pitr, .Machine$double.eps))
    }
    gkappa <- lapply(gkappa, function(x) x * Qir)
    gkappa <- lapply(gkappa, function(x) apply(x, 1, sum))
    gkappa <- Reduce(cbind, gkappa) %*% jacobian(alpha)
    
    dUdb <- matrix(NA, N, Ka)
    if (hier) dUdphi <- matrix(NA, N, length(phi))
    dUds <- if (correlation) matrix(NA, N, (0.5 * Ka * (Ka + 1))) else matrix(NA, N, Ka)  
    for (i in 1:nind){
      if (panel){
        anid <- theIds[i]
        theRows <- which(id == anid)
      }
      else theRows <- i
      if (hier) selecS <- lapply(Slist, function(x) x[i, , drop = FALSE]) else selecS <- NULL
      beta.r   <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                             omega = Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE], 
                             correlation = correlation, Pi = phi, Slist = selecS, mvar = mvar)
      dUdb[theRows, ] <- tcrossprod(eta[theRows, ], beta.r$d.mu)  
      dUds[theRows, ] <- tcrossprod(eta[theRows, ], beta.r$d.sigma)
      if (hier) dUdphi[theRows, ] <- tcrossprod(eta[theRows, , drop = FALSE], beta.r$d.pis)
    }
    if (correlation){
      vecX <- c()
      for (i in 1:Ka){
        vecX <- c(vecX, i:Ka)
      }
      Xac <- Xa[,vecX]
    } else {
      Xac <- Xa  
    }
    gbarfi <- if (fixed) Xc  * rowSums(eta) else numeric()
    gbarmi <- Xa  * dUdb
    gbarvi <- Xac * dUds
    gbarphi <- if (hier)  XaS * dUdphi else numeric()
    gbari  <- cbind(gkappa, gbarfi, gbarmi, gbarphi, gbarvi)
    colnames(gbari) <- names(theta)
    attr(lls, 'gradient') <- weights * gbari
  }
  if (make.estb){
    attr(lls,'bi') <- Br
    attr(lls,'Qir') <- Qir
    attr(lls,'probabilities') <- rowSums(Pitr) / R
  } 
  lls
}
