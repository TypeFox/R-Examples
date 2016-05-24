#########################################
# Maximum Likelihood Functions for gmnl
#########################################

## suml function from mlogit (Croissant, 2013)
suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))){
    d <- dim(x[[1]])
    s <- matrix(0,d[1], d[2])
    for (i in 1:n){
      x[[i]][is.na(x[[i]])] <- 0
      s <- s + x[[i]]
    }
  }
  else{
    s <- rep(0, length(x[[n]]))
    for (i in 1:n){
      x[[i]][is.na(x[[i]])] <- 0
      s <- s + x[[i]]
    }
  }
  s
}

## suml function for array
suml.array <- function(x){
  n <- length(x)
  d <- dim(x[[1]])
  s <- array(0, dim = c(d[1], d[2], d[3]))
  for (i in 1:n){
    x[[i]][is.na(x[[i]])] <- 0
    s <- s + x[[i]]
  }
  s
}

## Multinomial Model based on mlogit (Croissant, 2013)
ll.mlogit <- function(theta, y, X, gradient = TRUE, 
                      hessian =  TRUE, weights = NULL, 
                      get.bi = FALSE){
  if (is.null(weights)) weights <- 1
  exb  <- lapply(X, function(x) exp(crossprod(t(x), theta)))
  sexb <- suml(exb)
  Pni  <- lapply(exb, function(x){v <- x / sexb; 
                                  v[is.na(v)] <- 0;
                                  as.vector(v)})
  Pn <- Reduce("+", mapply("*", Pni, y, SIMPLIFY = FALSE))
  ll <- sum(log(Pn) * weights) 
  
  # Gradient
  if (gradient | hessian ){
    Px <- suml(mapply("*", X, Pni, SIMPLIFY = FALSE)) 
    yx <- suml(mapply("*", X, y, SIMPLIFY = FALSE))
    gradi <- (yx - Px) * weights ; colnames(gradi) <- names(theta)
    gradi[is.na(gradi)] <- 0
    attr(ll, "gradient") <- gradi
  }
  # Hessian
  if (hessian){
    dxpx <- lapply(X, function(x){
      d <- x - Px ;
      d[is.na(d)] <- 0 ;
      d})
    hess <- - suml(mapply(function(x, y) crossprod(x * y, y), Pni, dxpx , SIMPLIFY = FALSE))
    attr(ll, "hessian")
  }
  if (get.bi){
    Pni <- Reduce("cbind", Pni)
    colnames(Pni) <- names(y)
    attr(ll, "prob.alt") <- Pni
    attr(ll, "prob.ind") <- Pn
  }
  ll
}

## Simulated Maximum Likelihood for SMNL model
ll.smlogit <- function(theta, y, X, H = NULL, id = NULL, R,
                       weights = NULL, seed = 12345, notscale = NULL,
                       bound.err, gradient = TRUE, typeR = TRUE,
                       get.bi =  FALSE){
  ## Get globals
  K <- ncol(X[[1]])
  J <- length(X)
  N <- nrow(X[[1]])
  panel <- !is.null(id)
  if (panel){
    n <- length(unique(id))
    if (length(weights) == 1) weights <- rep(weights, N)
  } 
  het <- !is.null(H)
  
  ## Get parameters
  beta   <- theta[1L:K]
  tau    <- theta[K + 1L]
  if (het) delta <- theta[-c(1L:(K + 1L))]
  
  ## Make random draws
  set.seed(seed)
  epsilon <- if (typeR) truncnorm::rtruncnorm(R * ifelse(panel, n, N), a = - bound.err, b = bound.err) 
             else Make.epsi(R * ifelse(panel, n, N), bound.err)
  sigmaM  <- exp(tau * matrix(epsilon, ncol = ifelse(panel, n, N))) 
  sum.se  <- apply(sigmaM * matrix(epsilon, ncol = ifelse(panel, n, N)), 1, sum) 
  p.se    <- sum.se / apply(sigmaM, 1, sum) 
  sigbar  <- -log(apply(sigmaM , 1, mean))
  
  ## Random Part
  XB <- vector(mode = 'list', length = J)
  for (j in 1:J) XB[[j]] <- matrix(NA, N, R) 
  nind <- ifelse(panel, n, N)
  if (panel) theIds <- unique(id)
  if (get.bi) bi <- array(NA, dim = c(nind, R, K), 
                          dimnames = list(NULL, NULL, names(beta))) 
  for (i in 1:nind){
    if (panel){
      anid <- theIds[i]
      theRows <- which(id == anid)
    }
    else theRows <- i 
    err <- epsilon[((i - 1) * R + 1):(i * R)]
    if (het){ 
      sigma.nr <- exp(sigbar +  drop(tcrossprod(H[i, , drop = FALSE], t(delta))) + tau * err) 
    } else sigma.nr <- exp(sigbar + tau * err)
    br <- matrix(NA, K, R)
    for(k in 1:K){
      ns <- notscale[k]
      if (ns == 0) br[k, ] <- beta[k] * sigma.nr else br[k, ] <- beta[k]
    } 
    for (j in 1:J) XB[[j]][theRows, ] <- crossprod(t(X[[j]][theRows, , drop = FALSE]), br)
    if (get.bi) bi[i, , ] <- t(br)
  }
  ## Make log L
  EXB   <- lapply(XB, function(x) exp(x))
  SEXB  <- suml(EXB)
  Pntir <- lapply(EXB, function(x) x / SEXB)
  Pnr   <- suml(mapply("*", Pntir, y, SIMPLIFY =  FALSE))
  if (panel) Pnr <- apply(Pnr, 2, tapply, id, prod)
  Pn <- apply(Pnr, 1, mean)
  if (get.bi)  Qnr  <-  Pnr / (Pn * R)
  lnL <- if (panel) sum(log(Pn) * weights[!duplicated(id)]) else sum(log(Pn) * weights)
  
  ## Make Gradient
  if (gradient){
    lambda <- mapply(function(y, p) y - p, y, Pntir, SIMPLIFY =  FALSE) # J-list NT * R
    Qnr    <-  Pnr / (Pn * R) # N * R
    if (panel) Qnr <-  Qnr[as.character(id), ] # NT * R
    eta  <- lapply(lambda, function(x) x * Qnr)  # J- list NT * R
      
    dUdb <- dUdt <- vector(mode = 'list', length = J)
    if (het) dUdd <- vector(mode = 'list', length = J)
    d.mu <- d.tau <- matrix(NA, K, R)
    if (het) d.delta <- matrix(NA, K, R)
    for (j in 1:J){
      dUdt[[j]] <- dUdb[[j]] <- matrix(NA, N, K)
      if (het) dUdd[[j]] <- matrix(NA, N, K)
    }
    for (i in 1:nind){
      if (panel){
        anid <- theIds[i]
        theRows <- which(id == anid)
      }
      else theRows <- i
      err <- epsilon[((i - 1) * R + 1):(i * R)]
      if (het){ 
        sigma.nr <- exp(sigbar +  drop(tcrossprod(H[i, , drop = FALSE], t(delta))) + tau * err) 
      } else sigma.nr <- exp(sigbar + tau * err)
      for(k in 1:K){
        ns <- notscale[k]
        if(ns == 0){
          d.mu[k, ]  <- sigma.nr
          d.tau[k, ] <- sigma.nr * (- p.se + err) * beta[k]
          if (het) d.delta[k, ] <- sigma.nr * beta[k]
        } else {
          d.mu[k, ] <- 1 
          d.tau[k, ] <- 0 
          if (het) d.delta[k, ] <- 0
        }
      }  
      for (j in 1:J) {
        dUdt[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], d.tau)
        dUdb[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], d.mu)
        if (het) dUdd[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], d.delta)
      }   
    }
    tempt  <- mapply("*", dUdt, X, SIMPLIFY =  FALSE)
    gbart  <- suml(lapply(tempt, function(x) rowSums(x)))
    gbarmi <- suml(mapply("*", X, dUdb, SIMPLIFY =  FALSE))
    if (het){
      tempd  <- mapply("*", dUdd, X, SIMPLIFY =  FALSE)
      if(panel) H <- H[id, ]
      gbarde <- suml(lapply(tempd, function(x) rowSums(x))) * H
      gbari  <- cbind(gbarmi, gbart, gbarde)
    } else {
      gbari  <- cbind(gbarmi, gbart)
    }    
    colnames(gbari) <- names(theta)
    attr(lnL, 'gradient') <- gbari * weights
  }
  if (get.bi) {
    attr(lnL, "prob.alt") <- sapply(Pntir, function(x) apply(x, 1, mean))
    attr(lnL, "prob.ind") <- Pn
    attr(lnL, "bi") <- bi
    attr(lnL, 'Qir') <- Qnr
  }
  lnL 
}

## Simulated Maximum Likelihood for MIXL model
ll.mixlog <- function(theta, y, X, Z = NULL, id = NULL, ranp, R,
                      correlation, weights = NULL,
                      haltons =  NULL, seed = 12345,
                      mvar, gradient =  TRUE, get.bi = FALSE){
  ## Get globals
  K <- ncol(X[[1]])
  J <- length(X)
  N <- nrow(X[[1]])
  panel <- !is.null(id)
  if (panel){
    n <- length(unique(id))
    if (length(weights) == 1) weights <- rep(weights, N)
  } 
  hier <- !is.null(Z)
  
  ## Get variables
  Vara <- sort(match(names(ranp), colnames(X[[1]])))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  fixed <- !(Kc == 0)
  Xa   <- lapply(X, function(x) x[, Vara, drop = F])                        
  Xc   <- lapply(X, function(x) x[, Varc, drop = F]) 
  K <- Kc + Ka
  
  ## Get Parameters
  if (fixed) gamma <- theta[colnames(Xc[[1]])]
  beta.bar <- theta[colnames(Xa[[1]])]
  if (hier){
    # Parameters
    hname  <- unlist(lapply(names(mvar), function(x) outer(x, mvar[[x]], FUN = paste, sep = ".")))
    phi <- theta[hname]
    sigma <- theta[-c(1:(K + length(hname)))] 
    
    # Variables
    thecolumns <- lapply(mvar, function (x) {sort(match(x, colnames(Z)))}) 
    Zlist  <- lapply(thecolumns, function(x) Z[, x, drop = FALSE]) # Z variables in a list
  } else {
    sigma <- theta[-c(1:K)]
    phi  <- NULL
  }
  
  ## Make random draws
  set.seed(seed)
  Omega <- make.draws(R * ifelse(panel, n, N), Ka, haltons)
  
  # Fixed part of utility
  if (fixed) XBf <- lapply(Xc, function(x) as.vector(crossprod(t(as.matrix(x)), gamma)))
  # Random part of utility
  XBr <- vector(mode = 'list', length = J)
  if (hier) XaZ <- vector(mode = 'list', length = J)
  for (j in 1:J) {
    XBr[[j]] <- matrix(NA, N, R)
    if (hier) XaZ[[j]] <- matrix(NA, N, length(phi))
  }  
  nind <- ifelse(panel, n, N)
  if (panel) theIds <- unique(id)
  if (get.bi) bi <- array(NA, dim = c(nind, R, Ka), 
                          dimnames = list(NULL, NULL, names(beta.bar))) 
  for (i in 1:nind){
    if (panel){
      anid <- theIds[i]
      theRows <- which(id == anid)
    }
    else theRows <- i
    if (hier) selecZ <- lapply(Zlist, function(x) x[i, , drop = FALSE]) else selecZ <- NULL
    b <- Makeh.rcoef(beta.bar, sigma, ranp, Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE], 
                     correlation, Pi = phi, Slist = selecZ, mvar = mvar)
    for (j in 1:J) {
      XBr[[j]][theRows, ] <- crossprod(t(Xa[[j]][theRows, , drop = FALSE]), b$br) 
      if (hier && gradient) XaZ[[j]][theRows, ] <- Reduce(cbind, lapply(names(selecZ), function(x) kronecker(selecZ[[x]], rep(1, length(theRows))) * Xa[[j]][theRows, x]))
    }
    if (get.bi) bi[i, , ] <- t(b$br)
  }
  ## Make log L
  if (fixed){
    EXB <- mapply(function(x, y) exp(x + y), XBf, XBr, SIMPLIFY = FALSE)
  } else EXB <- lapply(XBr, function(x) exp(x))
  
  SEXB  <- suml(EXB)
  Pntir <- lapply(EXB, function(x) x / SEXB)
  Pnr   <- suml(mapply("*", Pntir, y, SIMPLIFY =  FALSE))
  if (panel) Pnr <- apply(Pnr, 2, tapply, id, prod)
  Pn <- apply(Pnr, 1, mean)
  if (get.bi)  Qnr  <-  Pnr / (Pn * R)
  lnL <- if (panel) sum(log(Pn) * weights[!duplicated(id)]) else sum(log(Pn) * weights)
  
  ## Make Gradient
  if (gradient){
    lambda <- mapply(function(y, p) y - p, y, Pntir, SIMPLIFY =  FALSE) 
    Qnr  <-  Pnr / (Pn * R) 
    if(panel) Qnr <- Qnr[as.character(id), ] 
    eta  <- lapply(lambda, function(x) x * Qnr)
    
    dUdb <- dUds <- vector(mode = 'list', length = J)
    if (hier) dUdphi <- vector(mode = 'list', length = J)
    for (j in 1:J){
      dUdb[[j]] <- matrix(NA, N, Ka)
      dUds[[j]] <- matrix(NA, N, length(sigma))
      if (hier) dUdphi[[j]] <- matrix(NA, N, length(phi))
    }
    for (i in 1:nind){
      if (panel){
        anid <- theIds[i]
        theRows <- which(id == anid)
      }
      else theRows <- i
      if (hier) selecZ <- lapply(Zlist, function(x) x[i, , drop = FALSE]) else selecZ <- NULL
      b <- Makeh.rcoef(beta.bar, sigma, ranp, Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE], 
                       correlation, Pi = phi, Slist = selecZ, mvar = mvar)
      for (j in 1:J) {
        dUdb[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.mu)
        dUds[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.sigma)
        if (hier) dUdphi[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.pis)
      }   
    }
    if (correlation){
      vecX <- c()
      for (i in 1:Ka){
        vecX <- c(vecX, i:Ka)
      }
      Xac <- lapply(Xa,  function(x) x[, vecX])
    } else{
      Xac <- Xa  
    }
    
    if (fixed){
      seta <- lapply(eta, function(x) rowSums(x))
      gbarfi <- suml(mapply("*", Xc, seta, SIMPLIFY =  FALSE))
    } else gbarfi <- numeric()
    gbarmi <- suml(mapply("*", Xa, dUdb, SIMPLIFY =  FALSE))
    gbarvi <- suml(mapply("*", Xac, dUds, SIMPLIFY =  FALSE))
    if (hier) {
      gbarphi <- suml(mapply("*", XaZ, dUdphi, SIMPLIFY =  FALSE))
      gbari  <- cbind(gbarfi, gbarmi, gbarphi, gbarvi)
    } else {
      gbari  <- cbind(gbarfi, gbarmi, gbarvi)
    }
    colnames(gbari) <- names(theta)
    attr(lnL, 'gradient') <- gbari * weights
  }
  
  if (get.bi) {
    attr(lnL, "prob.alt") <- sapply(Pntir, function(x) apply(x, 1, mean))
    attr(lnL, "prob.ind") <- Pn
    attr(lnL, "bi") <- bi
    attr(lnL, 'Qir') <- Qnr
  }
  lnL
}

## Simulated Maximum Likelihood for GMNL model
#' @import utils
ll.gmlogit <- function(theta, y, X, H = NULL, id = NULL, ranp, R,
                       correlation, weights = NULL, hgamma = "direct", 
                       haltons =  NULL, seed = 12345, bound.err,
                       notscale = NULL, gradient = TRUE, typeR = TRUE,
                       get.bi = FALSE){
  ## Get globals
  K <- ncol(X[[1]])
  J <- length(X)
  N <- nrow(X[[1]])
  panel <- !is.null(id)
  if (panel) {
    n <- length(unique(id))
    if (length(weights) == 1) weights <- rep(weights, N)
  } 
  het <- !is.null(H)
  
  ## Get variables
  Vara <- sort(match(names(ranp), colnames(X[[1]])))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  fixed <- !(Kc == 0)
  Xa   <- lapply(X, function(x) x[, Vara, drop = F])                       
  Xc   <- lapply(X, function(x) x[, Varc, drop = F])  
  K <- Kc + Ka
  
  ## Get parameters
  if (fixed) delta <- theta[colnames(Xc[[1]])]
  beta.bar <- theta[colnames(Xa[[1]])]
  if (het){
    names.het <- paste("het", colnames(H), sep = ".")
    het.par <- theta[names.het]
    stds <- head(theta[-c(1:(K + length(names.het)))], -2L) 
  } else stds    <- head(theta[-c(1:K)], -2L)
  tau    <- theta["tau"]
  gamma  <- theta["gamma"]
  
  ## Make random draws
  set.seed(seed)
  Omega <- make.draws(R * ifelse(panel, n, N), Ka, haltons)
  epsilon <- if (typeR) truncnorm::rtruncnorm(R * ifelse(panel, n, N), a = - bound.err, b = bound.err) else Make.epsi(R * ifelse(panel, n, N), bound.err)
  sigmaM  <- exp(tau * matrix(epsilon, ncol = ifelse(panel, n, N)))
  sum.se  <- apply(sigmaM * matrix(epsilon, ncol = ifelse(panel, n, N)), 1, sum) 
  p.se    <- sum.se / apply(sigmaM, 1, sum) 
  sigbar  <- - log(apply(sigmaM , 1, mean))
  
  ## Make utility for each individual
  XBr <- XBf <- vector(mode = 'list', length = J)
  for (j in 1:J) XBr[[j]] <- XBf[[j]] <- matrix(NA, N, R)
  nind <- ifelse(panel, n, N)
  if (panel) theIds <- unique(id)
  if (get.bi) bi <- array(NA, dim = c(nind, R, Ka + Kc), 
                          dimnames = list(NULL, NULL, colnames(X[[1]]))) 
  for (i in 1:nind){
    if (panel){
      anid <- theIds[i]
      theRows <- which(id == anid)
    }
    else theRows <- i
    err <- epsilon[((i - 1) * R + 1):(i * R)]
    if (het){
      Hb <- drop(tcrossprod(H[i, , drop = FALSE], t(het.par)))
      sigma.nr <- exp(sigbar + Hb + tau * err) 
    } else {
      sigma.nr <- exp(sigbar + tau * err)
      Hb <- NULL
    }    
    if (fixed){
      br <- matrix(NA, Kc, R)
      rownames(br) <- names(delta)
      for(k in 1:Kc){
        ns <- notscale[names(delta[k])]
        if (ns == 0) br[k, ] <- delta[k] * sigma.nr else br[k, ] <- delta[k]
      } 
    }
    b <- MakeGcoef(beta = beta.bar, stds = stds, ranp = ranp, 
                   omega = Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE],
                   correlation =  correlation, gamma = gamma, tau =  tau,
                   epsilon = err, sigbar = sigbar, p.se = p.se , 
                   hgamma = hgamma, notscale = notscale, Hb = Hb)
    for (j in 1:J) {
      XBr[[j]][theRows, ] <- crossprod(t(Xa[[j]][theRows, , drop = FALSE]), b$br) 
      if (fixed) XBf[[j]][theRows, ] <- crossprod(t(Xc[[j]][theRows, , drop = FALSE]), br) 
    }
    if (get.bi) bi[i, , ] <- if(fixed) t(rbind(br, b$br)) else t(b$br)
  }
  ## Make log L
  if (fixed){
    EXB <- mapply(function(x, y) exp(x + y), XBf, XBr, SIMPLIFY = FALSE)
  } else {
    EXB <- lapply(XBr, function(x) exp(x))
  }
  SEXB  <- suml(EXB)
  Pntir <- lapply(EXB, function(x) x / SEXB)
  Pnr   <- suml(mapply("*", Pntir, y, SIMPLIFY =  FALSE))
  if (panel) Pnr <- apply(Pnr, 2, tapply, id, prod)
  Pn <- apply(Pnr, 1, mean)
  if (get.bi)  Qnr  <-  Pnr / (Pn * R)
  lnL <- if (panel) sum(log(Pn) * weights[!duplicated(id)]) else sum(log(Pn) * weights)
  
  ## Make Gradient
  if (gradient){
    lambda <- mapply(function(y, p) y - p, y, Pntir, SIMPLIFY =  FALSE) 
    Qnr  <-  Pnr / (Pn * R) 
    if(panel) Qnr <-  Qnr[as.character(id), ] 
    eta  <- lapply(lambda, function(x) x * Qnr)
    
    dUdf <- dUdb <- dUds <- dUdt <- dUdg <- dUdh <- vector(mode = 'list', length = J)
    for (j in 1:J){
      dUdt[[j]] <- matrix(NA, N, K) 
      dUdg[[j]] <- matrix(NA, N, Ka)  
      dUdb[[j]] <- matrix(NA, N, Ka)
      if (fixed) dUdf[[j]] <- matrix(NA, N, Kc)
      if (het)   dUdh[[j]] <- matrix(NA, N, K)
      dUds[[j]] <- matrix(NA, N, length(stds))
    }
    for (i in 1:nind){
      if (panel){
        anid <- theIds[i]
        theRows <- which(id == anid)
      }
      else theRows <- i
      err <- epsilon[((i - 1) * R + 1):(i * R)]
      if (het){
        Hb <- drop(tcrossprod(H[i, , drop = FALSE], t(het.par)))
        sigma.nr <- exp(sigbar + Hb + tau * err) 
      } else {
        sigma.nr <- exp(sigbar + tau * err)
        Hb <- NULL
      }  
      if (fixed){
        d.muf <- d.tauf <- d.hetf <- matrix(NA, Kc, R)
        for(k in 1:Kc){
          ns <- notscale[names(delta[k])] # Look for name
          if (ns == 0) {
            d.muf[k, ]  <- sigma.nr
            d.tauf[k, ] <- sigma.nr * (- p.se + err) * delta[k]
            if (het) d.hetf[k, ] <- sigma.nr * delta[k]
          } else {
            d.muf[k, ]  <- 1
            d.tauf[k, ] <- 0
            if (het) d.hetf[k, ] <- 0
          }
        } 
      }   
      b <- MakeGcoef(beta = beta.bar, stds = stds, ranp = ranp, 
                     omega = Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE],
                     correlation =  correlation, gamma = gamma, tau =  tau,
                     epsilon = err, sigbar = sigbar, p.se = p.se , 
                     hgamma = hgamma, notscale =  notscale, Hb = Hb)
      for (j in 1:J) {
        if (fixed){
          dUdf[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], d.muf)
          dUdt[[j]][theRows, ] <- cbind(tcrossprod(eta[[j]][theRows, , drop = FALSE], d.tauf), 
                                        tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.tau))
          if (het) dUdh[[j]][theRows, ] <- cbind(tcrossprod(eta[[j]][theRows, , drop = FALSE], d.hetf), 
                                                 tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.het))
        } else {
          dUdt[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.tau)
          if (het) dUdh[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.het)
        } 
        dUdg[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.gamma)
        dUdb[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.mu)
        dUds[[j]][theRows, ] <- tcrossprod(eta[[j]][theRows, , drop = FALSE], b$d.stds)
      }   
    }
    if (correlation){
      vecX <- c()
      for (i in 1:Ka){
        vecX <- c(vecX, i:Ka)
      }
      Xac <- lapply(Xa,  function(x) x[, vecX])
    } else{
      Xac <- Xa  
    }
    
    if (fixed){
      tempt <- mapply(function(x, y, z) x * cbind(y, z), dUdt, Xc, Xa, SIMPLIFY =  FALSE)
    } else {
      tempt <- mapply(function(x, y) x * y, dUdt, Xa, SIMPLIFY =  FALSE)
    }
    gbart <-  suml(lapply(tempt, function(x) rowSums(x))) 
    if(het){ 
      if (fixed){
        temph <- mapply(function(x, y, z) x * cbind(y, z), dUdh, Xc, Xa, SIMPLIFY =  FALSE)
      } else {
        temph <- mapply(function(x, y) x * y, dUdh, Xa, SIMPLIFY =  FALSE)
      }
      if(panel) H <- H[id, ]
      ghet <- suml(lapply(temph, function(x) rowSums(x))) * H
    } else ghet <- numeric()
    tempg <- mapply("*", dUdg, Xa, SIMPLIFY =  FALSE)
    gbarg <- suml(lapply(tempg, function(x) rowSums(x))) 
    gbarfi <- if (fixed) suml(mapply("*", Xc, dUdf, SIMPLIFY =  FALSE)) else numeric()
    gbarmi <- suml(mapply("*", Xa, dUdb, SIMPLIFY =  FALSE))
    gbarvi <- suml(mapply("*", Xac, dUds, SIMPLIFY =  FALSE))
    gbari  <- cbind(gbarfi, gbarmi, ghet, gbarvi, gbart, gbarg)
    colnames(gbari) <- names(theta)
    attr(lnL, 'gradient') <- gbari * weights
  }
  if (get.bi) {
    attr(lnL, "prob.alt") <- sapply(Pntir, function(x) apply(x, 1, mean))
    attr(lnL, "prob.ind") <- Pn
    attr(lnL, "bi") <- bi
    attr(lnL, 'Qir') <- Qnr
  }
  lnL 
}

## Likelihood function for Latent class model
ll.mlogitlc <- function(theta, y, X, H, Q, id = NULL, weights = NULL,
                        gradient = TRUE, get.bi =  FALSE){
  # Parameters and globals
  K     <- ncol(X[[1]])
  N     <- nrow(X[[1]])
  J     <- length(X)
  panel <- !is.null(id)
  if (panel) if (length(weights) == 1) weights <- rep(weights, N)
  
  beta  <- matrix(theta[1L:(K * Q)], nrow = K, ncol = Q)
  rownames(beta) <- colnames(X[[1]]); colnames(beta) <- paste("class", 1:Q, sep = ':')
  gamma <- theta[-c(1L:(K * Q))]
    
  if (get.bi) bi <- t(beta)
  
  # Make weigths
  ew <- lapply(H, function(x) exp(crossprod(t(x), gamma)))
  sew <- suml(ew)
  Wnq <- lapply(ew, function(x){v <- x / sew; 
                                v[is.na(v)] <- 0;
                                as.vector(v)})
  Wnq <- Reduce(cbind, Wnq) # N*Q matrix: Probability for individual n in segement q
  
  # Make Multinomial Probability
  ep <- vector(mode = "list", length = J)
  for (j in 1:J) ep[[j]] <- matrix(NA, N, Q)
  for (j in 1:J) ep[[j]] <- exp(tcrossprod(X[[j]], t(beta))) 
  
  sep  <- suml(ep)
  Pnjq <- lapply(ep, function(x) x / sep) # list of N * Q
  Pnq  <- suml(mapply("*", Pnjq, y, SIMPLIFY = FALSE)) # Selected Probability
  if (panel) Pnq <- apply(Pnq, 2, tapply, id, prod)
  WPnq <- Wnq * Pnq # N * Q
  Ln   <- apply(WPnq, 1, sum)
  if (get.bi)  Qnr  <-  WPnq / Ln
  lnL <- if (panel) sum(log(Ln) * weights[!duplicated(id)]) else sum(log(Ln) * weights)
  
  
  ## Gradient
  if (gradient){
    lambda <- mapply(function(y, p) y - p, y, Pnjq, SIMPLIFY =  FALSE) # J list of N * Q
    Qnr  <-  WPnq / Ln # N * Q
    if (panel) Qnr <- Qnr[id, ] 
    eta  <- lapply(lambda, function(x) x * Qnr) # J list of N * Q
    etar <- lapply(eta,  function(x) x[, rep(1:Q, each = K)])
    Xg   <- lapply(X,  function(x) x[, rep(1:K, Q)]) # J list of N * (K*Q)
    grad.beta <- suml(mapply("*", Xg, etar, SIMPLIFY = FALSE))
    
    if (panel) {
      Wnq <- Wnq[id, ]
      H   <- lapply(H, function(x) x[id, ])
    }
    Wg <- vector(mode = "list", length = Q)
    IQ <- diag(Q)
    #for(q in 1:Q) Wg[[q]] <- matrix(NA, N, 1)
    for(q in 1:Q) Wg[[q]]<-  rowSums(Qnr * (repRows(IQ[q, ], N) - repCols(Wnq[, q], Q)))
    grad.gamma <- suml(mapply("*", H, Wg, SIMPLIFY = FALSE)) 
    gari <- cbind(grad.beta, grad.gamma)
    colnames(gari) <- names(theta)
    attr(lnL, "gradient") <- gari * weights
  }
  if (get.bi) {
    if (panel) Wnq <- Wnq[id, ]
    Pw <- lapply(Pnjq, function(x) x * Wnq)
    attr(lnL, "prob.alt") <- sapply(Pw, function(x) apply(x, 1, sum))
    attr(lnL, "prob.ind") <- Ln
    attr(lnL, "bi") <- bi
    attr(lnL, 'Qir') <- Qnr # WPnq / Ln
  }
  lnL 
}

## Likelihoood function for mixture of normals
ll.mnlogit <- function(theta, y, X, H, Q,
                       id =  NULL, ranp, R, correlation, weights = NULL,
                       haltons =  NULL, seed = 12345, 
                       gradient =  TRUE, get.bi = FALSE){
  # Get globals
  K <- ncol(X[[1]])
  J <- length(X)
  N <- nrow(X[[1]])
  panel <- !is.null(id)
  if (panel){
    n <- length(unique(id))
    if (length(weights) == 1) weights <- rep(weights, N)
  }  
  # Get parameters
  beta  <- matrix(theta[1L:(K * Q)], nrow = K, ncol = Q) # Matrix K * Q
  nstds <- if (!correlation) K * Q else (0.5 * K * (K + 1)) * Q
  stds  <- matrix(theta[(K * Q + 1):(K * Q + nstds)], ncol = Q)
  rownames(beta) <- colnames(X[[1]])
  colnames(beta) <- colnames(stds) <- paste("class", 1:Q, sep = ':')
  gamma <- theta[-c(1L:(K * Q + nstds))]
  
  # Get discrete weights
  ew <- lapply(H, function(x) exp(crossprod(t(x), gamma)))
  sew <- suml(ew)
  Wnq <- lapply(ew, function(x){v <- x / sew; 
                                v[is.na(v)] <- 0;
                                as.vector(v)})
  Wnq <- Reduce(cbind, Wnq) # N*Q matrix: Probability for individual n in segement q
  
  ## Make random draws
  set.seed(seed)
  Omega <- make.draws(R * ifelse(panel, n, N), K, haltons)
  
  # Make multinomial probability
  XBr <- vector(mode = 'list', length = J)
  for (j in 1:J) XBr[[j]] <- array(NA, dim = c(N, R, Q))
  nind <- ifelse(panel, n, N)
  if (panel) theIds <- unique(id)
  if (get.bi) bi <- array(NA, dim = c(nind, R, Q, K), 
                          dimnames = list(NULL, NULL, NULL, colnames(X[[1]]))) 
  for (i in 1:nind){
    if (panel){
      anid <- theIds[i]
      theRows <- which(id == anid)
    }
    else theRows <- i
    for (q in 1:Q){
      bq <- Makeh.rcoef(beta[, q], stds[, q], ranp, Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE], 
                       correlation, Pi = NULL, Slist = NULL, mvar = NULL)
      for (j in 1:J) {
        XBr[[j]][theRows, , q] <- crossprod(t(X[[j]][theRows, , drop = FALSE]), bq$br) 
      }
      if (get.bi) bi[i,, q,] <- t(bq$br)
    } 
  }
  
  EXB <- lapply(XBr, function(x) exp(x))
  SEXB  <- suml.array(EXB)
  Pntirq <- lapply(EXB, function(x) x / SEXB) 
  Pnrq   <- suml.array(mapply("*", Pntirq, y, SIMPLIFY =  FALSE))
  if (panel) Pnrq <- apply(Pnrq, c(2, 3), tapply, id, prod) 
  Pnq <- apply(Pnrq, c(1, 3), mean)
  WPnq <- Wnq * Pnq 
  Ln   <- apply(WPnq, 1, sum)
  if (get.bi)  Qir <- list(wnq = Wnq, Ln = Ln, Pnrq = Pnrq)
  lnL <- if (panel) sum(log(Ln) * weights[!duplicated(id)]) else sum(log(Ln) * weights)
  
  ## Gradient
  if (gradient){
    lambda <- mapply(function(y, p) y - p, y, Pntirq, SIMPLIFY =  FALSE) # J list of N * R * Q
    Wnq.mod  <- aperm(repmat(Wnq / Ln, dimen = c(1, 1, R)), c(1, 3, 2)) # n * R * Q
    Qnq.mod  <-  Wnq.mod * Pnrq # n * R * Q
    if (panel) Qnq.mod <- Qnq.mod[id,,] 
    eta  <- lapply(lambda, function(x) x * Qnq.mod) # J list of N * R * Q
    
    dUdb <- dUds <- vector(mode = 'list', length = J)
    for (j in 1:J){
      dUdb[[j]] <- array(NA, dim = c(N, K, Q))
      dUds[[j]] <- array(NA, dim = c(N, nrow(stds), Q))
    }
    
    for (i in 1:nind){
      if (panel){
        anid <- theIds[i]
        theRows <- which(id == anid)
      }
      else theRows <- i
      for (q in 1:Q){
        bq <- Makeh.rcoef(beta[, q], stds[, q], ranp, Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE], 
                          correlation, Pi = NULL, Slist = NULL, mvar = NULL)
        for (j in 1:J) {
          dUdb[[j]][theRows,, q] <- tcrossprod(eta[[j]][theRows,, q, drop = TRUE], bq$d.mu)
          dUds[[j]][theRows,, q] <- tcrossprod(eta[[j]][theRows,, q, drop = TRUE], bq$d.sigma)
        }
      } 
    }
    
    if (correlation){
      vecX <- c()
      for (i in 1:K){
        vecX <- c(vecX, i:K)
      }
      Xac <- lapply(X,  function(x) x[, vecX])
    } else{
      Xac <- X  
    }
    Xr   <- lapply(X,  function(x) x[, rep(1:K, Q)]) # J list of N * (K*Q)
    Xacr <- lapply(Xac,  function(x) x[, rep(1:ncol(Xac[[1]]), Q)]) # J list of N * (K*Q)
    dUdb <- lapply(dUdb,  function(x) matrix(x, nrow = N))
    dUds <- lapply(dUds,  function(x) matrix(x, nrow = N))
    grad.beta <- suml(mapply("*", Xr, dUdb, SIMPLIFY =  FALSE)) / R
    grad.stds <- suml(mapply("*", Xacr, dUds, SIMPLIFY =  FALSE)) / R
    
    # weight gradient
    Qnq <- WPnq / Ln
    if (panel) {
      Wnq <- Wnq[id, ]
      H   <- lapply(H, function(x) x[id, ])
      Qnq <- Qnq[id, ]
    }
    Wg <- vector(mode = "list", length = Q)
    IQ <- diag(Q)
    #for(q in 1:Q) Wg[[q]] <- matrix(NA, N, 1)
    for(q in 1:Q) Wg[[q]]<-  rowSums(Qnq * (repRows(IQ[q, ], N) - repCols(Wnq[, q], Q)))
    grad.gamma <- suml(mapply("*", H, Wg, SIMPLIFY = FALSE)) 
    
    gari <- cbind(grad.beta, grad.stds, grad.gamma)
    colnames(gari) <- names(theta)
    attr(lnL, "gradient") <- gari * weights
  }
  if (get.bi) {
    Pnjq <- lapply(Pntirq, function(x) apply(x, c(1, 3), mean))
    if (panel) Wnq <- Wnq[id, ]
    Pw <- lapply(Pnjq, function(x) x * Wnq)
    attr(lnL, "prob.alt") <- sapply(Pw, function(x) apply(x, 1, sum))
    attr(lnL, "prob.ind") <- Ln
    attr(lnL, "bi") <- bi
    attr(lnL, 'Qir') <- Qir
    attr(lnL, 'Wnq') <- Wnq
  }
  lnL
}


