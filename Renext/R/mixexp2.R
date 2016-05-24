
##=====================================
## probability (distribution) function
##=====================================

pmixexp2 <- function(q,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta,
                     log = FALSE) {
  
  res <- prob1 * pexp(q, rate = rate1, log.p = FALSE) +
      (1-prob1) * pexp(q, rate = rate2, log.p = FALSE)

  if (log) res <- log(res)
  res
  
}

##==================
## density function
##==================

dmixexp2 <- function(x,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta,
                     log = FALSE) {
  
  res <- prob1 * dexp(x, rate = rate1, log = FALSE) +
      (1 - prob1) * dexp(x, rate = rate2, log = FALSE) 
  
  if (log) res <- log(res)
  res
  
  
}

##============
## simulation
##============

rmixexp2 <- function(n,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta) {
  
  g.true <- rbinom(n, size = 1, prob = prob1)
  
  x <- rep(NA, n)
  ind <- g.true == 1
  n1 <- sum(ind)
  if (n1) x[ind] <- rexp(n1, rate = rate1)
  if (n1 < n) x[!ind] <-  rexp(n-n1, rate = rate2)

  x
  
}

##=============
## hazard rate
##=============

hmixexp2 <- function(x,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta) {
                       
  f <- prob1 * dexp(x, rate = rate1, log = FALSE) +
      (1 - prob1) * dexp(x, rate = rate2, log = FALSE)
  F <- prob1 * pexp(x, rate = rate1, log.p = FALSE) +
      (1 - prob1) * pexp(x, rate = rate2, log.p = FALSE)
  f / (1 - F)
  
}

##=======================
## cumulater hazard rate
##=======================

Hmixexp2 <- function(x,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta) {
  
  F <- prob1 * pexp(x, rate = rate1, log.p = FALSE) +
      (1 - prob1) * pexp(x, rate = rate2, log.p = FALSE)
  
  -log(1-F)
  
}

##=======================
## quantile function
##=======================

qmixexp2 <- function(p,
                     prob1,
                     rate1 = 1.0,
                     rate2 = rate1 + delta,
                     delta) {

  tol <- 1e-9
  n <- length(p) 

  if (rate1 > rate2) {
    ## warning("using the 'qmixep2' function with rate2 < rate1") 
    prob1 <- 1-prob1
    ratep <- rate1
    rate1 <- rate2
    rate2 <- ratep
  }

  ## set the smallest rate to 1.0 in order to improve
  ## conditionning
  
  rate1b <- rate1
  rate1 <- 1.0
  rate2 <- rate2 / rate1b

  Hs <- - log(1.0 - p)
  
  ## Added in versions > 1.00 for p close to 1
  ind <- (p > 0.99)
  
  if (any(ind)) {
    Hs[ind] <- - log1p(-p[ind])
  }

  lprob1 <- log(prob1)

  xs <- rep(NA, n)
  nit <- rep(NA, n)
  rate.bar <- prob1 + (1-prob1)*rate2

  ## cat(sprintf("XXX rate1 = %16.10f,  rate1 = %16.10f, rate.bar = %16.10f, lprob1 = %16.10f\n",
  ##             rate1, rate2, rate.bar, lprob1))
  
  for (i in 1L:n) {
    
    H.star <- Hs[i]

    x.L <-  max(c(H.star / rate.bar, H.star + lprob1))

    ## cat(sprintf("INI H.star = %e, x.prov = %ef\n", H.star, x.L))
    
    H.L <- Hmixexp2(x.L, prob1 = prob1, rate1 = 1.0, rate2 = rate2) 
    h.L <- hmixexp2(x.L, prob1 = prob1, rate1 = 1.0, rate2 = rate2)
    
    cvg <- FALSE
    iter <- 1
    
    while ( !cvg && (iter < 30) ) {
    
      x.prov <- x.L + (H.star - H.L) / h.L 

      ## cat(sprintf("iter = %3d x.prov = %16.10f\n", iter, x.prov)) 
      
      H.prov <- Hmixexp2(x.prov,
                         prob1 = prob1,
                         rate1 = 1.0,
                         rate2 = rate2)
      
      if ( abs(H.prov - H.star) < tol) {
        x.star <- x.prov
        cvg <- TRUE
      } else {
        
        ## For convexity reasons, H.prov should always be <= H.star
        if (H.prov <= H.star) {
          x.L <- x.prov
          h.L <- hmixexp2(x.prov, prob1 = prob1, rate1 = 1.0, rate2 = rate2)
          H.L <- H.prov
        } else {
          h.L <- 1.1*h.L
          ## cat(sprintf("** H.star - H.prov = %12.8f\n", H.star-H.prov))        
        } 
       
      }
      
      iter <- iter + 1

    }
    
    if (cvg) xs[i] <- x.star
    else {
      cat(sprintf("i = %d  it. = %d\n x.L = %e\n
                  H.L= %e H.star = %e\n",
                  i, iter, x.L, H.L, H.star))
      cat(sprintf("prob1 = %e  rate1 = %e rate2 = %e\n",
             prob1, rate1, rate2))
      cat(sprintf("p = %e\n",
                  p[i]))
      stop("diverged i =", i, "\n")
    }
    nit[i] <- iter
  }

  xs <- xs / rate1b
  attr(xs, "nit") <- nit
  ## print(nit)
  xs
  
}
