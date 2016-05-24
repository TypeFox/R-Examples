##############################
# Make Random Draws
##############################

## Halton function from mlogit package
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
#' @import stats
make.draws <- function(R, Ka, haltons){
  # Create the matrix of random numbers
  if (!is.null(haltons)) {
    length.haltons <- rep(R,Ka)
    prime <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
               47, 53, 59, 61, 71, 73, 79, 83, 89, 97, 101, 103,
               107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
               173, 179, 181, 191, 193, 197, 199)
    drop.haltons <- rep(16, Ka)
    if (!is.na(haltons) && !is.null(haltons$prime)) {
      if (length(haltons$prime) != Ka) {
        stop("wrong number of prime numbers indicated")
      }
      else{
        prime <- haltons$prime
      }
      if (!is.na(haltons) && !is.null(haltons$drop)) {
        if (!length(haltons$drop) %in% c(1,Ka)) stop("wrong number of drop indicated")
        if (length(haltons$drop) == 1) {
          drop.haltons <- rep(haltons$drop, Ka)
        }
        else{
          drop.haltons <- haltons$drop
        }
      }
    }
    random.nb <- numeric(0)
    i <- 0
    for (i in 1:Ka) {
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
  K <- (-1 + sqrt(1 + 8 * length(x))) / 2
  mat <- matrix(0, K, K)
  mat[lower.tri(mat, diag = TRUE)] <- x
  mat
}

## Make epsilon draws
#' @import stats
Make.epsi <- function(N, bound.err){
  epsilon <- qnorm(pnorm(-bound.err) * (1 - runif(N)) + pnorm(bound.err) * runif(N))
  epsilon
}

## Make Random Coefficients
#' @import stats
Makeh.rcoef <- function(beta, sigma, ranp, omega, correlation, Pi = NULL, Slist = NULL, mvar){
  names.r <- names(beta) 
  Ka    <- nrow(omega) 
  R     <- ncol(omega)  
  br    <- matrix(NA, Ka, R) 
  d.mu  <- d.sigma <- br   
  
  beta  <- drop(beta)
  sigma <- drop(sigma)
  names(beta) <- names(sigma) <- names(beta)
  rownames(br) <- rownames(d.mu) <- rownames(d.sigma) <- rownames(omega) <- names(beta)
  
  if (!is.null(Pi)) {
    d.pis           <- matrix(NA, length(Pi), R)
    names.pi        <- rep(names(mvar), sapply(mvar, length))
    names(Pi) <- rownames(d.pis) <- names.pi
  }
  
  if (correlation) {
    rownames(omega) <- NULL
    L <- makeL(sigma)
    L.omega <- tcrossprod(L, t(omega))
    rownames(L.omega) <- names(beta)
    br <- beta + L.omega
    d.mu <- matrix(1, Ka, R)
    rownames(d.mu) <- names.r
    d.sigma <- omega[rep(1:Ka, Ka:1), ]   #d U/ d omega
    if (!is.null(Pi)) d.pis[,] <- 1
    for (i in 1:Ka) {
      var  <- names(beta)[i]
      distr <- ranp[var]
      sigi <- i + cumsum(c(0, (Ka - 1):1))[1:i]
      if (any(names(Slist) == var)) {
        pic <- as.vector(Pi[names(Pi) == var])
        br[var, ] <- br[var, ] + drop(tcrossprod(Slist[[var]], t(pic)))
      }
      if (distr == "cn") {
        br[var, ]  <- pmax(br[var, ], 0)
        d.mu[var, ] <- as.numeric(br[var, ] > 0)
        d.sigma[sigi, ] <- repRows(as.numeric(br[var, ] > 0), length(sigi))  * d.sigma[sigi, ]
        if (any(names(Slist) == var)) {
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <- repRows(as.numeric(br[var, ] > 0), hmt)
        } 
      }
      if (distr == "ln") {
        br[var, ] <- exp(br[i, ])
        d.mu[var, ] <- br[var, ]
        d.sigma[sigi, ] <- repRows(br[var, ], length(sigi)) * d.sigma[sigi, ]
        if (any(names(Slist) == var)) {
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <-  repRows(br[var, ], hmt)
        } 
      }    
    }
  } else {
    for (j in 1:Ka) {
      var  <- names(beta)[j]
      distr <- ranp[var]
      if (distr == "n") {
        if (any(names(Slist) == var)) {
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * omega[var,, drop = F]
        } else {
          br[var, ] <- beta[var] + sigma[var] * omega[var,, drop = F]
        }
        d.mu[var, ]    <- 1
        d.sigma[var, ] <- omega[var, ]
        if (any(names(Slist) == var)) d.pis[rownames(d.pis) == var, ] <- 1   
      }
      if (distr == "ln") {
        if (any(names(Slist) == var)) {
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- exp(beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * omega[var, , drop = F])
        } else {
          br[var, ] <- exp(beta[var] + sigma[var] * omega[var, , drop = F])
        } 
        d.mu[var, ]    <- br[var, , drop = F]
        d.sigma[var,]  <- d.mu[var, ] * omega[var, ]
        if (any(names(Slist) == var)) {
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <- repRows(d.mu[var, ], hmt)
        } 
      }
      if (distr == "cn") {
        if (any(names(Slist) == var)) {
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- pmax(beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * omega[var, , drop = F], 0)
        }else{
          br[var, ] <- pmax(beta[var] + sigma[var] * omega[var, , drop = F], 0)
        }
        d.mu[var, ] <- as.numeric(br[var, ] > 0)
        d.sigma[var, ] <- d.mu[var, ] * omega[var, ]
        if (any(names(Slist) == var)) {
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <-  repRows(as.numeric(br[var, ] > 0), hmt)
        } 
      }
      if (distr == "u") {
        etauni <- 2 * pnorm(omega[var,, drop = F]) - 1
        if(any(names(Slist) == var)) {
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- beta[var] + drop(tcrossprod(Slist[[var]], t(pic)))  + etauni * sigma[var]
        }else{
          br[var, ] <- beta[var] + etauni * sigma[var]
        }
        d.mu[var, ] <- 1
        d.sigma[var, ] <- etauni
        if (any(names(Slist) == var)) d.pis[rownames(d.pis) == var, ] <- 1
      }  
      if (distr == "t") {
        etauni <- pnorm(omega[var,, drop = F])
        eta05  <- etauni < 0.5
        d.mu[var, ] <- 1
        d.sigma[var, ] <- eta05 * (sqrt(2 * etauni) - 1) + 
          !eta05 * (1 - sqrt(2 * (1 - etauni)))
        if (any(names(Slist) == var)) {
          pic <- as.vector(Pi[names(Pi) == var])
          br[var, ] <- beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * d.sigma[var, ] 
        } else {
          br[var, ] <- beta[var] + sigma[var] * d.sigma[var, ]
        }
        if (any(names(Slist) == var)) d.pis[rownames(d.pis) == var, ] <- 1 
      }
      if (distr == "sb") {
        if (any(names(Slist) == var)) {
          pic <- as.vector(Pi[names(Pi) == var])
          btemp <- beta[var] + drop(tcrossprod(Slist[[var]], t(pic))) + sigma[var] * omega[var, , drop = F]
        } else {
          btemp <- beta[var] + sigma[var] * omega[var, , drop = F]
        } 
        br[var, ] <- exp(btemp) / (1 + exp(btemp))
        d.mu[var, ]    <- br[var, ] - br[var, ] ^ 2
        d.sigma[var,]  <- d.mu[var, ] * omega[var, ]
        if (any(names(Slist) == var)) {
          hmt <- nrow(d.pis[rownames(d.pis) == var, , drop = F])
          d.pis[rownames(d.pis) == var, ] <- repRows(d.mu[var, ], hmt)
        } 
      }
    }
  }
  if (!is.null(Pi)) {
    list(br = br, d.mu = d.mu, d.sigma = d.sigma, d.pis = d.pis)
  } else {
    list(br = br, d.mu = d.mu, d.sigma = d.sigma)
  }
}

## Make Generalized Random Coefficients
#' @import stats
MakeGcoef <- function(beta, stds = NULL, ranp, omega = NULL, correlation = FALSE, 
                      gamma = NULL, tau, epsilon, sigbar, p.se,
                      hgamma = c("direct", "indirect"), notscale = NULL, Hb = NULL)
{
  names.r <- names(beta)  
  Ka    <- nrow(omega) 
  R     <- ncol(omega)  
  br    <- d.mu <- matrix(NA, Ka, R)
  rownames(br) <- rownames(d.mu) <- rownames(omega) <-  names(stds) <- names.r
  if (is.null(Hb)) {
    sigma.nr <- exp(sigbar + tau * epsilon)
    d.het <- NULL
  } else {
    sigma.nr <- exp(sigbar + Hb + tau * epsilon)
    d.het <- matrix(NA, Ka, R)
    rownames(d.het) <- names.r
  }
  gamexp <- if (hgamma == "direct") gamma else exp(gamma) / (1 + exp(gamma))
  d.tau <- d.gamma <- matrix(NA, Ka, R)
  rownames(d.tau) <- rownames(d.gamma) <- names.r
  if (correlation) {
    L <- makeL(stds)
    Lomega <- tcrossprod(L, t(omega))
    d.stds <- omega[rep(1:Ka, Ka:1), ] 
    name.ds <- c()
    for (i in 1:Ka) name.ds <- c(name.ds, rep(i:Ka, 1))
    rownames(d.stds) <- names.r[name.ds]
    for (j in 1:Ka) {
      ns <- notscale[names.r[j]]
      picksig <- j + cumsum(c(0, (Ka - 1):1))[1:j]
      if (ns == 0) {
        br[j, ]      <- beta[j] * sigma.nr + gamexp * Lomega[j, ] + (1 - gamexp) * sigma.nr * Lomega[j, ]
        d.mu[j, ]    <- sigma.nr 
        d.tau[j, ]   <- sigma.nr * (-p.se + epsilon) * (beta[j] + ((1 - gamexp) * Lomega[j, ]))
        if (hgamma == "direct") {
          d.gamma[j, ] <- Lomega[j, ] * (1 - sigma.nr)
        } else {
          d.gamma[j, ] <- (exp(gamma) / (1 + exp(gamma)) ^ 2) * Lomega[j, ] * (1 - sigma.nr)
        }
        if (!is.null(Hb)) d.het[j, ] <- beta[j] * sigma.nr + (1 - gamexp) * sigma.nr * Lomega[j, ]
        d.stds[picksig, ] <- d.stds[picksig, ] * (gamexp + ((1 - gamexp) * repRows(sigma.nr, length(picksig))))
      } else {
        br[j, ] <- beta[j] + Lomega[j, ]
        d.mu[j, ]    <- 1 
        d.tau[j, ]   <- 0
        d.gamma[j, ] <- 0
        if (!is.null(Hb)) d.het[j, ] <- 0
      }
    } 
  } else {
    d.stds <- matrix(NA, Ka, R)
    rownames(d.stds) <- names.r
    for (j in 1:Ka) {
      var <- names.r[j]
      ns <- notscale[var]
      distr <- ranp[var]
      draws <- switch(distr,
                      "n" = omega[var, ],
                      "ln" = omega[var, ],
                      "t" = {etauni <- pnorm(omega[var,, drop = F])
                              eta05  <- etauni < 0.5
                              eta05 * (sqrt(2 * etauni) - 1) + !eta05 * (1 - sqrt(2 * (1 - etauni)))
                       },
                      "u" = 2 * pnorm(omega[var,, drop = F]) - 1
        )
      if (ns == 0) {
        br[var, ] <- beta[var] * sigma.nr + gamexp * stds[var] * draws + (1 - gamexp) * sigma.nr * stds[var] * draws
        d.mu[var, ]    <- sigma.nr
        d.tau[var, ]   <- sigma.nr * (-p.se + epsilon) * (beta[var] + ((1 - gamexp) * stds[var] * draws))
        if (hgamma == "direct") {
          d.gamma[var, ] <- stds[var] * draws * (1 - sigma.nr)
        } else {
          d.gamma[var, ] <- (exp(gamma) / (1 + exp(gamma)) ^ 2) * stds[var] * draws * (1 - sigma.nr)
        }
        #if (!is.null(Hb)) d.het[var, ] <- beta[var] * sigma.nr
        if (!is.null(Hb)) d.het[var, ] <- beta[var] * sigma.nr + (1 - gamexp) * sigma.nr * stds[var] * draws
        d.stds[var, ]  <-  draws * (gamexp + ((1 - gamexp) * sigma.nr))
        if (distr == "ln") {
         br[var, ] <- exp(br[var, , drop = FALSE])
         d.mu[var, ] <- br[var, ] * d.mu[var, ]
         d.tau[var, ] <- br[var, ] * d.tau[var, ] 
         d.gamma[var, ] <- br[var, ] * d.gamma[var, ] 
         d.stds[var, ] <-  br[var, ] * d.stds[var, ]
         if (!is.null(Hb)) d.het[var, ] <- br[var, ] * d.het[var, ]
        }
      } else {
        br[var, ] <- beta[var] + stds[var] * draws
        d.mu[var, ]    <- 1
        d.tau[var, ]   <- 0
        d.gamma[var, ] <- 0
        if (!is.null(Hb)) d.het[var, ] <- 0
        d.stds[var, ]  <- draws
        if (distr == "ln") {
          br[var, ] <- exp(br[var, , drop = FALSE])
          d.mu[var, ] <- br[var, ]
          d.stds[var, ] <- br[var, ] * d.stds[var, ]
        }
      }
    }
  }
  list(br = br, d.mu = d.mu, d.tau = d.tau, d.gamma = d.gamma, d.stds = d.stds, d.het = d.het) 
}  



