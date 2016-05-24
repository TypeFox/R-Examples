#' @import expm
## #' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom stats dnbinom dpois pbeta pgamma rbeta rgamma rnbinom rpois runif
NULL

library(expm)
library(reshape2)

#' repmat
#'
#' This function replicates the matlab function repmat
#' @param X Target matrix
#' @param m New row dimension
#' @param n New column dimension
## #' @export
## #' @examples
## #' repmat(X, m, n)
repmat <- function(X, m, n) {
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  return(matrix(t(matrix(X, mx, nx*n)), mx*m, nx*n, byrow=T))
}

#' sensorMMPP
#'
#' This function provides the main MCMC inference engine
#' @param N Matrix of count data; axis 0 is the number of time intervals per day and axis 1 is the number of days in the data.
#' @param priors List with parameter values of prior distributions
#' @param ITERS Iteration controls: total # of iterations and # used for burn-in
#' @param EQUIV Parameter sharing controls <- c(S1, S2):
#' S1 <- force sharing of delta (day effect) among days,
#' S2 <- force sharing of eta (time of day) among days,
#' Values: 1 (all days share), 2 (weekdays/weekends), 3 (none)
#' @export
## #' @examples
## #' sensorMMPP(N, priors, c(50, 10), c(3, 3))
sensorMMPP <- function(N,
                       priors=list(aL=1,
                                   bL=1,
                                   aD=matrix(0, 1, 7)+5,
                                   aH=matrix(0, nrow=48, ncol=7)+1,
                                   z00=.99*10000,
                                   z01=.01*10000,
                                   z10=.25*10000,
                                   z11=.75*10000,
                                   aE=5,
                                   bE=1/3,
                                   MODE=0),
                       ITERS=c(50, 10),
                       EQUIV=c(3, 3)) {
  Nt <- replace(N, N == -1, NA)
  lenN <- dim(N)[2]
  priors$aL <- mean(N[N>-1])
  priors$bL <- 1
  # priors$aD <- c(mean(apply(Nt[, seq(1, lenN, 7)], 2, mean)),
  #                mean(apply(Nt[, seq(2, lenN, 7)], 2, mean)),
  #                mean(apply(Nt[, seq(3, lenN, 7)], 2, mean)),
  #                mean(apply(Nt[, seq(4, lenN, 7)], 2, mean)),
  #                mean(apply(Nt[, seq(5, lenN, 7)], 2, mean)),
  #                mean(apply(Nt[, seq(6, lenN, 7)], 2, mean)), 
  #                mean(apply(Nt[, seq(7, lenN, 7)], 2, mean)))
  # priors$aH <- matrix(cbind(apply(Nt[, seq(1, lenN, 7)], 1, mean),
  #                           apply(Nt[, seq(2, lenN, 7)], 1, mean),
  #                           apply(Nt[, seq(3, lenN, 7)], 1, mean),
  #                           apply(Nt[, seq(4, lenN, 7)], 1, mean),
  #                           apply(Nt[, seq(5, lenN, 7)], 1, mean),
  #                           apply(Nt[, seq(6, lenN, 7)], 1, mean),
  #                           apply(Nt[, seq(7, lenN, 7)], 1, mean)),
  #                     nrow=48)
  # priors$aD <- c(mean(N[, seq(1, lenN, 7)][N[, seq(1, lenN, 7)] > -1]),
  #                mean(N[, seq(2, lenN, 7)][N[, seq(2, lenN, 7)] > -1]),
  #                mean(N[, seq(3, lenN, 7)][N[, seq(3, lenN, 7)] > -1]),
  #                mean(N[, seq(4, lenN, 7)][N[, seq(4, lenN, 7)] > -1]),
  #                mean(N[, seq(5, lenN, 7)][N[, seq(5, lenN, 7)] > -1]),
  #                mean(N[, seq(6, lenN, 7)][N[, seq(6, lenN, 7)] > -1]),
  #                mean(N[, seq(7, lenN, 7)][N[, seq(7, lenN, 7)] > -1]))

  Niter <- ITERS[1]
  Nburn <- ITERS[2]
  Nplot <- ITERS[3]

  Z <- matrix(0, dim(N)[1], dim(N)[2])
  N0 <- pmax(N, 1)
  NE <- matrix(0, dim(N)[1], dim(N)[2])
  L <- (N+5)/2
  M <- matrix(c(.999, .001, .5, .5), nrow=2)
  xs <- seq(0, 1, 80)
  Nd <- 7
  Nh <- dim(N)[1]
  samples <- list()
  samples$L <- vector("list", Niter)
  samples$Z <- vector("list", Niter)
  samples$M <- vector("list", Niter)
  samples$N0 <- vector("list", Niter)
  samples$NE <- vector("list", Niter)
  samples$logp_NgLM <- vector("list", Niter)
  samples$logp_NgLZ <- vector("list", Niter)
  samples$L <- array(0, dim=c(dim(L)[1], dim(L)[2], Niter))
  samples$Z <- array(0, dim=c(dim(Z)[1], dim(Z)[2], Niter))
  samples$M <- array(0, dim=c(dim(M)[1], dim(M)[2], Niter))
  samples$N0 <- array(0, dim=c(dim(N0)[1], dim(N0)[2], Niter))
  samples$NE <- array(0, dim=c(dim(NE)[1], dim(NE)[2], Niter))
  samples$logp_NgLM <- matrix(0, 1, 50)
  samples$logp_NgLZ <- matrix(0, 1, 50)

  for (iter in 1:Niter+Nburn) {
    print(iter)
    L <- draw.L.given.N0(N0, priors, EQUIV)
    c(Z, N0, NE) := draw.Z.given.NLM(N, L, M, priors)
    M <- draw.M.given.Z(Z, priors)

    if (iter > Nburn) {
      samples$L[, ,iter-Nburn] <- L
      samples$Z[, ,iter-Nburn] <- Z
      samples$M[, ,iter-Nburn] <- M
      samples$N0[, ,iter-Nburn] <- N0
      samples$NE[, ,iter-Nburn] <- NE
      samples$logp_NgLM[iter-Nburn] <- prob.N.given.LM(N, L, M, priors)
      samples$logp_NgLZ[iter-Nburn] <- prob.N.given.LZ(N, L, Z, priors)
    }

    c(logpC, logpGD, logpGDz) := logp(N, samples, priors, iter-Nburn, EQUIV)
    logpC <- logpC/log(2)
    logpGD <- logpGD/log(2)
    logpGDz <- logpGDz/log(2)
    samples$logpC <- logpC
    samples$logpGD <- logpGD
  }
  return(
    list(L=melt(apply(samples$L, c(1, 2), mean)[, 1:7])$value,
    Z=melt(apply(samples$Z, c(1, 2), mean))$value))
}

#' dirichlet.log.pdf
#'
#' The log of the probability density function for the Dirichlet distribution.
#' Returns the belief that the probabilities of K rival events are x_i
#' given that each event has been observed A_i - 1 times.
#' @param K.probs vector of probabilities
#' @param A vector of concentration parameters.
#' @export
## #' @examples
## #' dirichlet.log.pdf(K.probs, A)
dirichlet.log.pdf <- function(K.probs, A) {
  if (length(K.probs) == 1) {
    # It _seems_ like this should be 0, since ln(1) = 0. However, Graham has it as one.
    return(1)
  }

  return(sum((A-1)*log(K.probs+.00000001))-sum(lgamma(A))+lgamma(sum(A)))
}

#' dirichlet.pdf
#'
#' The probability density function for the Dirichlet distribution.
#' Returns the belief that the probabilities of K rival events are x_i
#' given that each event has been observed A_i - 1 times.
#' @param K.probs Vector of probabilities
#' @param A Vector of concentration parameters.
#' @export
## #' @examples
## #' dirichlet.pdf(K.probs, A)
dirichlet.pdf <- function(K.probs, A) {
  if (length(K.probs) == 1) {
    return(1)
  }

  return(exp(dirichlet.log.pdf(K.probs, A)))
}

#' draw.Z.given.NLM
#'
#' Sample the given N, L, M
#' @param N Matrix of count data; axis 0 is the number of time intervals per day and axis 1 is the number of days in the data.
#' @param L Matrix containing the rate functions at every time slice
#' @param M Matrix containing the estimated transition probabilities for each iteration
#' @param priors List with parameter values of prior distributions
#' @export
## #' @examples
## #' draw.Z.given.NLM(N, L, M, priors)
#'
draw.Z.given.NLM <- function(N, L, M, priors) {
  N0 <- N
  NE <- 0*N
  Z <- 0*N
  ep <- 1e-50
  
  # First sample Z, N0, and NE:
  PRIOR <- M%^%100%*%as.vector(c(1, 0))
  po <- matrix(0, 2, length(N))
  p  <- matrix(0, 2, length(N))
  for (t in 1:length(N)) {
    if (N[t] != -1) {
      po[1, t] <- dpois(N[t], L[t])+ep
      po[2, t] <- sum(dpois(0:N[t], L[t])*
                        dnbinom(rev(0:N[t]), priors$aE,priors$bE/(1+priors$bE))
                      ) + ep
    }
    else {
      po[1, t] <- 1
      po[2, t] <- 1
    }
  }

  # Compute forward (filtering) posterior marginals
  p[, 1] <- PRIOR*po[, 1]
  p[, 1] <- p[, 1]/sum(p[, 1])
  for (t in 2:length(N)) {
    p[, t] <- (M%*%p[, t-1])*po[, t]
    p[, t] <- p[, t]/sum(p[, t])
  }

  # Do backward sampling
  for (t in rev(1:length(N))) {
    if (runif(1) > p[1, t]) {  # if event at time t
      if (N[t] != -1) {
        Z[t] <- 1
        # Likelihood of all possible event/normal combinations
        # (all possible values of N(E))
        ptmp <- log(dpois(0:N[t], L[t])) + log(dnbinom(rev(seq(0, N[t], 1)), priors$aE, priors$bE/(1+priors$bE)))
        ptmp <- ptmp-max(ptmp)
        ptmp <- exp(ptmp)
        ptmp <- ptmp/sum(ptmp)

        # Draw sample of N0 and compute NE
        N0[t] <- min(which(cumsum(ptmp) >= runif(1)))-1
        NE[t] <- N[t]-N0[t]
    }
    else {
        Z[t] <- 1
        N0[t] <- rpois(1, L[t])
        NE[t] <- rnbinom(1, priors$aE, priors$bE/(1+priors$bE))
      }
    }
    else {
      if (N[t] != -1) {
        Z[t] <- 0
        N0[t] <- N[t]
        NE[t] <- 0  # No event at time t
      }
      else {
        Z[t] <- 0
        N0[t] <- rpois(1, L[t])
        NE[t] <- 0
      }
    }

    ptmp <- matrix(0, 2, 1)
    ptmp[Z[t]+1] <- 1  # compute backward influence
    if (t>1) {
      p[, t-1] <- p[, t-1]*(t(M)%*%ptmp)
      p[, t-1] <- p[, t-1]/sum(p[, t-1])
    }
  }

  out <- list()
  out$Z <- Z
  out$N0 <- N0
  out$NE <- NE
  return(out)
}

#' draw.M.given.Z
#'
#' Sample the M given Z
#' @param Z Binary vector indicationg the presence (1) or absence (0) of an event at every time slice
#' @param prior Parameter values of a particular prior distribution
#' @export
## #' @examples
## #' draw.M.given.Z(Z, prior)
#'
draw.M.given.Z <- function(Z, prior) {
  n01 <- length(which(Z[1:length(Z)-1] == 0 & Z[2:length(Z)] == 1))
  n0 <- length(which(Z[1:length(Z)-1] == 0))
  n10 <- length(which(Z[1:length(Z)-1] == 1 & Z[2:length(Z) == 0]))
  n1 <- length(which(Z[1:length(Z)-1] == 1))
  z0 <- rbeta(1, n01+prior$z01, n0-n01+prior$z00)
  z1 <- rbeta(1, n10+prior$z10, n1-n10+prior$z11)
  M <- t(matrix(c(1-z0, z1, z0, 1-z1), nrow=2, ncol=2))
  return(M)
}

#' draw.L.given.N0
#'
#' Sample the L given N0
#' @param N0 Matrix containing the estimated baseline Poisson distribution for activity 
#' @param prior Parameter values of a particular prior distribution
#' @param EQUIV Parameter sharing controls <- c(S1, S2):  S1 <- force sharing of delta (day effect) among days, S2 <- force sharing of eta (time of day) among days, Values: 1 (all days share), 2 (weekdays/weekends), 3 (none)
#' @export
## #' @examples
## #' draw.L.given.N0(N0, prior, EQUIV)
#'
draw.L.given.N0 <- function(N0, prior, EQUIV) {
  Nd <- 7
  Nh <- dim(N0)[1]

  # First: Overall Average Rate
  if (prior$MODE) {
    L0 <- (sum(N0)+prior$aL)/(length(N0)+prior$bL)
  }
  else {
    L0 <- rgamma(1, shape=sum(N0)+prior$aL, scale=1/(length(N0)+prior$bL))
  }
  L <- matrix(0, dim(N0)[1], dim(N0)[2]) + L0

  # Second: Day Effect
  D <- matrix(0, 1, Nd)
  for(i in 1:length(D)) {
    alpha <- sum(N0[, seq(i, dim(N0)[2], 7)])+prior$aD[i]
    if (prior$MODE) {
      D[i] <- (alpha-1)  # mode of Gamma(a, 1) distribution
    }
    else {
      D[i] <- rgamma(1, alpha, scale=1)
    }
  }

  # Third: Time-of-day Effect
  A <- matrix(0, Nh, Nd)
  for (tau in 1:(dim(A)[2])) {
    for (i in 1:dim(A)[1]) {
      alpha <- sum(N0[i, seq(tau, dim(N0)[2], 7)])+prior$aH[i]
      if (prior$MODE) {
        A[i, tau] <- (alpha-1)  # mode of Gamma(a, 1) distribution
      }
      else {
        A[i, tau] <- rgamma(1, alpha, scale=1)
      }
    }
  }

  # Enforce parameter sharing
  if (EQUIV[1] == 1) {
    D[1:7] <- 1
  }
  else if (EQUIV[1] == 2) {
    D[c(1, 7)] <- mean(D[c(1, 7)])
    D[2:6] <- mean(D[2:6])
    D <- D/mean(D)
  }
  else if (EQUIV[1] == 3) {
    D <- D/mean(D)
  }

  ### FIX THIS
  # tau(t)
  if (EQUIV[2] == 1) {
    A[, 1:7] <- repmat(matrix(rowMeans(A)), 1, dim(A)[2])
  }
  else if (EQUIV[2] == 2) {
    A[, c(1, 7)] <- repmat(matrix(rowMeans(A[, c(1, 7)])), 1, 2)
    A[, 2:6] <- repmat(matrix(rowMeans(A[, 2:6])), 1, 5)
  }
  else if (EQUIV[2] == 3) {
    A <- A
  }

  for (tau in 1:dim(A)[2]) {
    A[, tau] = A[, tau]/mean(A[, tau])
  }

  # Compute L(t)
  for (d in 1:dim(L)[2]) {
    for (t in 1:dim(L)[1]) {
      dd <- (d-1)%%7+1
      L[t, d] <- L0*D[dd]*A[t, dd] #fix this line
    }
  }

  return(L)
}

#' logp
#'
#' Estimates the marginal likelihood of the data using the samples
#' @param N Matrix of count data; axis 0 is the number of time intervals per day and axis 1 is the number of days in the data.
#' @param samples List of different samples at all time periods
#' @param priors List with parameter values of prior distributions
#' @param iter Number of iterations over which to calcuate likelihood.
#' @param EQUIV Parameter sharing controls <- c(S1, S2):  S1 <- force sharing of delta (day effect) among days, S2 <- force sharing of eta (time of day) among days, Values: 1 (all days share), 2 (weekdays/weekends), 3 (none)
#' @export
## #' @examples
## #' logp(N, samples, priors, iter, EQUIV)
#'
logp <- function(N, samples, priors, iter, EQUIV) {
  tmp <- samples$logp_NgLZ[1:iter]
  tmpm <- mean(tmp)
  temp <- tmp - tmpm
  logpGDz <- log(1/mean(1/exp(tmp))) + tmpm # Gelfand-Dey estimate
  logpGD  <- log(1/mean(1/exp(tmp))) + tmpm # Gelfand-Dey estimate, marginalizing over Z

  Lstar <- apply(samples$L, c(1, 2), mean)
  Mstar <- apply(samples$M, c(1, 2), mean)
  logp_LMgN <- matrix(0, 1, iter)
  logp_LM <- prob.L.given.N0(Lstar, vector(), priors, EQUIV) + prob.M.given.Z(Mstar, 0, priors)
  logp_NgLM <- prob.N.given.LM(N, Lstar, Mstar, priors)
  for (ii in 1:iter) {
    logp_LMgN[ii] <- prob.L.given.N0(Lstar, samples$N0[, ,ii], priors, EQUIV)+prob.M.given.Z(Mstar, samples$Z[, ,ii], priors)
  }

  tmpm <- mean(exp(logp_LMgN))+tmpm
  logpC <- logp_NgLM + logp_LM - logp_LMgN  # Chib estimate
}

#' prob.M.given.Z
#'
#' This function evaluates p(M|Z)
#' @param M Matrix containing the estimated transition probabilities for each iteration
#' @param Z Binary vector indicationg the presence (1) or absence (0) of an event at every time slice
#' @param prior Parameter values of a particular prior distribution
#' @export
## #' @examples
## #' prob.M.given.Z(M, Z, prior)
#'
prob.M.given.Z <- function(M, Z, prior) {
  z1 <- M[1, 2]
  z0 <- M[2, 1]

  if (length(Z) != 0) {
    n01 <- sum(Z[1:length(Z)-1] == 0 & Z[2:length(Z)] == 1)
    n0  <- sum(Z[1:length(Z)-1] == 0)
    n10 <- sum(Z[1:length(Z)-1] == 1 & Z[2:length(Z)] == 0)
    n1  <- sum(Z[1:length(Z)-1] == 1)
  }
  else {
    n01 <- 0
    n0  <- 0
    n10 <- 0
    n1  <- 0
  }
  logp <- log(pbeta(z0, n01+prior$z01, n0-n01+prior$z00) + .001) +
          log(pbeta(z1, n10+prior$z10, n1-n10+prior$z11) + .001)

  return(logp)
}

#' prob.L.given.N0
#'
#' This function evaluates p(L|N0)
#' @param L Matrix containing the rate functions at every time slice
#' @param N0 Matrix containing the estimated baseline Poisson distribution for activity 
#' @param prior Parameter values of a particular prior distribution
#' @param EQUIV Parameter sharing controls <- c(S1, S2):  S1 <- force sharing of delta (day effect) among days, S2 <- force sharing of eta (time of day) among days, Values: 1 (all days share), 2 (weekdays/weekends), 3 (none)
#' @export
## #' @examples
## #' prob.L.given.N0(L, N0, prior, EQUIV)
#'
prob.L.given.N0 <- function(L, N0, prior, EQUIV) {
  L0 <- mean(L)
  Nd <- 7
  Nh <- dim(L)[1]
  A <- matrix(0, Nh, Nd)
  D <- rep(NA, Nd)
  for (i in 1:Nd) {
    D[i] <- mean(L[, i]/L0)
  }
  for (i in 1:Nd) {
    for (j in 1:Nh) {
      A[j, i] <- L[j, i]/L0/D[i]
    }
  }
  logp <- 0
  
  # Enforce parameter sharing
  paD <- prior$aD
  aD <- matrix(0, 1, Nd)
  paH <- prior$aH
  aH <- matrix(0, Nh, Nd)
  if (length(N0) != 0) {
    for (i in 1:Nd) {
      aD[i] <- sum(N0[, seq(i, dim(N0)[2], Nd)])  # Fix this line
    }
    for (i in 1:Nd) {
      for (j in 1:Nh) {
        aH[j, i] <- sum(N0[j, seq(i, dim(N0)[2], Nd)])
      }
    }
  }

  # d(t)
  if (EQUIV[1] == 1) {
    D <- sum(D)
    paD <- sum(paD)
    aD <- sum(aD)
  }
  else if (EQUIV[1] == 2) {
    D <- c(D[1] + D[7], sum(D[2:6]))
    paD <- c(paD[1] + paD[7], sum(paD[2:6]))
    aD <- c(aD[1] + aD[7], sum(aD[2:6]))
  }
  else if (EQUIV[1] == 3) {
    D <- D
    paD <- paD
    paH <- paH
  }

  # tau(t)
  if (EQUIV[2] == 1) {
    A <- matrix(rowSums(A)/Nd)
    aH <- matrix(rowSums(aH))
    paH <- matrix(rowSums(paH))
  }
  else if (EQUIV[2] == 2) {
    A <- matrix(c((A[, 1] + A[, 7])/2, rowSums(A[, 2:6])/5))
    aH <- matrix(c(aH[, 1] + aH[, 7], rowSums(aH[, 2:6])))
    paH <- matrix(c(paH[, 1] + paH[, 7], rowSums(paH[, 2:6])), nrow=1)
  }
  else if (EQUIV[2] == 3) {
    A <- A
    aH <- aH
    paH <- paH
  }

  logp <- logp + 
          log(pgamma(L0,
                     sum(sum(N0) + .00000001) + prior$aL,
                     1/(length(N0)+prior$bL))
              + .000000001)
  logp <- logp +
          log(dirichlet.pdf(D/Nd, aD + paD)
              + .000000001)

  for (i in 1:dim(A)[2]) {
    logp <- logp +
            log(dirichlet.pdf(A[, i]/Nh, aH[, i]+paH[, i])
                + .0000000001)
  }

  return(logp)
}

#' prob.N.given.LM
#'
#' This function evaluates p(N|L, M)
#' @param N Matrix of count data; axis 0 is the number of time intervals per day and axis 1 is the number of days in the data.
#' @param L Matrix containing the rate functions at every time slice
#' @param M Matrix containing the estimated transition probabilities for each iteration
#' @param prior Parameter values of a particular prior distribution
#' @export
## #' @examples
## #' prob.N.given.LM(N, L, M, prior)
#'
prob.N.given.LM <- function(N, L, M, prior) {
  PRIOR <- M%^%100%*%as.vector(c(1, 0))
  po <- matrix(0, 2, length(N))
  p  <- matrix(0, 2, length(N))

  for (t in 1:length(N)) {
    if (N[t] != -1) {
      po[1, t] <- dpois(N[t], L[t])
      po[2, t] <- sum(dpois(0:N[t], L[t])*dnbinom(rev(0:N[t]), prior$aE, prior$bE/(1+prior$bE)))
    }
    else {
      po[1, t] <- 1
      po[2, t] <- 1
    }
  }

  p[, 1] <- PRIOR*po[, 1]
  sp <- sum(p[, 1])
  logp <- log(sp)
  p[, 1] <- p[, 1]/sp
  for (t in 2:length(N)) {
    p[, t] <- (M%*%p[, t-1])*po[, t]
    sp <- sum(p[, t])
    logp <- logp + log(sp)
    p[, t] <- p[, t]/sp
  }

  return(logp)
}

#' prob.N.given.LZ
#'
#' This function evaluates p(N|L, Z)
#' @param N Matrix of count data; axis 0 is the number of time intervals per day and axis 1 is the number of days in the data.
#' @param L Matrix containing the rate functions at every time slice
#' @param Z Binary vector indicationg the presence (1) or absence (0) of an event at every time slice
#' @param prior Parameter values of a particular prior distribution
#' @export
## #' @examples
## #' prob.N.given.LZ(N, L, Z, prior)
#'
prob.N.given.LZ <- function(N, L, Z, prior) {
  logp <- 0
  for (t in 1:length(N)) {
    if (N[t] != -1) {
      if (Z[t] == 0) {
        logp <- logp + log(dpois(N[t], L[t]))
      }
      else {
        logp <- logp + log(sum(dpois(0:N[t], L[t])*dnbinom(rev(0:N[t]), prior$aE, prior$bE/(1+prior$bE))))
      }
    }
  }

  return(logp)
}

#' :=
#'
#' This function allows multiple assignments
#' @param lhs The value to be assigned to the left-hand argument (on the left)
#' @param rhs The value to be assigned to the right-hand argument (on the left)
#' @export
#' @examples
#' c(a, b):=c(3, 4)
#'
':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1) {
    lhs <- lhs[-1]
  }
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL))
  }

  if (is.function(rhs) || is(rhs, 'formula')) {
    rhs <- list(rhs)
  }
  if (length(lhs) > length(rhs)) {
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  }
  for (i in 1:length(lhs)) {
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  }

  return(invisible(NULL))
}
