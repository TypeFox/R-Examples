### Estimate the price of different European options on multiple stocks and
### compute risk measures based on aggregate losses.

require(qrng)
require(copula)


### 1) Functions ###############################################################

##' @title Produces Samples of a d-dimensional Clayton Copula
##' @param u (n, d+1)-matrix of samples in [0,1]
##' @param theta parameter of Clayton copula
##' @return (n,d)-matrix of Clayton copula samples
##' @author Mathieu Cambou, Christiane Lemieux, Marius Hofert
rClaytonMO <- function(u, theta)
{
  if(!is.matrix(u)) u <- rbind(u)

  dim. <- dim(u)
  n <- dim.[1]
  d <- dim.[2] - 1
  V <- qgamma(u[,d+1], shape=1/theta) # n-vector, frailty component
  E <- qexp(u[,seq_len(d)]) # (n,d)-matrix

  (1 + E/matrix(rep(V, d), ncol=d))^(-1/theta) # (n,d)-matrix
}

##' @title Produces Samples of a d-dimensional Geometric Brownian Motion (GBM)
##' @param u (n,d)-matrix of samples in [0,1]
##' @param S0 d-vector with initial GBM levels
##' @param mu d-vector containing the drifts of the GBMs
##' @param sigma d-vector containing the (daily) volatilities of the GBMs
##' @param T time horizon of the samples (in number of days)
##' @return (n,d)-matrix of GBM samples (d independent BMs)
##' @author Mathieu Cambou, Christiane Lemieux, Marius Hofert
##' @note To add dependency between the GBMs, a copula sample can be passed via the
##'       first argument. Further generalizations are possible, e.g.,
##'       sampling a full path for pricing path-dependent options
rGeoBM <- function(u, S0, mu, sigma, T)
{
    stopifnot(0 < u, u < 1, length(mu) == (d <- length(S0)), mu >= 0,
              length(sigma) == d, sigma >= 0, T > 0)
    log.diff <- qnorm(u) * matrix(rep(sigma, each=n), ncol=d) # (n,d)-matrix; or t(t(qnorm(u))*sigma)
    log.drft <- (mu - sigma^2 / 2) * T # d-vector
    log. <- matrix(rep(log.drft, n), ncol=d, byrow=TRUE) + log.diff # (n,d)-matrix
    matrix(rep(S0, n), ncol=d, byrow=TRUE) * exp(log.) # S_t, t in 1,..,T; (n,d)-matrix
}

##' @title Risk measures Value-at-Risk (VaR), Expected Shortfall (ES)
##'        and Capital Allocated on Aggregate of Losses
##' @param x (n, d)-matrix of losses
##' @param alpha confidence level
##' @return 5-vector of VaR, ES, Alloc(X.first), Alloc(X.mid), Alloc(X.last)
##' @author Mathieu Cambou, Christiane Lemieux, Marius Hofert
risk.measures <- function(x, alpha)
{
  if(!is.matrix(x)) x <- rbind(x)
  n <- nrow(x)
  d <- ncol(x)

  aloss  <- rowSums(x) # n-vector of aggregated losses
  VaR <- quantile(aloss, probs=alpha, names=FALSE) # VaR estimate
  l <- sum(xcd <- aloss > VaR)
  ES <- mean(aloss[xcd]) / (1-alpha) # ES estimate

  Alloc.first <- x[aloss>=VaR, 1] %*% rep(1/n, l) / (l/n) # capital allocated to X_1
  Alloc.mid   <- x[aloss>=VaR, floor(d/2)] %*% rep(1/n, l) / (l/n) # capital allocated to X_{floor(d/2)}
  Alloc.last  <- x[aloss>=VaR,d] %*% rep(1/n, l) / (l/n) # Capital Allocated to X_d

  ## return estimated risk measures
  c(VaR=VaR, ES=ES, Alloc.first=Alloc.first, Alloc.mid=Alloc.mid,
    Alloc.last=Alloc.last)
}

##' @title Payoff Function for Multiple Stocks Options
##' @param K strike of the option
##' @param N notional of the option
##' @param S0 a d-vector with stocks' current prices
##' @param S a (n, d)-matrix with stocks' final prices
##' @param type string "call" or "put"
##' @param method option type; available are "basket", "worst-of" and "best-of"
##'        option
##' @return n-vector values of the basket option's payoff function
##' @author Mathieu Cambou, Marius Hofert, Christiane Lemieux
##' @note This is the undiscounted value of the payoff function
##'       evaluated at the stock prices at maturity.
payoff <- function(K, N, S0, S, type = c("call", "put"),
                   method = c("basket", "worst.of", "best.of"))
{
  stopifnot(K >= 0, N >= 0, S0 >= 0, S >= 0, length(S0) == ncol(S))
  type <- match.arg(type)
  method <- match.arg(method)
  perf <- switch(method,
                 "basket" = {
                   rowMeans(t(t(S)/S0))
                 },
                 "worst.of" = {
                   apply(t(t(S)/S0), 1, min)
                 },
                 "best.of" = {
                   apply(t(t(S)/S0), 1, max)
                 },
                 stop("Wrong 'method'"))
  N * pmax(0, if(type=="call") perf - K else K - perf)
}


### 2) Case Study ##############################################################

### 2.1) Define parameters #####################################################

n <- 1e5 # Monte Carlo sample size
d <- 4 # dimension

## Stochastic process parameters
sigma <- rep(0.2, d) # volatilities
r <- 0.0001 # continuously compounded short rate
S0 <- rep(100, d) # initial stocks' levels
K <- 1.1 # option strike
N <- 1000 # option notional
T <- 1 # time horizon
alpha <- 0.99 # confidence level for VaR, ES

## Copulas
tau <- 0.5 # Kendall's tau
## Clayton
family.C <- "Clayton"
th.C <- iTau(getAcop(family.C), tau) # corresponding parameter
clayton.cop <- onacopulaL(family.C, nacList=list(th.C, 1:d)) # Clayton copula
## t_3
family.t <- "t"
nu <- 3 # degrees of freedom
th.t <- iTau(ellipCopula(family.t, df=nu), tau) # corresponding parameter
t.cop <- ellipCopula(family.t, param=th.t, dim=d, df=nu) # define copula object


### 2.2) Sampling ##############################################################

## Uniform samples for CDM
set.seed(271)
U.CDM  <- matrix(runif(n*d), ncol=d) # pseudo
set.seed(271)
U.CDM. <- ghalton(n, d=d) # quasi

## Uniform samples for MO
set.seed(271)
U.MO  <- matrix(runif(n*(d+1)), ncol=d+1) # pseudo
set.seed(271)
U.MO. <- ghalton(n, d=d+1) # quasi


## t samples via CDM
U.t.CDM  <- rtrafo(U.CDM,  cop=t.cop, inverse=TRUE) # pseudo
U.t.CDM. <- rtrafo(U.CDM., cop=t.cop, inverse=TRUE) # quasi

## Clayton samples via CDM
U.C.CDM  <- rtrafo(U.CDM,  cop=clayton.cop, inverse=TRUE) # pseudo
U.C.CDM. <- rtrafo(U.CDM., cop=clayton.cop, inverse=TRUE) # quasi

## Clayton samples via MO
U.C.MO  <- rClaytonMO(U.MO,  theta=th.C) # pseudo
U.C.MO. <- rClaytonMO(U.MO., theta=th.C) # quasi


## Geometric Brownian Motion samples
S.t.CDM <- rGeoBM(U.t.CDM, S0=S0, mu=rep(r, d), sigma=sigma, T=T)
S.C.CDM <- rGeoBM(U.C.CDM, S0=S0, mu=rep(r, d), sigma=sigma, T=T)
S.C.MO  <- rGeoBM(U.C.MO,  S0=S0, mu=rep(r, d), sigma=sigma, T=T)

## Quasi-geometric Brownian Motion samples
S.t.CDM. <- rGeoBM(U.t.CDM., S0=S0, mu=rep(r, d), sigma=sigma, T=T)
S.C.CDM. <- rGeoBM(U.C.CDM., S0=S0, mu=rep(r, d), sigma=sigma, T=T)
S.C.MO.  <- rGeoBM(U.C.MO.,  S0=S0, mu=rep(r, d), sigma=sigma, T=T)


### 2.3) Functional Calculation ################################################

erT <- exp(-r*T)

## Using pseudo-random samples

## Basket call
basket.t.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM))
basket.C.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM))
basket.C.MO  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO))

## Worst of call
worst.of.t.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM, method="worst.of"))
worst.of.C.CDM <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM, method="worst.of"))
worst.of.C.MO  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO,  method="worst.of"))

## Risk measures
rm.t.CDM <- risk.measures(S.t.CDM, alpha)
rm.C.CDM <- risk.measures(S.C.CDM, alpha)
rm.C.MO  <- risk.measures(S.C.MO,  alpha)


## Using quasi-random samples

## Basket call
basket.t.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM.))
basket.C.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM.))
basket.C.MO.  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO.))

## Worst of call
worst.of.t.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.t.CDM., method="worst.of"))
worst.of.C.CDM. <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.CDM., method="worst.of"))
worst.of.C.MO.  <- erT * mean(payoff(K, N=N, S0=S0, S=S.C.MO.,  method="worst.of"))

## Risk measures
rm.t.CDM. <- risk.measures(S.t.CDM., alpha)
rm.C.CDM. <- risk.measures(S.C.CDM., alpha)
rm.C.MO.  <- risk.measures(S.C.MO.,  alpha)


### 2.4) Results ###############################################################

res <- array(, dim=c(4,2,2), dimnames=list(type=c("basket", "worst.of",
                                           paste0("VaR.", alpha), paste0("ES.", alpha)),
                             copula=c("Clayton", paste0("t", nu)),
                             method=c("CDM", "MO")))
res["basket",,]     <- matrix(c(basket.C.CDM., NA, basket.C.MO., basket.t.CDM.), ncol=2)
res["worst.of",,]   <- matrix(c(worst.of.C.CDM., NA, worst.of.C.MO., worst.of.t.CDM.), ncol=2)
res[paste0("VaR.", alpha),,] <- matrix(c(rm.C.CDM.[1], NA, rm.C.MO.[1], rm.t.CDM.[1]), ncol=2)
res[paste0("ES.", alpha),,]  <- matrix(c(rm.C.CDM.[2], NA, rm.C.MO.[2], rm.t.CDM.[2]), ncol=2)
res