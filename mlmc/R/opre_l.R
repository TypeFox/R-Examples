# This code is derived and adapted from the original GPL-2 Matlab version by
# Mike Giles.  See http://people.maths.ox.ac.uk/~gilesm/mlmc/

sig_dW <- function(x, dW, h) {
  dW[2,] <- -0.5*dW[1,] + sqrt(0.75)*dW[2,]

  c(sqrt(pmax(0,x[2,]))*x[1,]*dW[1,],
    exp(-5*h)*0.25*sqrt(pmax(0,x[2,]))*dW[2,]);
}

mu <- function(x, h) {
  m <- c(0.05*x[1,],
         ((1-exp(-5*h))/h)*(0.04-x[2,]))
}

#' Financial options using an Euler-Maruyama discretisation
#'
#' Financial options based on scalar geometric Brownian motion and Heston models,
#' similar to Mike Giles' original 2008 Operations Research paper, using an
#' Euler-Maruyama discretisation
#'
#' This function is based on GPL-2 'Matlab' code by Mike Giles.
#'
#' @param l the level to be simulated.
#' @param N the number of samples to be computed.
#' @param option the option type, between 1 and 5.  The options are: \describe{
#'   \item{1 = European call;}{}
#'   \item{2 = Asian call;}{}
#'   \item{3 = lookback call;}{}
#'   \item{4 = digital call;}{}
#'   \item{5 = Heston model.}{}
#' }
#'
#' @author Louis Aslett <aslett@stats.ox.ac.uk>
#' @author Mike Giles <Mike.Giles@maths.ox.ac.uk>
#' @author Tigran Nagapetyan <nagapetyan@stats.ox.ac.uk>
#'
#' @references
#' M.B. Giles. Multilevel Monte Carlo path simulation. \emph{Operations Research}, 56(3):607-617, 2008.
#'
#' @examples
#' \dontrun{
#' # These are similar to the MLMC tests for the original
#' # 2008 Operations Research paper, using an Euler-Maruyama
#' # discretisation with 4^l timesteps on level l.
#' #
#' # The differences are:
#' # -- the plots do not have the extrapolation results
#' # -- two plots are log_2 rather than log_4
#' # -- the new MLMC driver is a little different
#' # -- switch to X_0=100 instead of X_0=1
#'
#' M    <- 4 # refinement cost factor
#' N0   <- 1000 # initial samples on coarse levels
#' Lmin <- 2 # minimum refinement level
#' Lmax <- 6 # maximum refinement level
#'
#' test.res <- list()
#' for(option in 1:5) {
#'   if(option==1) {
#'     cat("\n ---- Computing European call ---- \n")
#'     N      <- 2000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
#'   } else if(option==2) {
#'     cat("\n ---- Computing Asian call ---- \n")
#'     N      <- 2000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
#'   } else if(option==3) {
#'     cat("\n ---- Computing lookback call ---- \n")
#'     N      <- 2000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.01, 0.02, 0.05, 0.1, 0.2)
#'   } else if(option==4) {
#'     cat("\n ---- Computing digital call ---- \n")
#'     N      <- 3000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.02, 0.05, 0.1, 0.2, 0.5)
#'   } else if(option==5) {
#'     cat("\n ---- Computing Heston model ---- \n")
#'     N      <- 2000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
#'   }
#'
#'   test.res[[option]] <- mlmc.test(opre_l, M, N, L, N0, Eps, Lmin, Lmax, option=option)
#'
#'   # print exact analytic value, based on S0=K
#'   T   <- 1
#'   r   <- 0.05
#'   sig <- 0.2
#'   K   <- 100
#'
#'   d1  <- (r+0.5*sig^2)*T / (sig*sqrt(T))
#'   d2  <- (r-0.5*sig^2)*T / (sig*sqrt(T))
#'
#'   if(option==1) {
#'     val <- K*( pnorm(d1) - exp(-r*T)*pnorm(d2) )
#'     cat(sprintf("\n Exact value: %f, MLMC value: %f \n", val, test.res[[option]]$P[1]))
#'   } else if(option==3) {
#'     k   <- 0.5*sig^2/r
#'     val <- K*( pnorm(d1) - pnorm(-d1)*k - exp(-r*T)*(pnorm(d2) - pnorm(d2)*k) )
#'     cat(sprintf("\n Exact value: %f, MLMC value: %f \n", val, test.res[[option]]$P[1]))
#'   } else if(option==4) {
#'     val <- K*exp(-r*T)*pnorm(d2)
#'     cat(sprintf("\n Exact value: %f, MLMC value: %f \n", val, test.res[[option]]$P[1]))
#'   }
#'
#'   # plot results
#'   plot(test.res[[option]])
#' }
#' }
#'
#' # The level sampler can be called directly to retrieve the relevant level sums:
#' opre_l(l=7, N=10, option=1)
#'
#' @importFrom stats rnorm
#' @export
opre_l <- function(l, N, option) {
  M <- 4

  T   <- 1
  r   <- 0.05
  sig <- 0.2
  K   <- 100

  nf <- M^l
  nc <- nf/M

  hf <- T/nf
  hc <- T/nc

  sums <- rep(0, 6)

  for(N1 in seq(1, N, by=10000)) {
    N2 <- min(10000, N-N1+1)

    #
    # GBM model
    #
    if(option<5) {
      X0 <- K

      Xf <- rep(X0, N2)
      Xc <- Xf

      Af <- 0.5*hf*Xf
      Ac <- 0.5*hc*Xc

      Mf <- Xf
      Mc <- Xc

      if(l==0) {
        dWf <- sqrt(hf)*rnorm(N2)
        Xf  <- Xf + r*Xf*hf + sig*Xf*dWf
        Af <- Af + 0.5*hf*Xf
        Mf <- min(Mf,Xf)
      } else {
        for (n in 1:nc){
          dWc <- rep(0, N2)
          for (m in 1:M){
            dWf <- sqrt(hf)*rnorm(N2)
            dWc <- dWc + dWf
            Xf  <- Xf + r*Xf*hf + sig*Xf*dWf
            Af  <- Af + hf*Xf
            Mf  <- pmin(Mf,Xf)
          }
          Xc <- Xc + r*Xc*hc + sig*Xc*dWc
          Ac <- Ac + hc*Xc
          Mc <- pmin(Mc,Xc)
        }
        Af <- Af - 0.5*hf*Xf
        Ac <- Ac - 0.5*hc*Xc
      }

      if(option==1) {
        Pf <- pmax(0,Xf-K)
        Pc <- pmax(0,Xc-K)
      } else if(option==2) {
        Pf <- pmax(0,Af-K)
        Pc <- pmax(0,Ac-K)
      } else if(option==3) {
        beta <- 0.5826 # special factor for offset correction
        Pf <- Xf - Mf*(1-beta*sig*sqrt(hf))
        Pc <- Xc - Mc*(1-beta*sig*sqrt(hc))
      } else if(option==4) {
        Pf <- K * 0.5*(sign(Xf-K)+1)
        Pc <- K * 0.5*(sign(Xc-K)+1)
      }
      #
      # Heston model
      #
    } else {
      Xf <- matrix(c(K, 0.04), nrow=2, ncol=N2)
      Xc <- Xf

      if(l==0) {
        dWf <- sqrt(hf)*rnorm(N2)
        Xf  <- Xf + mu(Xf,hf)*hf + sig_dW(Xf,dWf,hf)
      } else {
        for(n in 1:nc) {
          dWc <- matrix(0, nrow=2, ncol=N2)
          for(m in 1:M) {
            dWf <- sqrt(hf)*rnorm(N2)
            dWc <- dWc + dWf
            Xf  <- Xf + mu(Xf,hf)*hf + sig_dW(Xf,dWf,hf)
          }
          Xc <- Xc + mu(Xc,hc)*hc + sig_dW(Xc,dWc,hc)
        }
      }

      Pf <- pmax(0,Xf[1,]-K)
      Pc <- pmax(0,Xc[1,]-K)
    }

    Pf <- exp(-r*T)*Pf
    Pc <- exp(-r*T)*Pc

    if(l==0) {
      Pc <- 0
    }

    sums[1] <- sums[1] + sum(Pf-Pc)
    sums[2] <- sums[2] + sum((Pf-Pc)^2)
    sums[3] <- sums[3] + sum((Pf-Pc)^3)
    sums[4] <- sums[4] + sum((Pf-Pc)^4)
    sums[5] <- sums[5] + sum(Pf)
    sums[6] <- sums[6] + sum(Pf^2)
  }
  sums
}
