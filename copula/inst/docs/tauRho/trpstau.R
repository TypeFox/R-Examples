library(copula) ##, lib.loc="../../copula.Rcheck")
# library(splines)
library(pspline)   ## avoids manually trying knots or df
# library(cobs)    ## allows pointwise constraint

source("gridsetup.R")


######################################################
## test with frank
######################################################
nsim <- 500  ## sample kendall's tau takes much longer to compute

thetaGrid <- seq(-.999, .999, by = .001)
frankTrFuns <- list(forwardTransf = function(x, ss) tanh(x / ss),
                    backwardTransf = function(x, ss) ss * atanh(x),
                    forwardDer = function(x, ss) (1 - tanh(x / ss)^2 ) / ss
                    )

## as alpha increases, tau goes to 1 much slower than rho; so ss is chosen to be bigger
.frankTau <- getAssoMeasFun(frankCopula(0), thetaGrid, nsim, "kendall", c(-1, 0, 1), c(-1, 0, 1), frankTrFuns, ss = 80, symmetrize = TRUE)


frankTauFun <- function(alpha) {
  ss <- .frankTau$ss
  forwardTransf <- .frankTau$trFuns$forwardTransf
  theta <- forwardTransf(alpha, ss)
  c(.frankTau$assoMeasFun$valFun(theta))
}

frankdTau <- function(alpha) {
  ss <- .frankTau$ss
  forwardTransf <- .frankTau$trFuns$forwardTransf
  forwardDer <- .frankTau$trFuns$forwardDer
  theta <- forwardTransf(alpha, ss)
  c(.frankTau$assoMeasFun$valFun(theta, 1)) * forwardDer(alpha, ss)
}

save(.frankTau, frankTauFun, frankdTau,
     file = "frankTau.rda")

########################################################
## try with plackett
## deal alpha <= 1 only;
## results for alpha > 1 can be obtained from 1 / alpha
########################################################
nsim <- 10000

plackettTrFuns <- list(forwardTransf = function(x, ss) x^(ss),
                       backwardTransf = function(x, ss) x^(1 / ss),
                       forwardDer = function(x, ss) ss * x^(ss - 1)
                       )

thetaGrid <- seq(.001, 1.5, by = .001)
.plackettTau <- getAssoMeasFun(plackettCopula(1), thetaGrid, nsim, "kendall", c(0, 1), c(-1, 0), plackettTrFuns, ss = 1/2, symmetrize = FALSE)

plackettTauFun <- function(alpha) {
  ss <- .plackettTau$ss
  forwardTransf <- .plackettTau$trFuns$forwardTransf
  valFun <- .plackettTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  idx <- theta <= 1
  val <- alpha
  val[idx] <- valFun(theta[idx])
  val[!idx] <- - valFun(1 / theta[!idx])
  val
  ## c(ifelse(theta <= 1, valFun(theta), -valFun(1/theta)))
}

plackettdTau <- function(alpha) {
  ss <- .plackettTau$ss
  forwardTransf <- .plackettTau$trFuns$forwardTransf
  forwardDer <- .plackettTau$trFuns$forwardDer
  valFun <- .plackettTau$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  idx <- theta <= 1
  val <- alpha
  val[idx] <- valFun(theta[idx], 1) * forwardDer(alpha[idx], ss)
  val[!idx] <-  valFun(1/theta[!idx], 1) * forwardDer(alpha[!idx], ss) / theta[!idx]^2
  val
  ## c(ifelse(alpha <= 1, valFun(theta, 1) * forwardDer(alpha, ss),
  ##         valFun(1/theta, 1) * forwardDer(alpha, ss) / theta^2))
}

save(.plackettTau, plackettTauFun, plackettdTau,
     file = "plackettTau.rda")
