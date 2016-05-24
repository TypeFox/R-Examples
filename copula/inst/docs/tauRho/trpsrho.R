library(copula) ##, lib.loc="../../copula.Rcheck")
# library(splines)
library(pspline)   ## avoids manually trying knots or df
# library(cobs)    ## allows pointwise constraint

source("gridsetup.R")


nsim <- 10000
######################################################
## test with frank
######################################################
thetaGrid <- seq(-.999, .999, by = .001)
frankTrFuns <- list(forwardTransf = function(x, ss) tanh(x / ss),
                    backwardTransf = function(x, ss) ss * atanh(x),
                    forwardDer = function(x, ss) (1 - tanh(x / ss)^2 ) / ss
                    )


.frankRho <- getAssoMeasFun(frankCopula(0), thetaGrid, nsim, "spearman", c(-1, 0, 1), c(-1, 0, 1), frankTrFuns, ss = 40, symmetrize = TRUE)


## using pointwise constraints does not neccesarily increase precision
## except at those constrained points
## frankRhoFun.cobs <- function(alpha) {
##   ss <- .frankRho$ss
##   forwardTransf <- .frankRho$trFuns$forwardTransf
##   theta <- forwardTransf(alpha, ss)
##   predict(.frankRho$assoMeasFun$fc, z = theta)[, "fit"]
## }

frankRhoFun <- function(alpha) {
  ss <- .frankRho$ss
  forwardTransf <- .frankRho$trFuns$forwardTransf
  theta <- forwardTransf(alpha, ss)
  c(.frankRho$assoMeasFun$valFun(theta))
}

frankdRho <- function(alpha) {
  ss <- .frankRho$ss
  forwardTransf <- .frankRho$trFuns$forwardTransf
  forwardDer <- .frankRho$trFuns$forwardDer
  theta <- forwardTransf(alpha, ss)
  c(.frankRho$assoMeasFun$valFun(theta, 1)) * forwardDer(alpha, ss)
}

## rhoTrue <- sapply(alphaGrid, function(x) rho(frankCopula(x)))
## dRhoTrue <- sapply(alphaGrid, function(x) dRho(frankCopula(x)))

save(.frankRho, frankRhoFun, frankdRho,
     file = "frankRho.rda")

########################################################
## try with t-4
########################################################
tTrFuns <- list(forwardTransf = function(x, ss) x,
                backwardTransf = function(x, ss) x,
                forwardDer = function(x, ss) rep(1, length(x))
                )
thetaGrid <- seq(-.999, .999, by = .001)

.t4Rho <- getAssoMeasFun(tCopula(0, df=4), thetaGrid, nsim, "spearman", c(-1, 0, 1), c(-1, 0, 1), tTrFuns, ss = 0, symmetrize = TRUE)

t4RhoFun <- function(alpha) {
  ss <- .t4Rho$ss
  forwardTransf <- .t4Rho$trFuns$forwardTransf
  theta <- forwardTransf(alpha, ss)
  c(.t4Rho$assoMeasFun$valFun(theta))
}

t4dRho <- function(alpha) {
  ss <- .t4Rho$ss
  forwardTransf <- .t4Rho$trFuns$forwardTransf
  forwardDer <- .t4Rho$trFuns$forwardDer
  theta <- forwardTransf(alpha, ss)
  c(.t4Rho$assoMeasFun$valFun(theta, 1)) * forwardDer(alpha, ss)
}

alphaGrid <- tTrFuns$backwardTransf(thetaGrid)
rhoTrue <- sapply(alphaGrid, function(x) rho(tCopula(x)))
dRhoTrue <-  6 / (pi * sqrt(4 - alphaGrid^2)) ## true for normal

save(.t4Rho, t4RhoFun, t4dRho,
     file = "t4Rho.rda")

## pdf("t4-rho.pdf", height=3, width=6, pointsize=9)
## par(mfrow=c(1,2), mgp=c(1.5, 0.5, 0), mar=c(3,3,0,0.5))
## plot(alphaGrid, rhoTrue, type="l", xlab=expression(theta), ylab=expression(rho(theta)))
## curve(t4RhoFun(x), add=TRUE, col="blue")

## curve(t4dRho(x), col="blue", xlab=expression(theta), ylab=bquote(partialdiff~ rho / partialdiff ~ theta))
## lines(alphaGrid, dRhoTrue)
## dev.off()

########################################################
## try with clayton
## deal with positive and negative alpha separately
########################################################
claytonTrFuns <- list(forwardTransf = function(x, ss)
                      ifelse(x <= 0, x, tanh(x / ss)),
                      backwardTransf = function(x, ss)
                      ifelse(x <= 0, x, ss * atanh(x)),
                      forwardDer = function(x, ss)
                      ifelse(x <= 0, 1, (1 - tanh(x / ss)^2) / ss)
                      )

## negative alpha
thetaGridNeg <- seq(-.999, -0.001, by = .001)
.claytonRhoNeg <- getAssoMeasFun(claytonCopula(0), thetaGridNeg, nsim, "spearman", c(-1, 0), c(-1, 0), claytonTrFuns, ss = 10, symmetrize = FALSE)

thetaGridPos <- seq(.001, .999, by = .001)
.claytonRhoPos <- getAssoMeasFun(claytonCopula(0), thetaGridPos, nsim, "spearman", c(0, 1), c(0, 1), claytonTrFuns, ss = 10, symmetrize = FALSE)

claytonRhoFun <- function(alpha) {
  ss <- .claytonRhoNeg$ss
  forwardTransf <- .claytonRhoNeg$trFuns$forwardTransf
  valFunNeg <- .claytonRhoNeg$assoMeasFun$valFun
  valFunPos <- .claytonRhoPos$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  c(ifelse(alpha <= 0, valFunNeg(theta), valFunPos(theta)))
}

claytondRho <- function(alpha) {
  ss <- .claytonRhoNeg$ss
  forwardTransf <- .claytonRhoNeg$trFuns$forwardTransf
  forwardDer <- .claytonRhoNeg$trFuns$forwardDer
  valFunNeg <- .claytonRhoNeg$assoMeasFun$valFun
  valFunPos <- .claytonRhoPos$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  c(ifelse(alpha <= 0, valFunNeg(theta, 1), valFunPos(theta, 1))) * forwardDer(alpha, ss)
}

save(.claytonRhoNeg, .claytonRhoPos, claytonRhoFun, claytondRho,
     file = "claytonRho.rda")

########################################################
## try with gumbel
## transform from (1, \infty) to (0, 1)
########################################################
n <- 10000

gumbelTrFuns <- list(forwardTransf = function(x, ss) tanh((x - 1) / ss),
                      backwardTransf = function(x, ss) ss * atanh(x) + 1,
                      forwardDer = function(x, ss) (1 - tanh( (x - 1)/ ss)^2) / ss
                      )

## negative alpha
thetaGrid <- seq(.001, .999, by = .001)
.gumbelRho <- getAssoMeasFun(gumbelCopula(1), thetaGrid, nsim, "spearman", c(0, 1), c(0, 1), gumbelTrFuns, ss = 20, symmetrize = FALSE)

gumbelRhoFun <- function(alpha) {
  ss <- .gumbelRho$ss
  forwardTransf <- .gumbelRho$trFuns$forwardTransf
  valFun <- .gumbelRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  c(valFun(theta))
}

gumbeldRho <- function(alpha) {
  ss <- .gumbelRho$ss
  forwardTransf <- .gumbelRho$trFuns$forwardTransf
  forwardDer <- .gumbelRho$trFuns$forwardDer
  valFun <- .gumbelRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  c(valFun(theta, 1)) * forwardDer(alpha, ss)
}

save(.gumbelRho, gumbelRhoFun, gumbeldRho,
     file = "gumbelRho.rda")
