library(copula)  ##, lib.loc="../../copula.Rcheck")
library(pspline) ## avoids manually trying knots or df

source("gridsetup.R")

######################################################
## test with galambos
######################################################
thetaGrid <- seq(0.01, .999, by = .001)
galambosTrFuns <- list(forwardTransf = function(x, ss) (tanh(x / ss))^0.33,
                    backwardTransf = function(x, ss) ss * atanh(x^(1/0.33)),
                    forwardDer = function(x, ss) 0.33 * (tanh(x / ss))^(0.33 - 1) * (1 - tanh(x / ss)^2 ) / ss
                    )

nsim <- 0
.galambosTau <- getAssoMeasFun(galambosCopula(0), thetaGrid, nsim, "kendall", c(0, 1), c(0, 1), galambosTrFuns, ss = 15, symmetrize = FALSE)

nsim <- 0
.galambosRho <- getAssoMeasFun(galambosCopula(0), thetaGrid, nsim, "spearman", c(0, 1), c(0, 1), galambosTrFuns, ss = 5, symmetrize = FALSE)


save(.galambosTau, .galambosRho,
     file = "galambos.rda")

######################################################
## test with huslerReiss
######################################################
thetaGrid <- seq(0.001, .999, by = .001)
huslerReissTrFuns <- list(forwardTransf = function(x, ss) (tanh(x / ss))^0.33,
                    backwardTransf = function(x, ss) ss * atanh(x^(1/0.33)),
                    forwardDer = function(x, ss) 0.33 * (tanh(x / ss))^(0.33 - 1) * (1 - tanh(x / ss)^2 ) / ss
                    )

nsim <- 0
.huslerReissTau <- getAssoMeasFun(huslerReissCopula(0), thetaGrid, nsim, "kendall", c(0, 1), c(0, 1), huslerReissTrFuns, ss = 15, symmetrize = FALSE)

nsim <- 0
.huslerReissRho <- getAssoMeasFun(huslerReissCopula(0), thetaGrid, nsim, "spearman", c(0, 1), c(0, 1), huslerReissTrFuns, ss = 5, symmetrize = FALSE)

save(.huslerReissTau, .huslerReissRho,
     file = "huslerReiss.rda")

######################################################
## test with tev, df = 4
######################################################
thetaGrid <- seq(0.001, .999, by = .001)
tevTrFuns <-  list(forwardTransf = function(x, ss) x^(ss),
                       backwardTransf = function(x, ss) x^(1 / ss),
                       forwardDer = function(x, ss) ss * x^(ss - 1)
                       )

nsim <- 0
.tevTau <- getAssoMeasFun(tevCopula(0, df=4), thetaGrid, nsim, "kendall", c(0, 1), c(0, 1), tevTrFuns, ss = 2, symmetrize = FALSE)

nsim <- 0
.tevRho <- getAssoMeasFun(tevCopula(0, df=4), thetaGrid, nsim, "spearman", c(0, 1), c(0, 1), tevTrFuns, ss = 1, symmetrize = FALSE)

save(.tevTau, .tevRho,
     file = "tev.rda")
