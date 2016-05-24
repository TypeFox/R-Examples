## sample from a given copula and compute sample tau or rho
TauRhoSample <- function(copula, nsim = 10000,
                         method = c("spearman", "kendall")) {
  u <- rCopula(nsim, copula)
  u <- na.omit(u)
  cor(u, method = method)[1,2]
}

## for evCopula, tau is expressed as a 1-dim integral with A and A''
tauEvCopula <- function(copula) {
  integrand <- function(x) x * (1 - x) / A(copula, x) * dAdu(copula, x)$der2
  integrate(integrand, 0, 1)$value
}

## for evCopula, rho is expressed as a 1-dim integral with A
rhoEvCopula <- function(copula) {
  integrand <- function(x) 1 / (A(copula, x) + 1)^2
  12 * integrate(integrand, 0, 1)$value - 3
}

## This function returns a matching grid between copula param theta and
## rho (or tau) that will be fed to iPsi for pspline fitting
getGrid <- function(copula, paramGrid, nsim = 10000,
                    method = c("spearman", "kendall"),
                    backwardTransf, ss) {
  alphaGrid <- backwardTransf(paramGrid, ss)
  ngrid <- length(alphaGrid)
  rhoGrid <- rep(NA, ngrid)
  for (i in 1:ngrid) {
    ## cat("alpha = ", alphaGrid[i], "\n")
    copula@parameters[1] <- alphaGrid[i] ## note: it wouldn't work for tCopula without [1] because df will be removed
    if (nsim == 0) {
      ## for evCopula, do integral instead of sampling
      if (method == "spearman") rhoGrid[i] <- rhoEvCopula(copula)
      else rhoGrid[i] <- tauEvCopula(copula)
    }
    else rhoGrid[i] <- TauRhoSample(copula, nsim, method)
  }
  list(param = paramGrid, val = rhoGrid)
}

## This function takes in the matching grid and
## returns a function that can give rho (or tau) and its derivative
## heavy: the weight given to known points (cannot be too big, it
##        breaks sm.spline
iPsi <- function(grid, norder = 3,
                   paramKnown, valKnown, heavy=999,...) {
  good <- !(grid$param %in% paramKnown)
  param <- c(paramKnown, grid$param[good])
  ## give known points weight heavy and other points weight 1
  val <- c(valKnown, grid$val[good])
  weight <- c(rep(heavy, length(valKnown)), rep(1, sum(good)))

  ## fit from pspline
  fm <- sm.spline(param, val, weight, norder = norder)# , method=3) # default method is 3 if cv == FALSE and df == NULL and spar == NULL
  ## fit from cobs
  ## con <- cbind(rep(0, length(paramKnown)), paramKnown, valKnown)
  ## fc <- cobs(grid$param, grid$val, "increase", lambda = -1, pointwise = con)
  rm(good, param, val, weight)
  valFun <- function(x, nderiv = 0) {
    c(predict(fm, x, nderiv = nderiv )) ## c is added on Jan. 20, 2009
    ## otherwise, it returns a matrix. Why?
  }
  ## list(fm = fm, fc = fc, valFun = valFun)
  list(fm = fm, valFun = valFun)
}

#############d##############################################################
## thetaGrid: a grid on [-1, 1] or [0, 1] or a valid range of rho or tau
## method: can be "spearman" or "kendall", which is passed to getGrid
## paramKnown and valKnwon: take in known info about the map
## trFuns: a list containing backward and forward transform function
## ss: a tuning parameter for the transform
## symmetrize: if TRUE, make the function symmetric about zero.
############################################################################
## the function returns a list of
##   assocmeasure:  from iPsi, a list containing the fitted object
##                  from sm.spline and a function using prediction to
##                  approximate rhoFun or tauFun and its derivative
##   trFuns and ss
############################################################################
getAssoMeasFun <- function(copula,
                           thetaGrid,
                           nsim,
                           method = "spearman",
                           paramKnown, valKnown,
                           trFuns,
                           ss,
                           symmetrize = FALSE) {
  ## trFuns is a list of functions
  forwardTransf <- trFuns$forwardTransf
  backwardTransf <- trFuns$backwardTransf
  forwardDer <- trFuns$forwardDer

  alphaGrid <- backwardTransf(thetaGrid, ss)
  mygrid <- getGrid(copula, thetaGrid, nsim, method, backwardTransf, ss)
  mygrid$param <- thetaGrid
  if (symmetrize) mygrid$val <- sign(mygrid$param) * (abs(mygrid$val) + abs(mygrid$val[rev(1:length(alphaGrid))])) / 2

  assoMeasFun <- iPsi(mygrid, norder = 3, paramKnown = paramKnown, valKnown = valKnown)

  list(assoMeasFun = assoMeasFun, trFuns = trFuns, ss = ss)
}
