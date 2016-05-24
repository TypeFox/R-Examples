#### A "fast" and simplified version of demo(TGforecasts)
#### see ../demo/TGforecasts.R
####     ~~~~~~~~~~~~~~~~~~~~~

require(simsalapar)
## also  require(fGarch)

(hasRmpi <- .Platform$OS.type!="windows" # <- as Rmpi produces errors on the winbuilder
 && require("Rmpi"))

### Variable list ##############################################################

## Create varlist (default type is 'frozen'; default expressions are generated)
## note: no n.sim here
vList <-
    varlist(forecaster = list(type="grid", expr = quote(italic(forecaster)),
                              value = c("statistician", "optimist", "pessimist")),
	    scoringfun = list(type="inner", expr = quote(S),
			      value = list(SE = function(x, y) mean((x-y)^2),
					   AE = function(x, y) mean(abs(x-y)))),
            ## GARCH(1,1): X_t = sigma_t * eps_t,
            ## sigma_t^2 = alpha_0 + alpha_1*X_{t-1}^2 + beta_1*sigma_{t-1}^2
            alpha0     = list(expr = quote(alpha[0]), value = 0.05),
            alpha1     = list(expr = quote(alpha[1]), value = 0.2),
            beta1      = list(expr = quote(beta[1]), value = 0.75),
            ## constant prediction of optimist:
            pred.optimist  = list(expr = "Optimist", value = 5),
            ## constant prediction of pessimist:
            pred.pessimist = list(expr = "Pessimist", value = 0.05),
            ## n
            n = list(value = 64),
            ## simulation method
            method = list(value = "fGarch"))

### doOne() including ingredient functions #####################################
##  -------
##  instead of 'local', use a construction like
##  doOne <- function(x, <nonGrids>, simGARCH11, predictTG)

doOne <- local({

    ## function for simulating a squared GARCH(1,1)
    simGARCH11 <- function(n, a0, a1, b1, method=c("fGarch", "Gneiting"))
    {
        method <- match.arg(method)
        switch(method,
               "fGarch" = {
                   ## specify the GARCH(1,1) process
		   spec <- fGarch::garchSpec(model =
					     list(omega = a0, alpha = a1, beta = b1))
		   ## simulate; extended=TRUE => data.frame
		   ##		with columns garch, sigma, eps  time series
                   fGarch::garchSim(spec, n=n, extended=TRUE)
               },
               "Gneiting" = { # Gneiting's code (with minor adjustments)
                   alpha <- a1
                   beta <- b1
                   gamma <- a0
                   z <- rep(0, n+1)
                   y <- rep(0, n+1)
                   sigmasq <- rep(1, n+1)
                   z[1] <- rnorm(1)
                   for (t in 2:(n+1)) {
                       sigmasq[t] <- alpha*z[t-1]^2 + beta*sigmasq[t-1] + gamma
                       z[t] <- rnorm(1, sd=sqrt(sigmasq[t]))
                   }
                   ts(cbind(garch=z[-1], sigma=sqrt(sigmasq[-1])))
               },
               stop("wrong 'method'"))
    }

    ## prediction function
    predictTG <- function(x, n, forecaster, pred.opt, pred.pes) {
        stopifnot(is.character(forecaster))
        switch(forecaster,
               "optimist"  = rep(pred.opt, n),
               "pessimist" = rep(pred.pes, n),
               ## statistician predicts hat{X}_t = E[Y_t | sigma_t^2] = sigma_t^2
               "statistician" = x[,"sigma"]^2,
               stop("forecaster not supported"))
    }

    ##' Statistic (use it like this to include simulation in time measurement)

    ##' @title Function to Compute the Simulation Results for One Grid Line
    ##' @param x one-line data.frame containing a combination of grid variables
    ##' @param nonGrids values of non-"grid"-variables
    ##' @return return value of doCallWE()
    ##' @author Marius Hofert
    function(n, alpha0, alpha1, beta1, forecaster,
             pred.optimist, pred.pessimist, method, scoringfun)
    {
        ## simulate squared GARCH(1,1)
        X <- simGARCH11(n, a0=alpha0, a1=alpha1, b1=beta1, method=method)
        Y <- X[,"garch"]^2              # see (2) in Gneiting (2011)
        ## predict
        pred <- predictTG(X, n=n, forecaster=forecaster,
                          pred.opt=pred.optimist, pred.pes=pred.pessimist)

        ## basic check
        stopifnot(length(pred)==length(Y))

        ## compute scoring functions (for all 'inner' variables simultaneously)
        c(SE = scoringfun[["SE"]](pred, y=Y),
          AE = scoringfun[["AE"]](pred, y=Y))
    }
})

## build grid and non-grid variables
stopifnot(dim(pGrid <- mkGrid(vList)) == c(3,1), # only the forecaster here
	  get.nonGrids(vList)$n.sim == 1) # the "automatic n.sim" if there is none

tex.vl <- toLatex(vList)

### Computation ################################################################

S.T <- system.time

## 'reproducing' Table 4 in Gneiting (2011)

## sequentially
S.T(res <- doLapply(vList, doOne=doOne, monitor=TRUE))

## our own monitoring function:
myMoni <- function(i.sim, j, pGrid, res4, n.sim) {
    cat(SysI(),": ",formG(pGrid, j), "; SE=",
        format(res4$value[["SE"]], digits=9, width=15), "\n", sep = "")
}
# Trick, so we can use the printInfo utility functions:
environment(myMoni) <- environment(printInfo[["default"]])


## in parallel

## due to 'R CMD check --as-cran' allowing only <= 2 cores
## note: if doExtras, the check with '--as-cran' fails
(doExtras <- simsalapar:::doExtras())
(nc <- simsalapar:::nCores4test())
(nc.win <- if(.Platform$OS.type=="windows") 1 else nc) # otherwise win-builder fails

## if simsalapar no longer *depends* on parallel:
makeCluster <- parallel::makeCluster

S.T(resM.<- doMclapply    (vList, cores  =nc.win, doOne=doOne, monitor=TRUE))
S.T(resM <- doMclapply    (vList, cores  =nc.win, doOne=doOne, monitor=myMoni))
S.T(resF <- doForeach     (vList, cluster=makeCluster(nc, type="PSOCK"),
                           doOne=doOne, monitor=printInfo[["fileEach"]]))
## For this problem, these are much slower on a typical 4-core desktop:
S.T(resC <- doClusterApply(vList, cluster=makeCluster(nc, type="PSOCK"),
                           doOne=doOne, monitor=printInfo[["gfile"]]))
if(hasRmpi)
S.T(resR <- doRmpi        (vList, nslaves=nc, doOne=doOne, monitor=TRUE))
if(doExtras) { # use all available cores
  print(S.T(resC.<- doClusterApply(vList, cluster=makeCluster(nc, type="PSOCK"),
                                   doOne=doOne, monitor=myMoni)))
  if(hasRmpi)
  print(S.T(resR.<- doRmpi	  (vList, nslaves=nc, doOne=doOne, monitor=myMoni)))
}

## check that we have the same results
stopifnot(doRes.equal(resM, res),
	  doRes.equal(resM, resM.),
	  doRes.equal(resC, res),
          if(hasRmpi) doRes.equal(resR, res) else TRUE,
	  if(doExtras)
	  doRes.equal(resC, resC.) &&
	  (if(hasRmpi) doRes.equal(resR, resR.) else TRUE) else TRUE,
	  doRes.equal(resF, res))


### Analysis ###################################################################

## extract results
val  <- getArray(res) # array of values
err  <- getArray(res, "error") # array of error indicators
warn <- getArray(res, "warning") # array of warning indicators
time <- getArray(res, "time") # array of user times in ms

## warnings, errors
if(any(err > 0))
ftable(100* err, col.vars="forecaster") # percentage of errors
if(any(warn > 0))
ftable(100*warn, col.vars="forecaster") # percentage of warnings

## 'reproduce' Table 4 of Gneiting (2011)
tab1 <- apply(val, c("forecaster", "scoringfun"), mean)
print(apply(tab1, 2, format, digits=3, scipen=-1), quote=FALSE)

names(attributes(tab1)[["dim"]]) <- NULL # so that check below works
## not ok names(dim(tab1)) <- NULL # so that check below works

## *without* mkAL (as check):
res0 <- doLapply(vList, doAL=FALSE, doOne=doOne, monitor=FALSE)
tab0 <- t(sapply(res0, `[[`, "value")) # values for the 3 forecasters: statistician, optimist, pessimist (see monitoring output)
dimnames(tab0)[[1]] <- c("statistician", "optimist", "pessimist")
names(dimnames(tab0)) <- c("forecaster", "scoringfun")
stopifnot(identical(tab0, tab1)) # check


### Testing other demo()s  etc

if(doExtras) {
    .checking <- TRUE # to run a smaller demo:
    demo(robcovMCD)
}
