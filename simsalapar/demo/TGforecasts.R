### Squared GARCH(1,1) forecast example from Gneiting (2011),
### Making and Evaluating Point Forecasts,
### Journal of the American Statistical Association 106, 746--762


### Setup ######################################################################

require(simsalapar)
## we also use require(fGarch)


### Variable list ##############################################################

## Create varlist (default type is 'frozen'; default expressions are generated)
## note: no n.sim here
varList <-
    varlist(forecaster = list(type="grid", expr = quote(italic(forecaster)),
                              value = c("statistician", "optimist", "pessimist")),
            scoringfun = list(type="inner", expr = quote(S),
                              value = list(
                              ## squared error
                              SE  = function(x, y) mean((x-y)^2),
                              ## absolute error
                              AE  = function(x, y) mean(abs(x-y)),
                              ## absolute percentage error
                              APE = function(x, y) mean(abs((x-y)/y)),
                              ## relative error
                              RE  = function(x, y) mean(abs((x-y)/x)))),
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
            n = list(value = 1e5),
            ## simulation method
            method = list(value = "fGarch"))

## build grid and non-grid variables
(pGrid <- mkGrid(varList)) # only the forecaster here
get.n.sim(varList) # just 1 here
str(ng <- get.nonGrids(varList))


### doOne() including ingredient functions #####################################
##  -------
##  instead of 'local', you may use a (simpler, but less elegant) construction
##  like doOne <- function(x, nonGrids, simGARCH11, predictTG) or pass auxiliary
##  functions as 'frozen' variables

doOne <- local({

    ## function for simulating a squared GARCH(1,1)
    simGARCH11 <- function(n, a0, a1, b1, method=c("fGarch", "Gneiting"))
    {
        method <- match.arg(method)
        switch(method,
               "fGarch" = {
                   ## specify the GARCH(1,1) process
                   spec <- fGarch::garchSpec(model = list(omega = a0, alpha = a1, beta = b1))
                   ## simulate (extended => data.frame with columns garch, sigma, eps time series)
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

    ##' Statistic

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
        c(SE  = scoringfun[["SE" ]](pred, y=Y),
          AE  = scoringfun[["AE" ]](pred, y=Y),
          APE = scoringfun[["APE"]](pred, y=Y),
          RE  = scoringfun[["RE" ]](pred, y=Y))
    }
})



### Computation ################################################################

## 'reproducing' Table 4 in Gneiting (2011)

## sequentially
system.time(res <- doLapply(varList, #sfile="TGforecasts_res_l.rds",
			    doOne=doOne, monitor=TRUE))

## parallel
require(parallel) # e.g. for  makeCluster() below:
canFork <- (.Platform$OS.type != "windows")
(hasRmpi <- require("Rmpi"))# not on Windows

useSock <- TRUE
system.time(resM <- doMclapply(varList, sfile="TGforecasts_res_M.rds",
			       doOne=doOne, monitor=TRUE))
system.time(resC <- doClusterApply(varList, sfile="TGforecasts_res_C.rds",
				   doOne=doOne, monitor=printInfo[["gfile"]]))
if(useSock)
system.time(resCP<- doClusterApply(varList, cluster=makeCluster(detectCores(),
                                   type="PSOCK"), sfile="TGforecasts_res_CP.rds",
				   doOne=doOne, monitor=printInfo[["gfile"]]))
if(canFork)
system.time(resCF<- doClusterApply(varList, cluster=makeCluster(detectCores(),
                                   type="FORK"), sfile="TGforecasts_res_CF.rds",
				   doOne=doOne, monitor=printInfo[["gfile"]]))
if(hasRmpi)
system.time(resR <- doRmpi(varList, sfile="TGforecasts_res_R.rds",
			   doOne=doOne, monitor=printInfo[["gfile"]]))
system.time(resF <- doForeach(varList, sfile="TGforecasts_res_F.rds",
			      doOne=doOne, monitor=TRUE))

## check that we have the same results
stopifnot(doRes.equal(resM, res),
	  doRes.equal(resC, res),
	  if(useSock) doRes.equal(resCP,res) else TRUE,
	  if(canFork) doRes.equal(resCF,res) else TRUE,
	  doRes.equal(resR, res),
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
res0 <- doLapply(varList, sfile="TGforecasts_res_l_doAL=FALSE.rds", doAL=FALSE,
                 doOne=doOne, monitor=TRUE)
tab0 <- t(sapply(res0, `[[`, "value")) # values for the 3 forecasters: statistician, optimist, pessimist (see monitoring output)
dimnames(tab0)[[1]] <- c("statistician", "optimist", "pessimist")
names(dimnames(tab0)) <- c("forecaster", "scoringfun")
stopifnot(identical(tab0, tab1)) # check

## plot ("separate" the 4 scoring functions)
mayplot(val, varList, row.vars="scoringfun", xvar="forecaster",
        ylim="local", log="y")



### (improved) Figure 1 Gneiting (2011) ########################################

## extract variables (except 'n' -- T = 200 here)
ng. <- ng$nonGrids
alpha0 <- ng.[["alpha0"]]
alpha1 <- ng.[["alpha1"]]
beta1  <- ng.[["beta1"]]
opt.pr <- ng.[["pred.optimist"]]
pes.pr <- ng.[["pred.pessimist"]]
method <- ng.[["method"]]

## reproducing Figure 1 in Gneiting (2011)
stopifnot(require(fGarch))
set.seed(0)
T <- 200
simG11 <- environment(doOne)$simGARCH11
X <- simG11(T, a0=alpha0, a1=alpha1, b1=beta1, method=method)
predTG <- environment(doOne)$predictTG
pred.stat <- predTG(X, n=T, forecaster="statistician")
pred.opt  <- predTG(X, n=T, forecaster="optimist",  pred.opt=opt.pr)
pred.pes  <- predTG(X, n=T, forecaster="pessimist", pred.pes=pes.pr)
col3 <- c("blue", "orange", "red")
plot(cbind(X[,"garch"]^2, stat=pred.stat, optm=pred.opt, pess=pred.pes),
     plot.type="single", col=c("black", col3), lwd=c(1,2,2,2),
           type="l", xlab="Trading day", ylab="Asset price",
     main="Improved Figure 1 from Gneiting (2011)")
legend("topleft", inset=0.04, lty=1, lwd=2, col=col3,
       bty="n", legend=c("Statistician", "Optimist", "Pessimist"))
