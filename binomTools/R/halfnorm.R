### Halfnorm points and optional envelope returned as a list
halfnormList <- function(object, resType, env, nsim)
{
    res <- Residuals(object, type=resType)

    ## Absolute values of residuals in ascending order plotted against Phi^(-1){(i+n-1/8)/(2n+1/2)}
    no.res <- length(object$resid)
    expected <- qnorm((1:no.res + no.res - 1/8)/(2*no.res + 1/2))
    sorted <- sort(abs(res))

    List <- list(resType=resType, residuals=sorted, expected=expected) 
    if(env) List$EnvList <- halfnormEnv(object, resType, nsim)
    List
}

### Simulated envelope
halfnormEnv <- function(object, resType, nsim) 
{
  
    mf <- model.frame(object)
    mr <- model.response(mf)
    groupSize <- apply(mr, 1, sum)
    if(!is.null(model.weights(mf))) mw <- model.weights(mf)
    else mw <- rep(1, nrow(mf))
    start <- coef(object)
    call <- update(object, formula=.YY ~ ., data = simData, weights =
    mw, evaluate = F, start=start)

    ## Simulation of data
    simData <- object$data
    simData$mw <- mw
    envList <- vector('list', 3)
    for(i in 1:nsim) {
        suc <- rbinom(dim(object$data)[1], groupSize, object$fitted)
        simData$.YY <- cbind(suc, groupSize-suc)
        simObj <- eval(call)
        envRes <- sort(abs(Residuals(simObj, type=resType)))
        names(envRes) <- NULL
        if(is.null(unlist(envList))) envList <- list(minValues=envRes,
            meanValues=envRes/nsim, maxValues=envRes)  
        else {envList$minValues[envRes < envList$minValues] <- envRes[envRes < envList$minValues]
            envList$meanValues <- envList$meanValues + envRes/nsim
              envList$maxValues[envRes > envList$maxValues] <- envRes[envRes > envList$maxValues]
        }
    }

    envList
}

### Plot of points and optional envelope
halfplot <- function(list, env, type=NULL, ...)
{
    xy <- xy.coords(list$expected, list$residuals); x <- xy$x; y <- xy$y; 
    minEnv <- list$EnvList$minValues
    meanEnv <- list$EnvList$meanValues
    maxEnv <- list$EnvList$maxValues
    resType <- list$resType

    resType <- switch(resType, approx.deletion =
                      "approx. deletion", exact.deletion =
                      "exact deletion", standard.deviance =
                      "standardized. deviance", standard.pearson =
                      "standardized Pearson", deviance = "deviance", pearson =
                      "Pearson", working = "working", response =
                      "response", partial = "partial") 

    if(env) ylim <- c(0, max(c(maxEnv, y)))
    else ylim <- c(0, max(y))

    oldPar <- par(mar=c(5, 5.5, 4, 2) + 0.1)
    on.exit(par(oldPar))
    plot(x, y, pch=20, cex.axis=0.8, las=1, tcl=-0.3, mgp=c(3, 0.6, 0),
        ylim=ylim,
        xlab='Expected value of half-normal order statistic', 
        ylab=paste('Absolute value of \n', resType, 'residuals'), 
        main='Half-normal plot of residuals', type=type, ...)
    if(env) {
        lines(x, minEnv, lty=3)
        lines(x, meanEnv)
        lines(x, maxEnv, lty=3)
    }
}

### Function to identify potential outliers
idHalfplot <- function(list, n, env, pch=NULL, tolerance=0.25, ...)
{
    xy <- xy.coords(list$expected, list$residuals); x <- xy$x; y <- xy$y
    if(is.null(pch)) pch <- names(list$residuals)
    ind <- rep(FALSE, length(x)); obsList <- integer(0)

    while(sum(ind) < n) {
        obs <- identify(x[!ind], y[!ind], n=1, plot=FALSE, tolerance=tolerance, ...)
        obs <- which(!ind)[obs]
        obsList <- c(obsList, obs)
        ind[obs] <- TRUE
        halfplot(list, env, type='n')
        points(x[which(!ind)], y[which(!ind)], pch=20)
        text(x[obsList], y[obsList], labels=pch[obsList], cex=0.75)
    }
}

### Final halfnorm function
halfnorm <- function(object, resType=c("approx.deletion",
                     "exact.deletion", "standard.deviance",
                     "standard.pearson", "deviance", "pearson",
                     "working", "response", "partial"), env=T, nsim =
                     20, plot=T, identify=F, n=2)  {  

    resType <- match.arg(resType)
    List <- halfnormList(object, resType, env, nsim)

    if(plot) {
        halfplot(List, env)
        if(identify)
            idHalfplot(List, n, env)}
    else List
}


