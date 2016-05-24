evm <- function(y, data, family=gpd, ...){
  theCall <- match.call()
  if (!missing(data)) {
      y <- ifelse(deparse(substitute(y))== "substitute(y)", deparse(y),deparse(substitute(y)))
      y <- formula(paste(y, "~ 1"))
      y <- model.response(model.frame(y, data=data))
  }
  UseMethod("evm", y)
}

evm.default <-
function (y, data, family=gpd, th= -Inf, qu,
          ..., # arguments specific to family such as phi = ~ 1
          penalty = NULL, prior = "gaussian",
          method = "optimize", cov="observed",
          start = NULL, priorParameters = NULL,
          maxit = 10000, trace=NULL,
          iter = 40500, burn=500, thin = 4,
          proposal.dist = c("gaussian", "cauchy"),
          jump.cov, jump.const=NULL, R=100, verbose=TRUE) {

    modelParameters <- texmexParameters(theCall, family,...)

    ##################### Sort out method, penalty/prior, trace...

    method <- texmexMethod(method)
    prior <- texmexPrior(prior, penalty, method, priorParameters)

    trace <- texmexTrace(trace, method)
    otrace <- trace[1]; trace <- trace[2]

    ############################## Construct data to use...

    if (missing(data)){ data <- NULL }
    else { y <- deparse(substitute(y)) }

    # Get list containing response (y) and design matrix for each parameter
    modelData <- texmexPrepareData(y, data, modelParameters)

    if (missing(th) & !missing(qu)) {
        th <- quantile(modelData$y, qu)
    }

    if (!is.finite(th)){ rate <- 1 }
    else { rate <- mean(modelData$y > th) }

    modelData <- texmexThresholdData(th, modelData)

    ###################### If family does not give info matrix...
    if (is.null(family$info)){ cov <- "numeric" }

    ###################### Check and sort out prior parameters...
    priorParameters <- texmexPriorParameters(prior, priorParameters, modelData)

    ################################## Do the optimization....

    o <- evmFit(data = modelData, family=family, th=th,
                 prior=prior,
                 start=start, hessian = cov == "numeric",
                 priorParameters = priorParameters,
                 maxit = maxit, trace = otrace)

    if (o$convergence != 0 | o$value == 10^6) {
        warning("Non-convergence in evm.default")
    }

    ##### Construct object containing the penalized likelihood estimates

    o <- constructEVM(o, family, th, rate, prior, modelParameters, theCall,
                      modelData, data, priorParameters, cov)

    #### Simulate from posteriors....
    if (method == "s"){
        proposal.dist <- match.arg(proposal.dist)
        o <- evmSim(o, priorParameters=priorParameters,
                    prop.dist=proposal.dist,
                    jump.const=jump.const, jump.cov=jump.cov,
                    iter=iter, start=start, verbose=verbose,
                    thin=thin, burn=burn,
                    trace=trace, theCall)
    } # Close else
    else if (method == "b"){
        o <- evmBoot(o, R=R)
    }

    o
}

