svyrc <- function(formula, design, nd=1, symmetric=FALSE, diagonal=FALSE,
                  weighting=c("marginal", "uniform", "none"), rowsup=NULL, colsup=NULL,
                  Ntotal=nrow(design), exclude=c(NA, NaN),
                  se=c("none", "replicate"), ncpus=getOption("boot.ncpus"),
                  family=quasipoisson, weights=NULL, start=NULL, etastart=NULL,
                  tolerance=1e-8, iterMax=5000, trace=FALSE, verbose=TRUE, ...) {
  model <- svyassocmod("rc", formula=formula, design=design, se=se,
                       nd=nd, symmetric=symmetric, diagonal=diagonal,
                       weighting=weighting, rowsup=rowsup, colsup=colsup, Ntotal=Ntotal, exclude=exclude,
                       ncpus=ncpus, family=family, weights=weights, start=start, etastart=etastart,
                       tolerance=tolerance, iterMax=iterMax, trace=trace, verbose=verbose, ...)

  # We want the actual call the user made: without this, we have the line above,
  # which fails when running replicates because arguments are not in scope
  model$call <- match.call()

  model
}

svyhmskew <- function(formula, design, nd.symm=NA, diagonal=FALSE,
                      weighting=c("marginal", "uniform", "none"),  rowsup=NULL, colsup=NULL,
                      Ntotal=nrow(design), exclude=c(NA, NaN),
                      se=c("none", "replicate"), ncpus=getOption("boot.ncpus"),
                      family=quasipoisson, weights=NULL, start=NULL, etastart=NULL,
                      tolerance=1e-8, iterMax=5000, trace=FALSE, verbose=TRUE, ...) {
  model <- svyassocmod("hmskew", formula=formula, design=design, se=se,
                       nd.symm=nd.symm, diagonal=diagonal,
                       weighting=weighting, rowsup=rowsup, colsup=colsup, Ntotal=Ntotal, exclude=exclude,
                       ncpus=ncpus, family=family, weights=weights, start=start, etastart=etastart,
                       tolerance=tolerance, iterMax=iterMax, trace=trace, verbose=verbose, ...)

  # We want the actual call the user made: without this, we have the line above,
  # which fails when running replicates because arguments are not in scope
  model$call <- match.call()

  model
}

svyyrcskew <- function(formula, design, nd.symm=NA, nd.skew=1, diagonal=FALSE,
                       weighting=c("marginal", "uniform", "none"), rowsup=NULL, colsup=NULL,
                       Ntotal=nrow(design), exclude=c(NA, NaN),
                       se=c("none", "replicate"), ncpus=getOption("boot.ncpus"),
                       family=quasipoisson, weights=NULL, start=NA, etastart=NULL,
                       tolerance=1e-8, iterMax=15000, trace=FALSE, verbose=TRUE, ...) {
  model <- svyassocmod("yrcskew", formula=formula, design=design, se=se,
                       nd.symm=nd.symm, nd.skew=nd.skew, diagonal=diagonal,
                       weighting=weighting, rowsup=rowsup, colsup=colsup, Ntotal=Ntotal, exclude=exclude,
                       ncpus=ncpus, family=family, weights=weights, start=start, etastart=etastart,
                       tolerance=tolerance, iterMax=iterMax, trace=trace, verbose=verbose, ...)

  # We want the actual call the user made: without this, we have the line above,
  # which fails when running replicates because arguments are not in scope
  model$call <- match.call()

  model
}


svyrcL <- function(formula, design, nd=1,
                   layer.effect=c("homogeneous.scores", "heterogeneous", "none"),
                   symmetric=FALSE, diagonal=c("none", "heterogeneous", "homogeneous"),
                   weighting=c("marginal", "uniform", "none"),
                   Ntotal=nrow(design), exclude=c(NA, NaN),
                   se=c("none", "replicate"), ncpus=getOption("boot.ncpus"),
                   family=quasipoisson, weights=NULL, start=NULL, etastart=NULL,
                   tolerance=1e-8, iterMax=5000, trace=FALSE, verbose=TRUE, ...) {
  model <- svyassocmod("rcL", formula=formula, design=design, se=se,
                       nd=nd, layer.effect=layer.effect, symmetric=symmetric, diagonal=diagonal,
                       weighting=weighting, Ntotal=Ntotal, exclude=exclude,
                       ncpus=ncpus, family=family, weights=weights, start=start, etastart=etastart,
                       tolerance=tolerance, iterMax=iterMax, trace=trace, verbose=verbose, ...)

  # We want the actual call the user made: without this, we have the line above,
  # which fails when running replicates because arguments are not in scope
  model$call <- match.call()

  model
}

svyrcL.trans <- function(formula, design, nd=1,
                         symmetric=FALSE, diagonal=c("none", "heterogeneous", "homogeneous"),
                         weighting=c("marginal", "uniform", "none"),
                         Ntotal=nrow(design), exclude=c(NA, NaN),
                         se=c("none", "replicate"), ncpus=getOption("boot.ncpus"),
                         family=quasipoisson, weights=NULL, start=NULL, etastart=NULL,
                         tolerance=1e-8, iterMax=5000, trace=FALSE, verbose=TRUE, ...) {
  model <- svyassocmod("rcL.trans", formula=formula, design=design, se=se,
                       nd=nd, diagonal=diagonal,
                       weighting=weighting, Ntotal=Ntotal, exclude=exclude,
                       ncpus=ncpus, family=family, weights=weights, start=start, etastart=etastart,
                       tolerance=tolerance, iterMax=iterMax, trace=trace, verbose=verbose, ...)

  # We want the actual call the user made: without this, we have the line above,
  # which fails when running replicates because arguments are not in scope
  model$call <- match.call()

  model
}

svyhmskewL <- function(formula, design, nd.symm=NA,
                       layer.effect.skew=c("homogeneous.scores", "heterogeneous", "none"),
                       layer.effect.symm=c("heterogeneous", "uniform", "homogeneous.scores", "none"),
                       diagonal=c("none", "heterogeneous", "homogeneous"),
                       weighting=c("marginal", "uniform", "none"),
                       Ntotal=nrow(design), exclude=c(NA, NaN),
                       se=c("none", "replicate"), ncpus=getOption("boot.ncpus"),
                       family=quasipoisson, weights=NULL, start=NULL, etastart=NULL,
                       tolerance=1e-8, iterMax=5000, trace=FALSE, verbose=TRUE, ...) {
  model <- svyassocmod("hmskewL", formula=formula, design=design, se=se,
                       nd.symm=nd.symm, layer.effect.symm=layer.effect.symm,
                       layer.effect.skew=layer.effect.skew, diagonal=diagonal,
                       weighting=weighting, Ntotal=Ntotal, exclude=exclude,
                       ncpus=ncpus, family=family, weights=weights, start=start, etastart=etastart,
                       tolerance=tolerance, iterMax=iterMax, trace=trace, verbose=verbose, ...)

  # We want the actual call the user made: without this, we have the line above,
  # which fails when running replicates because arguments are not in scope
  model$call <- match.call()

  model
}

svyassocmod <- function(model.function, formula, design,
                        Ntotal=nrow(design), exclude=c(NA, NaN),
                        se=c("none", "replicate"), ncpus=getOption("boot.ncpus"), ...) {
  se <- match.arg(se)

  if(!inherits(design, c("survey.design", "svyrep.design")))
      stop("'design' must be an object of class \"survey.design\" or \"svyrep.design\"")

  if(se == "replicate" && !inherits(design, "svyrep.design"))
      stop("only objects of class \"svyrep.design\" are supported with se=\"replicate\"; see ?as.svrepdesign")

  tab <- survey::svytable(formula, design, Ntotal=Ntotal, exclude=exclude)

  model <- do.call(model.function, list(tab=tab, se="none", ...))

  if(is.null(model))
      return(NULL)

  if(se == "replicate") {
      jb <- jackboot(se, ncpus, NA, tab, model,
                     assoc1=getS3method("assoc", class(model)), assoc2=NULL,
                     model.function=model.function, formula=formula, design=design,
                     Ntotal=Ntotal, exclude=exclude, ...)
      model$assoc$covmat <- jb$covmat
      model$assoc$adj.covmats <- jb$adj.covmats
      model$assoc$boot.results <- jb$boot.results
      model$assoc$jack.results <- jb$jack.results
      model$assoc$svyrep.results <- jb$svyrep.results
  }
  else {
      model$assoc$covmat <- numeric(0)
      model$assoc$adj.covmats <- numeric(0)
      model$assoc$boot.results <- numeric(0)
      model$assoc$jack.results <- numeric(0)
      model$assoc$svyrep.results <- numeric(0)
  }

  model$assoc$covtype <- se

  model
}

svyrep <- function(formula, design, theta, ...,
                   Ntotal=NULL, exclude=c(NA, NaN), cl=NULL, scale.weights=FALSE,
                   load.balancing=FALSE)
{
    # Adapted from withReplicates()
    wts <- design$repweights
    scale <- design$scale
    rscales <- design$rscales

    if(scale.weights)
        pwts <- design$pweights/sum(design$pweights)
    else
        pwts <- design$pweights

    if (inherits(wts, "repweights_compressed")) {
        if(scale.weights)
            wts$weights <- sweep(wts$weights, 2, drop(colSums(wts$weights)), "/")
    }
    else {
        if(scale.weights)
            wts <- sweep(wts, 2, drop(colSums(wts)), "/")
    }

    rpwts <- if(design$combined.weights) 1
             else pwts

    wts <- as.matrix(wts)

    tabl <- lapply(1:ncol(wts), function(i) {
        # Trick svytable() into using current replicate weights instead of sampling weights
        design$pweights <- rpwts * wts[, i]
        survey::svytable(formula, design, Ntotal=Ntotal, exclude=exclude)
    })

    if(!is.null(cl) && requireNamespace("parallel")) {
        if(load.balancing)
            u <- parallel::parLapplyLB(cl, 1:length(tabl), function(i, ...) theta(tab=tabl[[i]], ...), ...)
        else
            u <- parallel::parLapply(cl, 1:length(tabl), function(i, ...) theta(tab=tabl[[i]], ...), ...)
    }
    else {
        u <- lapply(1:length(tabl), function(i) theta(tab=tabl[[i]], ...))
    }

    do.call(rbind, u)
}
