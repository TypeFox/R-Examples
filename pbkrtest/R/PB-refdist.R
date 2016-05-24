### ###########################################################
###
### Parallel computing of reference distribution
###
### ###########################################################

PBrefdist <- function(largeModel, smallModel, nsim=1000, seed=NULL, cl=NULL, details=0){
    UseMethod("PBrefdist")
}

PBrefdist.lm <- function(largeModel, smallModel, nsim=1000, seed=NULL, cl=NULL, details=0){

  ##cat(".....PBrefdist.lm\n")
  t0 <- proc.time()
  .is.cluster <- !is.null(cl) && inherits(cl, "cluster")

  if (!.is.cluster){
    ref <- .lm_refDist(largeModel, smallModel, nsim, seed=seed)
  } else {
    nsim2 <- round(nsim/length(cl))
    if (details>=1)
      cat(sprintf("* Using %i clusters and %i samples per cluster\n", length(cl), nsim2))
    clusterExport(cl, ls(envir=.GlobalEnv), envir = .GlobalEnv)
    clusterSetRNGStream(cl)
    ref <- unlist(clusterCall(cl, .lm_refDist, largeModel, smallModel, nsim2))
  }

  ref <- ref[ref>0]
  ctime <- (proc.time()-t0)[3]
  attr(ref,"ctime") <- ctime
  if (details>0)
    cat(sprintf("Reference distribution with %i samples; computing time: %5.2f secs. \n",
                length(ref), ctime))

  ref
}

.lm_refDist <- function(lg, sm, nsim=20, seed=NULL, simdata=simulate(sm, nsim=nsim, seed=seed)){
    ##simdata <- simulate(sm, nsim, seed=seed)
    ee  <- new.env()
    ee$simdata <- simdata

    ff.lg <- update.formula(formula(lg),simdata[,ii]~.)
    ff.sm <- update.formula(formula(sm),simdata[,ii]~.)
    environment(ff.lg) <- environment(ff.sm) <- ee

    cl.lg <- getCall(lg)
    cl.sm <- getCall(sm)

    cl.lg$formula <- ff.lg
    cl.sm$formula <- ff.sm

    ref <- rep.int(NA, nsim)
    for (ii in 1:nsim){
        ref[ii] <- 2*(logLik(eval(cl.lg))-logLik(eval(cl.sm)))
    }
    ref
}

.merMod_refDist <- function(lg, sm, nsim=20, seed=NULL, simdata=simulate(sm, nsim=nsim, seed=seed)){
  #simdata <- simulate(sm, nsim=nsim, seed=seed)
  unname(unlist(lapply(simdata, function(yyy){
    sm2     <- refit(sm, newresp=yyy)
    lg2     <- refit(lg, newresp=yyy)
    2*(logLik( lg2, REML=FALSE ) - logLik( sm2, REML=FALSE ))
  })))
}

PBrefdist.mer <-
PBrefdist.merMod <- function(largeModel, smallModel, nsim=1000, seed=NULL, cl=NULL, details=0){

    t0 <- proc.time()
    if (getME(smallModel, "is_REML"))
        smallModel <- update( smallModel, REML=FALSE )
    if (getME(largeModel, "is_REML"))
        largeModel <- update( largeModel, REML=FALSE )
    
    .is.cluster <- !is.null(cl) && inherits(cl, "cluster")
    if (!.is.cluster){
        ref <- .merMod_refDist(largeModel, smallModel, nsim=nsim, seed=seed)
    } else {
        nsim.cl <- nsim %/% length(cl)
        clusterSetRNGStream(cl)
        ref <- unlist( clusterCall(cl, fun=.merMod_refDist, largeModel, smallModel, nsim=nsim.cl) )
    }
    
    ctime <- (proc.time()-t0)[3]
    attr(ref,"ctime")   <- ctime
    attr(ref,"samples") <- c(nsim=nsim, npos=sum(ref>0))
    if (details>0)
        cat(sprintf("Reference distribution with %5i samples; computing time: %5.2f secs. \n",
                    length(ref), ctime))
    
    ref
}


















