


sigma.hat.lm <- function(object,...){
    sigma <- summary(object)$sigma
    return (sigma)
}

sigma.hat.glm <- function(object,...){
   dispersion <- if (is.null(object$dispersion)){
                    summary(object)$dispersion
                  }
                  else{
                    object$dispersion
                  }
    if (object$family$family == "gaussian") {
      sigma <- sqrt(dispersion)
    }
    else {
      sigma <- summary(object, correlation = TRUE)$sigma
      #sigma <- sqrt(deviance(object)/df.residual(object))
    }
    return(sigma)
}


sigma.hat.sim <- function(object,...){
    sigma <- object@sigma
    return (sigma)
}

sigma.hat.merMod <- function(object,...){
    #object <- summary (object)
    fcoef <- fixef(object)
    #useScale <- attr (VarCorr (object), "sc")  # =sc?
    #useScale <- object@dims["useSc"]
    useScale <- getME(object, "devcomp")$dims["useSc"]
    #ngrps <- lapply(object@flist, function(x) length(levels(x)))
    #n.groupings <- length (ngrps)
    varc <- VarCorr (object)
    sc <- attr(varc, "sc")  # =useScale
    recorr <- lapply(varc, function(el) attr(el, "correlation"))
    reStdDev <- c(lapply(varc, function(el) attr(el, "stddev")), list(Residual = sc))
    n.groupings <- length(recorr)
    sigmas <- as.list (rep (NA, n.groupings+1))
    sigmas[1] <- ifelse (useScale, sc, 1) #####if NA, sd=1
    cors <- as.list (rep (NA, n.groupings+1))
    names (sigmas) <- names (cors) <- c ("data", names (varc))
    for (k in 1:n.groupings){
      sigmas[[k+1]] <- reStdDev[[k]]
      cors[[k+1]] <- as.matrix (recorr[[k]])
      if (length (cors[[k+1]]) == 1) cors[[k+1]] <- NA
    }
    return (list (sigma=sigmas, cors=cors))
}


sigma.hat.sim.merMod <- function(object,...)
{
    sigma <- object@sigma
    return (sigma)
}
