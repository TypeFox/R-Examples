# S3 methods

print.maxlikeFit <- function(x, ...) {
    converge <- x$optim$converge
    cat("\nCall:", paste(deparse(x$call)), "\n\n")
    cat("Coefficients:\n")
    print(x$Est, ...)
    cat("\nAIC:", x$AIC, "\n\n")
    if(converge != 0)
        warning("Model did not converge")
    }



coef.maxlikeFit <- function(object, ...) object$Est[,"Est"]


vcov.maxlikeFit <- function(object, ...) object$vcov





summary.maxlikeFit <- function(object, digits=3, ...) {
    s <- object$Est
    z <- s[,"Est"] / s[,"SE"]
    p <- 2*pnorm(abs(z), lower.tail = FALSE)
    s <- cbind(s, z=z, "P(>|z|)"=p)
    converge <- object$optim$converge
    cat("\nCall:", paste(deparse(object$call)), "\n", fill=TRUE)
    cat("Coefficients:\n")
    print(s, digits=digits, ...)
    cat("\noptim convergence code:", converge, "\n")
    cat("\nAIC:", object$AIC, "\n\n")
    if(converge != 0)
        warning("Model did not converge")
    invisible(s)
    }




logLik.maxlikeFit <- function(object, ...) {
    -object$optim$value
}


AIC.maxlikeFit <- function(object, ..., k=2) {
    2*object$optim$value + k*length(coef(object)[object$not.fixed])
}




predict.maxlikeFit <- function(object, ...) {
    e <- coef(object)
    rasters <- object$rasters
    if(is.null(rasters)) {
        rasters <- try(get(as.character(object$call$rasters)))
        if(identical(class(rasters)[1],  "try-error"))
            stop("could not find the raster data")
        warning("raster data were not saved with object, using the data found in the workspace instead.")
    }
    link <- object$link
    cd.names <- names(rasters)
    npix <- prod(dim(rasters)[1:2])
    z <- as.data.frame(matrix(getValues(rasters), npix))
    names(z) <- cd.names
    formula <- object$call$formula
    varnames <- all.vars(formula)
    if(!all(varnames %in% cd.names))
        stop("at least 1 covariate in the formula is not in rasters.")
    Z.mf <- model.frame(formula, z, na.action="na.pass")
    Z.terms <- attr(Z.mf, "terms")
    Z <- model.matrix(Z.terms, Z.mf)
    eta <- drop(Z %*% coef(object))
    if(identical(link, "logit"))
        psi.hat <- plogis(eta)
    else if(identical(link, "cloglog"))
        psi.hat <- 1-exp(-exp(eta))
    else
        stop("link function should be either 'logit' or 'cloglog'")
    psi.mat <- matrix(psi.hat, dim(rasters)[1], dim(rasters)[2],
                      byrow=TRUE)
    psi.raster <- raster(psi.mat)
    extent(psi.raster) <- extent(rasters)
    psi.raster
}









# Garbage
chisq <- function(object, ...) UseMethod("chisq")

chisq.maxlikeFit <- function(object, fact, ...) {
    stop("This is wrong")
    E <- predict(object)
    if(!missing(fact)) {
        sum2 <- function(x, na.rm=TRUE) {
            if(all(is.na(x)))
                return(NA)
            else
                return(sum(x, na.rm=TRUE))
        }
        E <- aggregate(E, fact=fact, fun=sum2, expand=FALSE)
    }
    xy <- object$points.retained
    # Need area of each new pixel
    Ep <- E / cellStats(E, sum)
    npix <- ncell(Ep)
    cellID <- cellFromXY(Ep, xy)
    cellID <- factor(cellID, levels=1:npix)
    n <- table(cellID)
    p <- values(Ep)
    observed <- expected <- Ep
    values(observed) <- as.numeric(n)
    values(expected) <- as.numeric(nrow(xy)*p)
    observed[is.na(expected)] <- NA
    rast.stack <- stack(observed, expected)
    names(rast.stack) <- c("observed", "expected")
    keep <- !is.na(p)
    nr <- n[keep]
    pr <- p[keep]
    if(any(pr==0)) {
        warning("Zero probs were converted to .Machine$double.xmin")
        pr[pr==0] <- .Machine$double.xmin
}
    out <- list(test=chisq.test(nr, p=pr, rescale.p=TRUE, ...),
                obsexp=rast.stack)
    return(out)
}

