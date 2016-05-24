## #################################################################
##
## Calculates K %*% beta and related quantities for various
## type of model objects.
##
## #################################################################

linest <- function(object, K=NULL, level=0.95, ...){
    UseMethod("linest")
}

linest.lm <- function(object, K=NULL, level=0.95, ...){
    bhat <- coef(object)
    if (is.null(K))
        K <- .defineK( bhat )
    if (!is.matrix(K))
        K <- matrix(K, nrow=1)

    is.est <- is_estimable(K, null_basis( object ))

    VV0  <- vcov(object)
    ddf  <- object$df.residual
    ddf.vec <- rep(ddf, nrow( K ))
    res    <- .getKb( K, bhat, VV0, ddf.vec, is.est)

    p.value <- 2*pt(abs(res[,"t.stat"]), df=res[,"df"], lower.tail=FALSE)

    qq  <- qt(1-(1-level)/2, df=res[,"df"])
    lwr <- res[,"estimate"] - qq * res[,"se"]
    upr <- res[,"estimate"] + qq * res[,"se"]
    #res <- cbind(res, p.value, lwr, upr)
    res <- cbind( res, p.value )

    .finalize_linest(res, K)
}

linest.glm <- function(object, K=NULL, level=0.95, type=c("link","response"), ...){
    type <- match.arg(type)
    bhat <- coef(object)
    if (is.null(K))
        K <- .defineK( bhat )
    if (!is.matrix(K))
        K <- matrix(K, nrow=1)

    is.est <- is_estimable(K, null_basis( object ))

    VV0  <- vcov(object)
    ddf.vec <- rep(object$df.residual, nrow( K ))
    res     <- .getKb( K, bhat, VV0, ddf.vec, is.est)

    if(family(object)[1] %in% c("poisson","binomial","quasipoisson","quasibinomial")){
        p.value <- 2*pnorm(abs(res[,"t.stat"]), lower.tail=FALSE)
        qq <- qnorm(1-(1-level)/2)
        colnames(res)[4] <- "z.stat"
        res <- res[,-3] # NO df's
    } else {
        p.value <- 2*pt(abs(res[,"t.stat"]), df=res[,"df"], lower.tail=FALSE)
        qq <- qt(1-(1-level)/2, df=res[,"df"])
    }

    lwr <- res[,"estimate"] - qq * res[,"se"]
    upr <- res[,"estimate"] + qq * res[,"se"]

    if (type=="response"){
        fit    <- family(object)$linkinv(res[,"estimate"])
        se.fit <- res[,"se"] * abs(family(object)$mu.eta(res[,"estimate"]))
        res[,"estimate"]  <- fit
        res[,"se"] <- se.fit
        lwr <- family(object)$linkinv(lwr)
        upr <- family(object)$linkinv(upr)
    }

    ##res <- cbind(res, p.value, lwr, upr)
    res <- cbind( res, p.value )
    .finalize_linest(res, K)

}


linest.geeglm <- function(object, K=NULL, level=0.95, type=c("link","response"), ...){
    type <- match.arg(type)
    bhat <- coef(object)
    if (is.null(K))
        K <- .defineK( bhat )
    if (!is.matrix(K))
        K <- matrix(K, nrow=1)

    is.est <- is_estimable(K, null_basis( object ))

    VV0  <- summary(object)$cov.scaled
    ddf.vec <- rep(1, nrow(K))
    res     <- .getKb( K, bhat, VV0, ddf.vec, is.est)

    p.value <- 2*pnorm(abs(res[,"t.stat"]), lower.tail=FALSE)
    qq <- qnorm(1-(1-level)/2)
    colnames(res)[4] <- "z.stat"
    res <- res[,-3, drop=FALSE] # NO df's
    ## print("00000000000000000")
    lwr <- res[,"estimate"] - qq * res[,"se"]
    upr <- res[,"estimate"] + qq * res[,"se"]


    if (type=="response"){
        fit    <- family(object)$linkinv(res[,"estimate"])
        se.fit <- res[,"se"] * abs(family(object)$mu.eta(res[,"estimate"]))
        res[,"estimate"]  <- fit
        res[,"se"] <- se.fit
        lwr <- family(object)$linkinv(lwr)
        upr <- family(object)$linkinv(upr)
    }

    ## res <- cbind(res, p.value, lwr, upr)
    res <- cbind( res, p.value )
    .finalize_linest(res, K)
}

linest.lmerMod <- function(object, K=NULL, level=0.95, adjust.df=TRUE, ...){

    bhat <- lme4::fixef(object)

    if (is.null(K))
        K <- .defineK( bhat )
    if (!is.matrix(K))
        K <- matrix(K, nrow=1)

    is.est <- is_estimable(K, null_basis( object ))

    if (adjust.df){
        if (requireNamespace("pbkrtest", quietly=TRUE)){
            VVu  <- vcov(object)
            VV   <- pbkrtest::vcovAdj(object)
            ddf.vec <- unlist(lapply(1:nrow(K),
                                     function(ii) pbkrtest::ddf_Lb(VV , K[ii,], VVu)))
        } else {
            stop("adjustment of degrees of freedom requires that 'pbkrtest' is installed")
        }
    } else {
        a   <- logLik(object)
        ddf <- attributes(a)$nall - attributes(a)$df
        ddf.vec <- rep(ddf, length(bhat))
        VV  <- vcov(object)
    }


    res     <- .getKb( K, bhat, VV, ddf.vec, is.est)
    p.value <- 2*pt(abs(res[,"t.stat"]), df=res[,"df"], lower.tail=FALSE)
    qq  <- qt(1-(1-level)/2, df=res[,"df"])
    lwr <- res[,"estimate"] - qq * res[,"se"]
    upr <- res[,"estimate"] + qq * res[,"se"]
    ##res <- cbind(res, p.value, lwr, upr)
    res <- cbind( res, p.value )
    .finalize_linest(res, K)

}

linest.merMod <- function(object, K=NULL, level=0.95, ...){
    cl <- match.call()
    cl[[1]] <- as.name("linest.lmerMod")
    cl$adjust.df <- FALSE
    eval(cl)
}



### UTILITIES ###

.getKb <- function(K, bhat, VV, ddf.vec, is.est, level=0.95){
    #' cat(".getKb")
    #' print(attributes(K))
    off <- attr(K, "offset")
    #' print(off)

    used       <- which(!is.na(bhat))
    bhat.used  <- bhat[used]
    K   <- K[, used, drop=FALSE]
    res <- matrix(NA, nrow=nrow(K), ncol=3)
    for (ii in 1:nrow(res)){
        kk <- K[ii,]
        if (is.est[ii]){
            est  <- sum(kk * bhat.used)
            se   <- sqrt(sum(kk * (VV %*% kk)))
            df2  <- ddf.vec[ii]
            res[ii,] <- c(est, se, df2)
        }
    }

    if (!is.null(off))
        res[,1] <- res[,1] + off[[1]]

    colnames(res) <- c("estimate","se","df")
    t.stat        <- res[,"estimate"] / res[,"se"]
    cbind(res, t.stat)
}

.defineK <- function( bhat ){
    K <- diag( 1, length( bhat ) )
    rownames( K ) <- names( bhat )
    ## print(K)
    K
}

.finalize_linest <- function(.coef, K){
    ## print(rownames(K))
    if (!is.null(rownames(K)))
        rownames(.coef) <- rownames(K)
    res  <- list(coef=.coef, grid=attr(K,"grid"), K=K)
    class(res) <- "linearEstimate"
    res
}

print.linearEstimate <- function(x, ...){
    print(cbind(x$coef, x$grid))
    invisible(x)
}

