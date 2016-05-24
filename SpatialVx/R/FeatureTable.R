FeatureTable <- function(x, fudge = 1e-8, hits.random = NULL, correct.negatives = NULL, fA = 0.05) {

    if(class(x) != "matched") stop("FeatureTable: invalid type for argument x.")

    a <- attributes(x)
    a$names <- NULL
    out <- list()
    attributes(out) <- a

    xdim <- dim(x$X.labeled)

    m <- x$matches

    hits <- dim(m)[1]
    miss <- length(x$unmatched$X)
    fa   <- length(x$unmatched$Xhat)

    No <- max(x$X.labeled, na.rm = TRUE)
    Nf <- max(x$Y.labeled, na.rm = TRUE)

    if(is.null(correct.negatives)) {

            bigD <- (1 - fA) * No / fA - Nf

	    if(bigD < 0) {

		warning("FeatureTable: attempted to estimate a value for correct.negatives, but got a negative answer.  Setting to NA.")
		bigD <- NA

	    } else if(bigD == 0) bigD <- fudge

    } else bigD <- correct.negatives

    if(is.null(hits.random)) {

	hits.random <- (Nf * No) / (Nf + miss + bigD)

    } # end of if must calculate 'hits.random' stmt.

    denom <- Nf + miss - hits.random
    if(denom == 0) denom <- fudge

    GSS <- (hits - hits.random) / denom

    s <- (hits + miss) / bigD

    POD <- hits / (hits + miss + fudge)
    SH2 <- POD * (1 - POD) / (hits + miss + fudge)
    POD.se <- sqrt(SH2)

    FArate <- fa / (fa + bigD)
    SF2 <- FArate * (1 - FArate) / (fa + bigD)
    FArate.se <- sqrt(SF2)

    FAR <- fa / (hits + fa + fudge)
    FAR.se <- sqrt( (FAR^4) * ( (1 - POD) / (hits + fudge) + (1 - POD) / (fa + fudge) ) * ( hits^2 / (fa^2 + fudge) ) )

    HSS <- 2 * (hits * bigD - fa * miss) / ( (hits + miss ) * (miss + bigD) + (hits + fa) * (fa + bigD) + fudge )
    SHSS2 <- SF2 * (HSS^2) * ( 1/(POD - FArate + fudge) + (1 - s) * (1 - 2 * s))^2 + SH2 * (HSS^2) * (1/ (POD - FArate + fudge) - s * (1 - 2 * s))^2
    HSS.se <- sqrt(SHSS2)

    if(2 - HSS == 0) GSS.se <- sqrt(4 * SHSS2/fudge)
    else GSS.se <- sqrt(4 * SHSS2 / ((2 - HSS)^4))

    nomen <- c("GSS", "POD", "false alarm rate", "FAR", "HSS")
    theta <- c(GSS, POD, FArate, FAR, HSS)
    theta.se <- c(GSS.se, POD.se, FArate.se, FAR.se, HSS.se)

    names(theta) <- nomen
    names(theta.se) <- nomen

    out$estimates <- theta
    out$se <- theta.se

    tab <- c(hits, miss, fa, bigD)
    names(tab) <- c("hits", "misses", "fales alarms", "correct negatives")

    out$feature.contingency.table <- tab

    class(out) <- "FeatureTable"
    return(out)

} # end of 'FeatureTable' function.


ci.FeatureTable <- function(x, alpha = 0.05, ...) {

    a <- attributes(x)
    if(!is.null(a$data.name)) dname <- a$data.name
    else dname <- deparse(substitute(x))

    z.alpha <- qnorm(alpha / 2, lower.tail = FALSE)

    val <- x$estimates
    se <- x$se

    res <- cbind( val - z.alpha * se, val, val + z.alpha * se)

    conf.level <- paste(round((1 - alpha) * 100, digits = 2), "%", sep = "")

    colnames(res) <- c(paste(conf.level, " lower CI", sep = ""), "Estimate", paste(conf.level, " upper CI", sep = ""))

    attr(res, "data.name") <- dname
    attr(res, "method") <- "Normal Approximation"
    attr(res, "conf.level") <- (1 - alpha) * 100

    class(res) <- "ci"
    return(res)

} # end of 'ci.FeatureTable' function.

print.FeatureTable <- function(x, ...) {

    print(summary(object = x, ...))
    invisible()

} # end of 'print.FeatureTable' function.

summary.FeatureTable <- function(object, ...) {

    cat("\n\nFeature-based Contingency Table\n")
    print(object$feature.contingency.table)

    out <- ci(object, ...)
    print(out)

    invisible(out)

} # end of 'summary.FeatureTable' function.


