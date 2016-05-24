## This is part of the QCA3 package
## By ronggui HUANG 2010.

directCalibration <- function(x,fullin, fullout, crossover, infz=0.953, outfz=0.047,details=FALSE) {
    ## Ragin, C. 2008. "Measurement versus Calibration: A Set-theoretic Approach." in The Oxford handbook of political methodology, Oxford University Press, : 174-198.
    ## Also, Ragin, Charles C. 2008. "Redesigning social inquiry: fuzzy sets and beyond." Chapter 5. Chicago: University of Chicago Press.
    if (fullin < fullout || crossover < fullout || fullin < crossover)
        stop("It should be that: fullin > crossover > fullout.")
    var.dev <- x - crossover
    inLogOdd <- log(infz/(1-infz))
    outLogOdd <- log(outfz/(1-outfz))
    inScalar <- inLogOdd/(fullin - crossover)
    outScalar <- outLogOdd/(fullout - crossover)
    scalars <-  rep(NA,length(x))
    scalars[var.dev > 0] <- inScalar
    scalars[var.dev < 0] <- outScalar
    product <-  scalars * var.dev
    fz <- exp(product )/(1+exp(product))
    if (details) {
        ans <- data.frame(x=x,deviations=var.dev,scalars=scalars,logOdd=product,membership=fz)
    } else ans <- fz
    ans
}

## indirectCalibration <- function(x, qualcode, ...){
##     ## for the indirect method, may use the mpf package
##     require(mfp)
##     md <- mfp(qualcode~fp(x,df=2),family=binomial(logit))
##     ans <- predict(md,type="response")
##     ans
## }
