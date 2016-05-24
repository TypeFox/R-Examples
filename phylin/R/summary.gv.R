summary.gv <- 
function(object, ...) {
    mtest <- mtest.gv(object)
    nl <- "\n\r"
    cat(nl)
    cat(paste("Summary for", deparse(substitute(object)), nl))
    cat(nl)
    cat(paste("Number of samples: ", ncol(object$x), nl))
    cat(paste("Variogram parameters", nl))
    cat(paste("   lag:             ", round(object$param$lag, 3), nl))
    cat(paste("   lag tolerance:   ", round(object$param$tol, 3), nl))
    cat(paste("   maximum distance:", round(object$param$lmax, 3), nl))
    cat(paste("   lags with data:  ", 
              length(object$lag[!is.na(object$lag)]), nl))
    cat(paste("   semi-variance summary:", nl))
    print(summary(object$gamma))
    cat(nl)

    if(mtest) {
        # calculate R squared for the model
        ssres <- sum((object$residuals)**2)
        sstot <- sum((object$gamma - mean(object$gamma))**2)
        R2 <- 1 - (ssres/sstot)

        models <- c("gaussian", "exponential", "spherical", "linear")
        cat(paste("Semivariogram with ", models[object$model$type], 
                  " model (R squared = ", round(R2, 4), ")", nl, sep=""))
        cat(paste("Model parameters:", nl))
        cat(paste("   sill:  ", round(object$model$sill, 3), nl))
        cat(paste("   range: ", round(object$model$range, 3), nl))
        cat(paste("   nugget:", round(object$model$nugget, 3), nl))
        if (object$model$type == 4) 
            cat(paste("   slope:", round(object$model$slope, 3), nl))
    } else {
        cat(paste("No model found.", nl))
    }
}
