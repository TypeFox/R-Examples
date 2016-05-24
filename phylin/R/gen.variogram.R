gen.variogram <-
function(x, y, lag = mean(x)/sqrt(nrow(x)), tol=lag/2, lmax = NA,
         bootstraps = 999, verbose = FALSE) {

    ## Check x class status and attempt coercion to matrix
    if (!(class(x) == "matrix")) {
        x <- as.matrix(x)
    }

    ## Check if there are multiple genetic distances and
    ## attempt coercion to matrix
    multi <- FALSE    
    if (class(y) == "list") {
        multi <- TRUE
        for (i in 1:length(y)) {
            y[[i]] <- as.matrix(y[[i]])
        }
    } else {
        y <- as.matrix(y)
    }
            
    if (multi) {
        gv <- .gen.variogram.multi(x, y, lag = lag, tol = tol,
                                  lmax = lmax,
                                  bootstraps = bootstraps,
                                  verbose = verbose)
    } else {
        gv <- .gen.variogram.single(x, y, lag = lag, tol = tol,
                                   lmax = lmax)
    }

    ## Throw a warning if lag used generates NAs
    if (sum(is.na(gv$gamma)) > 0) {
        warning("The variogram contains NAs.",
                "Consider adjusting the lag.")
    }    

    gv
}


