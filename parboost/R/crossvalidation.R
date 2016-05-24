##' Internal cross-validation function for determing mstop
##' in each mboost model.
##'
##' Takes an mboost object, a matrix of folds and the grid
##' of values to be used during cross-validation. Returns
##' the optimal mstop value. Uses mclapply for
##' parallelization. The implementation is based
##' on \code{\link{mboost}{cvrisk.mboost}}
##' @title Cross-validation for mboost models
##' @param object mboost object
##' @param folds A matrix of folds
##' @param grid The grid of mstop values used during cross-validation
##' @return The optimal mstop value
##' @author Ronert Obst
##' @keywords internal
#' @references Torsten Hothorn, Peter Buehlmann, Thomas Kneib,
#' Matthias Schmid and Benjamin Hofner (2012).
#' Model-Based Boosting. R package version 2.1-2.
cv_subsample <- function(object, folds, grid = 1:mstop(object), cores_cv = detectCores()) {
    weights <- model.weights(object)
    if (any(weights == 0)) warning("zero weights")

    fitfct <- object$update
    oobrisk <- matrix(0, nrow = ncol(folds), ncol = length(grid))

    fam_name <- object$family@name
    call <- deparse(object$call)

    dummyfct <- function(weights, oobweights) {
        mod <- fitfct(weights = weights, oobweights = oobweights)
        mod[max(grid)]
        mod$risk()[grid]
    }

    OOBweights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
    OOBweights[folds > 0] <- 0 # inverse of folds
    if (Sys.info()[1] != "Windows") {
        oobrisk <- mclapply(1:ncol(folds),
                            function(i) dummyfct(weights = folds[, i],
                                                 oobweights = OOBweights[,
                                                     i]), mc.cores = cores_cv)
    } else {
            oobrisk <- lapply(1:ncol(folds),
                            function(i) dummyfct(weights = folds[, i],
                                                 oobweights = OOBweights[, i]))
    }

    ## get errors if mclapply is used
    if (any(idx <- sapply(oobrisk, is.character)))
        stop(sapply(oobrisk[idx], function(x) x))

    oobrisk <- t(as.data.frame(oobrisk))
    oobrisk <- oobrisk / colSums(OOBweights)
    colnames(oobrisk) <- grid
    rownames(oobrisk) <- 1:nrow(oobrisk)

    as.numeric(names(colSums(oobrisk))[which.min(colSums(oobrisk))])
}
