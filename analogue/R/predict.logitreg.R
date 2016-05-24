##' Predict the posterior probability of analogue-ness between fossil
##' and training set samples based on logistic regression fits.
##'
##' @title Posterior probability of analogue-ness for fossil samples
##' @param object object of class \code{"logitreg"}
##' @param newdata matrix of dissimilarities between training and
##' fossil samples. Should be an object of class \code{"distance"}.
##' @param group The group to plot the logit model for. Can be one or
##' more of the group labels or \code{"Combined"} to draw the individual
##' logit models. Alternatively, and the default, is to use \code{"all"},
##' which divides the plotting region into the required number of
##' plotting regions draws all the fitted curves.
##' @param k numeric; the number of close modern analogues to consider.
##' Currently not to be used!
##' @param ... additional arguments passed to other methods.
##' @return A matrix of posterior probabilities is returned.
##' @author Gavin Simpson
##' @method predict logitreg
##' @S3method predict logitreg
##' @keywords methods
`predict.logitreg` <- function(object, newdata, group = "all",
                               k = 1, ...) {
    if(missing(newdata)) {
        return(fitted(object)) ## fitted.logitreg defined
    }
    ## want the close modern analogues - so this all comes from]
    ## newdata. Code is based on cma() - must refactor!
    nsamp <- ncol(newdata)
    nams <- colnames(newdata)
    close <- lapply(seq_len(nsamp),
                    function(i, newdata) {
                        sort(newdata[, i])
                    }, newdata = newdata)
    names(close) <- nams

    ## minimum distance to a group
    mindist <- t(sapply(close, .minDijGroup,
                        groups = object$groups, k = k))

    ## posterior probability of analogue-ness per group
    posterior <- .postProbGroup(mindist, object$models)

    ## return object -- just return posterior for now
    ## is a matrix, cols = groups; rows = fossil samples
    posterior
}

## Function to extract the min Dij for each biome type
.minDijGroup <- function(x, groups, k = 1) {
    lev <- levels(groups)
    out <- numeric(length = length(lev))
    names(out) <- lev
    for(i in seq_along(lev)) {
        out[i] <- unname(x[which(groups[names(x)] == lev[i])[k]])
    }
    out
}

## Function to compute the probability of analogue for the minimum Dij
.postProbGroup <- function(x, model) {
    foo <- function(i, x, model) {
        dat <- data.frame(Dij = x[,i])
        predict(model[[i]], newdata = dat, type = "response")
    }
    out <- sapply(colnames(x), FUN = foo, x = x, model = model)
    out
}
