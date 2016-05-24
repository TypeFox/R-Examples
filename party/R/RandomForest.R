
# $Id: RandomForest.R 540 2014-06-27 07:41:46Z thothorn $

### the fitting procedure
cforestfit <- function(object, controls, weights = NULL, fitmem = NULL, ...) {

    if (!extends(class(object), "LearningSample"))
        stop(sQuote("object"), " is not of class ", sQuote("LearningSample"))
    if (!extends(class(controls), "ForestControl"))
        stop(sQuote("controls"), " is not of class ", sQuote("ForestControl"))

    if (is.null(fitmem)) 
        fitmem <- ctree_memory(object, TRUE)
    if (!extends(class(fitmem), "TreeFitMemory"))
        stop(sQuote("fitmem"), " is not of class ", sQuote("TreeFitMemory"))

    if (is.null(weights))
        weights <- object@weights
    storage.mode(weights) <- "double"
    USER_WEIGHTS <- is.matrix(weights)
    if (USER_WEIGHTS) {
        controls@ntree <- ncol(weights)
        weights <- as.data.frame(weights)
        if (nrow(weights) != object@nobs)
            stop(sQuote("weights"), " are not a double matrix of ", 
                 object@nobs, " rows")
        bweights <- weights
        bwhere <- vector(mode = "list", length = controls@ntree)
        ### grow the tree
        ensemble <- .Call("R_Ensemble_weights", object, bwhere, bweights, fitmem, controls,
                          PACKAGE = "party")
    } else {
        if (length(weights) != object@nobs || storage.mode(weights) != "double")
            stop(sQuote("weights"), " are not a double vector of ", 
                 object@nobs, " elements")
        bweights <- vector(mode = "list", length = controls@ntree)
        bwhere <- vector(mode = "list", length = controls@ntree)
        ### grow the tree
        ensemble <- .Call("R_Ensemble", object, weights, bwhere, bweights, fitmem, controls,
                          PACKAGE = "party")
    }

    ### prepare the returned object
    RET <- new("RandomForest")
    RET@ensemble <- ensemble
    RET@where <- bwhere
    RET@weights <- bweights
    if (USER_WEIGHTS) {
        RET@initweights <- as.double(rep(1.0, object@nobs)) ### <FIXME>
    } else {
        RET@initweights <- as.double(weights)
    }
    RET@responses <- object@responses
    if (inherits(object, "LearningSampleFormula"))
        RET@data <- object@menv

    RET@update <- function(weights = NULL) {
        cforestfit(object = object, controls = controls,
                   weights = weights, fitmem = fitmem, ...)
    }


    ### (estimated) conditional distribution of the response given the
    ### covariates
    RET@cond_distr_response <- function(newdata = NULL, mincriterion = 0, ...) { 
        
        pw <- RET@prediction_weights(newdata = newdata, mincriterion =
                                     mincriterion, ...)

        response <- object@responses

        ### survival: estimated Kaplan-Meier
        if (any(response@is_censored)) {
            resp <- response@variables[[1]]
            RET <- lapply(pw, function(w) 
                mysurvfit(resp, weights = w))
            return(RET)
        }

        ### classification: estimated class probabilities
        ### regression: the means, not really a distribution
        RET <- lapply(pw, function(w) w %*% response@predict_trafo / sum(w))
        return(RET)
    }

    ### predict in the response space, always!
    RET@predict_response <- function(newdata = NULL, mincriterion = 0, 
        type = c("response", "prob"), ...) { 

        type <- match.arg(type)
        cdresp <- RET@cond_distr_response(newdata = newdata, 
                                          mincriterion = mincriterion, ...)
        if (type == "prob" || object@responses@ninputs > 1)
            return(cdresp)

        response <- object@responses
        ### classification: classes
        if (all(response@is_nominal || response@is_ordinal)) {
            lev <- levels(response@variables[[1]])
            RET <- factor(lev[unlist(lapply(cdresp, which.max))],
                          levels = levels(response@variables[[1]]))
            return(RET)
        }

        ### survival: median survival time
        if (any(response@is_censored)) {
            RET <- sapply(cdresp, mst)
            return(RET)
        }

        ### regression: mean (median would be possible)
        RET <- matrix(unlist(cdresp),
                      nrow = length(cdresp), byrow = TRUE)
        ### <FIXME> what about multivariate responses?
        colnames(RET) <- names(response@variables)
        ### </FIXME>
        return(RET)
    }

    ### get terminal node numbers
    RET@get_where <- function(newdata = NULL, mincriterion = 0, ...) {

        if (is.null(newdata) && mincriterion == 0) {
            if (all(where > 0)) return(RET@where)
        }

        newinp <- newinputs(object, newdata)

        lapply(ensemble, function(e) 
            R_get_nodeID(e, newinp, mincriterion))
    }

    RET@prediction_weights <- function(newdata = NULL, 
                                       mincriterion = 0, OOB = FALSE) {

        newinp <- newinputs(object, newdata)

        RET <- .Call("R_predictRF_weights", ensemble, bwhere, bweights, 
                     newinp, mincriterion, OOB && is.null(newdata),
                     PACKAGE = "party")
        names(RET) <- rownames(newinp@variables)
        RET
    }
    return(RET)
}

### the unfitted forest, an object of class `StatModel'
### see package `modeltools'
RandomForest <- new("StatModel",
                    capabilities = new("StatModelCapabilities"),
                    name = "random forest",
                    dpp = ctreedpp,
                    fit = cforestfit,
                    predict = function(object, ...) 
                        object@predict_response(...))

cforest_control <- function(teststat = "max", 
                            testtype = "Teststatistic",
                            mincriterion = qnorm(0.9),
                            savesplitstats = FALSE,
                            ntree = 500, mtry = 5, replace = TRUE, 
                            fraction = 0.632, 
                            trace = FALSE, ...) {

    if (is.null(mtry)) mtry <- 0
    tmp <- ctree_control(teststat = teststat, testtype = testtype,
                         mincriterion = mincriterion, 
                         savesplitstats = savesplitstats, 
                         mtry = mtry, ...)
    RET <- new("ForestControl")
    RET@ntree <- as.integer(ntree)
    RET@replace <- as.logical(replace)
    RET@fraction <- as.double(fraction)
    RET@trace <- as.logical(trace)
    RET <- copyslots(tmp, RET)
    if (!validObject(RET))
        stop("RET is not a valid object of class", class(RET))
    RET
}

cforest_classical <- function(...) cforest_control(teststat = "max",
                            testtype = "Teststatistic",
                            mincriterion = qnorm(0.9), 
                            replace = TRUE, ...)

cforest_unbiased <- function(...) cforest_control(teststat = "quad", 
                            testtype = "Univ",
                            mincriterion = 0,
                            replace = FALSE, 
                            fraction = 0.632, ...) 
    
### the top-level convenience function
cforest <- function(formula, data = list(), subset = NULL, weights = NULL, 
                    controls = cforest_unbiased(),
                    xtrafo = ptrafo, ytrafo = ptrafo, scores = NULL) {

    ### setup learning sample
    ls <- dpp(RandomForest, formula, data, subset, xtrafo = xtrafo, 
              ytrafo = ytrafo, scores = scores)

    ### setup memory
    fitmem <- ctree_memory(ls, TRUE)

    ### fit and return a conditional tree
    fit(RandomForest, ls, controls = controls, weights = weights, 
        fitmem = fitmem)
}

###
### variable importance for `cforest'
###
### see ?importance (in `randomForest'), too
###
###

### extract ID of _all_ variables the tree uses for splitting
varIDs <- function(node) {

    v <- c()
    foo <- function(node) {
        if (node[[4]]) return(NULL)
        v <<- c(v, node[[5]][[1]])
        foo(node[[8]])
        foo(node[[9]])
    }
    foo(node)
    return(v)
}

### calculate proximity matrix: p[i,j] = number of times obs i and j are 
### in the same terminal node
proximity <- function(object, newdata = NULL) {

    if (is.null(newdata)) {
        wh <- object@where
        rn <- rownames(object@data@get("response"))
    } else {
        wh <- object@get_where(newdata = newdata)
        rn <- rownames(newdata)
    }
    ### extract prediction weights
    prox <- .Call("R_proximity", wh, package = "party")
    prox <- matrix(unlist(prox), ncol = length(prox))
    rownames(prox) <- rn
    colnames(prox) <- rn
    prox
}


### FIXME: newdata may be missing, reuse weights
### partialPlot.BinaryTree?
partialPlot.party <-
    function (x, newdata, x.var, which.class, weights, plot = TRUE, add = FALSE,
              n.pt = min(length(unique(newdata[, xname])), 51), rug = TRUE,
              xlab = deparse(substitute(x.var)), ylab = "",
              main = paste("Partial Dependence on", deparse(substitute(x.var))),
              ...) 
{
    classRF <- all(x@responses@is_nominal || x@responses@is_ordinal)

    x.var <- substitute(x.var)
    xname <- if (is.character(x.var)) x.var else {
        if (is.name(x.var)) deparse(x.var) else {
            eval(x.var)
        }
    }

    if (!xname %in% names(newdata))
        stop("variable", " ", xname, " ", "not known")
    xv <- newdata[, xname]
    n <- nrow(newdata)

    if (missing(weights)) weights <- rep(1, n)

    if (classRF) {
        if (missing(which.class)) {
            focus <- 1
        }
        else {
            focus <- charmatch(which.class, levels(x@responses@variables[[1]]))
            if (is.na(focus)) 
                stop(which.class, " ", "is not one of the class labels.")
        }
    }
    if (is.factor(xv)) { ### includes ordered
        x.pt <- levels(xv)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- newdata
            x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt, ordered = is.ordered(xv))
            if (classRF) {
                pr <- treeresponse(x, newdata = x.data)
                pr <- matrix(unlist(pr), nrow = length(pr), byrow = TRUE)
                y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] > 0,
                                                    pr[, focus], 1)) -
                                         rowMeans(log(ifelse(pr > 0, pr, 1))),
                                         weights, na.rm=TRUE)
            } else y.pt[i] <- weighted.mean(predict(x, newdata = x.data), weights, na.rm = TRUE)

        }
        if (add) {
            points(1:length(x.pt), y.pt, type="h", lwd=2, ...)
        } else {
            if (plot) barplot(y.pt, width=rep(1, length(y.pt)), col="blue",
                              xlab = xlab, ylab = ylab, main=main,
                              names.arg=x.pt, ...)
        }
    } else {
        x.pt <- seq(min(xv), max(xv), length = n.pt)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- newdata
            x.data[, xname] <- rep(x.pt[i], n)
            class(x.data[, xname]) <- class(newdata[, xname])
            storage.mode(x.data[, xname]) <- storage.mode(newdata[, xname])
            if (classRF) {
                pr <- treeresponse(x, newdata = x.data)
                pr <- matrix(unlist(pr), nrow = length(pr), byrow = TRUE)
                y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] == 0, 1, pr[, focus]))
                                         - rowMeans(log(ifelse(pr == 0, 1, pr))),
                                         weights, na.rm=TRUE)
            } else {
                y.pt[i] <- weighted.mean(predict(x, newdata = x.data), weights, na.rm=TRUE)
            }
        }
        if (add) {
            lines(x.pt, y.pt, ...)
        } else {
            if (plot) plot(x.pt, y.pt, type = "l", xlab=xlab, ylab=ylab,
                           main = main, ...)
        }
        if (rug && plot) {
            if (n.pt > 10) {
                rug(quantile(xv, seq(0.1, 0.9, by = 0.1)), side = 1)
            } else {
                rug(unique(xv, side = 1))
            }
        }
    }
    invisible(list(x = x.pt, y = y.pt))
}
