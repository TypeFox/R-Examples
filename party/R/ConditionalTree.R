
# $Id: ConditionalTree.R 536 2014-06-26 12:08:54Z thothorn $

### the fitting procedure
ctreefit <- function(object, controls, weights = NULL, fitmem = NULL, ...) {

    if (!extends(class(object), "LearningSample"))
        stop(sQuote("object"), " is not of class ", sQuote("LearningSample"))
    if (!extends(class(controls), "TreeControl"))
        stop(sQuote("controls"), " is not of class ", sQuote("TreeControl"))

    if (is.null(fitmem)) 
        fitmem <- ctree_memory(object, TRUE)
    if (!extends(class(fitmem), "TreeFitMemory"))
        stop(sQuote("fitmem"), " is not of class ", sQuote("TreeFitMemory"))

    if (is.null(weights))
        weights <- object@weights
    storage.mode(weights) <- "double"
    if (length(weights) != object@nobs || storage.mode(weights) != "double")
        stop(sQuote("weights"), " are not a double vector of ", 
             object@nobs, " elements")
    if (max(abs(floor(weights) -  weights)) > sqrt(.Machine$double.eps))
        stop(sQuote("weights"), " contains real valued elements; currently
             only integer values are allowed") 

    where <- rep(0, object@nobs)
    storage.mode(where) <- "integer"

    ### grow the tree
    tree <- .Call("R_TreeGrow", object, weights, fitmem, controls, where,
                  PACKAGE = "party")

    ### create S3 classes and put names on lists
    tree <- prettytree(tree, names(object@inputs@variables), 
                       object@inputs@levels)

    ### prepare the returned object
    RET <- new("BinaryTree")
    RET@tree <- tree
    RET@where <- where
    RET@weights <- weights
    RET@responses <- object@responses
    if (inherits(object, "LearningSampleFormula"))
        RET@data <- object@menv

    RET@update <- function(weights = NULL) {
        ctreefit(object = object, controls = controls, 
                 weights = weights, fitmem = fitmem, ...)
    }

    ### get terminal node numbers
    RET@get_where <- function(newdata = NULL, mincriterion = 0, ...) {

        if (is.null(newdata) && mincriterion == 0) {
            if (all(where > 0)) return(where)
        }

        newinp <- newinputs(object, newdata)

        R_get_nodeID(tree, newinp, mincriterion)
    }

    ### (estimated) conditional distribution of the response given the
    ### covariates
    RET@cond_distr_response <- function(newdata = NULL, mincriterion = 0, ...) { 
        
        wh <- RET@get_where(newdata = newdata, mincriterion = mincriterion)

        response <- object@responses

        ### survival: estimated Kaplan-Meier
        if (any(response@is_censored)) {
            swh <- sort(unique(wh))
#            w <- .Call("R_getweights", tree, swh,
#                       PACKAGE = "party")
            RET <- vector(mode = "list", length = length(wh))
            resp <- response@variables[[1]]
            for (i in 1:length(swh)) {
                w <- weights * (where == swh[i])
                RET[wh == swh[i]] <- list(mysurvfit(resp, weights = w))
            }
            return(RET)
        }

        ### classification: estimated class probabilities
        ### regression: the means, not really a distribution
        RET <- .Call("R_getpredictions", tree, wh, PACKAGE = "party")
        return(RET)
    }

    ### predict in the response space, always!
    RET@predict_response <- function(newdata = NULL, mincriterion = 0, 
        type = c("response", "node", "prob"), ...) { 

        type <- match.arg(type)
        if (type == "node")
            return(RET@get_where(newdata = newdata, 
                                 mincriterion = mincriterion, ...))

        cdresp <- RET@cond_distr_response(newdata = newdata, 
                                          mincriterion = mincriterion, ...)
        if (type == "prob")
            return(cdresp)

        ### <FIXME> multivariate responses, we might want to return
        ###         a data.frame
        ### </FIXME>
        if (object@responses@ninputs > 1) 
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
        RET <- unlist(cdresp)
        RET <- matrix(unlist(RET),
                      nrow = length(RET), byrow = TRUE)
        ### <FIXME> what about multivariate responses?
        if (response@ninputs == 1)
            colnames(RET) <- names(response@variables)
        ### </FIXME>
        return(RET)
    }

    RET@prediction_weights <- function(newdata = NULL, 
                                       mincriterion = 0, ...) {

        wh <- RET@get_where(newdata = newdata, mincriterion = mincriterion)

        swh <- sort(unique(wh))
#        w <- .Call("R_getweights", tree, swh,
#                   PACKAGE = "party")
        RET <- vector(mode = "list", length = length(wh))   
        
        for (i in 1:length(swh))
            RET[wh == swh[i]] <- list(weights * (where == swh[i]))
        return(RET)
    }
    return(RET)
}

### data pre-processing (ordering, computing transformations etc)
ctreedpp <- function(formula, data = list(), subset = NULL, 
                     na.action = NULL, xtrafo = ptrafo, ytrafo = ptrafo, 
                     scores = NULL, ...) {

    dat <- ModelEnvFormula(formula = formula, data = data, 
                           subset = subset, designMatrix = FALSE, 
                           responseMatrix = FALSE, ...)
    inp <- initVariableFrame(dat@get("input"), trafo = xtrafo, 
                             scores = scores)

    response <- dat@get("response")

    if (any(is.na(response)))
        stop("missing values in response variable not allowed")

    resp <- initVariableFrame(response, trafo = ytrafo, response = TRUE,
                              scores = scores)

    RET <- new("LearningSampleFormula", inputs = inp, responses = resp,
               weights = rep(1, inp@nobs), nobs = inp@nobs,
               ninputs = inp@ninputs, menv = dat)
    RET
}

### the unfitted conditional tree, an object of class `StatModel'
### see package `modeltools'
conditionalTree <- new("StatModel",
                   capabilities = new("StatModelCapabilities"),
                   name = "unbiased conditional recursive partitioning",
                   dpp = ctreedpp,
                   fit = ctreefit,
                   predict = function(object, ...) 
                                 object@predict_response(...) )

### we need a `fit' method for data = LearningSample
setMethod("fit", signature = signature(model = "StatModel",
                                       data = "LearningSample"),
    definition = function(model, data, ...)
        model@fit(data, ...)
)

### control the hyper parameters
ctree_control <- function(teststat = c("quad", "max"),
                          testtype = c("Bonferroni", "MonteCarlo", "Univariate", "Teststatistic"),
                          mincriterion = 0.95, minsplit = 20, minbucket = 7, stump = FALSE,
                          nresample = 9999, maxsurrogate = 0, mtry = 0, 
                          savesplitstats = TRUE, maxdepth = 0) {

    teststat <- match.arg(teststat)
    testtype <- match.arg(testtype)
    RET <- new("TreeControl")
    if (teststat %in% levels(RET@varctrl@teststat)) {
        RET@varctrl@teststat <- factor(teststat, 
            levels = levels(RET@varctrl@teststat))
    } else {
        stop(sQuote("teststat"), teststat, " not defined")
    }

    if (testtype %in% levels(RET@gtctrl@testtype))
        RET@gtctrl@testtype <- factor(testtype, 
            levels = levels(RET@gtctrl@testtype))
    else
        stop(testtype, " not defined")

    if (RET@gtctrl@testtype == "Teststatistic")
        RET@varctrl@pvalue <- as.logical(FALSE)

    RET@gtctrl@nresample <- as.integer(nresample)
    RET@gtctrl@mincriterion <- as.double(mincriterion)
    if (all(mtry > 0)) {
        RET@gtctrl@randomsplits <- as.logical(TRUE)
        RET@gtctrl@mtry <- as.integer(mtry)
    }
    RET@tgctrl@savesplitstats <- as.logical(savesplitstats)
    RET@splitctrl@minsplit <- as.double(minsplit)
    RET@splitctrl@maxsurrogate <- as.integer(maxsurrogate)
    RET@splitctrl@minbucket <- as.double(minbucket)
    RET@tgctrl@stump <- as.logical(stump)
    RET@tgctrl@maxdepth <- as.integer(maxdepth)
    RET@tgctrl@savesplitstats <- as.logical(savesplitstats)
    if (!validObject(RET))
        stop("RET is not a valid object of class", class(RET))
    RET
}

### the top-level convenience function
ctree <- function(formula, data = list(), subset = NULL, weights = NULL, 
                  controls = ctree_control(), xtrafo = ptrafo, 
                  ytrafo = ptrafo, scores = NULL) {

    ### setup learning sample
    ls <- dpp(conditionalTree, formula, data, subset, xtrafo = xtrafo, 
              ytrafo = ytrafo, scores = scores)

    ### setup memory
    fitmem <- ctree_memory(ls, TRUE)

    ### fit and return a conditional tree
    fit(conditionalTree, ls, controls = controls, weights = weights, 
        fitmem = fitmem)
}
