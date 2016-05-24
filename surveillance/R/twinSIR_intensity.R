################################################################################
# Authors: Sebastian Meyer, with contributions by Michael Hoehle
# Date: 02 June 2009, modified 25 Mar 2011, 27 Jun 2012
#
# This file contains functions related to calculating and plotting intensities.
################################################################################


################################################################################
# Calculate the two components of the intensity lambda(t|H_t) for each row
# of the event history.
# Be aware that the function assumes atRiskY == 1 in all rows!
#
# ARGS:
#  theta - parameter vector c(alpha,beta), where
#          beta also contains the baseline coefficients in the first place
#  X     - covariate matrix related to alpha, i.e. the epidemic component
#  Z     - covariate matrix related to beta, i.e. the Cox-like endemic component
#
# RETURNS: a numeric matrix with two columns e and h and nrow(X)==nrow(Z) rows
################################################################################

.eh <- function(theta, X, Z)
{
    # Extracting params from theta
    dimX <- dim(X)
    nRows <- dimX[1] # = nrow(Z)
    px <- dimX[2]
    pz <- ncol(Z)
    alpha <- theta[seq_len(px)]
    beta <- theta[px + seq_len(pz)]
    
    # Calculate the epidemic component e(t|H_t) and the endemic component h(t)
    e <- if (px > 0L) drop(X %*% alpha) else numeric(nRows)
    h <- if (pz > 0L) drop(exp(Z %*% beta)) else numeric(nRows)
    
    # Return the two components of the infection intensity related to the
    # rows of the event history in a two column matrix
    eh <- cbind(e = e, h = h)
    return(eh)
}


################################################################################
# Cumulative hazard function
#
#      \Lambda(t) = \int_{timeRange[1]}^t \lambda(s) ds,
#
#  where \lambda(s) = \sum_{i=1}^n \lambda_i(s)
#
# Be aware that the function assumes atRiskY == 1 for all rows of X/Z/survs !!!
#
# ARGS:
#  t     - scalar time point until we want to integrate, must be non-negative
#  theta - parameter vector c(alpha,beta), where
#          beta also contains the baseline coefficients in the first place
#  X     - covariate matrix related to alpha, i.e. the epidemic component
#  Z     - covariate matrix related to beta, i.e. the Cox-like endemic component
#  survs - data.frame with columns id, start, stop, event; "timeRange" attribute
#  weights - vector of length nrow(X) indicating the number of individuals
#            with the same covariates. weights are allowed to change over time.
#            Note: it is assumed that none of the individuals covered by
#            "weights"  can have an actual event, if so they need to have their
#            own row
#
# RETURNS: value of the cumulative hazard function at time t
################################################################################

Lambda <- function(t, theta, X, Z, survs, weights)
{
    timeRange <- attr(survs, "timeRange")
    eh <- if (!isScalar(t) || t < timeRange[1L]) {
        stop("invalid argument 't': must be a scalar >= ", timeRange[1L],
             " (beginning of observation period)")
    } else if (t == timeRange[1L]) {
        return(0)
    } else if (t < timeRange[2L]) {
        # We have to extract the relevant intervals
        sortedStop <- sort(unique(survs$stop))
        # Find first stop time beyond t
        idx <- match(TRUE, sortedStop >= t)
        firstBeyondt <- sortedStop[idx]
        includeSurvsRow <- survs$stop <= firstBeyondt
        # If t between start and stop of an interval we need to chop...
        if (firstBeyondt != t) {
            survs$stop[survs$stop == firstBeyondt] <- t
        }
        # Extract relevant parts
        survs <- survs[includeSurvsRow,]
        weights <- weights[includeSurvsRow]
        .eh(theta, X[includeSurvsRow,], Z[includeSurvsRow,])
    } else { # if t >= attr(survs, "timeRange")[2], we take all rows
        .eh(theta, X, Z)
    }
    
    lambda <- rowSums(eh)
    dt <- survs$stop - survs$start
    intlambda <- sum(weights * lambda * dt)   # no individual sums as in loglik
    return(intlambda)
}



################################################################################
# Function to plot the path of the infection intensity or the proportions of
# the endemic or epidemic component, either on an individual basis or related
# to the total intensity at each event (=infection) time.
# The function works with objects of class "simEpidata"
# as well as with objects of class "twinSIR".
################################################################################

# 'model' is the result of getModel(x)
# if x is of class "twinSIR": theta = (alpha, beta) = (alpha, (h0coefs, betarest)) 
# if x is of class "simEpidata": theta = (alpha, 1, betarest)
# per default, the function uses the fitted or true parameters, respectively
intensityplot_twinSIR <- function(model,
    which = c("epidemic proportion", "endemic proportion", "total intensity"),
    aggregate = TRUE, theta = NULL, plot = TRUE, add = FALSE, rug.opts = list(), ...)
{
    which <- match.arg(which)
    
    ## model components
    survs <- model$survs
    start <- attr(survs, "timeRange")[1L]
    end <- attr(survs, "timeRange")[2L]
    timeIntervals <- unique(survs[c("start", "stop")])
    timepoints <- unique(c(timeIntervals$stop, end))
    # need 'end' here, because model does only contain rows with atRiskY == 1,
    # otherwise would terminate in advance if all individuals have been infected
    nTimes <- length(timepoints)
    idlevels <- levels(survs$id)
    
    ## helper function for use with by()
    intensity <- function(iddata, what) {
        # 'iddata' will be a subset of survs, 'what' will be "wlambda" or "we"
        y <- numeric(nTimes)            # zeroes
        y[match(iddata$stop, timepoints)] <- iddata[[what]]
        y
    }
    
    ## Calculate epidemic (e) and endemic (h) component in each row of the model
    eh <- do.call(".eh", args = c(list(theta = theta), model[c("X", "Z")]))
    
    ## Calculate individual _total intensity_ paths
    lambda <- rowSums(eh)
    survs$wlambda <- as.vector(model$weights * lambda)
    
    ## put individual intensity paths into a matrix [nTimes x n]
    wlambdaID <- by(data = survs, INDICES = survs["id"],
                    FUN = intensity, what = "wlambda", simplify = FALSE)
    # initially infectious individuals (without re-infection) don't appear in
    # survs, since they are never atRiskY => wlambdaID[[i]] is NULL for such an
    # individual i but should be a 0-vector of length nTimes
    initiallyInfected <- names(which(sapply(wlambdaID, is.null)))
    #if (length(initiallyInfected) > 0L)   # not necessary
    wlambdaID[initiallyInfected] <- rep(list(numeric(nTimes)), length(initiallyInfected))
    wlambdaIDmatrix <- as.matrix(as.data.frame(c(wlambdaID), optional = TRUE))

    ## alternative way but slower:
    ## wlambdaIDmatrix <- matrix(0, nrow = nTimes, ncol = length(idlevels),
    ##                           dimnames = list(NULL, idlevels))
    ## for (ID in idlevels) {
    ##     iddata <- survs[survs$id == ID,]
    ##     wlambdaIDmatrix[match(iddata$stop, timepoints), ID] <- iddata$wlambda
    ## }
    
    if (which != "total intensity") {
        ## Calculate individual _epidemic intensity_ paths
        survs$we <- {
            px <- ncol(model$X)
            if (px == 0L) {
                stop("nothing to do, model does not contain both components")
            }
            as.vector(model$weights * eh[,1])
        }
        ## put individual epidemic intensity paths into a matrix [nTimes x n]
        weID <- by(data = survs, INDICES = list(id = survs$id),
                   FUN = intensity, what = "we", simplify = FALSE)
        # we have to replace NULL entries by numeric(nTimes) (cf. wlambdaID)
        weID[initiallyInfected] <- rep(list(numeric(nTimes)), length(initiallyInfected))
        weIDmatrix <- as.matrix(as.data.frame(c(weID), optional = TRUE))

        ## alternative code which is slower:
        ## weIDmatrix <- matrix(0, nrow = nTimes, ncol = length(idlevels),
        ##                       dimnames = list(NULL, idlevels))
        ## for (ID in idlevels) {
        ##     iddata <- survs[survs$id == ID,]
        ##     weIDmatrix[match(iddata$stop, timepoints), ID] <- iddata$we
        ## }
    }
    
    ## Generate matrix with data for 'matplot'
    ydata2plot <-
        if (which == "total intensity") {
            if (aggregate) {
                rowSums(wlambdaIDmatrix)
            } else {
                wlambdaIDmatrix
            }
        } else {   # calculate epidemic proportion
            if (aggregate) {
                rowSums(weIDmatrix) / rowSums(wlambdaIDmatrix)
            } else {
                weIDmatrix / wlambdaIDmatrix
            }
        }
    if (which == "endemic proportion") {
        ydata2plot <- 1 - ydata2plot
    }
    ydata2plot <- as.matrix(ydata2plot)
    colnames(ydata2plot) <- if (aggregate) which else idlevels
    
    if (which != "total intensity") {
        # there may be NAs in data2plot where the total intensity equals 0
        # => when calculating proportions we get 0 / 0 = NA
        # we redefine those values to 0. (0-intensity => 0-proportion)
        ydata2plot[is.na(ydata2plot)] <- 0
    }
    
    # prepend time (x) column
    data2plot <- cbind(stop = timepoints, ydata2plot)
    
    # if the epidemic is SIRS or SIS (re-susceptibility), there may be time
    # blocks during the observation period, where no individual is susceptible:
    # Problem: those time blocks are not included in the model component,
    #          which only contains rows with atRiskY == 1
    # Solution: fill the missing time periods with 0 intensity (or proportion)
    innerStart <- timeIntervals[-1L, "start"]
    innerStop <- timeIntervals[-nrow(timeIntervals), "stop"]
    noSusceptiblesStopTimes <- innerStart[innerStop != innerStart]
    if (length(noSusceptiblesStopTimes) > 0L) {
        data2plot <- rbind(data2plot,
            cbind(noSusceptiblesStopTimes,
                  matrix(0, nrow = length(noSusceptiblesStopTimes),
                            ncol = ncol(ydata2plot))
            )
        )
        data2plot <- data2plot[order(data2plot[,1L]),]
    }
    
    ## Plot and return data
    if (plot) {
        dotargs <- list(...)
        nms <- names(dotargs)
        if(! "xlab" %in% nms) dotargs$xlab <- "time"
        if(! "ylab" %in% nms) dotargs$ylab <- which
        if(! "lty" %in% nms) dotargs$lty <- 1
        do.call("matplot",
                args = c(list(x = c(start, data2plot[,1L]),
                              y = rbind(data2plot[1L, -1L, drop = FALSE],
                                        data2plot[  , -1L, drop = FALSE]),
                              type = "S", add = add),
                         dotargs))
        if (is.list(rug.opts)) {
            if (is.null(rug.opts$ticksize)) rug.opts$ticksize <- 0.02
            if (is.null(rug.opts$quiet)) rug.opts$quiet <- TRUE
            do.call("rug", args = c(list(x = attr(survs, "eventTimes")),
                                    rug.opts))
        }
        invisible(data2plot)
    } else {
        data2plot
    }
}


### intensityplot-methods for objects of classes "twinSIR" and "simEpidata"

intensityplot.twinSIR <- function ()
{
    cl <- match.call()
    cl[[1]] <- as.name("intensityplot_twinSIR")
    names(cl)[names(cl) == "x"] <- "model"
    cl$model <- quote(getModel(x))

    if (is.null(theta)) {
        cl$theta <- quote(coef(x))
    }

    eval(cl)
}

intensityplot.simEpidata <- function ()
{
    cl <- match.call()
    cl[[1]] <- as.name("intensityplot_twinSIR")
    names(cl)[names(cl) == "x"] <- "model"
    cl$model <- quote(getModel(x))

    if (is.null(theta)) {
        config <- attr(x, "config")
        cl$theta <- quote(c(config$alpha, 1, config$beta))   # 1 is for true h0
    }

    message("Note: the (true) baseline hazard is only evaluated",
            " at the beginning of the time intervals")
    eval(cl)
}

formals(intensityplot.twinSIR) <- formals(intensityplot.simEpidata) <-
    c(alist(x=), formals(intensityplot_twinSIR)[-1])


