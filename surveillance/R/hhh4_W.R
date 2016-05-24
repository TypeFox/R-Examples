################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Helper functions for neighbourhood weight matrices in hhh4()
###
### Copyright (C) 2012-2016 Sebastian Meyer
### $Revision: 1687 $
### $Date: 2016-04-01 21:40:25 +0200 (Fre, 01. Apr 2016) $
################################################################################


checkNeighbourhood <- function (neighbourhood)
{
    ## setValidity() in sts.R only guarantees correct 'dim' and 'dimnames'
    ## we also assert numeric or logical matrix with non-NA entries
    ## FIXME: However, we currently don't check for symmetry and for zeros on
    ## the diagonal...
    stopifnot(is.matrix(neighbourhood),
              nrow(neighbourhood) == ncol(neighbourhood),
              is.numeric(neighbourhood) | is.logical(neighbourhood),
              is.finite(neighbourhood))
    invisible(TRUE)
}


### calculate the weighted sum of counts of adjacent (or all other) regions
### i.e. the nTime x nUnit matrix with elements ne_ti = sum_j w_jit * y_jt
## W is either a nUnits x nUnits matrix of time-constant weights w_ji
## or a nUnits x nUnits x nTime array of time-varying weights

weightedSumNE <- function (observed, weights, lag)
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  tY <- t(observed)                     # -> nUnits x nTime
  
  res <- apply(weights, 2L, function (wi)
               ## if dim(weights)==2 (time-constant weights), length(wi)=nUnits,
               ## if dim(weights)==3, wi is a matrix of size nUnits x nTime
               .colSums(tY * wi, nUnits, nTime, na.rm=TRUE))
  
  rbind(matrix(NA_real_, lag, nUnits),
        res[seq_len(nTime-lag),,drop=FALSE])
}


### normalize weight matrix such that each row sums to 1 (at each time point)

normalizeW <- function (W)
{
    dimW <- dim(W)
    if (length(dimW) == 2L) {
        W / .rowSums(W, dimW[1L], dimW[2L])
    } else { # time-varying weights
        res <- apply(W, 3L, normalizeW)
        dim(res) <- dimW
        res
    }
}


### scale and/or normalize a weight matrix/array

scaleNEweights.default <- function (weights, scale = NULL, normalize = FALSE)
{
    if (!is.null(scale))
        weights <- scale * weights
    if (normalize)
        weights <- normalizeW(weights)
    weights
}

## update parametric weights functions w, dw, d2w
scaleNEweights.list <- function (weights, scale = NULL, normalize = FALSE)
{
    if (is.null(scale) && !normalize)
        return(weights)
    
    if (normalize) {
        dprod <- function (u, v, du, dv) du * v + u * dv
        dfrac <- function (u, v, du, dv) (du * v - u * dv) / v^2
        w <- function (...)
            scaleNEweights.default(weights$w(...), scale, TRUE)
        dw <- function (...) {
            W <- scaleNEweights.default(weights$w(...), scale)
            dW <- clapply(X = weights$dw(...), # matrix or list thereof
                          FUN = scaleNEweights.default,
                          scale = scale) # always returns a list
            dimW <- dim(W)
            normW <- .rowSums(W, dimW[1L], dimW[2L])
            normdW <- lapply(X = dW, FUN = .rowSums, m = dimW[1L], n = dimW[2L])
            mapply(FUN = dfrac, du = dW, dv = normdW,
                   MoreArgs = list(u = W, v = normW),
                   SIMPLIFY = FALSE, USE.NAMES = FALSE)
        }
        ## for d2w() we need all the stuff from dw() -> substitute
        d2w <- as.function(c(alist(...=), substitute({
            dWnorm <- DWBODY
            d2W <- clapply(X = weights$d2w(...), # matrix or list thereof
                           FUN = scaleNEweights.default,
                           scale = scale) # always returns a list
            normd2W <- lapply(X = d2W, FUN = .rowSums, m = dimW[1L], n = dimW[2L])
            ## order of d2w is upper triangle BY ROW
            dimd <- length(dW)
            ri <- rep.int(seq_len(dimd), rep.int(dimd, dimd)) # row index
            ci <- rep.int(seq_len(dimd), dimd)                # column index
            uppertri <- ci >= ri
            mapply(FUN = function (k, l, d2W, normd2W) {
                       dfrac(dW[[k]], normW, d2W, normdW[[l]]) -
                           dprod(W/normW, normdW[[k]]/normW, dWnorm[[l]],
                                 dfrac(normdW[[k]], normW, normd2W, normdW[[l]]))
                   },
                   k = ri[uppertri], l = ci[uppertri],
                   d2W = d2W, normd2W = normd2W,
                   SIMPLIFY = FALSE, USE.NAMES = FALSE)
        }, list(DWBODY = body(dw)))))
    } else {
        w <- function (...)
            scaleNEweights.default(weights$w(...), scale)
        dw <- function (...)
            clapply(X = weights$dw(...),
                    FUN = scaleNEweights.default,
                    scale = scale)
        d2w <- function (...)
            clapply(X = weights$d2w(...),
                    FUN = scaleNEweights.default,
                    scale = scale)
    }
    
    ## return list with updated functions
    list(w = w, dw = dw, d2w = d2w, initial = weights$initial)
}


##################################
### check ne$weights specification
##################################


### checks for a fixed matrix/array

checkWeightsArray <- function (W, nUnits, nTime, name = deparse(substitute(W)),
                               check0diag = FALSE, islands = FALSE)
{
    if (!is.array(W) || !(length(dim(W)) %in% 2:3))
        stop("'", name, "' must return a matrix or 3-dim array")
    if (any(dim(W)[1:2] != nUnits) || isTRUE(dim(W)[3] != nTime))
        stop("'", name, "' must conform to dimensions ",
             nUnits, " x ", nUnits, " (x ", nTime, ")")
    if (any(is.na(W))) {
        if (islands) # normalization of parametric weights yields division by 0
            warning("neighbourhood structure contains islands")
        stop("missing values in '", name, "' are not allowed")
    }
    if (check0diag) {
        diags <- if (is.matrix(W)) diag(W) else apply(W, 3, diag)
        if (any(diags != 0))
            warning("'", name, "' has nonzeros on the diagonal",
                    if (!is.matrix(W)) "s")
    }
}


### check parametric weights specification consisting of a list of:
## - three functions: w, dw, and d2w
## - a vector of initial parameter values

checkWeightsFUN <- function (object)
{
    fnames <- paste0(c("","d","d2"), "w")
    if (any(!sapply(object[fnames], is.function)))
        stop("parametric weights require functions ",
             paste0("'", fnames, "'", collapse=", "))
    if (any(!sapply(object[fnames], function(FUN) length(formals(FUN)) >= 3L)))
        stop("parametric weights functions must accept (not necessarily use)",
             "\n  at least 3 arguments (parameter vector, ",
             "neighbourhood order matrix, data)")
    if (!is.vector(object$initial, mode="numeric") ||
        length(object$initial) == 0L)
        stop("parametric weights require initial parameter values")
    TRUE
}


### entry function for checks in hhh4()

checkWeights <- function (weights, nUnits, nTime,
                          nbmat, data,  # only used for parametric weights
                          check0diag = FALSE)
{
    name <- deparse(substitute(weights)) # "control$ne$weights"

    ## check specification
    testweights <- if (is.array(weights)) weights else {
        if (is.list(weights) && checkWeightsFUN(weights)
            && checkNeighbourhood(nbmat)) {
            if (all(nbmat %in% 0:1))
                warning("'", deparse(substitute(nbmat)),
                        "' is binary (should contain",
                        " general neighbourhood orders)")
            weights$w(weights$initial, nbmat, data)
        } else {
            stop("'", name, "' must be a matrix/array or a list of functions")
        }
    }
    
    ## apply matrix/array checks
    if (is.list(weights)) { # parametric weights
        if (length(dim(testweights)) > 2L)
            warning("time-varying parametric weights are not fully supported")
        checkWeightsArray(testweights, nUnits, nTime, name = paste0(name, "$w"),
                          check0diag = check0diag,
                          islands = any(.rowSums(nbmat, nUnits, nUnits) == 0))
        dim.d <- length(weights$initial)
        dw <- weights$dw(weights$initial, nbmat, data)
        d2w <- weights$d2w(weights$initial, nbmat, data)
        if (dim.d == 1L && !is.list(dw) && !is.list(d2w)) {
            checkWeightsArray(dw, nUnits, nTime, name=paste0(name, "$dw"))
            checkWeightsArray(d2w, nUnits, nTime, name=paste0(name, "$d2w"))
        } else {
            if (!is.list(dw) || length(dw) != dim.d)
                stop("'", name, "$dw' must return a list (of matrices/arrays)",
                     " of length ", dim.d)
            if (!is.list(d2w) || length(d2w) != dim.d*(dim.d+1)/2)
                stop("'", name, "$d2w' must return a list (of matrices/arrays)",
                     " of length ", dim.d*(dim.d+1)/2)
            lapply(dw, checkWeightsArray, nUnits, nTime,
                   name=paste0(name, "$dw[[i]]"))
            lapply(d2w, checkWeightsArray, nUnits, nTime,
                   name=paste0(name, "$d2w[[i]]"))
        }
    } else checkWeightsArray(testweights, nUnits, nTime, name = name,
                             check0diag = check0diag)
    
    ## Done
    invisible(TRUE)
}



#############################################
### Utility functions for fitted hhh4-objects
#############################################


### extract the (final) weight matrix/array from a fitted hhh4 object

getNEweights <- function (object, pars = coefW(object),
                          scale = ne$scale, normalize = ne$normalize)
{
    ne <- object$control$ne
    weights <- if (is.list(ne$weights)) { # parametric weights
        nd <- length(ne$weights$initial)
        if (length(pars) != nd) stop("'pars' must be of length ", nd)
        ne$weights$w(pars, neighbourhood(object$stsObj), object$control$data)
    } else { # NULL or fixed weight structure
        ne$weights
    }
    if (is.null(normalize)) normalize <- FALSE  # backward compatibility < 1.9-0
    scaleNEweights.default(weights, scale, normalize)
}


### extract parameters of neighbourhood weights from hhh4-object or coef vector

coefW <- function (object)
{
    coefs <- if (inherits(object, "hhh4")) object$coefficients else object
    coefW <- coefs[grep("^neweights", names(coefs))]
    names(coefW) <- sub("^neweights\\.", "", names(coefW))
    coefW
}
