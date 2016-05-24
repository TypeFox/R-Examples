# This file contains copies of methods from the pls package adapted to methods not included
# in the pls package itself. Copies were made of pls 2.4-3 in march 2015.
#' @title Multivariate regression function
#' 
#' @description Adaptation of \code{mvr} from package \code{pls} v 2.4.3.
#' 
#' @param formula	a model formula. Most of the lm formula constructs are supported. See below.
#' @param ncomp	the number of components to include in the model (see below).
#' @param Y.add	a vector or matrix of additional responses containing relevant information about the observations. Only used for cppls.
#' @param data an optional data frame with the data to fit the model from.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action	a function which indicates what should happen when the data contain missing values. The default is set by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful. See na.omit for other alternatives.
#' @param shrink optional shrinkage parameter for \code{stpls}.
#' @param method the multivariate regression method to be used. If "model.frame", the model frame is returned.
#' @param scale	numeric vector, or logical. If numeric vector, X is scaled by dividing each variable with the corresponding element of scale. If scale is TRUE, X is scaled by dividing each variable by its sample standard deviation. If cross-validation is selected, scaling by the standard deviation is done for every segment.
#' @param validation character. What kind of (internal) validation to use. See below.
#' @param model	a logical. If TRUE, the model frame is returned.
#' @param x	a logical. If TRUE, the model matrix is returned.
#' @param y	a logical. If TRUE, the response is returned.
#' @param ...	additional arguments, passed to the underlying fit functions, and mvrCv.
#' 
#' @import parallel
#' @seealso \code{\link[pls]{mvr}}
#' @export
mvrV <- function(formula, ncomp, Y.add, data, subset, na.action, shrink,
                method = c("truncation","stpls","model.frame"), #, "ridgepls"
                scale = FALSE, validation = c("none", "CV", "LOO"),
                model = TRUE, x = FALSE, y = FALSE, ...)
{
  ret.x <- x                          # More useful names
  ret.y <- y
  
  ## Get the model frame
  mf <- match.call(expand.dots = FALSE)
  if (!missing(Y.add)) {
    ## Temporarily add Y.add to the formula
    Y.addname <- as.character(substitute(Y.add))
    mf$formula <- update(formula, paste("~ . +", Y.addname))
  }
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]                # Retain only the named arguments
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  method <- match.arg(method, c("stpls","truncation", "model.frame"))
  if (method == "model.frame") return(mf)
  ## Get the terms
  mt <- attr(mf, "terms")        # This is to include the `predvars'
  # attribute of the terms
  ## Get the data matrices
  Y <- model.response(mf, "numeric")
  if (is.matrix(Y)) {
    if (is.null(colnames(Y)))
      colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
  } else {
    Y <- as.matrix(Y)
    colnames(Y) <- deparse(formula[[2]])
  }
  if (missing(Y.add)) {
    Y.add <- NULL
  } else {
    Y.add <- mf[,Y.addname]
    ## Remove Y.add from the formula again
    mt <- drop.terms(mt, which(attr(mt, "term.labels") == Y.addname),
                     keep.response = TRUE)
  }
  X <- delete.intercept(model.matrix(mt, mf))
  
  nobj <- dim(X)[1]
  npred <- dim(X)[2]
  
  ## model.matrix prepends the term name to the colnames of matrices.
  ## If there is only one predictor term, and the corresponding matrix
  ## has colnames, remove the prepended term name:
  if (length(attr(mt, "term.labels")) == 1 &&
      !is.null(colnames(mf[[attr(mt, "term.labels")]])))
    colnames(X) <- sub(attr(mt, "term.labels"), "", colnames(X))
  
  ## Set or check the number of components:
  if (missing(ncomp)) {
    ncomp <- min(nobj - 1, npred)
    ncompWarn <- FALSE              # Don't warn about changed `ncomp'
  } else {
    if (ncomp < 1 || ncomp > min(nobj - 1, npred))
      stop("Invalid number of components, ncomp")
    ncompWarn <- TRUE
  }
  
  ## Handle any fixed scaling before the the validation
  sdscale <- identical(TRUE, scale)   # Signals scaling by sd
  if (is.numeric(scale))
    if (length(scale) == npred)
      X <- X / rep(scale, each = nobj)
    else stop("length of 'scale' must equal the number of x variables")
    
    ## Optionally, perform validation:
    switch(match.arg(validation),
           CV = {
             val <- mvrCvV(X, Y, ncomp, Y.add = Y.add, method = method, scale = sdscale, shrink = shrink, ...)
           },
           LOO = {
             segments <- as.list(1:nobj)
             attr(segments, "type") <- "leave-one-out"
             val <- mvrCvV(X, Y, ncomp, Y.add = Y.add, method = method, scale = sdscale, shrink = shrink,
                          segments = segments, ...)
           },
           none = {
             val <- NULL
           }
    )
    ## Check and possibly adjust ncomp:
    if (identical(TRUE, ncomp > val$ncomp)) {
      ncomp <- val$ncomp
      if (ncompWarn) warning("`ncomp' reduced to ", ncomp,
                             " due to cross-validation")
    }
    
    ## Select fit function:
    fitFunc <- switch(method,
                      truncation = trunc_fit,
                      stpls = stpls.fit)
#                      ridgepls = ridgepls.fit)
    
    ## Perform any scaling by sd:
    if (sdscale) {
      ## This is faster than sd(X), but cannot handle missing values:
      scale <- sqrt(colSums((X - rep(colMeans(X), each = nobj))^2) /
                      (nobj - 1))
      if (any(abs(scale) < .Machine$double.eps^0.5))
        warning("Scaling with (near) zero standard deviation")
      X <- X / rep(scale, each = nobj)
    }
    
    ## Fit the model:
    start.time <- proc.time()[3]
    z <- fitFunc(X, Y, ncomp, Y.add = Y.add, shrink = shrink, ...)
    z$fit.time <- proc.time()[3] - start.time
    
    ## Build and return the object:
    class(z) <- c("mvrV","mvr")
    z$na.action <- attr(mf, "na.action")
    z$ncomp <- ncomp
    if(!missing(shrink))
      z$shrink <- shrink
    z$method <- method
    if (is.numeric(scale)) z$scale <- scale
    z$validation <- val
    z$call <- match.call()
    z$terms <- mt
    if (model) z$model <- mf
    if (ret.x) z$x <- X
    if (ret.y) z$y <- Y
    z
}


## The basic cross-validation function
mvrCvV <- function(X, Y, ncomp, Y.add = NULL, weights = NULL, shrink,
                   method = c("truncation","stpls","model.frame"), #, "ridgepls"
                   scale = FALSE, segments = 10,
                   segment.type = c("random", "consecutive", "interleaved"),
                   length.seg, jackknife = FALSE, trace = FALSE, ...)
{
  ## Initialise:
  Y <- as.matrix(Y)
  if(!(missing(Y.add) || is.null(Y.add)))
    Y.add <- as.matrix(Y.add)
  
  ## Save dimnames:
  dnX <- dimnames(X)
  dnY <- dimnames(Y)
  
  ## Remove dimnames for performance (doesn't seem to matter; in fact,
  ## as far as it has any effect, it hurts a tiny bit in most situations).
  ## dimnames(X) <- dimnames(Y) <- NULL
  
  ## Save dimensions:
  nobj <- dim(X)[1]
  npred <- dim(X)[2]
  nresp <- dim(Y)[2]
  if(!missing(shrink)){
    nshrink <- length(shrink)
  } else {
    shrink <- FALSE
  }
  
  ## Check the `scale' parameter:
  if (!is.logical(scale) || length(scale) != 1)
    stop("'scale' must be 'TRUE' or 'FALSE'")
  
  ## Set up segments:
  if (is.list(segments)) {
    if (is.null(attr(segments, "type")))
      attr(segments, "type") <- "user supplied"
  } else {
    if (missing(length.seg)) {
      segments <- cvsegments(nobj, k = segments, type = segment.type)
    } else {
      segments <- cvsegments(nobj, length.seg = length.seg,
                             type = segment.type)
    }
  }
  
  ## Reduce ncomp, if neccessary:
  ncomp <- min(ncomp, nobj - max(sapply(segments, length)) - 1)
  
  ## Select fit function:
  method <- match.arg(method,c("truncation","stpls")) #, "ridgepls"
  fitFunc <- switch(method,
                    truncation = trunc_fit,
                    stpls = stpls.fit)
#                    ridgepls = ridgepls.fit)
  
  ## Helper function to perform the cross-validatoin for one segment.
  ## Defined inside mvrCv to be able to access local variables:
  mvrCvSeg <- function(n.seg) {
    if (trace) cat(n.seg, "")
    
    ## Set up train data:
    seg <- segments[[n.seg]]
    Xtrain <- X[-seg,, drop=FALSE]
    if (scale) {
      ntrain <- nrow(Xtrain)
      ## This is faster than sd(X), but cannot handle missing values:
      sdtrain <-
        sqrt(colSums((Xtrain - rep(colMeans(Xtrain), each = ntrain))^2) /
               (ntrain - 1))
      if (any(abs(sdtrain) < .Machine$double.eps^0.5))
        warning("Scaling with (near) zero standard deviation")
      Xtrain <- Xtrain / rep(sdtrain, each = ntrain)
    }
    
    ## Fit the model:
    fit <- fitFunc(Xtrain, Y[-seg,, drop=FALSE], ncomp, shrink = shrink,
                   Y.add = Y.add[-seg,, drop=FALSE], stripped = TRUE,
                   weights = weights[-seg], ...)
    
    ## Set up test data:
    Xtest <- X
    if (scale) Xtest <- Xtest / rep(sdtrain, each = nobj)
    Xtest <- Xtest - rep(fit$Xmeans, each = nobj)
    
    ## Predict test data:
    if(!is.logical(shrink)){ # ST-PLS shrinkage
      pred <- array(dim = c(nobj, nresp, ncomp, nshrink))
      for(aa in 1:nshrink){
        for (a in 1:ncomp){ 
          pred[, , a, aa] <- sweep(Xtest %*% fit$coefficients[, , a, aa], 2, fit$Ymeans, "+")
        }
      }    
      cvPred <- pred[seg,,,, drop = FALSE]
    } else { # No shrinkage
      pred <- array(0, dim = c(nobj, nresp, ncomp))
      Ymeansrep <- rep(fit$Ymeans, each = nobj)
      for (a in 1:ncomp)
        pred[,,a] <- Xtest %*% fit$coefficients[,,a] + Ymeansrep
      cvPred <- pred[seg,,, drop=FALSE]
    }
    if(!is.logical(shrink)){ # ST-PLS shrinkage
      adjCur <- length(seg) * apply((pred - c(Y))^2,2:4,sum)
    } else {
      adjCur <- length(seg) * colSums((pred - c(Y))^2)
    }
    return(list(adj = adjCur,
                cvPred = cvPred, #pred[seg,,, drop=FALSE],
                gammas = if (method == "cppls") fit$gammas else NULL,
                cvCoef = if (jackknife) fit$coefficients else NULL
    ))
  }
  
  ## Perform the cross-validation, optionally in parallel:
  if (trace) cat("Segment: ")
  results <- lapplyFunc(pls.options()$parallel, seq_along(segments), mvrCvSeg)
  if (trace) cat("\n")
  
  ## Variables to save CV results in:
  if(!is.logical(shrink)){
    adj <- array(0,dim=c(nresp, ncomp, nshrink))
    cvPred <- array(dim = c(nobj, nresp, ncomp, nshrink))
    if (jackknife)
      cvCoef <- array(dim = c(npred, nresp, ncomp, nshrink, length(segments)))
    ## Collect the results:
    for (n.seg in seq_along(segments)) {
      res <- results[[n.seg]]
      adj <- adj + res$adj
      cvPred[segments[[n.seg]],,,] <- res$cvPred
      if (jackknife) cvCoef[,,,,n.seg] <- res$cvCoef
    }
    
    ## Calculate validation statistics:
    PRESS0 <- apply(Y, 2, var) * nobj^2 / (nobj - 1) # FIXME: Only correct for loocv!
    PRESS <- apply((cvPred - c(Y))^2, 2:4, sum)

    ## Add dimnames:
    shrinknames <- paste("Shrink", shrink)
    objnames <- dnX[[1]]
    if (is.null(objnames)) objnames <- dnY[[1]]
    respnames <- dnY[[2]]
    nCompnames <- paste(1:ncomp, "comps")
    names(PRESS0) <- respnames
    dimnames(adj) <- dimnames(PRESS) <-
      list(respnames, nCompnames, shrinknames)
    dimnames(cvPred) <- list(objnames, respnames, nCompnames, shrinknames)
    if (jackknife)
      dimnames(cvCoef) <- list(dnX[[2]], respnames, nCompnames, shrinknames,
                               paste("Seg", seq_along(segments)))
  } else {
    adj <- matrix(0, nrow = nresp, ncol = ncomp)
    cvPred <- array(0, dim = c(nobj, nresp, ncomp))
    if (jackknife)
      cvCoef <- array(dim = c(npred, nresp, ncomp, length(segments)))
    if (method == "truncation") gammas <- list()
    
    ## Collect the results:
    for (n.seg in seq_along(segments)) {
      res <- results[[n.seg]]
      adj <- adj + res$adj
      cvPred[segments[[n.seg]],,] <- res$cvPred
      if (jackknife) cvCoef[,,,n.seg] <- res$cvCoef
      if (method == "truncation") gammas[[n.seg]] <- res$gammas
    }
    
    ## Calculate validation statistics:
    PRESS0 <- apply(Y, 2, var) * nobj^2 / (nobj - 1) # FIXME: Only correct for loocv!
    PRESS <- colSums((cvPred - c(Y))^2)
    
    ## Add dimnames:
    objnames <- dnX[[1]]
    if (is.null(objnames)) objnames <- dnY[[1]]
    respnames <- dnY[[2]]
    nCompnames <- paste(1:ncomp, "comps")
    names(PRESS0) <- respnames
    dimnames(adj) <- dimnames(PRESS) <-
      list(respnames, nCompnames)
    dimnames(cvPred) <- list(objnames, respnames, nCompnames)
    if (jackknife)
      dimnames(cvCoef) <- list(dnX[[2]], respnames, nCompnames,
                               paste("Seg", seq_along(segments)))
  }
  
  list(method = "CV", pred = cvPred, coefficients = if (jackknife) cvCoef,
       gammas = if (method == "cppls") gammas,
       PRESS0 = PRESS0, PRESS = PRESS, adj = adj / nobj^2,
       segments = segments, ncomp = ncomp)
}

## delete.intercept: utilitiy function that deletes the response coloumn from
## a model matrix, and adjusts the "assign" attribute:
delete.intercept <- function(mm) {
  ## Save the attributes prior to removing the intercept coloumn:
  saveattr <- attributes(mm)
  ## Find the intercept coloumn:
  intercept <- which(saveattr$assign == 0)
  ## Return if there was no intercept coloumn:
  if (!length(intercept)) return(mm)
  ## Remove the intercept coloumn:
  mm <- mm[,-intercept, drop=FALSE]
  ## Update the attributes with the new dimensions:
  saveattr$dim <- dim(mm)
  saveattr$dimnames <- dimnames(mm)
  ## Remove the assignment of the intercept from the attributes:
  saveattr$assign <- saveattr$assign[-intercept]
  ## Restore the (modified) attributes:
  attributes(mm) <- saveattr
  ## Return the model matrix:
  mm
}

print.mvrV <- function(x, ...) {
  switch(x$method,
         truncation = {
           ana = "Partial least squares regression"
           alg = "cppls"
         },
         stpls = {
           ana = "Partial least squares regression"
           alg = "kernel"
         },
         # ridgepls = {
         #   ana = "Partial least squares regression"
         #   alg = "ridge pls"
         # },
         stop("Unknown fit method.")
  )
  cat(ana, ", fitted with the", alg, "algorithm.")
  if (!is.null(x$validation))
    cat("\nCross-validated using", length(x$validation$segments),
        attr(x$validation$segments, "type"), "segments.")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

## The explained X variance:
.explvarV <- function(object)
  switch(class(object)[1],
         mvrV = 100 * object$Xvar / object$Xtotvar,
         mvr = 100 * object$Xvar / object$Xtotvar,
         princomp =,
         prcomp = 100 * object$sdev^2 / sum(object$sdev^2),
         scores =,
         loadings = attr(object, "explvar")
  )

## Summary method for mvrV objects
#' @title Summary method for stpls and trunc
#' @description Adaptation of \code{summary.mvr} from the \code{pls} package v 2.4.3.
#' 
#' @param object an mvrV object
#' @param what one of "all", "validation" or "training"
#' @param digits integer. Minimum number of significant digits in the output. Default is 4.
#' @param print.gap	Integer. Gap between coloumns of the printed tables.
#' @param ...	Other arguments sent to underlying methods.

#' @seealso \code{\link[pls]{summary.mvr}}
#' @export
summary.mvrV <- function(object, what = c("all", "validation", "training"),
                        digits = 4, print.gap = 2, ...)
{
  what <- match.arg(what)
  if (what == "all") what <- c("validation", "training")
  if (is.null(object$validation)) what <- "training"
  
  nobj <- nrow(object$scores)
  nresp <- length(object$Ymeans)
  yvarnames <- respnames(object)
  cat("Data: \tX dimension:", nobj, length(object$Xmeans),
      "\n\tY dimension:", nobj, nresp)
  cat("\nFit method:", object$method)
  cat("\nNumber of components considered:", object$ncomp)
  
  for (wh in what) {
    if (wh == "training") {
      cat("\nTRAINING: % variance explained\n")
      xve <- .explvarV(object)
      r2  <- R2(object, estimate = "train",
                intercept = FALSE)$val
      yve <- 100 * drop(r2)
      if(!is.null(object$shrink)){
        tbl <- t(cbind(apply(xve,2,cumsum), yve))
        dimnames(tbl) <- list(c(paste(colnames(xve),"X",sep=", "), paste(dimnames(r2)[[4]],yvarnames,sep=", ")),
                              paste(1:object$ncomp, "comps"))
      } else {
        tbl <- rbind(cumsum(xve), yve)
        dimnames(tbl) <- list(c("X", yvarnames),
                              paste(1:object$ncomp, "comps"))
      }
      print(tbl, digits = digits, print.gap = print.gap, ...)
    } else {
      cat("\n\nVALIDATION: RMSEP")
      cat("\nCross-validated using", length(object$validation$segments),
          attr(object$validation$segments, "type"), "segments.\n")
      print(RMSEP(object), digits = digits, print.gap = print.gap, ...)
    }
  }
}

## Internal function to apply FUN over X, optionally in parallel:
lapplyFunc <- function(parSpec, X, FUN, nonForkInit) {
  if (is.null(parSpec) || (is.numeric(parSpec) && parSpec == 1)) {
    ## Serially
    results <- lapply(X, FUN)
  } else {
    ## Parallel
#    requireNamespace(parallel, quietly = TRUE)
    stop_cluster <- FALSE           # Whether to kill the workers afterwards
    
    if (is.numeric(parSpec) && parSpec > 1) {
      ## Number => number of workers with mclapply
      results <- mclapply(X, FUN, mc.cores = parSpec)
    } else {
      if (is.call(parSpec)) {
        ## Unevaluated call => evaluate it to create the cluster:
        parSpec <- eval(parSpec)
        stop_cluster <- TRUE
      }
      
      if (inherits(parSpec, "cluster")) {
        ## Run library(pls) on cluster if type != FORK
        if (!inherits(parSpec[[1]], "forknode")
            && !missing(nonForkInit)) {
          eval(nonForkInit)
        }
        results <- parLapply(parSpec, X, FUN)
        
        if (stop_cluster) {
          stopCluster(parSpec)
        }
      } else {
        stop("Unknown parallelity specification: '", parSpec, "'")
      }
    }
  }
  
  return(results)
}

## Calculate the validation statistics needed for (R)MSEP and R^2.
## Note that it accepts any values for `estimate', but only calculates
## statistics for "train", "test" and "CV".
mvrValstats <- function(object, estimate,
                        newdata, ncomp = 1:object$ncomp, comps,
                        intercept = cumulative, se = FALSE, ...)
{
  if(is.null(object$shrink)){
    shrink <- FALSE
  } else {
    shrink <- object$shrink
  }
  ## Makes the code slightly simpler:
  cumulative <- missing(comps) || is.null(comps)
  
  if (any(estimate == "CV")) {
    ## Check that cross-validation is possible:
    if (!cumulative)
      stop("Cross-validation is not supported when `comps' is specified")
    if (is.null(object$validation))
      stop("`object' has no `validation' component")
  }
  
  ## The calculated stuff:
  nestimates <- length(estimate)
  nresp <- dim(fitted(object))[2]
  respnames <- dimnames(fitted(object))[[2]]
  if(!is.logical(shrink)){ # With shrinkage
    SSE <- array(dim = c(nestimates, nresp,
                         if(cumulative) 1 + length(ncomp) else 2, length(shrink)),
                 dimnames = list(estimate = estimate,
                                 response = respnames,
                                 model = if (cumulative) {
                                   c("(Intercept)", paste(ncomp, "comps"))
                                 } else {
                                   c("(Intercept)", paste("(Intercept), Comp",
                                                          paste(comps, collapse = ", ")))
                                 },
                                 shrinkage = paste("Shrink",shrink)
                 ))
  } else {
    SSE <- array(dim = c(nestimates, nresp,
                         if(cumulative) 1 + length(ncomp) else 2),
                 dimnames = list(estimate = estimate,
                                 response = respnames,
                                 model = if (cumulative) {
                                   c("(Intercept)", paste(ncomp, "comps"))
                                 } else {
                                   c("(Intercept)", paste("(Intercept), Comp",
                                                          paste(comps, collapse = ", ")))
                                 }
                 ))
  }
  SST <- array(dim = c(nestimates, nresp),
               dimnames = list(estimate = estimate, response = respnames))
  nobj <- numeric(nestimates)
  names(nobj) <- estimate
  
  ## Calculate the statistics:
  for (i in seq(along = estimate)) {
    switch(estimate[i],
           train = {
             resp <- as.matrix(model.response(model.frame(object)))
             nobj[i] <- nrow(resp)
             if (inherits(object$na.action, "exclude")) {
               resp <- napredict(object$na.action, resp) # inserts NAs
             }
             if(!is.logical(shrink)){
               res <- if (cumulative)
                 residuals(object, ...)[,,ncomp,, drop=FALSE]
               else
                 resp - predict(object, comps = comps, ...)
               SST[i,] <- apply(resp, 2, var, na.rm = TRUE) *
                 (nobj[i] - 1)
               SSE[i,,1,] <- SST[i,]
               SSE[i,,-1,] <- apply(res^2,2:4,sum, na.rm = TRUE)
             } else {
               res <- if (cumulative)
                 residuals(object, ...)[,,ncomp, drop=FALSE]
               else
                 resp - predict(object, comps = comps, ...)
               SST[i,] <- apply(resp, 2, var, na.rm = TRUE) *
                 (nobj[i] - 1)
               SSE[i,,] <- cbind(SST[i,], colSums(res^2, na.rm = TRUE))
             }
           },
           test = {
             if (missing(newdata)) stop("Missing `newdata'.")
             ## Remove any observations with NAs:
             newdata <- model.frame(formula(object), data = newdata)
             resp <- as.matrix(model.response(newdata))
             pred <- if (cumulative)
               predict(object, ncomp = ncomp, newdata = newdata,...)
             else
               predict(object, comps = comps, newdata = newdata,...)
             nobj[i] <- nrow(newdata)
             SST[i,] <- apply(resp, 2, var) * (nobj[i] - 1)
             if(!is.logical(shrink)){
               SSE[i,,1,]  <- apply(sweep(resp, 2, object$Ymeans)^2, 2:4, sum) # FIXME: ikke testet
               SSE[i,,-1,] <- apply((pred - c(resp))^2, 2:4, sum)
             } else {
               SSE[i,,] <- cbind(colSums(sweep(resp, 2, object$Ymeans)^2),
                                 colSums((pred - c(resp))^2))
             }
           },
           CV = {
             resp <- as.matrix(model.response(model.frame(object)))
             nobj[i] <- nrow(resp)
             SST[i,] <- apply(resp, 2, var) * (nobj[i] - 1)
             if(!is.logical(shrink)){
               SSE[i,,1,]  <- object$validation$PRESS0
               SSE[i,,-1,] <- object$validation$PRESS[,ncomp,, drop=FALSE]
             } else {
               SSE[i,,] <-
                 cbind(object$validation$PRESS0,
                       object$validation$PRESS[,ncomp, drop=FALSE])
             }
           }
    )
  }
  
  if (cumulative) comps <- ncomp
  ## Either remove the intercept or add a "zeroth" component:
  if (intercept)
    comps <- c(0, comps)
  else
    if(!is.logical(shrink)){
      SSE <- SSE[,,-1,, drop=FALSE]
    } else {
      SSE <- SSE[,,-1, drop=FALSE]
    } 
  return(list(SSE = SSE, SST = SST, nobj = nobj, comps = comps,
              cumulative = cumulative))
}


## R2: Return R^2
R2 <- function(object, ...) UseMethod("R2")
R2.mvr <- function(object, estimate, newdata, ncomp = 1:object$ncomp, comps,
                   intercept = cumulative, se = FALSE, ...) {
  ## Makes the code slightly simpler:  FIXME: maybe remove
  cumulative <- missing(comps) || is.null(comps)
  
  ## Figure out which estimate(s) to calculate:
  allEstimates <- c("all", "train", "CV", "test")
  if (missing(estimate)) {
    ## Select the `best' available estimate
    if (!missing(newdata)) {
      estimate = "test"
    } else {
      if (!is.null(object$validation)) {
        estimate = "CV"
      } else {
        estimate = "train"
      }
    }
  } else {
    estimate <- allEstimates[pmatch(estimate, allEstimates)]
    if (any(is.na(estimate))) stop("`estimate' should be a subset of ",
                                   paste(allEstimates, collapse = ", "))
    if (any(estimate == "all")) {
      estimate <- allEstimates[-1] # Try all estimates (except "all")
      if(missing(newdata)) estimate <- setdiff(estimate, "test")
      if(is.null(object$validation) || !cumulative)
        estimate <- setdiff(estimate, "CV")
    }
  }
  
  ## Get the needed validation statistics:
  cl <- match.call(expand.dots = FALSE)
  cl$estimate <- estimate             # update estimate argument
  cl[[1]] <- as.name("mvrValstats")
  valstats <- eval(cl, parent.frame())
  
  ## Calculate the R^2s:
  R2 <- 1 - valstats$SSE / c(valstats$SST)
  
  return(structure(list(val = R2, type = "R2", comps = valstats$comps,
                        cumulative = valstats$cumulative, call = match.call()),
                   class = "mvrVal"))
}


## MSEP: Return MSEP
MSEP <- function(object, ...) UseMethod("MSEP")
MSEP.mvr <- function(object, estimate, newdata, ncomp = 1:object$ncomp, comps,
                     intercept = cumulative, se = FALSE, ...)
{
  ## Makes the code slightly simpler:
  cumulative <- missing(comps) || is.null(comps)
  
  ## Figure out which estimate(s) to calculate:
  allEstimates <- c("all", "train", "CV", "adjCV", "test")
  if (missing(estimate)) {
    ## Select the `best' available estimate
    if (!missing(newdata)) {
      estimate = "test"
    } else {
      if (!is.null(object$validation)) {
        estimate = c("CV", "adjCV")
      } else {
        estimate = "train"
      }
    }
  } else {
    estimate <- allEstimates[pmatch(estimate, allEstimates)]
    if (any(is.na(estimate))) stop("`estimate' should be a subset of ",
                                   paste(allEstimates, collapse = ", "))
    if (any(estimate == "all")) {
      estimate <- allEstimates[-1] # Try all estimates (except "all")
      if(missing(newdata)) estimate <- setdiff(estimate, "test")
      if(is.null(object$validation) || !cumulative)
        estimate <- setdiff(estimate, c("CV", "adjCV"))
    }
  }
  
  ## adjCV needs the statistics for CV and train, so we optionally
  ## have to add them:
  if (adjCV <- any(estimate == "adjCV")) {
    ## Note: this removes any duplicate elements
    calcestimates <- union(estimate, c("train", "CV"))
  } else {
    calcestimates <- estimate
  }
  ## Get the needed validation statistics:
  cl <- match.call(expand.dots = FALSE)
  cl$estimate <- calcestimates        # update estimate argument
  cl[[1]] <- as.name("mvrValstats")
  valstats <- eval(cl, parent.frame())
  
  ## Calculate the MSEPs:
  MSEP <- valstats$SSE / valstats$nobj
  if (adjCV) {
    if(!is.null(object$shrink)){
      ## Calculate the adjusted CV
      MSEP["adjCV",,,] <- MSEP["CV",,,]
      if (intercept) {
        MSEP["adjCV",,-1,] <- MSEP["adjCV",,-1,] + MSEP["train",,-1,] -
          object$validation$adj[,ncomp,]
      } else {
        MSEP["adjCV",,,] <- MSEP["adjCV",,,] + MSEP["train",,,] -
          object$validation$adj[,ncomp,]
      }
    } else {
      ## Calculate the adjusted CV
      MSEP["adjCV",,] <- MSEP["CV",,]
      if (intercept) {
        MSEP["adjCV",,-1] <- MSEP["adjCV",,-1] + MSEP["train",,-1] -
          object$validation$adj[,ncomp]
      } else {
        MSEP["adjCV",,] <- MSEP["adjCV",,] + MSEP["train",,] -
          object$validation$adj[,ncomp]
      }
    }
    ## Remove any specially added estimates (this also adds back any
    ## duplicate elements):
    if(!is.null(object$shrink)){
      MSEP <- MSEP[estimate,,,, drop=FALSE]
    } else {
      MSEP <- MSEP[estimate,,, drop=FALSE]
    }
  }
  
  return(structure(list(val = MSEP, type = "MSEP", comps = valstats$comps,
                        cumulative = valstats$cumulative, call = match.call()),
                   class = "mvrVal"))
}

# RMSEP: A wrapper around MSEP to calculate RMSEPs
RMSEP <- function(object, ...) UseMethod("RMSEP")
RMSEP.mvr <- function(object, ...) {
  cl <- match.call()
  cl[[1]] <- as.name("MSEP")
  z <- eval(cl, parent.frame())
  z$val <- sqrt(z$val)
  z$type <- "RMSEP"
  z$call[[1]] <- as.name("RMSEP")
  z
}

## Print method for mvrVal objects:
print.mvrVal <- function(x, digits = 4, print.gap = 2, ...) {
  nresp <- dim(x$val)[2]
  yvarnames <- dimnames(x$val)[[2]]
  names(dimnames(x$val)) <- NULL
  for (i in 1:nresp) {
    if (nresp > 1) cat("\nResponse:", yvarnames[i], "\n")
    if(length(dim(x$val)) == 4){
      print(x$val[,i,,], digits = digits, print.gap = print.gap, ...)
    } else {
      print(x$val[,i,], digits = digits, print.gap = print.gap, ...)
    }
  }
  invisible(x)
}
