### optimisation.R:  Classes, methods and functions for searching for
### optimal baselines.
### $Id: optimisation.R 179 2011-01-09 14:29:07Z bhm $

## The function mvrValstats has been copied from package pls (overlapping authors)
## because eval(parent.frame()) seems to break call to pls::msep.


###
### Common generic functions and default methods
###


## A function for performing the testing.
## It will have methods for the predictionTest subclasses, and the
## baselineAlgTest class.
setGeneric("runTest", function(object, X, y, ...) standardGeneric("runTest"))

## A function for extracting the quality measures from objects.
## It will have methods for predictionResult and baselineAlgResult.
setGeneric("qualMeas", function(object, ...) standardGeneric("qualMeas"))

## A function to extract parameter settings
## It will have methods for predictionResult, baselineAlgTest,
## baselineAlgResult and baselineAlg.
setGeneric("param", function(object) standardGeneric("param"))


###
### Classes and methods for prediction testing
###


## A virtual class for all prediction tests.  It will have subclasses for each
## type of predictor (PLSR, RR, ...).
setClass("predictionTest")

## A class for the prediction test results
setClass("predictionResult",
         representation(param = "numeric",
                        qualMeas = "numeric",
                        ind.min = "numeric", ##?
                        minQualMeas = "numeric", ##?
                        param.min = "numeric", ##?
                        qualMeasName = "character",
                        paramName = "character"
         ))

## Extractor generics and methods:
setMethod("param", "predictionResult", function(object) object@param)
setMethod("qualMeas", "predictionResult", function(object, ...) object@qualMeas)
setGeneric("ind.min", function(object) standardGeneric("ind.min"))
setMethod("ind.min", "predictionResult", function(object) object@ind.min)
setGeneric("minQualMeas", function(object) standardGeneric("minQualMeas"))
setMethod("minQualMeas", "predictionResult", function(object) object@minQualMeas)
setGeneric("param.min", function(object) standardGeneric("param.min"))
setMethod("param.min", "predictionResult", function(object) object@param.min)
setGeneric("qualMeasName", function(object) standardGeneric("qualMeasName"))
setMethod("qualMeasName", "predictionResult", function(object) object@qualMeasName)
setGeneric("paramName", function(object) standardGeneric("paramName"))
setMethod("paramName", "predictionResult", function(object) object@paramName)

##
## PLSR
##

## The class
setClass("PLSRTest", contains = "predictionTest",
         representation(ncomp = "numeric", cvsegments = "list"))

## The prediction test method
setMethod("runTest", "PLSRTest",
          function(object, X, y) {
            if(requireNamespace("pls", quietly = TRUE)){
              ncomp <- object@ncomp
              cvsegments <- object@cvsegments
              ## FIXME: Perhaps use mvrCv directly
              res <- pls::plsr(y ~ X, ncomp = ncomp, validation = "CV",
                               segments = cvsegments)
              msep <- drop(pls::MSEP(res, estimate = "adjCV")$val)
              if (NCOL(y) == 1) {
                rmsep <- sqrt(msep)
              } else {
                ## For multiple responses, use the square root of the
                ## mean (over the responses) relative (to 0 components)
                ## MSEP:
                rmsep <- sqrt(colMeans(msep / msep[,1]))
              }
              ind.min <- which.min(rmsep)
              return(new("predictionResult",
                         param = 0:ncomp,
                         qualMeas = rmsep, ind.min = ind.min, minQualMeas = min(rmsep),
                         param.min = ind.min - 1,
                         paramName = "ncomp", qualMeasName = "RMSEP"))
            } else {
              warning('Package pls not installed')
              return(list())
            }
          })

##
## Ridge regression
##

## The class
setClass("ridgeRegressionTest", contains = "predictionTest",
         representation(lambda = "numeric"))

## The prediction test method
setMethod("runTest", "ridgeRegressionTest",
          function(object, X, y) {
            if(requireNamespace("MASS", quietly = TRUE)){
              lambda <- object@lambda
              res <- MASS::lm.ridge(y ~ X, lambda = lambda)
              ind.min <- which.min(res$GCV)
              return(new("predictionResult",
                         param = lambda,
                         qualMeas = res$GCV, ind.min = ind.min, minQualMeas = min(res$GCV),
                         param.min = res$lambda[ind.min],
                         paramName = "lambda", qualMeasName = "GCV"))
            } else {
              warning('Package MASS not installed')
              return(list())
            }
          })



###
### Classes and methods for describing baseline algorithms
### FIXME: Put in separate file!


setClass("baselineAlg",
         representation(name = "character",
                        description = "character",
                        funcName = "character",
                        param = "data.frame"
         ),
         prototype(param = data.frame(
           name = NA, integer = NA, min = NA, incl.min = NA,
           default = NA, max = NA, incl.max = NA)[0,]
         ),
         validity = function(object) {
           if (!identical(names(object@param),
                          c("name","integer","min","incl.min","default","max","incl.max")))
             return("The param slot does not have the correct coloumn names")
           return(TRUE)
         }
)

## Accessor functions:
setGeneric("name", function(object) standardGeneric("name"))
setMethod("name", "baselineAlg", function(object) object@name)
setGeneric("description", function(object) standardGeneric("description"))
setMethod("description", "baselineAlg", function(object) object@description)
setGeneric("funcName", function(object) standardGeneric("funcName"))
setMethod("funcName", "baselineAlg", function(object) object@funcName)
setMethod("param", "baselineAlg", function(object) object@param)


## A list with objects for each implemented algorithm:
baselineAlgorithms <- list(
  
  als = new("baselineAlg",
            name = "als",
            description = "Asymmetric Least Squares",
            funcName = "baseline.als",
            param = data.frame(  # FIXME
              name = c("lambda","p"), # maxit
              integer = c(FALSE, FALSE),
              min = c(0, 0),
              incl.min = c(TRUE, TRUE),
              default = c(6, 0.05),
              max = c(Inf, 1),
              incl.max = c(FALSE, TRUE)
            )),
  fillPeaks = new("baselineAlg",
                  name = "fillPeaks",
                  description = "Iterateive baseline correction algorithm based on mean suppression",
                  funcName = "baseline.fillPeaks",
                  param = data.frame(  # FIXME
                    name = c("lambda","hwi","it", "int"),
                    integer = c(FALSE, TRUE, TRUE, TRUE),
                    min = c(0, 1, 1, 3),
                    incl.min = c(TRUE, TRUE, TRUE, TRUE),
                    default = c(NA, NA, NA, NA),
                    max = c(Inf, Inf, Inf, Inf),
                    incl.max = c(FALSE, FALSE, FALSE, FALSE)
                  )),
  irls = new("baselineAlg",
             name = "irls",
             description = "Iterative Restricted Least Squares",
             funcName = "baseline.irls",
             param = data.frame( # FIXME
               name = c("lambda1", "lambda2", "maxit", "wi"),
               integer = c(FALSE, FALSE, TRUE, TRUE),
               min = c(0, 0, 0, 0),
               incl.min = c(TRUE, TRUE, TRUE, TRUE),
               default = c(5, 9, 200, 0.05),
               max = c(Inf, Inf, Inf, 1),
               incl.max = c(FALSE, FALSE, FALSE, TRUE)
             )),
  ## FIXME: Not tested, because baseline.lowpass gives strange results!
  lowpass = new("baselineAlg",
                name = "lowpass",
                description = "Low-pass filter based on fast Fourier transform",
                funcName = "baseline.lowpass",
                param = data.frame(  # FIXME
                  name = c("steep","half"),
                  integer = c(FALSE, TRUE),
                  min = c(0, 1),
                  incl.min = c(TRUE, TRUE),
                  default = c(2, 5),
                  max = c(Inf, Inf),
                  incl.max = c(FALSE, FALSE)
                )),
  medianWindow = new("baselineAlg",
                     name = "medianWindow",
                     description = "Local medians",
                     funcName = "baseline.medianWindow",
                     param = data.frame(  # FIXME
                       name = c("hwm","hws"),         # end
                       integer = c(TRUE, TRUE),
                       min = c(0,1),
                       incl.min = c(TRUE, TRUE),
                       default = c(NA, NA),
                       max = c(Inf, Inf),
                       incl.max = c(FALSE, FALSE)
                     )),
  modpolyfit = new("baselineAlg",
                   name = "modpolyfit",
                   description = "Modified iterative polynomial fitting",
                   funcName = "baseline.modpolyfit",
                   param = data.frame(  # FIXME
                     name = c("degree", "tol", "rep"), # t
                     integer = c(TRUE, FALSE, TRUE),
                     min = c(1, 0, 0),
                     incl.min = c(TRUE, TRUE, TRUE),
                     default = c(4, 1e-3, 100),
                     max = c(Inf, Inf, Inf),
                     incl.max = c(FALSE, FALSE, FALSE)
                   )),
  ## FIXME: Transformed parameters??
  peakDetection = new("baselineAlg",
                      name = "peakDetection",
                      description = "Peak detection",
                      funcName = "baseline.peakDetection",
                      param = data.frame(  # FIXME
                        name = c("left","right","lwin", "rwin", "snminimum", "mono", "multiplier"),
                        integer = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE),
                        min = c(1,1,1,1,0,0,0),
                        incl.min = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
                        default = c(NA),
                        max = c(Inf,Inf,Inf,Inf,Inf,1,Inf),
                        incl.max = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
                      )),
  ## FIXME:  Perhaps easier to specify span than NoXP??
  rfbaseline = new("baselineAlg",
                   name = "rfbaseline",
                   description = "Robust Baseline Estimation by Ruckstuhl etal",
                   funcName = "baseline.rfbaseline",
                   param = data.frame(  # FIXME
                     name = c("NoXP","b"),        # Lots of arguments!
                     integer = c(TRUE, FALSE),
                     min = c(3, 0),
                     incl.min = c(TRUE, FALSE),
                     default = c(NA, 3.5),
                     max = c(Inf, Inf),
                     incl.max = c(FALSE, FALSE)
                   )),
  rollingBall = new("baselineAlg",
                    name = "rollingBall",
                    description = "Rolling Ball",
                    funcName = "baseline.rollingBall",
                    param = data.frame(  # FIXME
                      name = c("wm","ws"),
                      integer = c(TRUE, TRUE),
                      min = c(2, 2),
                      incl.min = c(TRUE, TRUE),
                      default = c(NA, NA),
                      max = c(Inf, Inf),
                      incl.max = c(FALSE, FALSE)
                    )),
  shirley = new("baselineAlg",
					name = "shirley",
					description = "A shirley baseline correction algorithm",
					funcName = "baseline.shirley",
					param = data.frame(
					  name = c("maxit", "err"), # maxit
					  integer = c(FALSE, FALSE),
					  min = c(1, 1e-8),
					  incl.min = c(TRUE, TRUE),
					  default = c(50, 1e-6),
					  max = c(Inf, 1e-3),
					  incl.max = c(FALSE, FALSE)
					))
  
) ## End of baselineAlgorithms <- list(



###
### Classes and methods for baseline algorithm tests
###


## One class for all algorithms
setClass("baselineAlgTest",
         representation(algorithm = "baselineAlg",
                        param = "list",
                        extraArgs = "list"),
         prototype(extraArgs = list()))

## A class for the baseline algorithm test results:
setClass("baselineAlgResult",
         representation(param = "list",
                        qualMeas = "matrix", ## or multidimensional array?
                        qualMeas.ind.min = "numeric",
                        minQualMeas = "numeric", ##?
                        param.ind.min = "numeric",
                        param.min = "list",
                        qualMeasName = "character"
         ))

## Accessor functions:
setGeneric("algorithm", function(object) standardGeneric("algorithm"))
setMethod("algorithm", "baselineAlgTest", function(object) object@algorithm)
setMethod("param", "baselineAlgTest", function(object) object@param)
setGeneric("extraArgs", function(object) standardGeneric("extraArgs"))
setMethod("extraArgs", "baselineAlgTest", function(object) object@extraArgs)
## An extractor method for extracting the func name from a baselineAlgTest
setMethod("funcName", "baselineAlgTest",
          function(object) funcName(object@algorithm))
setMethod("param", "baselineAlgResult", function(object) object@param)
setGeneric("qualMeas.ind.min", function(object) standardGeneric("qualMeas.ind.min"))
setMethod("qualMeas.ind.min", "baselineAlgResult", function(object) object@qualMeas.ind.min)
setMethod("minQualMeas", "baselineAlgResult", function(object) object@minQualMeas)
setGeneric("param.ind.min", function(object) standardGeneric("param.ind.min"))
setMethod("param.ind.min", "baselineAlgResult", function(object) object@param.ind.min)
setMethod("param.min", "baselineAlgResult", function(object) object@param.min)
setMethod("qualMeasName", "baselineAlgResult", function(object) object@qualMeasName)

## A method for doing the algorithm testing
setMethod("runTest", "baselineAlgTest",
          function(object, X, y, predictionTest, postproc, verbose = FALSE) {
            ## Build a data frame with all parameter combinations:
            params <- expand.grid(object@param, KEEP.OUT.ATTRS = FALSE)
            nparams <- nrow(params)
            ## Variables to accumulate results:
            qualMeas <- list()
            param.min <- minQualMeas <- regpar.ind.min <- numeric(nparams)
            if (verbose) cat("Looping through the", nparams,
                             "baseline parameter settings\n")
            for (i in seq_len(nparams)) {
              if (verbose) cat(i, ": correcting... ", sep = "")
              ## Do the baseline correction
              bl <- do.call("baseline",
                            c(list(spectra = X, method =  name(object@algorithm)),
                              params[i,, drop=FALSE], object@extraArgs))
              Xcorr <- getCorrected(bl)
              ## Perform any post processing:
              if (!missing(postproc) && !is.null(postproc)) {
                if (verbose) cat("postprocessing... ")
                Xcorr <- postproc(Xcorr)
              }
              ## Test the predictor on the baseline corrected data
              if (verbose) cat("prediction testing...")
              res <- runTest(predictionTest, Xcorr, y)
              qualMeas[[i]] <- qualMeas(res)
              regpar.ind.min[i] <- ind.min(res)
              minQualMeas[i] <- minQualMeas(res)
              param.min[i] <- param.min(res)
              if (verbose) cat("\n")
            }
            qualMeas <- do.call("cbind", qualMeas)
            ## Find and return the best results
            best.i <- which.min(minQualMeas)
            ## Join the regression parameter and the baseline parameters:
            param <- c(list(param(res)), object@param)
            names(param)[1] <- paramName(res)
            param.ind.min <-
              c(regpar.ind.min[best.i],
                unlist(expand.grid(lapply(object@param, seq_along))[best.i,, drop=FALSE]))
            names(param.ind.min)[1] <- paramName(res)
            param.min <- c(param.min[best.i], as.list(params[best.i,, drop=FALSE]))
            names(param.min)[1] <- paramName(res)
            return(new("baselineAlgResult",
                       param = param,
                       qualMeas = qualMeas,
                       qualMeas.ind.min = c(regpar.ind.min[best.i], best.i),
                       minQualMeas = minQualMeas[best.i],
                       param.ind.min = param.ind.min,
                       param.min = param.min,
                       qualMeasName = qualMeasName(res)
            ))
          })


###
### Method for extracting test result values from a baseline algorithm
### test
###
setMethod("qualMeas", "baselineAlgResult",
          function(object, ..., MIN, AVG, DEFAULT = c("all", "cond.min", "overall.min", "avg")) {
            DEFAULT <- match.arg(DEFAULT)
            ## Extract needed information from object:
            res <- object@qualMeas
            params <- object@param
            param.ind.min <- object@param.ind.min
            
            ## Collect the named parameters:
            setparams <- list(...)
            
            ## 0) Substitute for any setparams == "overall"
            atoverall <- names(setparams)[setparams == "overall"]
            if (length(atoverall) > 0) {
              setparams[atoverall] <- param.ind.min[atoverall]
            }
            ## 0b) Substitute for any setparams == "all"
            useall <- names(setparams)[setparams == "all"]
            if (length(useall) > 0) {
              setparams[useall] <- lapply(params[useall], seq_along)
            }
            
            ## Figure out which parameters to take (conditional) min over:
            MINns <- AVGns <- character()
            if (DEFAULT == "cond.min") {
              ## Take min over all remaining parameters
              if(missing(AVG)) {
                MINns <- setdiff(names(params), names(setparams))
              } else {
                MINns <- setdiff(names(params), c(names(setparams), AVG))
              }
              ## Figure out which parameters to take average over:
            } else {
              if (DEFAULT == "avg") {
                ## Take average over all remaining parameters
                if(missing(MIN)) {
                  AVGns <- setdiff(names(params), names(setparams))
                } else {
                  AVGns <- setdiff(names(params), c(names(setparams), MIN)) }
              }
            }
            ## Figure out which parameters to take (conditional) min over:
            if (!missing(MIN) && length(MINns)==0) {
              ## Take min over parameters specified in MIN
              MIN <- substitute(MIN)
              wasList <- is.call(MIN)
              MINns <- sapply(MIN, as.character) # Convert to char, if
              # needed; also make it a vector
              ## If MIN was a list, the first element will be "list":
              if (wasList) MINns <- MINns[-1]
              ## Sanity check:
              if (!all(MINns %in% names(params)))
                stop("Non-existing parameter(s) specified in 'MIN'.")
            }
            ## Figure out which parameters to take average over:
            if (!missing(AVG) && length(AVGns)==0) {
              ## Take average over parameters specified in AVG
              AVG <- substitute(AVG)
              wasList <- is.call(AVG)
              AVGns <- sapply(AVG, as.character) # Convert to char, if
              # needed; also make it a vector
              ## If AVG was a list, the first element will be "list":
              if (wasList) AVGns <- AVGns[-1]
              ## Sanity check:
              if (!all(AVGns %in% names(params)))
                stop("Non-existing parameter(s) specified in 'AVG'.")
            }
            
            ## If DEFAULT is "overall.min", set any remaining parameters to their
            ## value corresponding to the overall min:
            if (DEFAULT == "overall.min") {
              nms <- setdiff(names(params), c(names(setparams), MINns, AVGns))
              for (nm in nms) {
                setparams[[nm]] <- param.ind.min[nm]
              }
            }
            
            ## 1) Convert the results into an array
            newdim <- sapply(params, length)
            dim(res) <- newdim
            dimnames(res) <- params
            
            ## 2) Reorder the results into the order given in params
            
            ## Internal function for applying without dropping dimensions
            kapply <- function(x, MARGIN, FUN, lab){
              dims <- dim(x)
              command <- "x["
              if(MARGIN==1){
                command <- "x[1"
                dimnames(x)[[1]][1] <- lab
              } else {
                command <- "x["
              }
              for(i in 2:length(dims)){
                if(MARGIN==i){
                  command <- paste(command,", 1", sep="")
                  dimnames(x)[[i]][1] <- lab
                } else {
                  command <- paste(command,", ", sep="")
                }
              }
              command <- paste(command, ", drop=FALSE]", sep="")
              y <- eval(parse(text=command))
              y[] <- apply(x, setdiff(1:length(dims),MARGIN), FUN)
              y
            }
            
            ## Internal function for subsetting without collapsing
            ksubset <- function(x, setparams){
              dn <- names(dimnames(x))
              sn <- names(setparams)
              dims <- dim(x)
              command <- "x["
              if(sn[1]%in%dn){
                command <- paste("x[setparams$", sn[1], sep="")
              } else {
                command <- "x["
              }
              for(i in 2:length(dims)){
                if(sn[i]%in%dn){
                  command <- paste(command, ", setparams$", sn[i], sep="")
                } else {
                  command <- paste(command,", ", sep="")
                }
              }
              command <- paste(command, ", drop=FALSE]", sep="")
              eval(parse(text=command))
            }
            
            ## 3) Extract the desired sub-array:
            if(length(setparams)>0)
              res <- ksubset(res, setparams)
            res.names <- names(dimnames(res))
            
            ## 4) Take the min over any parameters in MINns
            if( length(MINns) > 0){
              for( i in 1:length(MINns)){
                res <- kapply(res, which(res.names==MINns[i]), min, 'min')
              }
            }
            
            ## 5) Take the mean over any parameters in AVGns
            if( length(AVGns) > 0){
              for( i in 1:length(AVGns)){
                res <- kapply(res, which(res.names==AVGns[i]), mean, 'avg')
              }
            }
            
            return(res)
          })

###
### Function for testing several baseline algorithms
###

doOptim <- function(baselineTests, X, y, predictionTest,
                    postproc = NULL, tmpfile = "tmp.baseline",
                    verbose = FALSE, cleanTmp = FALSE) {
  nalgs <- length(baselineTests)
  ## Variables to collect return values:
  results <- param.min <- list()
  savefiles <- character(nalgs)
  minQualMeas <- numeric(nalgs)
  
  if (verbose) cat("Looping through the", nalgs, "baseline algorithm tests\n")
  for (i in seq_len(nalgs)) {
    fname <- funcName(baselineTests[[i]])
    blname <- names(baselineTests)[i]
    if (is.null(blname) || !nzchar(blname))
      blname <- name(algorithm(baselineTests[[i]]))
    if (verbose) cat(i, ": ", blname, ":\n", sep = "")
    savefiles[i] <- paste(tmpfile, blname, "RData", sep = ".")
    if (file.exists(savefiles[i])) {
      load(savefiles[i])
      if (verbose) cat(" Loaded pre-calculated results from file", savefiles[i], "\n")
    } else {
      res <- runTest(baselineTests[[i]], X, y, predictionTest, postproc, verbose)
      save(res, file = savefiles[i])
    }
    results[[i]] <- res
    minQualMeas[i] <- minQualMeas(res)
    param.min[[i]] <- param.min(res)
  }
  ## Optionally, delete temporary files:
  if (isTRUE(cleanTmp)) unlink(savefiles)
  ## Find and return the best results
  best.i <- which.min(minQualMeas)
  minQualMeas <- minQualMeas[best.i]
  param.min <- param.min[[best.i]]
  names(results) <- names(baselineTests)
  return(list(baselineTests = baselineTests,
              results = results,
              minQualMeas = minQualMeas,
              baselineAlg.min = name(algorithm(baselineTests[[best.i]])),
              param.min = param.min))
}

## Function to extract minimum from optimisation:
overall.min <- function(results) {
  with(results, list(qualMeas = minQualMeas, algorithm = baselineAlg.min,
                     param = param.min))
}

## The function mvrValstats has been copied (2013-08-10) from package pls (overlapping authors)
## because eval(parent.frame()) seems to break call to pls::msep.
mvrValstats <- function (object, estimate, newdata, ncomp = 1:object$ncomp, 
                         comps, intercept = cumulative, se = FALSE, ...) 
{
  cumulative <- missing(comps) || is.null(comps)
  if (any(estimate == "CV")) {
    if (!cumulative) 
      stop("Cross-validation is not supported when `comps' is specified")
    if (is.null(object$validation)) 
      stop("`object' has no `validation' component")
  }
  nestimates <- length(estimate)
  nresp <- dim(fitted(object))[2]
  respnames <- dimnames(fitted(object))[[2]]
  SSE <- array(dim = c(nestimates, nresp, if (cumulative) 1 + 
                         length(ncomp) else 2), dimnames = list(estimate = estimate, 
                                                                response = respnames, model = if (cumulative) {
                                                                  c("(Intercept)", paste(ncomp, "comps"))
                                                                } else {
                                                                  c("(Intercept)", paste("(Intercept), Comp", paste(comps, 
                                                                                                                    collapse = ", ")))
                                                                }))
  SST <- array(dim = c(nestimates, nresp), dimnames = list(estimate = estimate, 
                                                           response = respnames))
  nobj <- numeric(nestimates)
  names(nobj) <- estimate
  for (i in seq(along = estimate)) {
    switch(estimate[i], train = {
      resp <- as.matrix(model.response(model.frame(object)))
      nobj[i] <- nrow(resp)
      if (inherits(object$na.action, "exclude")) {
        resp <- napredict(object$na.action, resp)
      }
      res <- if (cumulative) residuals(object, ...)[, , 
                                                    ncomp, drop = FALSE] else resp - predict(object, 
                                                                                             comps = comps, ...)
      SST[i, ] <- apply(resp, 2, var, na.rm = TRUE) * (nobj[i] - 
                                                         1)
      SSE[i, , ] <- cbind(SST[i, ], colSums(res^2, na.rm = TRUE))
    }, test = {
      if (missing(newdata)) stop("Missing `newdata'.")
      newdata <- model.frame(formula(object), data = newdata)
      resp <- as.matrix(model.response(newdata))
      pred <- if (cumulative) predict(object, ncomp = ncomp, 
                                      newdata = newdata, ...) else predict(object, 
                                                                           comps = comps, newdata = newdata, ...)
      nobj[i] <- nrow(newdata)
      SST[i, ] <- apply(resp, 2, var) * (nobj[i] - 1)
      SSE[i, , ] <- cbind(colSums(sweep(resp, 2, object$Ymeans)^2), 
                          colSums((pred - c(resp))^2))
    }, CV = {
      resp <- as.matrix(model.response(model.frame(object)))
      nobj[i] <- nrow(resp)
      SST[i, ] <- apply(resp, 2, var) * (nobj[i] - 1)
      SSE[i, , ] <- cbind(object$validation$PRESS0, object$validation$PRESS[, 
                                                                            ncomp, drop = FALSE])
    })
  }
  if (cumulative) 
    comps <- ncomp
  if (intercept) 
    comps <- c(0, comps)
  else SSE <- SSE[, , -1, drop = FALSE]
  return(list(SSE = SSE, SST = SST, nobj = nobj, comps = comps, 
              cumulative = cumulative))
}
