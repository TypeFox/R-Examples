################################
# sperrorest
# Spatial error estimation and variable importance
# Alexander Brenning
# University of Waterloo
################################

##################################################################
# History:
#
# 2009 - 2011 - general code development
#
# Oct-Dec 2011
# ------------
# - package project and documentation created
#
# 29 Dec 2011
# -----------
# - built internal release 0.1-1
#
# 29 Jan 2012
# -----------
# - internal release 0.1-2
# - some bug fixes, e.g. in err.* functions
# - improved support of pooled versus unpooled error estimation
# - changed some argument names
# - this version was used for Angie's analyses
#
# 15 Feb 2012
# -----------
# - updated help pages
# - supports a variety of block bootstrap types
# - some minor bug fixes
#
# 1 Mar 2012 (0.1-5)
# ------------------
# - made training set estimation optional
# - robustified code using try()
#
# 19 June 2012 (0.2-0)
# --------------------
# - last pre-release version
# - replaced Stoyan's data set with Jannes Muenchow's data, adapted examples
#
# 19 June 2012 (0.2-1) - FIRST RELEASE ON CRAN
# --------------------------------------------
#
##################################################################


#' Summarize error statistics obtained by \code{sperrorest}
#'
#' \code{summary.sperroresterror} calculates mean, standard deviation, median etc. of the calculated error measures at the specified level (overall, repetition, or fold).
#' \code{summary.sperrorestpoolederror} does the same with the pooled error, at the overall or repetition level.
#' @name summary.sperroresterror
#' @method summary sperroresterror
#' @param object \code{sperroresterror} resp. \code{sperrorestcombinederror} error object calculated by \code{\link{sperrorest}}
#' @param level Level at which errors are summarized: 0: overall; 1: repetition; 2: fold
#' @param pooled If \code{TRUE} (default), mean and standard deviation etc are calculated between fold-level error estimates. If \code{FALSE}, apply first a \code{\link{weighted.mean}} among folds before calculating mean, standard deviation etc among repetitions. See also Details.
#' @param na.rm Remove \code{NA} values? See \code{\link{mean}} etc.
#' @param ... additional arguments (currently ignored)
#' @return Depending on the level of aggregation, a \code{list} or \code{data.frame} with mean, and at level 0 also standard deviation, median and IQR of the error measures.
#' @details Let's use an example to explain the \code{pooled} argument. E.g., assume we are using 100-repeated 10-fold cross-validation. 
#' If \code{pooled=TRUE} (default), the mean and standard deviation calculated when summarizing at \code{level=0}
#' are calculated across the error estimates obtained for each of the \code{100*10 = 1000} folds.
#' If \code{pooled=FALSE}, mean and standard deviation are calculated across the \code{100} repetitions, using the weighted average of the fold-level errors to calculate an error value for the entire sample. This will essentially not affect the mean value but of course the standard deviation of the error. \code{pooled=FALSE} is not recommended, it is mainly for testing purposes; when the test sets are small (as in leave-one-out cross-validation, in the extreme case), consider running \code{\link{sperrorest}} with \code{err.pooled=TRUE} and examine only the \code{pooled.error} component of its result.
#' @seealso \code{\link{sperrorest}}
#' @export
summary.sperroresterror = function(object, level = 0, pooled = TRUE, na.rm = TRUE, ...)
{
    err = unclass(object)
    if (pooled) {
        if (level <= 2)
            err = lapply(err, function(x) t(sapply(x, function(y) data.frame(train = y$train, test = y$test, 
                distance = ifelse(any(names(y)=="distance"),y$distance,-1)))) )
        if (level <= 1) {
            errdf = err[[1]]
            if (length(err) > 1) {
                for (i in 2:length(err)) {
                    errdf = rbind(errdf, err[[i]])
                }
            }
            rownames(errdf) = NULL
            err = as.data.frame(errdf)
        }
        if (level <= 0) {
            err = data.frame(
                mean   = apply(err, 2, function(y) mean(unlist(y), na.rm = na.rm)),
                sd     = apply(err, 2, function(y) sd(unlist(y), na.rm = na.rm)),
                median = apply(err, 2, function(y) median(unlist(y), na.rm = na.rm)),
                IQR    = apply(err, 2, function(y) IQR(unlist(y), na.rm = na.rm)) )
        }
    } else {
        if (level <= 2)
            err = lapply(err, function(x) t(sapply(x, function(y) data.frame(train = y$train, test = y$test, 
                distance = ifelse(any(names(y)=="distance"),y$distance,-1)))) )
        if (level <= 1) {
            ###w = summary.partition(resampling) ?????
            err = lapply( err, function(x) apply(x, 2, function(y) weighted.mean(unlist(y), na.rm = na.rm)) )
            nms = names(err)
            err = as.data.frame(t(as.data.frame(err)))
            rownames(err) = nms
        }
        if (level <= 0) {
            err = data.frame(
                mean = sapply(err, mean), 
                sd = sapply(err, sd),
                median = sapply(err, median),
                IQR = sapply(err, IQR))
        }
    }
    return(err)
}

#' @rdname summary.sperrorest
#' @inheritParams summary.sperroresterror
#' @name summary.sperrorestpoolederror
#' @method summary sperrorestpoolederror
#' @export
summary.sperrorestpoolederror = function(object, level = 0, na.rm = TRUE, ...)
{
    class(object) = NULL
    object = as.data.frame(object)
    if (level <= 0) {
        object = data.frame(
            mean = sapply(object, mean), 
            sd = sapply(object, sd),
            median = sapply(object, median),
            IQR = sapply(object, IQR))
    }
    return(object)
}

#' Summarize variable importance statistics obtained by \code{sperrorest}
#'
#' \code{summary.sperrorestimportance} calculated mean, standard deviation, median etc. of the calculated error measures at the specified level (overall, repetition, or fold).
#' @name summary.sperrorestimportance
#' @method summary sperrorestimportance
#' @param object \code{\link{sperrorestimportance}} object calculated by \code{\link{sperrorest}} called with argument \code{importance=TRUE}
#' @inheritParams summary.sperroresterror
#' @param which optional character vector specifying selected variables for which the importances should be summarized (to do: check implementation)
#' @return a list or data.frame, depending on the \code{level} of aggregation
#' @export
summary.sperrorestimportance = function(object, level = 0, na.rm = TRUE, which = NULL, ...)
{
    arrdim = c( length(object), length(object[[1]]), dim(object[[1]][[1]]) )
    arrdimnames = list( names(object), names(object[[1]]),
            rownames(object[[1]][[1]]), colnames(object[[1]][[1]]) )
    arr = array( NA, dim = arrdim, dimnames = arrdimnames )
    for (i in 1:length(object))
        for (j in 1:length(object[[i]]))
            arr[i,j,,] = as.matrix(object[[i]][[j]])
    if (level <= 1)
        arr = apply(arr, c(1,3,4), mean, na.rm = na.rm)
    if (level <= 0) {
        if (is.null(which)) {
            arr = data.frame(
                mean = apply(arr, c(2,3), mean, na.rm = na.rm),
                sd = apply(arr, c(2,3), sd, na.rm = na.rm),
                median = apply(arr, c(2,3), median, na.rm = na.rm),
                IQR = apply(arr, c(2,3), IQR, na.rm = na.rm) )
        } else {
            arr = arr[ , , which ]
            arr = data.frame(
                mean = apply(arr, 2, mean, na.rm = na.rm),
                sd = apply(arr, 2, sd, na.rm = na.rm),
                median = apply(arr, 2, median, na.rm = na.rm),
                IQR = apply(arr, 2, IQR, na.rm = na.rm) )
        }
    }            
    return(arr)
}


#' Draw stratified random sample
#'
#' \code{resample.strat.uniform} draws a stratified random sample (with or without replacement) from the samples in \code{data}. Stratification is over the levels of \code{data[,param$response]}. The same number of samples is drawn within each level.
#' @param data a \code{data.frame}, rows represent samples
#' @param param a list with the following components: \code{strat} is either the name of a factor variable in \code{data} that defines the stratification levels, or a vector of type factor and length \code{nrow(data)}; \code{n} is a numeric value specifying the size of the subsample; \code{replace} determines if sampling is with or without replacement
#' @return a \code{data.frame} containing a subset of the rows of \code{data}.
#' @details If \code{param$replace=FALSE}, a subsample of size \code{min(param$n,nrow(data))} will be drawn from \code{data}. If \code{param$replace=TRUE}, the size of the subsample is \code{param$n}.
#' @seealso \code{\link{resample.uniform}}, \code{\link{sample}}
#' @examples
#' data(ecuador) # Muenchow et al. (2012), see ?ecuador
#' d = resample.strat.uniform(ecuador, param = list(strat = "slides", nstrat = 100))
#' nrow(d) # == 200
#' sum(d$slides == "TRUE") # == 100
#' @export
resample.strat.uniform = function(data, param = list( strat = "class", nstrat = Inf, replace = FALSE ))
{
    # Old version:
    if (!is.null(param$response)) {
        warning("'param$response' argument in 'resample.strat.uniform' renamed to 'strat';\n modify your code accordingly")
        if (is.null(param$strat)) param$strat = param$response
    }

    # Use defaults if not specified:
    if (is.null(param$strat))    param$strat = "class"
    if (is.null(param$nstrat))   param$nstrat = Inf
    if (is.null(param$replace))  param$replace = FALSE
    
    stopifnot( (length(param$strat) == 1) | (length(param$strat) == nrow(data)) )
    if (length(param$strat == 1)) {
        strat = data[,param$strat]
    } else strat = param$strat
    if (!is.factor(strat))
        stop("'strat' must either be a vector of factor type, or the name of a factor variable in 'data'")
    # Each factor level must have at least one sample, otherwise sampling within this level is impossible:
    minstrat = min(tapply(strat, strat, length))
    stopifnot(minstrat >= 1)
    # might want to change this to a warning.???

    if (!param$replace) param$nstrat = min(param$nstrat, minstrat)

    # Uniform sampling within each stratum:
    sel = c()
    for (lev in levels(strat)) {
        wh = sample( which(strat == lev), size = param$nstrat, replace = param$replace )
        sel = c(sel, wh)
    }
    return(data[sel,])
}
# To do: allow nstrat to be a named vector


#' Draw uniform random (sub)sample
#'
#' \code{resample.uniform} draws a random (sub)sample (with or without replacement) from the samples in \code{data}.
#' @param data a \code{data.frame}, rows represent samples
#' @param param a list with the following components: \code{n} is a numeric value specifying the size of the subsample; \code{replace} determines if sampling is with or without replacement
#' @return a \code{data.frame} containing a subset of the rows of \code{data}.
#' @details If \code{param$replace=FALSE}, a subsample of size \code{min(param$n,nrow(data))} will be drawn from \code{data}. If \code{param$replace=TRUE}, the size of the subsample is \code{param$n}.
#' @seealso \code{\link{resample.strat.uniform}}, \code{\link{sample}}
#' @examples
#' data(ecuador) # Muenchow et al. (2012), see ?ecuador
#' d = resample.uniform(ecuador, param = list(strat = "slides", n = 200))
#' nrow(d) # == 200
#' sum(d$slides == "TRUE")
#' @export
resample.uniform = function(data, param = list(n = Inf, replace = FALSE))
{
    # Apply defaults if missing from parameter list:
    if (is.null(param$n)) param$n = Inf
    if (is.null(param$replace)) param$replace = FALSE
    
    if (!param$replace) param$n = min( param$n, nrow(data) )

    # Uniform sampling with or without replacement:
    sel = sample(nrow(data), size = param$n, replace = param$replace)
    
    return( data[ sel, ] )
}


#' Perform spatial error estimation and variable importance assessment
#'
#' \code{sperrorest} is a flexible interface for multiple types of spatial and non-spatial cross-validation 
#' and bootstrap error estimation and permutation-based assessment of spatial variable importance.
#' @inheritParams partition.cv
#' @param data a \code{data.frame} with predictor and response variables. Training and test samples 
#' will be drawn from this data set by \code{train.fun} and \code{test.fun}, respectively.
#' @param formula A formula specifying the variables used by the \code{model}. Only simple formulas 
#' without interactions or nonlinear terms should be used, e.g. \code{y~x1+x2+x3} but not 
#' \code{y~x1*x2+log(x3)}. Formulas involving interaction and nonlinear terms may possibly work 
#' for error estimation but not for variable importance assessment, but should be used with caution.
#' @param coords vector of length 2 defining the variables in \code{data} that contain the x and y coordinates of sample locations
#' @param model.fun Function that fits a predictive model, such as \code{glm} or \code{rpart}. The 
#' function must accept at least two arguments, the first one being a formula and the second a 
#' data.frame with the learning sample.
#' @param model.args Arguments to be passed to \code{model.fun} (in addition to the \code{formula} and \code{data} argument, which are provided by \code{sperrorest})
#' @param pred.fun Prediction function for a fitted model object created by \code{model}. 
#' Must accept at least two arguments: the fitted \code{object} and a \code{data.frame} \code{newdata} with data 
#' on which to predict the outcome.
#' @param pred.args (optional) Arguments to \code{pred.fun} (in addition to the fitted model object and the \code{newdata} argument, which are provided by \code{sperrorest})
#' @param smp.fun A function for sampling training and test sets from \code{data}. E.g., \code{\link{partition.kmeans}} for spatial cross-validation using spatial \emph{k}-means clustering.
#' @param smp.args (optional) Arguments to be passed to \code{est.fun}
#' @param train.fun (optional) A function for resampling or subsampling the training sample in order to achieve, e.g., uniform sample sizes on all training sets, or maintaining a certain ratio of positives and negatives in training sets. E.g., \code{\link{resample.uniform}} or \code{\link{resample.strat.uniform}}
#' @param train.param (optional) Arguments to be passed to \code{resample.fun}
#' @param test.fun (optional) Like \code{train.fun} but for the test set.
#' @param test.param (optional) Arguments to be passed to \code{test.fun}
#' @param err.fun A function that calculates selected error measures from the known responses in \code{data} and the model predictions delivered by \code{pred.fun}. E.g., \code{\link{err.default}} (the default). See example and details below.
#' @param err.unpooled logical (default: \code{TRUE}): calculate error measures on each fold within a resampling repetition
#' @param err.pooled logical (default: \code{FALSE}): calculate error measures based on the pooled predictions of all folds within a resampling repetition
#' @param err.train logical (default: \code{TRUE}): calculate error measures on the training set (in addition to the test set estimation)
#' @param imp.variables (optional; used if \code{importance=TRUE}) Variables for which permutation-based variable importance assessment is performed. If \code{importance=TRUE} and \code{imp.variables} is \code{NULL}, all variables in \code{formula} will be used.
#' @param imp.permutations (optional; used if \code{importance=TRUE}) Number of permutations used for variable importance assessment.
#' @param importance logical: perform permutation-based variable importance assessment?
#' @param ... currently not used
#' @param distance logical (default: \code{FALSE}): if \code{TRUE}, calculate mean nearest-neighbour distances from test samples to training samples using \code{\link{add.distance.represampling}}
#' @param do.gc numeric (default: 1): defines frequency of memory garbage collection by calling \code{\link{gc}}; if \code{<1}, no garbage collection; if \code{>=1}, run a \code{gc()} after each repetition; if \code{>=2}, after each fold
#' @param do.try logical (default: \code{FALSE}): if \code{TRUE} [untested!!], use \code{\link{try}} to robustify calls to \code{model.fun} and \code{err.fun}; use with caution!
#' @param silent If \code{TRUE}, show progress on console (in Windows Rgui, disable 'Buffered output' in 'Misc' menu)
#' @return A list (object of class \code{sperrorest}) with (up to) four components:
#' \item{error}{a \code{sperroresterror} object containing predictive performances at the fold level}
#' \item{represampling}{a \code{\link{represampling}} object}
#' \item{pooled.error}{a \code{sperrorestpoolederror} object containing predictive performances at the repetition level}
#' \item{importance}{a \code{sperrorestimportance} object containing permutation-based variable importances at the fold level}
#' @return An object of class \code{sperrorest}, i.e. a list with components \code{error} (of class \code{sperroresterror}), \code{represampling} (of class \code{represampling}), \code{pooled.error} (of class \code{sperrorestpoolederror}) and \code{importance} (of class \code{sperrorestimportance}).
#' @note To do: (1) Parallelize the code; (2) Optionally save fitted models, training and test samples in the results object; (3) Optionally save intermediate results in some file, and enable the function to continue an interrupted sperrorest call where it was interrupted. (3) Optionally have sperrorest dump the result of each repetition into a file, and to skip repetitions for which a file already exists. (4) Save sperrorest version number in results object.
#' @references Brenning, A. 2012. Spatial cross-validation and bootstrap for the assessment of 
#' prediction rules in remote sensing: the R package 'sperrorest'. 
#' IEEE International Symposium on Geoscience and Remote Sensing IGARSS, in press.
#'
#' Brenning, A. 2005. Spatial prediction models for landslide hazards: review, comparison and evaluation. Natural Hazards and Earth System Sciences, 5(6): 853-862.
#'
#' Brenning, A., S. Long & P. Fieguth. Forthcoming. Detecting rock glacier flow structures using Gabor filters and IKONOS imagery. Submitted to Remote Sensing of Environment.
#'
#' Russ, G. & A. Brenning. 2010a. Data mining in precision agriculture: Management of spatial information. In 13th International Conference on Information Processing and Management of Uncertainty, IPMU 2010; Dortmund; 28 June - 2 July 2010.  Lecture Notes in Computer Science, 6178 LNAI: 350-359.
#'
#' Russ, G. & A. Brenning. 2010b. Spatial variable importance assessment for yield prediction in Precision Agriculture. In Advances in Intelligent Data Analysis IX, Proceedings, 9th International Symposium, IDA 2010, Tucson, AZ, USA, 19-21 May 2010.  Lecture Notes in Computer Science, 6065 LNCS: 184-195.
#' @seealso \pkg{ipred}
#' @aliases sperroresterror sperrorestimportance
#' @examples
#' data(ecuador) # Muenchow et al. (2012), see ?ecuador
#' fo = slides ~ dem + slope + hcurv + vcurv + 
#'      log.carea + cslope
#' 
#' # Example of a classification tree fitted to this data:
#' library(rpart)
#' ctrl = rpart.control(cp = 0.005) # show the effects of overfitting
#' fit = rpart(fo, data = ecuador, control = ctrl)
#' par(xpd = TRUE)
#' plot(fit, compress = TRUE, main = "Stoyan's landslide data set")
#' text(fit, use.n = TRUE)
#'
#' # Non-spatial 5-repeated 10-fold cross-validation:
#' mypred.rpart = function(object, newdata) predict(object, newdata)[,2]
#' nspres = sperrorest(data = ecuador, formula = fo,
#'     model.fun = rpart, model.args = list(control = ctrl),
#'     pred.fun = mypred.rpart,
#'     smp.fun = partition.cv, smp.args = list(repetition=1:5, nfold=10))
#' summary(nspres$error)
#' summary(nspres$represampling)
#' plot(nspres$represampling, ecuador)
#'
#' # Spatial 5-repeated 10-fold spatial cross-validation:
#' spres = sperrorest(data = ecuador, formula = fo,
#'     model.fun = rpart, model.args = list(control = ctrl),
#'     pred.fun = mypred.rpart,
#'     smp.fun = partition.kmeans, smp.args = list(repetition=1:5, nfold=10))
#' summary(spres$error)
#' summary(spres$represampling)
#' plot(spres$represampling, ecuador)
#' 
#'  smry = data.frame(
#'      nonspat.training = unlist(summary(nspres$error,level=1)$train.auroc),
#'      nonspat.test     = unlist(summary(nspres$error,level=1)$test.auroc),
#'      spatial.training = unlist(summary(spres$error,level=1)$train.auroc),
#'      spatial.test     = unlist(summary(spres$error,level=1)$test.auroc))
#' boxplot(smry, col = c("red","red","red","green"), 
#'      main = "Training vs. test, nonspatial vs. spatial",
#'      ylab = "Area under the ROC curve")
#' @export
sperrorest = function(formula, data, coords = c("x", "y"),
    model.fun, model.args = list(),
    pred.fun = NULL, pred.args = list(),
    smp.fun = partition.loo, smp.args = list(),
    train.fun = NULL, train.param = NULL,
    test.fun = NULL, test.param = NULL,
    err.fun = err.default,
    err.unpooled = TRUE,
    err.pooled = FALSE,
    err.train = TRUE,
    imp.variables = NULL,
    imp.permutations = 1000,
    importance = !is.null(imp.variables),
    distance = FALSE,
    do.gc = 1,
    do.try = FALSE,
    silent = FALSE, ...)
{
    # Some checks:
    if (missing(model.fun)) stop("'model.fun' is a required argument")
    if (as.character(attr(terms(formula),"variables"))[3] == "...")
        stop("formula of the form lhs ~ ... not accepted by 'sperrorest'\nspecify all predictor variables explicitly")
    stopifnot(is.function(model.fun))
    stopifnot(is.function(smp.fun))
    if (!is.null(train.fun)) stopifnot(is.function(train.fun))
    if (!is.null(test.fun)) stopifnot(is.function(test.fun))
    stopifnot(is.function(err.fun))
    if (importance) {
        if (!err.unpooled) {
            warning("'importance=TRUE' currently only supported with 'err.unpooled=TRUE'.\nUsing 'importance=FALSE'")
            importance = FALSE
        }
        stopifnot(is.numeric(imp.permutations))
        if (!is.null(imp.variables)) stopifnot(is.character(imp.variables))
    }
    stopifnot(is.character(coords))
    stopifnot(length(coords) == 2)
    if (importance & !err.unpooled)
        stop("variable importance assessment currently only supported at the unpooled level")

    # Check if user is trying to bypass the normal mechanism for generating training and test data sets 
    # and for passing formulas:
    if (any(names(model.args) == "formula")) stop("'model.args' cannot have a 'formula' element")
    if (any(names(model.args) == "data")) stop("'model.args' cannot have a 'data' element")
    if (any(names(pred.args) == "object")) stop("'pred.args' cannot have an 'object' element:\nthis will be generated by 'sperrorest'")
    if (any(names(pred.args) == "newdata")) stop("'pred.args' cannot have a 'newdata' element:\nthis will be generated by 'sperrorest'")

    # Some checks related to recent changes in argument names:
    dots.args = list(...)
    if (length(dots.args) > 0) {
        if (any(names(dots.args) == "predfun")) stop("sorry: argument names have changed; 'predfun' is now 'pred.fun'")
        if (any(names(dots.args) == "model")) stop("sorry: argument names have changed; 'model' is now 'model.fun'")
        if (any(names(dots.args) == "err.combined")) stop("sorry: argument names have changed; 'err.combined' is now 'err.pooled'")
        if (any(names(dots.args) == "err.uncombined")) stop("sorry: argument names have changed; 'err.uncombined' is now 'err.unpooled'")
        warning("'...' arguments currently not supported:\nuse 'model.args' to pass list of additional arguments to 'model.fun'")
    }

    # Name of response variable:
    response = as.character(attr(terms(formula),"variables"))[2]
    
    smp.args$data = data
    smp.args$coords = coords

    resamp = do.call(smp.fun, args = smp.args)
    if (distance)
        # Parallelize this function???
        resamp = add.distance(resamp, data, coords = coords, fun = mean)

    if (err.unpooled) {
        res = lapply(resamp, unclass)
        class(res) = "sperroresterror"
    } else res = NULL
    pooled.err = NULL
    # required to be able to assign levels to predictions if appropriate:
    is.factor.prediction = NULL

    ### Permutation-based variable importance assessment (optional):
    impo = NULL
    if (importance) {
        # Importance of which variables:
        if (is.null(imp.variables))
            imp.variables = strsplit(as.character(formula)[3]," + ",fixed=TRUE)[[1]]
        # Dummy data structure that will later be populated with the results:
        impo = resamp
        # Create a template that will contain results of variable importance assessment:
        imp.one.rep = as.list( rep(NA, length(imp.variables)) )
        names(imp.one.rep) = imp.variables
        tmp = as.list(rep(NA, imp.permutations))
        names(tmp) = as.character(1:imp.permutations)
        for (vnm in imp.variables) imp.one.rep[[vnm]] = tmp
        rm(tmp)
    }

    # For each repetition:
    for (i in 1:length(resamp)) {
        if (!silent) cat(date(), "Repetition", names(resamp)[i], "\n")

        # Collect pooled results in these data structures:
        if (err.train) pooled.obs.train = pooled.pred.train = c()
        pooled.obs.test = pooled.pred.test = c()

        # Parallelize this???
        # For each fold:
        for (j in 1:length(resamp[[i]]))
        {
            if (!silent) cat(date(), "- Fold", j, "\n")
            
            # Create training sample:
            nd = data[ resamp[[i]][[j]]$train , ]
            if (!is.null(train.fun))
                nd = train.fun(data = nd, param = train.param)

            # Train model on training sample:
            margs = c( list(formula = formula, data = nd), model.args )
            
            if (do.try) {
                fit = try(do.call(model.fun, args = margs), silent = silent)
            
                # Error handling:
                if (class(fit) == "try-error") {
                    fit = NULL
                    if (err.unpooled) {
                        if (err.train) res[[i]][[j]]$train = NULL
                        res[[i]][[j]]$test = NULL
                        if (importance) impo[[i]][[j]] = c() # ???
                    }
                    if (do.gc >= 2) gc()
                    next # skip this fold
                }
                
            } else {
                fit = do.call(model.fun, args = margs)
            }

            if (err.train) {
                # Apply model to training sample:
                pargs = c( list(object = fit, newdata = nd), pred.args )
                if (is.null(pred.fun)) {
                    pred.train = do.call(predict, args = pargs)
                } else {
                    pred.train = do.call(pred.fun, args = pargs)
                }
                rm(pargs)
            
                # Calculate error measures on training sample:
                if (err.unpooled)
                    if (do.try) {
                        err.try = try(err.fun(nd[,response], pred.train), silent = silent)
                        if (class(err.try) == "try-error") err.try = NULL # ???
                        res[[i]][[j]]$train = err.try
                    } else {
                        res[[i]][[j]]$train = err.fun(nd[,response], pred.train)
                    }
                if (err.pooled) {
                    pooled.obs.train = c( pooled.obs.train, nd[,response] )
                    pooled.pred.train = c( pooled.pred.train, pred.train )
                }
            } else {
                if (err.unpooled) res[[i]][[j]]$train = NULL
            }
            
            # Create test sample:
            nd = data[ resamp[[i]][[j]]$test , ]
            if (!is.null(test.fun))
                nd = test.fun(data = nd, param = test.param)
            # Create a 'backup' copy for variable importance assessment:
            if (importance) nd.bak = nd
            # Apply model to test sample:
            pargs = c( list(object = fit, newdata = nd), pred.args )
            if (is.null(pred.fun)) {
                pred.test  = do.call(predict, args = pargs)
            } else {
                pred.test  = do.call(pred.fun, args = pargs)
            }
            rm(pargs)
            
            # Calculate error measures on test sample:
            if (err.unpooled) {
                if (do.try) {
                    err.try = try(err.fun(nd[,response], pred.test), silent = silent)
                    if (class(err.try) == "try-error") err.try = NULL # ???
                    res[[i]][[j]]$test = err.try
                } else {
                    res[[i]][[j]]$test  = err.fun(nd[,response], pred.test)
                }
            }
            if (err.pooled) {
                pooled.obs.test = c( pooled.obs.test, nd[,response] )
                pooled.pred.test = c( pooled.pred.test, pred.test )
                is.factor.prediction = is.factor(pred.test)
            }
            
            ### Permutation-based variable importance assessment:
            if (importance & err.unpooled) {
            
                if (is.null(res[[i]][[j]]$test)) {
                    impo[[i]][[j]] = c()
                    if (!silent) cat(date(), "-- skipping variable importance\n")
                } else {

                    if (!silent) cat(date(), "-- Variable importance\n")
                    imp.temp = imp.one.rep
    
                    # Parallelize this: ???
                    for (cnt in 1:imp.permutations) {
                        # Some output on screen:
                        if (!silent & (cnt>1))
                            if (log10(cnt)==floor(log10(cnt))) 
                                cat(date(), "   ", cnt, "\n")
    
                        # Permutation indices:
                        permut = sample(1:nrow(nd), replace = FALSE)
    
                        # For each variable:
                        for (vnm in imp.variables) {
                            # Get undisturbed backup copy of test sample:
                            nd = nd.bak
                            # Permute variable vnm:
                            nd[,vnm] = nd[,vnm][permut]
                            # Apply model to perturbed test sample:
                            pargs = c( list(object = fit, newdata = nd), pred.args )
                            if (is.null(pred.fun)) {
                                pred.test  = do.call(predict, args = pargs)
                            } else {
                                pred.test  = do.call(pred.fun, args = pargs)
                            }
                            rm(pargs)
                            
                            # Calculate variable importance:
                            if (do.try) {
                                permut.err = try(err.fun(nd[,response], pred.test), silent = silent)
                                if (class(permut.err) == "try-error") {
                                    imp.temp[[vnm]][[cnt]] = c() # ???
                                } else {
                                    imp.temp[[vnm]][[cnt]] = 
                                        as.list( unlist(res[[i]][[j]]$test) - unlist(permut.err) )
                                        # (apply '-' to corresponding list elements; only works
                                        # if all list elements are scalars)
                                }
                            } else {
                                permut.err = err.fun(nd[,response], pred.test)
                                imp.temp[[vnm]][[cnt]] = 
                                    as.list( unlist(res[[i]][[j]]$test) - unlist(permut.err) )
                                    # (apply '-' to corresponding list elements; only works
                                    # if all list elements are scalars)
                            }
                        }
                    }
                    # average the results obtained in each permutation:
                    impo[[i]][[j]] = as.data.frame(
                        t(sapply(imp.temp, 
                            function(y) sapply(as.data.frame(t(
                                sapply( y, as.data.frame ))), 
                                    function(x) mean(unlist(x)) ))))
                    rm(nd.bak, nd) # better safe than sorry...
                } # end of else if (!is.null(res[[i]][[j]]$test))
            }
            
        }
        
        # Put the results from the pooled estimation into the pooled.err data structure:
        if (err.pooled) {
            if (is.factor(data[,response])) {
                lev = levels(data[,response])
                if (err.train) pooled.obs.train = factor(lev[pooled.obs.train], levels = lev)
                pooled.obs.test = factor(lev[pooled.obs.test], levels = lev)
                if (is.factor.prediction) {
                    if (err.train) pooled.pred.train = factor(lev[pooled.pred.train], levels = lev)
                    pooled.pred.test = factor(lev[pooled.pred.test], levels = lev)
                }
            }
            pooled.err.train = NULL
            if (err.train)
                pooled.err.train = err.fun( pooled.obs.train, pooled.pred.train )
            if (i == 1) {
                pooled.err = t(unlist( list( train = pooled.err.train,
                                           test  = err.fun( pooled.obs.test,  pooled.pred.test ) ) ))
            } else {
                pooled.err = rbind( pooled.err, 
                             unlist( list( train = pooled.err.train,
                                           test  = err.fun( pooled.obs.test,  pooled.pred.test ) ) ) )
            }
            
            if (do.gc >= 2) gc()
        } # end for each fold
        
        if ((do.gc >= 1) & (do.gc < 2)) gc()
    } # end for each repetition

    # convert matrix(?) to data.frame:
    if (err.pooled) {
        pooled.err = as.data.frame(pooled.err)
        rownames(pooled.err) = NULL
        class(pooled.err) = "sperrorestpoolederror"
    }
    
    if (!silent) cat(date(), "Done.\n")

    if (importance) class(impo) = "sperrorestimportance"
    
    RES = list(
        error = res, 
        represampling = resamp, 
        pooled.error = pooled.err,
        importance = impo )
    class(RES) = "sperrorest"
    
    return( RES )
}

#' Summary and print methods for sperrorest results
#'
#' Summary methods provide varying level of detail while print methods provide full details.
#' @name summary.sperrorest
#' @method summary sperrorest
#' @param object a \code{\link{sperrorest}} object
#' @param ... additional arguments for \code{\link{summary.sperroresterror}} or \code{\link{summary.sperrorestimportance}}
#' @param x Depending on method, a \code{\link{sperrorest}}, \code{\link{sperroresterror}} or \code{\link{sperrorestimportance}} object
#' @seealso \code{\link{sperrorest}}, \code{\link{sperroresterror}}, \code{\link{sperrorestimportance}}, 
#' \code{\link{summary.sperroresterror}}, \code{\link{summary.sperrorestimportance}}
#' @export
summary.sperrorest = function(object, ...) {
    list(
        error = summary(object$error, ...),
        represampling = summary(object$represampling, ...),
        pooled.error = summary(object$pooled.error, ...),
        importance = summary(object$importance, ...) )
}

#' @rdname summary.sperrorest
#' @name print.sperrorestimportance
#' @method print sperrorestimportance
#' @export
print.sperrorestimportance = function(x, ...) print(unclass(summary(x, level = Inf, ...)))

#' @rdname summary.sperrorest
#' @name print.sperroresterror
#' @method print sperroresterror
#' @export
print.sperroresterror = function(x, ...) print(unclass(summary(x, level = Inf, ...)))

#' @rdname summary.sperrorest
#' @name print.sperrorestpoolederror
#' @method print sperrorestpoolederror
#' @export
print.sperrorestpoolederror = function(x, ...) print(unclass(summary(x, level = Inf, ...)))

#' @rdname summary.sperrorest
#' @name print.sperrorest
#' @method print sperrorest
#' @export
print.sperrorest = function(x, ...) print(unclass(summary(x, level = Inf, ...)))
