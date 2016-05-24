#' @name gapfill-package
#' @aliases Gapfill-Package Gapfill-package gapfill-Package 
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Overview
#' @description
#' The package provides tools to fill-in missing values in satellite data.
#' It can be used to gap-fill, e.g., MODIS NDVI data and
#' is helpful when developing new gap-fill algorithms.
#' The methods are tailored to data (images) observed at equally-spaced points in time. 
#' This is typically the case for MODIS land surface products and AVHRR NDVI data, among others. \cr
#' The predictions of the missing values are based on a subset-predict procedure, i.e.,
#' each missing value is predicted separately by
#' (1) selecting subsets of the data that are in a neighborhood around the missing point and
#' (2) predicting the missing value based on the subset.\cr
#' The main function of the package is \code{\link{Gapfill}}.
#' 
#' @section Features:
#' \itemize{
#' \item Gap-filling can be executed in parallel. 
#'
#' \item Users may define new \code{\link{Subset}} and \code{\link{Predict}} functions
#' and run alternative prediction algorithms with little effort.
#' See \link{Extend} for more information and examples.  
#'
#' \item Visualization of space-time data are simplified through the \code{ggplot2}-based
#' function \code{\link{Image}}.
#' }
#'
#' @references F. Gerber, R. Furrer, G. Schaepman-Strub, R. de Jong, M. E. Schaepman, 2016,
#' Predicting missing values in spatio-temporal satellite data.
#' \url{http://arxiv.org/abs/1605.01038}. 
#' @seealso \code{\link{Gapfill}}, \code{\link{Subset-Predict}}, \code{\link{Extend}}, \code{\link{Image}}.
#' @docType package
#' @name gapfill-package
#' @keywords package
NULL

#' @name Gapfill
#' @aliases gapfill gap-fill Gap-fill
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Main Function for Gap-Filling
#' @description
#' The function fills (predicts) missing values in satellite data.
#' We illustrate it with MODIS NDVI data,
#' but it can also be applied to other data, that is recorded at equally spaced points in time.
#' Moreover, the function provides infrastructure for the development of new gap-fill algorithms.
#' The predictions of the missing values are based on a subset-predict procedure, i.e.,
#' each missing value is predicted separately by
#' (1) selcting a subset of the data to a neighborhood around the missing value and
#' (2) predicting the values based on that subset.
#' 
#' @param data Numeric array with four dimensions. The input (satellite) data to be gap-filled.
#' Missing values should be encoded as \code{NA}. When using the default \code{\link{Subset}} and \code{\link{Predict}}
#' functions, the data should have the dimensions: x coordinate, y coordinate, seasonal index (e.g., day of the year), and year.
#' See the \code{\link{ndvi}} dataset for an example. 
#' @param fnSubset Function to subset the \code{data} around a missing value.
#' See \code{\link{Subset}} and \link{Extend} for more information.
#' @param fnPredict Function to predict a missing value based on the return value of \code{fnSubset}.
#' See \code{\link{Predict}} and \link{Extend} for more information.
#' @param iMax Integer vector of length 1.
#' The maximum number of iterations until \code{NA} is returned as predicted value.  
#' @param nPredict Integer vector of length 1. Specifies the length of the vector returned from \code{fnPredict}.
#' Values larger than 1 may increase memory usage considerably. 
#' @param subset If \code{"missing"}, all missing values in \code{data} are filled.
#' If a logical array of the same dimensions as \code{data} or 
#' a vector with positive integers, only the missing elements of \code{data[subset]} are filled.
#' Note that this does not the same effect as selecting a subset of the input data, since
#' independent of the specified subset, all values in \code{data} are used to inform the predictions.
#' @param clipRange Numeric vector of length 2.
#' Specifies the lower and the upper bound of the filled data.
#' Values outside this range are clipped accordingly.
#' If \code{nPredict} is larger than 2, only the first return value of \code{fnPredict} will be clipped. 
#' @param dopar Logical vector of length 1.
#' If \code{TRUE}, the \code{\%dopar\%} construct from the R package foreach is used.
#' This allows the function to predict several missing values in parallel,
#' if a parallel back-end (e.g., from the R package doParallel or doMpi) is available.
#' See the example below and \code{\link[foreach]{foreach}} for more information. 
#' @param verbose Logical vector of length 1.
#' If \code{TRUE}, messages are printed.
#' @param ... Additional arguments passed to \code{fnSubset} and \code{fnPredict}. 
#'
#' @return List of length 4 with the entries:
#' \itemize{
#'  \item{\code{fill}}{ contains the gap-filled data.
#' if \code{nPredict = 1}, \code{fill} is an array of dimension \code{dim(data)},
#' otherwise the array is of dimension \code{c(dim(data), nPredict)}.}
#'  \item{\code{mps}}{ integer vector of length equaling the number of predicted values.
#' Contains the (1 dimensional) indices of the predicted values.}
#' \item{\code{time}}{ list of length 4 containing timing information.
#'                     \itemize{\item{\code{start}}{ start date and time.}
#'                              \item{\code{end}}{ end date and time.}
#'                              \item{\code{elapsedMins}}{ elapsed minutes.}
#'                              \item{\code{elapsedSecsPerNA}}{ elapsed seconds per predicted value.}
#'                              }
#'                    }
#'  \item{\code{call}}{ call used to produce the object.}
#' }
#'
#' @details
#' The predictions of the missing values are based on a subset-predict procedure, i.e.,
#' each missing value is predicted separately by
#' (1) selecting a subset of the data to a
#' neighborhood around it and (2) predicting the values based on
#' that subset. The following gives more information on this subset-predict strategy.\cr
#' Missing values are often unevenly distributed in \code{data}.
#' Therefore, the size of a reasonable subset may be different depending on the position of the considered missing value.  
#' The search strategy to find that subset is encoded in \code{fnSubset}.
#' The function returns different subsets depending on the argument \code{i}.
#' The decision whether or not a subset is suitable and the prediction itself is
#' implemented in \code{fnPredict}.
#' To be more specific, the subset-predict procedure loops over the following two steps to predict one missing value:
#' \describe{
#' \item{(1) }{The function \code{fnSubset} is provided with the argument \code{i = i} (where \code{i <- 0} in the first iteration) and
#' returns a subset around the missing value.}
#' \item{(2) }{The function \code{fnPredict} decides whether the subset contains enough information to predict the missing value.
#' If so, the predicted value is returned.
#' Otherwise, the function returns \code{NA} and the algorithm increases \code{i} by one (\code{i <- i + 1})
#' before continuing with step (1).}
#' }
#' The procedure stops if one of the following criteria is met:
#' \itemize{
#' \item \code{fnPredict} returns a non-\code{NA} value,
#' \item \code{iMax} tries have be completed,
#' \item \code{fnSubset} returns the same subset two times in a row. 
#' }
#' 
#' 
#' 
#' @note
#' The default \code{\link{Predict}} function implements the prediction of the missing value
#' and can also return lower and upper bounds of an approximated 90\% prediction interval.
#' See the help page of \code{\link{Predict}} for more information on the prediction interval.
#' The example section below shows how the prediction interval can be calculated and displayed.  
#' 
#' To tailor the procedure to a specific dataset, it might be necessary to
#' adapt the subset and/or the prediction strategy.
#' On the one hand, this can be done by changing the default arguments of \code{\link{Subset}} and
#' \code{\link{Predict}} through the argument \code{...} of \code{Gapfill}.
#' See the help of the corresponding functions for more information about their arguments.
#' On the other hand, the user can define a new subset and predict functions, and pass them to \code{Gapfill}
#' through the arguments \code{fnSubset} and \code{fnPredict}.
#' See \link{Extend} for more information. 
#' 
#' The current implementation of \code{\link{Subset}} does not take into account
#' that values at the boundaries of \code{data} can be neighboring to each other.
#' For example, if global data (entire sphere) are considered,
#' \code{data[1,1,,]} is a neighbor of \code{data[dim(data)[1], dim(data)[2],,]}.
#' Similar considerations apply when data are available for an entire year. 
#' To take this into account, the \code{Subset} function can be redefined accordingly or
#' the data can be augmented.
#' 
#' There are two strategies to run the gap-filling in parallel.
#' The first one is to set the argument \code{dopar} of \code{Gapfill} to \code{TRUE} and
#' to use an openMP or MPI parallel back-end.
#' The parallel back-end needs to be setup before the call to \code{Gapfill}.
#' An example using the R package \code{doParallel} is given below.
#' Note that there exist other parallel back-ends implemented in other packages; such as, e.g., the package \code{doMpi}.
#' Some parallel back-ends are platform dependent. 
#' While this approach shortens the process time by distributing the computational workload,
#' it does not reduce the memory footprint of the procedure.
#' The second strategy, which  also reduces memory usage, is to split the \code{data} into several independent chunks.
#' Whether data chunks are independent or not depends on the function provided to \code{fnSubset}. 
#' For example, the default \code{\link{Subset}} function never includes data that
#' is further apart from the missing value than 1 seasonal index.
#' Hence, \code{data[,,1:3,]} can be used to gap-fill \code{data[,,2,]}.\cr
#'
#' @references F. Gerber, R. Furrer, G. Schaepman-Strub, R. de Jong, M. E. Schaepman, 2016,
#' Predicting missing values in spatio-temporal satellite data.
#' \url{http://arxiv.org/abs/1605.01038}. 
#' @seealso \code{\link{Extend}}, \code{\link{Subset-Predict}}, \code{\link{Image}}.
#' @examples
#' \dontrun{
#' out <- Gapfill(ndvi, clipRange = c(0, 1))
#'
#' ## look at input and output
#' str(ndvi)
#' str(out)
#' Image(ndvi)
#' Image(out$fill)
#'
#' ## run on 2 cores in parallel
#' if(require(doParallel)){
#'   registerDoParallel(2)
#'   out <- Gapfill(ndvi, dopar = TRUE)
#' }
#'
#' ## return also the prediction interval
#' out <- Gapfill(ndvi, nPredict = 3, predictionInterval = TRUE)
#'
#' ## dimension has changed according to 'nPredict = 3'
#' dim(out$fill)
#'
#' ## clip values outside the valid parameter space [0,1].
#' out$fill[out$fill < 0] <- 0
#' out$fill[out$fill > 1] <- 1
#'
#' ## images of the output:
#' ## predicted NDVI
#' Image(out$fill[,,,,1])
#' ## lower bound of the prediction interval
#' Image(out$fill[,,,,2])
#' ## upper bound of the prediction interval
#' Image(out$fill[,,,,3])
#' ## prediction interval length
#' Image(out$fill[,,,,3] - out$fill[,,,,2])
#' 
#' }
#' @export
#' @importFrom foreach foreach %dopar%
Gapfill <- function (data,
                     fnSubset = Subset,
                     fnPredict = Predict,
                     iMax = Inf, 
                     nPredict = 1L,
                     subset = "missings",
                     clipRange = c(-Inf, Inf),
                     dopar = FALSE,
                     verbose = TRUE,
                     ...) 
{
    stopifnot(is.array(data))
    stopifnot(identical(length(dim(data)), 4L))
    stopifnot(is.function(fnSubset))
    fnSubset_formals <- formals(fnSubset)
    stopifnot(all(names(fnSubset_formals)[1:3] == c("data", "mp", 
        "i")))
    stopifnot(is.function(fnPredict))
    fnPredict_formals <- formals(fnPredict)
    stopifnot(all(names(fnPredict_formals)[1:2] == c("a", "i")))
    stopifnot(is.numeric(iMax))
    stopifnot(identical(length(iMax), 1L))
    stopifnot(iMax >= 0L)
    stopifnot(is.numeric(nPredict))
    stopifnot(identical(length(nPredict), 1L))
    stopifnot(nPredict >= 1L)
    stopifnot(subset[1] == "missings" || (identical(dim(subset), 
        dim(data)) && is.array(subset) && is.logical(subset)) ||
        (is.numeric(subset) && is.vector(subset)))
    if (subset[1] == "missings") 
        mps <- which(c(is.na(data)))
    else if (is.logical(subset))
        mps <- which(c(is.na(data) & subset))
    else {
        if (any(duplicated(subset)))
            stop("Argument \"subset\" contains duplicated values.")
        mps <- which(c(is.na(data)))
        mps <- mps[mps %in% subset]
        if (!all(subset %in% mps))
            stop("The values(s) ", paste(subset[!(subset %in% mps)], collapse = ", "),
                 " of argument \"subset\" do not match missing values of \"data\".")
    }
        
    stopifnot(is.numeric(clipRange))
    stopifnot(identical(length(clipRange), 2L))
    stopifnot(clipRange[1] < clipRange[2])
    stopifnot(is.logical(dopar))
    stopifnot(identical(length(dopar), 1L))
    stopifnot(is.logical(verbose))
    stopifnot(identical(length(verbose), 1L))
    dotArgs <- list(...)
    if (!all(names(dotArgs) %in% c(names(formals(fnSubset)), 
        names(formals(fnPredict))))) 
        stop("The \"...\" argument(s) ", paste(sQuote(names(dotArgs)[!(names(dotArgs) %in% 
            c(names(formals(fnSubset)), names(formals(fnPredict))))]), 
            collapse = ","), " is/are not used by fnSubset() or fnPredict()")
    
    dotArgs_fnSubset <- dotArgs[names(dotArgs) %in% names(formals(fnSubset))]
    dotArgs_fnPredict <- dotArgs[names(dotArgs) %in% names(formals(fnPredict))]
    if (nPredict == 1) 
        fill <- data
    else {
        fill <- array(NA, c(dim(data), nPredict))
        fill[,,,,1] <- data
    }
    time <- list(start =  Sys.time(), end = NULL,
                 elapsedMins = NULL, elapsedSecsPerNA = NULL)
    if (verbose) {
        verb_data <- length(data)
        verb_observed <- sum(!is.na(data))
        verb_observedP <- round(verb_observed/verb_data * 100, 1)
        verb_missing <- verb_data - verb_observed
        verb_missingP <- round(verb_missing/verb_data * 100, 1)
        verb_topred <- length(mps)
        verb_topredP <- round(verb_topred/verb_data * 100, 1)
        message(sprintf(gsub("X", nchar(verb_data),
                             "data has %Xi values: %Xi (%2.1f%%) observed,\n                      %Xi (%2.1f%%) missing,\n                      %Xi (%2.1f%%) to predict."), 
            verb_data, verb_observed, verb_observedP, verb_missing, 
            verb_missingP, verb_topred, verb_topredP))
        message("started at ", time$start, ".")
    }
    mcall <- match.call()
    if (!dopar) {
        for (i_pixel in mps) {
            mp <- IndexOneFour(i_pixel, dim(data))
            i <- 0L
            a <- do.call(fnSubset,
                         c(list(data = data, mp = mp, i = i),
                           dotArgs_fnSubset))
            z <- do.call(fnPredict, c(list(a = a, i = i), dotArgs_fnPredict))
            while (is.na(z[1]) && i < iMax) {
                i <- i + 1L
                aNew <- do.call(fnSubset, c(list(data = data, mp = mp, i = i),
                                            dotArgs_fnSubset))
                if (identical(dim(aNew), dim(a))) 
                  break
                a <- aNew
                z <- do.call(fnPredict, c(list(a = a, i = i),
                                          dotArgs_fnPredict))
            }
            if (nPredict == 1) 
                fill[i_pixel] <- z[1]
            else
                fill[mp[1], mp[2], mp[3], mp[4],] <- z[1:nPredict]
        }
    }
    else {
        foreach_fill <- foreach::`%dopar%`(foreach::foreach(i_pixel = mps, 
                                                            .combine = rbind), {
            mp <- IndexOneFour(i_pixel, dim(data))
            i <- 0L
            a <- do.call(fnSubset, c(list(data = data, mp = mp, i = i),
                                     dotArgs_fnSubset))
            z <- do.call(fnPredict, c(list(a = a, i = i), dotArgs_fnPredict))
            while (is.na(z[1]) && i < iMax) {
                i <- i + 1L
                aNew <- do.call(fnSubset, c(list(data = data, mp = mp, i = i),
                                            dotArgs_fnSubset))
                if (identical(dim(aNew), dim(a))) 
                  break
                a <- aNew
                z <- do.call(fnPredict, c(list(a = a, i = i), 
                                          dotArgs_fnPredict))
            }
            z[1:nPredict]
        })
        if (nPredict == 1) 
            fill[mps] <- foreach_fill[,1]
        else {
            id <- c(array(mps, c(length(mps), nPredict)) +
                  array(rep(seq(0, length(data) * (nPredict - 1), length(data)), 
                            each = length(mps)), c(length(mps), nPredict)))
            fill[id] <- foreach_fill[, 1:nPredict]
        }
    }
    if (!identical(clipRange, c(-Inf, Inf))) {
        if (nPredict == 1) {
            fill[fill < clipRange[1]] <- clipRange[1]
            fill[fill > clipRange[2]] <- clipRange[2]
        } else {
            fill[which(fill[,,,,1] < clipRange[1])] <- clipRange[1]
            fill[which(fill[,,,,1] > clipRange[2])] <- clipRange[2]
        }
    }
    time$end <- Sys.time()
    time$elapsedMins <- difftime(time$end, time$start, units = "mins")
    time$elapsedSecsPerNA <- difftime(time$end, time$start, units = "sec")/length(mps)
    if (verbose) 
        message("elapsed time is ", round(time$elapsedMins, 3), 
                " mins (", round(time$elapsedSecsPerNA, 5), " secs per NA).")
    list(fill = fill, mps = mps, time = time, call = mcall)
}


#' @name Extend
#' @aliases extend alternative
#' @rdname Extend
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Implement an Alternative Gap-fill Algorithm
#' @description
#' By default, the \code{\link{Gapfill}} function uses the \code{\link{Subset}}
#' and \code{\link{Predict}} functions to predict missing values.
#' To implement alternative gap-fill procedures, these functions can be replaced
#' by user defined ones and passed to the \code{\link{Gapfill}} function via the arguments
#' \code{fnSubset} and \code{fnPredict}.\cr 
#' The example section below gives two such extensions:
#' \describe{\item{Example 1: }{Illustration of the concept.
#' The prediction is the mean of the subset around a missing value.}
#' \item{Example 2: }{An algorithm using the \code{\link{Score}} and the \code{\link[stats]{lm}} functions.}
#' }
#' 
#' @details
#' To work properly the user-defined \code{Subset} function needs to have the arguments:
#' \describe{\item{\code{data}: }{The input data array.}
#' \item{\code{mp}: }{Numeric vector of length 4 specifying the index of the
#' currently treated missing value.}
#' \item{\code{i}: }{Integer vector of length 1. Number of non-successfully tried subsets.}
#' }
#' The function user-defined \code{\link{Predict}} function, needs to have the arguments:
#' \describe{
#' \item{\code{a}: }{Return value of the \code{Subset} function.}
#' \item{\code{i}: }{Integer vector of length 1. Number of non-successfully tried subsets.}
#' }
#' Both functions may take additional arguments.
#' The default values of these arguments can be changed via
#' the \code{...} arguments of \code{\link{Gapfill}}.
#'
#' @references F. Gerber, R. Furrer, G. Schaepman-Strub, R. de Jong, M. E. Schaepman, 2016,
#' Predicting missing values in spatio-temporal satellite data.
#' \url{http://arxiv.org/abs/1605.01038}.
#' 
#' @seealso \code{\link{Gapfill}}, \code{\link{Subset-Predict}}, \code{\link{Score}}, \code{\link[stats]{lm}}.
#' @examples
#' \dontrun{
#' ## Example 1: mean ----------------------------------
#' ## define a predict function
#' PredictMean <- function (a, i) mean(a, na.rm = TRUE)
#'
#' out1 <- Gapfill(data = ndvi, fnPredict = PredictMean)
#' Image(out1$fill)
#'
#' ## start with a smaller subset
#' args(Subset)
#' out2 <- Gapfill(data = ndvi, fnPredict = PredictMean,
#'                 initialSize = c(0, 0, 1, 6))
#' Image(out2$fill)
#'
#' ## require at least "nNotNA" non-NA values
#' ## return predicted value and number of iterations i
#' PredictMean2 <- function (a, i, nNotNA) {
#'     if (sum(!is.na(a)) < nNotNA)
#'         return (c(NA, NA))
#'     c(mean(a, na.rm = TRUE), i)
#' }
#' out3 <- Gapfill(data = ndvi, fnPredict = PredictMean2, nPredict = 2,
#'                 initialSize = c(0, 0, 1, 6), nNotNA = 0)
#' stopifnot(identical(c(out2$fill), c(out3$fill[,,,,1])))
#' Image(out3$fill[,,,,2])  # number of used iterations i
#'
#' out4 <- Gapfill(data = ndvi, fnPredict = PredictMean2, nPredict = 2,
#'                 initialSize = c(0, 0, 1, 6), nNotNA = 50)
#' Image(out4$fill[,,,,1])  # fill values
#' Image(out4$fill[,,,,2])  # number of used iterations i
#'
#'
#' ## Example 2: Score() and lm() ----------------------
#' PredictLm <- function (a, i, nNotNA = 50, minScores = 2){
#'     if (sum(!is.na(a)) < nNotNA)
#'         return (NA)
#'     am <- Array2Matrix(a)
#'     sx <- Score(t(am))
#'     lsx <- length(sx)
#'     if (lsx < minScores)
#'         return (NA)
#'     sy <- Score(am)
#'     lsy <- unique(length(sy))
#'     if (lsy < minScores)
#'         return (NA)
#'     df <- data.frame(z = c(am),
#'                      sx = rep(sx, ncol(am)),
#'                      sy = rep(sy, each = nrow(am)))
#'     newdata <- df[IndexTwoOne(attr(am, "mp"), dim(am)),]
#'     m <- lm(z ~ sx * sy, data = df)
#'     predict(m, newdata = newdata)
#' }
#' 
#' ## test PredictLm() by running it
#' ## manually for one missing value
#' mp <- IndexOneFour(which(is.na(ndvi))[1], dim(ndvi))
#' a <- Subset(data = ndvi, mp = mp, i = 0)
#' PredictLm(a = a, i = 0)
#'
#' ## run PredictLm() on ndvi data
#' out5 <- Gapfill(data = ndvi, fnPredict = PredictLm,
#'                 nNotNA = 50)
#' Image(out5$fill)
#' }
NULL

#' @name Subset-Predict
#' @aliases Subset Predict fnSubset fnPredict
#' @rdname Subset-Predict
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Subset and Predict Functions
#' @description
#' The \code{Subset} and \code{Predict} function used in the default configuration of \code{\link{Gapfill}}.
#' To predict a missing value, the two function are called sequentially as described the help page of \code{\link{Gapfill}}.
#' 
#' @param data Numeric array with four dimensions. The input (satellite) data to be gap-filled.
#' Missing values should be encoded as \code{NA}.
#' The data should have the dimensions: x coordinate, y coordinate, seasonal index (e.g., day of the year), and year.
#' See the \code{ndvi} dataset for an example. 
#' @param mp Integer vector of length 4 encoding the position of the missing value in \code{data} to predict.
#' @param i Integer vector of length 1. The number of tried subsets that lead to a \code{NA} return value from \code{Predict}.
#' @param initialSize Integer vector of length 4, that provides the size of the subset for \code{i = 0}.
#' @return \code{Subset} returns an array with 4 dimensions containing the missing value
#' at the position indicated by the attribute \code{mp}.
#' @export
Subset <- function(data, mp, i,
                   initialSize = c(10L, 10L, 1L, 5L)){
    ArrayAround(data = data, mp = mp,
                size = initialSize + c(i, i, 0L, 0L))
}

#' @rdname Subset-Predict
#' @param a Return value of \code{Subset()}.
#' @param nTargetImage Integer vector of length 1. Minimum number of non-NA values in the image containing the missing value.
#' If the criterion is not met, \code{NA} is returned. 
#' @param nImages Integer vector of length 1. Minimum number of non-empty images.
#' If the criterion is not met, \code{NA} is returned.
#' @param nQuant Integer vector of length 1. Parameter passed to \code{\link{EstimateQuantile}}.
#' @param predictionInterval Logical vector of length 1.
#' If \code{TRUE}, the predicted value together with the lower and upper bounds
#' of an approximated 90\% prediction interval are returned.
#' When \code{predictionInterval} \code{= TRUE}, the function returns 3 values, and hence,
#' the argument \code{nPredict} of \code{\link{gapfill}} has to be set to 3 in order to store all returned values. 
#' @return \code{Predict} returns a numeric vector containing the predicted value
#' (and if \code{predictionInterval} is \code{TRUE}, the lower and upper bounds of the prediction interval),
#' or \code{NA}, if no prediction was feasible. 
#' @seealso \code{\link{Gapfill}}, \code{\link{Extend}},
#' \code{\link{EstimateQuantile}}, \code{\link{Score}}, \code{\link{ndvi}}.
#' @details
#' The \code{Subset} function defines the search strategy to find a
#' relevant subset by calling the function \code{\link{ArrayAround}}.
#' The size of the initial subset is given by the argument \code{initialSize}.
#' Its default values is \code{c(5L, 5L, 1L, 5L)}, which corresponds to a spatial extend of 5 pixels
#' in each direction from the missing value and includes time points having the previous, the same or the next seasonal index and
#' are not further apart than 5 years.
#' With an increase of the argument \code{i}, the spatial extent of the subset increases.
#' 
#' The \code{Predict} function decides whether the subset \code{a} is suitable and
#' calculates the prediction (fill value) when a suitable subset is provided.
#' To formulate the conditions that are used to decide if a subset is suitable,
#' consider the subset \code{a} as a collection of images.
#' More precisely, if \code{dim(a)} \code{=} \code{c(d1, d2, d3, d4)},
#' it can be seen as a collection of \code{d3*d4} images with an extent of \code{d1} by \code{d2} pixels.    
#' Using this terminology, we require the following conditions to be fulfilled
#' in order to predict the missing value:
#' \itemize{
#' \item \code{a} contains at least \code{nTargetImage} non-NA values in the image containing the missing value,
#' \item \code{a} contains at least \code{nImages} non-empty images.
#' }
#' The prediction itself is based on sorting procedures (see \code{\link{Score}} and
#' \code{\link{EstimateQuantile}}) and the quantile regression function \code{\link[quantreg]{rq}}.
#' 
#' If the argument \code{predictionInterval} is \code{TRUE} the \code{Predict} functions returns
#' the predicted value together with the lower and upper bounds of an approximated 90\% prediction interval.
#' The interval combines the uncertainties introduced by \code{\link{Score}}
#' and \code{\link{EstimateQuantile}}.
#' @references F. Gerber, R. Furrer, G. Schaepman-Strub, R. de Jong, M. E. Schaepman, 2016,
#' Predicting missing values in spatio-temporal satellite data.
#' \url{http://arxiv.org/abs/1605.01038}. 
#' @note
#' The current implementation of \code{Subset} does not take into account
#' that locations at the boundary of \code{data} can be neighboring to each other.
#' For example, if global data (entire sphere) are considered, the location
#' \code{data[1,1,,]} is a neighbor of \code{data[dim(data)[1], dim(data)[2],,]}.
#' Similar considerations apply when data are available for an entire year. 
#' To take this into account, the \code{Subset} function can be redefined accordingly or
#' the data can be augmented.
#' @examples
#' ## Assume we choose c(5, 5, 1, 5) as initalSize of the subset
#' iS <- c(5, 5, 1, 5)
#' ## case 1: initial subset leads to prediction -------
#' i <- 0
#' a <- Subset(data = ndvi, mp = c(1, 3, 1, 2), i = i, initialSize = iS)
#' p <- Predict(a = a, i = i)
#' p
#' stopifnot(identical(a, ArrayAround(data = ndvi, mp = c(1, 3, 1, 2),
#'                                    size = c(5 + i, 5 + i, 1, 5))))
#' stopifnot(identical(p, Gapfill(data = ndvi, subset = 1807,
#'                                initialSize = iS, verbose = FALSE)$fill[1807]))
#'
#' ## case 2: two tries are necessary ------------------
#' i <- 0
#' a <- Subset(data = ndvi, mp = c(20, 1, 1, 2), i = i, initialSize = iS)
#' p <- Predict(a = a, i = i)
#' p
#'
#' ## Increase i and try again.
#' i <- i + 1
#' a <- Subset(data = ndvi, mp = c(20, 1, 1, 2), i = i, initialSize = iS)
#' p <- Predict(a = a, i = i)
#' p
#' stopifnot(identical(a, ArrayAround(data = ndvi, mp = c(20, 1, 1, 2),
#'                                    size = c(5 + i, 5 + i, 1, 6))))
#' stopifnot(identical(p, Gapfill(data = ndvi, subset = 1784,
#'                                initialSize = iS, verbose = FALSE)$fill[1784]))
#' @export
#' @importFrom quantreg rq coef.crq predict.rq
#' @importFrom stats quantile
Predict <- function(a,
                    i,
                    nTargetImage = 5,
                    nImages = 4,
                    nQuant = 2,
                    predictionInterval = FALSE){

    ## arrange data in matrix
    am <- Array2Matrix(a)
    mp <- attr(am, "mp")
     
    ## require 20 non missing pixel in target image
    if(sum(!is.na(am[,mp[2]])) < nTargetImage[1])
        return(NA)
    ## require 5 non missing images
    if(sum(apply(!is.na(am), 2, any)) < nImages[1])
       return(NA)

    tau <- EstimateQuantile(a = a, mp = attr(a, "mp"), nQuant = nQuant[1],
                            predictionInterval = predictionInterval[1])
    if(any(is.na(tau)))
        return(NA)

    s <- Score(mat = am)
    r <- rank(s, na.last = "keep")
    
    df <- data.frame(z = c(am),
                     rank = rep(rank(s, na.last = "keep"),
                                each = dim(am)[1]))
    m <- quantreg::rq(z ~ rank, tau[1], data = df)
    p <- quantreg::coef.crq(m)[1] +
         quantreg::coef.crq(m)[2] * r[mp[2]]
    p <- unname(p)

    if(!predictionInterval[1])
        return(p)

    mu <- quantreg::rq(z ~ rank, max(tau), data = df)
    ml <- quantreg::rq(z ~ rank, min(tau), data = df)

    if(r[mp[2]] == max(df$rank, na.rm = TRUE)){
        pu <- quantreg::predict.rq(mu, newdata = data.frame(rank = r[mp[2]]))
    } else {
        pu <- stats::quantile(quantreg::predict.rq(mu), .95, na.rm = TRUE)
    }
    if(r[mp[2]] == min(df$rank, na.rm = TRUE)){
        pl <- quantreg::predict.rq(ml, newdata = data.frame(rank = r[mp[2]]))
    } else {
        pl <- stats::quantile(quantreg::predict.rq(ml), .05, na.rm = TRUE)    
    }
    pu <- max(p, pu)
    pl <- min(p, pl)
    unname(c(p, sort(c(pl, pu))))
}

#' @name ArrayAround
#' @rdname ArrayAround
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Subset an Array with 4 dimensions
#' @description Given an array \code{data} with 4 dimensions,
#' a subset around the element with coordinates \code{mp} ("missing position") is extracted.
#' The size of the subset in all four directions
#' from \code{mp} is specified by \code{size}. \cr
#' \code{ArrayAroundRandom} returns a subset around a
#' random location in \code{data}.
#' @param data Array with 4 dimensions. 
#' @param mp Integer vector of length 4
#' indexing an element in \code{data}. 
#' @param size Integer vector of length 4, that provides
#' the size of the subset in all four dimensions
#' (around \code{mp}).
#' @return Array with 4 dimensions corresponding to the specified subset.
#' The attribute \code{mp} of the returned array is an integer vector
#' of length 4 giving \code{mp} relative to the
#' returned array. 
#' @note
#' When \code{size = c(0, 0, 0, 0)}, the returned subset consists of one value
#' (the value of \code{data} indexed with \code{mp}.) \cr
#' When \code{mp} is near the boundaries of \code{data},
#' the returned subset may be smaller than indicated by the argument \code{size}
#' and the attribute \code{mp} may indicate an element near the boundaries of the subset.
#' @examples
#' a <- array(1:16, c(2, 2, 2, 2))
#' ArrayAround(data = a, mp = c(1, 1, 1, 1), size = c(0, 0, 0, 0))
#' ## returns the first element a[1,1,1,1]
#'
#' ArrayAround(data = a, mp = c(2, 2, 2, 2), size = c(0, 0, 0, 0))
#' ## returns the last element a[2,2,2,2]
#'
#' ArrayAround(data = a, mp = c(1, 1, 1, 1), size = c(1, 0, 0, 0))
#' ## returns a[1:2,1,1,1]
#'
#' ArrayAround(data = a, mp = c(1, 1, 1, 1), size = c(1, 1, 1, 1))
#' ## returns a
#' 
#' @export
ArrayAround <- function(data, mp, size){
    stopifnot(is.array(data))
    stopifnot(identical(length(dim(data)), 4L))

    stopifnot(is.numeric(mp))
    stopifnot(identical(length(mp), 4L))
    stopifnot(all(1L <= mp) && all(mp <= dim(data)))

    stopifnot(is.numeric(size))
    stopifnot(identical(length(size), 4L))
    stopifnot(all(0L <= size))

    d1Lo <- max(1L, mp[1] - size[1])
    d1Up <- min(dim(data)[1], mp[1] + size[1])
    d2Lo <- max(1L, mp[2] - size[2])
    d2Up <- min(dim(data)[2], mp[2] + size[2])
    d3Lo <- max(1L, mp[3] - size[3])
    d3Up <- min(dim(data)[3], mp[3] + size[3])
    d4Lo <- max(1L, mp[4] - size[4])
    d4Up <- min(dim(data)[4], mp[4] + size[4])

    ## get box around missing point
    data <- data[d1Lo:d1Up, d2Lo:d2Up,
                 d3Lo:d3Up, d4Lo:d4Up, drop=FALSE]
    attr(data, "mp") <- c(min(size[1] + 1L, mp[1]),
                          min(size[2] + 1L, mp[2]),
                          min(size[3] + 1L, mp[3]),
                          min(size[4] + 1L, mp[4]))
    data
}

#' @name ArrayAroundRandom
#' @rdname ArrayAround
#' @param target One of 
#' \code{c("all", "missing", "observed")}.
#' Indicates from which subset of \code{data} a random
#' location is sampled.
#' @param verbose Logical vector of length 1.
#' If \code{TRUE}, messages are printed.
#' 
#' @examples
#' 
#' ArrayAroundRandom(a)
#' ArrayAroundRandom(a, size = c(1, 2, 1, 2))
#' @export
ArrayAroundRandom <- function(data, size = c(0L, 0L, 0L, 0L),
                              target = c("all", "missing", "observed"),
                              verbose = TRUE){

    stopifnot(is.array(data))
    stopifnot(identical(length(dim(data)), 4L))
    ## 'size' used and tested in ArrayAround()
    target <- match.arg(target)

    index <- switch(target,
                    all = 1:length(data),
                    missing = which(is.na(c(data))),
                    observed = which(!is.na(c(data))))
    if(length(index) < 1)
        stop("No \"", target, "\" values in \"data\".") 
    mp <- IndexOneFour(index[sample(1:length(index), 1)],
                       dim(data))
    if(verbose)
        message("mp = c(", mp[1],", ",
                mp[1],", ", mp[1],", ", mp[1],")")
    ArrayAround(data = data, mp = mp, size = size)
}

#' @name Array2Matrix
#' @rdname Array2Matrix
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Convert an Array with 4 Dimensions into a Matrix
#' @description Converts the array, \code{a}, with 4 dimensions, \code{c(d1, d2, d3, d4)},
#' into a matrix with \code{d1*d2} rows and \code{d3*d4} columns.
#'
#' @param a Array with 4 dimensions. 
#' @return A matrix. If \code{a} has the attribute \code{mp}, the transformed attribute is returned as well.
#' See \code{\link{ArrayAround}} for more information about \code{mp}.
#' @seealso \link{Index}, \code{\link{ArrayAround}}.
#' @examples
#' a <- array(data = 1:16, dim = c(2, 2, 2, 2))
#' Array2Matrix(a = a)
#' attr(a, "mp") <- c(1, 2, 2, 1)
#' Array2Matrix(a = a)
#'
#' Array2Matrix(ArrayAround(data = a, mp = c(1, 1, 1, 1),
#'                          size = c(1, 1, 2, 2))) 
#' @export
Array2Matrix <- function(a){
    stopifnot(length(dim(a)) == 4)

    if(!is.null(attr(a, "mp"))){
        mp <- attr(a, "mp")
        attr(a, "mp") <- c(IndexTwoOne(attr(a, "mp")[1:2],
                                       dim(a)[1:2]),
                           IndexTwoOne(attr(a, "mp")[3:4],
                                       dim(a)[3:4]))
    }
    dim(a) <- c(prod(dim(a)[1:2]), prod(dim(a)[3:4]))
    class(a) <- "matrix"
    a
}


#' @name Index
#' @aliases IndexTwoOne
#' @rdname Index
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Index Conversions
#' @description
#' Converts an index from the first length to the second length.
#' For example, assume that \code{c(2, 2)} indexes an element in a matrix with 2 rows and 5 columns.
#' If the matrix is transformed to a vector,
#' the same element can be accessed with the index \code{IndexTwoOne(c(2, 2), c(2, 5))} (=4). 
#' @param index Integer vector of length 2 for
#' \code{IndexTwoOne()} and length 1 for \code{IndexOneFour()}.
#' @param dimTwo Integer vector of length 2
#' indicating the dimension of the 2 dimensional array.
#' @return Vector of length 1 for \code{IndexTwoOne()} and
#' length 4 for \code{IndexOneFour()}. 
#' @examples
#' ## IndexTwoOne
#' IndexTwoOne(c(2, 2), c(2, 5))
#' v <- 1:10
#' dimTwo <- c(2, 5)
#' m <- array(v, dimTwo)
#' stopifnot(v[IndexTwoOne(c(2, 2), dimTwo)] == m[2,2])
#' 
#' @export
IndexTwoOne <- function(index, dimTwo){
    stopifnot(is.numeric(index))
    stopifnot(length(index) == 2)
    index <- as.integer(index)
        
    stopifnot(is.numeric(dimTwo))
    stopifnot(length(dimTwo) == 2)
    dimTwo <- as.integer(dimTwo)
    
    (index[2] - 1) * dimTwo[1] + index[1]
}

#' @name Index
#' @aliases IndexOneFour
#' @rdname Index
#' @param dimFour Integer vector of length 4
#' indicating the dimension of the 4 dimensional array. 
#' @examples
#'  
#' ## IndexOneFour
#' IndexOneFour(13, c(2, 2, 2, 2))
#' w <- 1:16
#' dimFour <- c(2, 2, 2, 2)
#' a <- array(w, dimFour)
#' stopifnot(a[1,1,2,2] == w[13])
#' 
#' @export
IndexOneFour <- function(index, dimFour){
    stopifnot(is.numeric(index))
    stopifnot(length(index) == 1)
    index <- as.integer(index)
    stopifnot(is.numeric(dimFour))
    stopifnot(length(dimFour) == 4)
    dimFour <- as.integer(dimFour)
    
    index <- index - 1L
    r1 <- index %/% dimFour[1]
    r2 <- r1 %/% dimFour[2]
    r3 <- r2 %/% dimFour[3]
    c(index %% dimFour[1], r1 %% dimFour[2],
      r2 %% dimFour[3], r3) + 1L
}

#' @name Score
#' @title Score Columns of a Matrix Containing NAs by its Values
#' @description Helper function for \code{\link{Predict}} used to
#' score the columns of a matrix according to their values.  
#' The scoring of a given column is done by pair-wise comparisons with all other columns.
#' The comparison of columns is done by pair-wise comparisons of the non-missing values.
#' This procedure is robust to missing values, if all columns of the matrix
#' have a similar (potentially shifted) distribution of values.
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @param mat Numeric matrix. May contain \code{NA} values.
#' @return Numeric vector of length \code{ncol(mat)}.
#' @note Interfaces a C++ function. The R package Rcpp is used.  
#' @references F. Gerber, R. Furrer, G. Schaepman-Strub, R. de Jong, M. E. Schaepman, 2016,
#' Predicting missing values in spatio-temporal satellite data.
#' \url{http://arxiv.org/abs/1605.01038}. 
#' @importFrom Rcpp evalCpp
#' @useDynLib gapfill
#' @examples
#' mat <- rbind(c( 1,  2, NA),
#'              c(NA, NA,  1),
#'              c( 2, NA,  3),
#'              c( 1,  5, NA),
#'              c(NA,  2,  5))
#' s <- Score(mat)
#'
#' ## manual calculation in R
#' Mean <- function(x) mean(x, na.rm = TRUE)
#' sByHand <- c(Mean(c(Mean(mat[,1] > mat[,2]),
#'                     Mean(mat[,1] > mat[,3]))),
#'              Mean(c(Mean(mat[,2] > mat[,1]),
#'                     Mean(mat[,2] > mat[,3]))),
#'              Mean(c(Mean(mat[,3] > mat[,1]),
#'                     Mean(mat[,3] > mat[,2]))))
#' stopifnot(identical(s, sByHand))
#' 
#' @export
Score <- function(mat) {
    .Call('gapfill_Score_cpp', PACKAGE = 'gapfill', mat)
}

#' @name EstimateQuantile
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Estimate the Quantile of a Missing Value
#' @description Helper function for \code{\link{Predict}}.
#' The function estimates the quantile of the missing value at position \code{mp} from the
#' data \code{a} relative to its image a[,,mp[3],mp[4]].
#' @param a Numeric array with 4 dimensions.
#' @param mp Integer vector of length 4 indexing the position of the
#' missing value to predict.
#' @param nQuant Integer vector of length 1. Minimum number of non-missing values in \code{a[mp[1], mp[2],,]} required to estimate the quantile.
#' If \code{a[mp[1], mp[2],,]} contains less non-missing values, the neighboring values of \code{a[mp[1], mp[2],,]} are also taken into account.
#' @param predictionInterval Logical vector of length 1.
#' If \code{TRUE}, the estimated quantile together with lower and upper bounds of an approximate 90\% uncertainty interval
#' is returned.  
#' @return If \code{predictionInterval} is \code{FALSE}, a numeric vector of length 1 being the estimated quantile of the missing value
#' \code{a[mp[1], mp[2], mp[3], mp[4]]} is returned.
#' Otherwise, a numeric vector of length 3 containing the estimated quantile and the lower and upper bounds of an
#' approximate 90\% uncertainty interval is returned.
#' @seealso \code{\link{Predict}}.
#' @references F. Gerber, R. Furrer, G. Schaepman-Strub, R. de Jong, M. E. Schaepman, 2016,
#' Predicting missing values in spatio-temporal satellite data.
#' \url{http://arxiv.org/abs/1605.01038}.
#' @examples
#' a <- Subset(data = ndvi, mp = c(1, 3, 1, 2), i = 0)
#' EstimateQuantile(a = a, mp = attr(a, "mp"), nQuant = 2)
#'
#' @importFrom stats ecdf
#' @export
EstimateQuantile <- function (a, mp, nQuant, predictionInterval = FALSE) {

    ref <- a[mp[1], mp[2],,,drop = FALSE]
    am <- Array2Matrix(a)
    i <- 0
    while (sum(!is.na(ref)) < nQuant[1] && i <= max(dim(a)[1:2])) {
        i <- i+1
        ref <- ArrayAround(a, mp, c(i, i, Inf, Inf))
    }
    if (all(is.na(ref)))
        return (NA)
    dim(ref) <- c(prod(dim(ref)[1:2]), prod(dim(ref)[3:4]))
    tmp <- array(NA, dim(ref))
    for(it in 1:dim(am)[2]) {
        if (all(is.na(am[,it])))
            tmp[,it] <- NA
        else
            tmp[,it] <- ecdf(am[,it])(ref[,it])
    }
    tt <- apply(tmp, 2, mean, na.rm = TRUE)
    mu <- mean(tt, na.rm = TRUE)
    if(!predictionInterval[1])
        return(mu)
    interval <- unname(quantile(tt, prob = c(0.05, 0.95), na.rm = TRUE))
    c("mu" = mu, "lower" = interval[1], "upper" = interval[2])
}

#' @name Validate
#' @aliases validate
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Validation with RMSE
#' @description The function summarizes the validation scenario and
#' returns the root mean squared error (RMSE) of the predictions.
#' The typical validation procedure is: start with the \code{trueData}.
#' Remove some validation points to obtain artificially generated \code{dataObserved}.
#' Predicting the validation points based on \code{dataObserved} leads to \code{dataFilled}.  
#' @param dataObserved Numeric vector containing the observed data.
#' @param dataFilled Numeric vector containing the filled (predicted) data.
#' Needs to have the same length as \code{dataObserved}. 
#' @param dataTrue Numeric vector containing the true data.
#' Needs to have the same length as \code{dataObserved}. 
#' @param include Logical vector indicating which element to include in the
#' calculation of the RMSE.
#' @return Numeric matrix with one 1 row and 6 columns having the entries:
#' \itemize{
#' \item \code{nNA}: number of missing values in \code{dataObserved},
#' \item \code{nFilled}: number of predicted values,
#' \item \code{nNotFilled}: number of not predicted missing values,
#' \item \code{ratioFilled}: ratio: \code{nFilled} / \code{nNA},
#' \item \code{nCrossvali}: number of values for validation,
#' \item \code{RMSE}: root mean squared error.
#' }
#' @examples
#' Validate(c(1, NA, 2, NA), c(1, 2, 2, NA), c(1, 1, 2, 2))
#'
#' ## validate gap-fill predictions: consider the ndvi data
#' Image(ndvi)
#'
#' ## define some validation points vp
#' ## in the image of the day 145 of the year 2004
#' vp <- 300 + c(5:10) + rep(21 * c(0:5), each = 6)
#'
#' ## remove the vp values from the data
#' nn <- ndvi
#' nn[vp] <- NA
#' Image(nn)
#'
#' ## predict the vp values 
#' out <- Gapfill(nn, subset = vp)
#' Validate(dataObserved = nn, dataFilled = out$fill,
#'          dataTrue = ndvi)
#' 
#' @export
Validate <- function (dataObserved,
                      dataFilled,
                      dataTrue,
                      include = rep(TRUE, length(dataObserved))) {
    
    stopifnot(identical(length(dataObserved), length(dataFilled)))
    stopifnot(identical(length(dataObserved), length(dataTrue)))
    stopifnot(identical(length(dataObserved), length(include)))
    stopifnot(is.numeric(dataObserved))
    stopifnot(is.numeric(dataFilled))
    stopifnot(is.numeric(dataTrue))
    stopifnot(is.logical(include))
    
    vp <-  is.na(dataObserved) & !is.na(dataFilled) &
           !is.na(dataTrue) & include
    diff <- dataFilled - dataTrue
    diff[!vp] <- NA

    nNA <- sum(is.na(dataObserved))
    nFilled <- sum(is.na(dataObserved) & !is.na(dataFilled))
    nNotFilled <- nNA - nFilled
    ratioFilled <- nFilled / nNA
    nCrossvali <- sum(vp)
 
    cbind(nNA         = nNA,
          nFilled     = nFilled,
          nNotFilled  = nNotFilled,
          ratioFilled = ratioFilled,
          nCrossvali  = nCrossvali,
          RMSE        = sqrt(mean(diff^2, na.rm = TRUE)))
}

#' @name Image 
#' @author Florian Gerber, \email{florian.gerber@@math.uzh.ch}.
#' @title Image Panels
#' @description
#' Creates an image panel to visualize data in 4, 3 or 2 dimensional arrays (e.g., space-time data). 
#' The function returns a \code{ggplot2} object, which 
#' can be modified using ggplot2 (and/or grid) syntax. 
#'
#' @param x Numeric array with 4, 3, or 2 dimensions
#' containing the data to be plotted.
#' @param zlim Numeric vector of length 2.
#' Gives the upper and lower bound of the plotted values. 
#' @param col Vector of colors. 
#' @param na.value Vector of length one.
#' The color to be used for NA values. 
#' @param byrow Logical vector of length one.
#' Indicates the ordering of the panels.
#' @param xlab Character vector (or expression) of length one giving the x-axis label.
#' @param ylab Character vector (or expression) of length one giving the y-axis label.
#' @param colbarTitle Character vector (or expression) of length one giving the colorbar label.
#' @param ... Additional arguments are passed to \code{ggplot}.
#' @return Object (plot) of class \code{c("gg", "ggplot2")}.
#' @seealso \code{\link{ndvi}}, \code{\link[ggplot2]{ggplot}}.
#' @examples
#' Image(ndvi)
#' 
#' p1 <- Image(ndvi, colbarTitle = "NDVI", byrow  = FALSE)
#' p1
#' 
#' p2 <- Image(ndvi[,,3,2], na.value = "white", colbarTitle = "NDVI") +
#'       theme(strip.text.x = element_blank(),
#'             strip.text.y = element_blank(),
#'             panel.border = element_rect(fill = NA, size = 1))
#' p2
#' 
#' ## place modified color bar left
#' p2 + guides(fill = guide_colorbar(title = "NDVI", 
#'                                   barwidth = 1,
#'                                   barheight = 20,
#'                                   label.position = "left", 
#'                                   legend.position = c(0, 0))) +
#'      theme(legend.position = "left")
#' 
#' ## place color bar at bottom
#' p2 + guides(fill = guide_colorbar(title = "NDVI", 
#'                                   barwidth = 7,
#'                                   barheight = .7,
#'                                   label.position = "bottom", 
#'                                   legend.position = c(0, 0)),
#'                                   direction = "horizontal") +
#'      theme(legend.position = "bottom")
#'
#' @importFrom fields tim.colors
#' @import ggplot2
#' @export
Image <- function(x = NULL,
                  zlim = range(x, na.rm = TRUE),
                  col = fields::tim.colors(1000),
                  na.value = "black",
                  byrow = TRUE,
                  xlab = "",
                  ylab = "",
                  colbarTitle = "",
                  ...){

    stopifnot(2 <= length(dim(x)) && length(dim(x)) <= 4)
    if(length(dim(x)) == 2)
        dim(x) <- c(dim(x), 1, 1)
    if(length(dim(x)) == 3)
        dim(x) <- c(dim(x), 1)
    if(is.logical(x))
        x <- array(as.numeric(x), dim(x))
    stopifnot(is.numeric(x))
    
    stopifnot(is.numeric(zlim))
    stopifnot(identical(length(zlim), 2L))

    na.value <- na.value[1]

    stopifnot(is.logical(byrow))
    byrow <- byrow[1]
    
    stopifnot(is.character(xlab))
    xlab <- xlab[1]
    
    stopifnot(is.character(ylab))
    ylab <- ylab[1]

    stopifnot(is.character(colbarTitle))
    colbarTitle <- colbarTitle[1]
              
    if(!byrow)
        x <- aperm(x, c(1,2,4,3))
    
    if(is.null(dimnames(x)))
        dimnames(x) <- list(1:dim(x)[1], 1:dim(x)[2],
                            1:dim(x)[3], 1:dim(x)[4])

    plotData <- expand.grid(x = 1:dim(x)[1],
                            y = 1:dim(x)[2],
                            doy = dimnames(x)[[3]],
                            year = dimnames(x)[[4]])
    plotData$z <- c(x)
    
    p <- ggplot(data = plotData,
                mapping = aes_string(x = "x", y = "y", fill = "z"),
                ...) +
                   geom_tile() +
        theme_bw() + xlab(xlab) + ylab(ylab) +
        scale_fill_gradientn(colors = col,
                             limits = zlim, na.value = na.value) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        coord_fixed() + facet_grid(year ~ doy) +
        theme(line = element_blank(),
              strip.background = element_rect(fill = "white",
                                              color = "white"),
              strip.text.x = element_text(),
              strip.text.y = element_text(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              line = element_blank(),
              legend.position = "right",
              panel.border = element_blank(),
              plot.margin = unit(c(1,1,1,1), "mm"),
              panel.margin = unit(2, "mm")) +
    guides(fill = guide_colorbar(title = colbarTitle,
                                 barwidth = 1,
                                 barheight = 10,
                                 label.position = "right",
                                 legend.position=c(0,0))) 
    p
}


#' @name ndvi
#' @title NDVI Data from Alaska
#' @description
#' The dataset was created to test gap-fill algorithms.
#' It mimics a subset of the MODIS NDVI data (product MOD13A1) in the region of Alaska.
#' The data product features one image per 16-day time interval, i.e., 24 images per year.  
#' The indicated images (see \code{Image(ndvi)}) were downloaded and stored as a 4 dimensional array.
#' Its dimensions correspond to longitude, latitude, day of the year, and year.
#' 
#' @docType data
#' @usage ndvi
#' @format Numeric array with 4 dimensions. As indicated by the dimnames of the array:
#' \itemize{
#' \item{dim 1:} {longitude,}
#' \item{dim 2:} {latitude,}
#' \item{dim 3:} {day of the year,}
#' \item{dim 4:} {year.}
#' }
#' The values are NDVI values, and hence, between 0 and 1. Missing values are encoded as \code{NA}. 
#' @source The actual MOD13A data product is available from NASA EOSDIS Land Processes DAAC,
#' USGS Earth Resources Observation and Science (EROS) Center, Sioux Falls, South Dakota \url{https://lpdaac.usgs.gov}.
#' MODIS data can be downloaded with the R package MODIS \url{https://r-forge.r-project.org/projects/modis/}.
#' @examples
#' str(ndvi)
#' Image(ndvi)
NULL



