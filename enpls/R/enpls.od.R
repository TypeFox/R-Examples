#' Ensemble Partial Least Squares for Outlier Detection
#'
#' This function performs outlier detection with ensemble partial least squares.
#'
#' This function performs outlier detection with ensemble partial least squares.
#'
#' @param x predictor matrix
#' @param y response vector
#' @param maxcomp Maximum number of components included within the models, 
#' if not specified, default is the variable (column) numbers in x.
#' @param MCtimes times of Monte-Carlo
#' @param method \code{"mc"} or \code{"bootstrap"}. Default is \code{"mc"}.
#' @param ratio sample ratio used when \code{method = "mc"}
#' @param parallel Integer. Number of parallel processes to use. 
#' Default is \code{1}, which means run serially.
#'
#' @return A list containing four components:
#' \itemize{
#' \item \code{error.mean} - error mean for all samples (absolute value)
#' \item \code{error.median} - error median for all samples
#' \item \code{error.sd} - error sd for all samples
#' \item \code{predict.error.matrix} - the original prediction error matrix
#' }
#'
#' @author Nan Xiao <\email{road2stat@@gmail.com}>
#'
#' @seealso See \code{\link{enpls.fs}} for feature selection with ensemble PLS. 
#' See \code{\link{enpls.en}} for ensemble PLS regression.
#'
#' @export enpls.od
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#'
#' @references
#' DongSheng Cao, Yizeng Liang, Qingsong Xu, Hongdong Li, and Xian Chen. 
#' "A new strategy of outlier detection for QSAR/QSPR."
#' \emph{Journal of computational chemistry} 31, no. 3 (2010): 592--602.
#'
#' Dongsheng Cao, Yizeng Liang, Qingsong Xu, Yifeng Yun, and Hongdong Li. 
#' "Toward better QSAR/QSPR modeling: simultaneous outlier detection and 
#' variable selection using distribution of model features." 
#' \emph{Journal of computer-aided molecular design} 25, no. 1 (2011): 67--80.
#'
#' @examples
#' data(alkanes)
#' x = alkanes$x
#' y = alkanes$y
#'
#' set.seed(42)
#' od = enpls.od(x, y, MCtimes = 100)
#' print(od)
#' plot(od)
#' plot(od, criterion = 'sd')

enpls.od = function(x, y, 
                    maxcomp = NULL, 
                    MCtimes = 500L, 
                    method = c('mc', 'bootstrap'), ratio = 0.8, 
                    parallel = 1L) {

  if (is.null(maxcomp)) maxcomp = ncol(x)

  method = match.arg(method)

  x.row = nrow(x)
  samp.idx = vector('list', MCtimes)
  samp.idx.remain = vector('list', MCtimes)

  if (method == 'mc') {
    for (i in 1L:MCtimes) {
      samp.idx[[i]] = sample(1L:x.row, floor(x.row * ratio))
      samp.idx.remain[[i]] = setdiff(1L:x.row, samp.idx[[i]])
    }
  }

  if (method == 'bootstrap') {
    for (i in 1L:MCtimes) {
      samp.idx[[i]] = sample(1L:x.row, x.row, replace = TRUE)
      samp.idx.remain[[i]] = setdiff(1L:x.row, unique(samp.idx[[i]]))
    }
  }

  plsdf = as.data.frame(cbind(x, y))

  if (parallel < 1.5) {

    errorlist = vector('list', MCtimes)
    for (i in 1L:MCtimes) {
      plsdf.sample = plsdf[samp.idx[[i]], ]
      plsdf.remain = plsdf[samp.idx.remain[[i]], ]
      errorlist[[i]] = suppressWarnings(enpls.od.core(plsdf.sample, plsdf.remain, maxcomp))
    }

  } else {

    registerDoParallel(parallel)
    errorlist = foreach(i = 1L:MCtimes) %dopar% {
      plsdf.sample = plsdf[samp.idx[[i]], ]
      plsdf.remain = plsdf[samp.idx.remain[[i]], ]
      enpls.od.core(plsdf.sample, plsdf.remain, maxcomp)
    }

  }

  prederrmat = matrix(NA, ncol = x.row, nrow = MCtimes)
  for (i in 1L:MCtimes) {
    for (j in 1L:length(samp.idx.remain[[i]])) {
      prederrmat[i, samp.idx.remain[[i]][j]] = errorlist[[i]][j]
    }
  }

  errmean   = abs(colMeans(prederrmat, na.rm = TRUE))
  errmedian = apply(prederrmat, 2L, median, na.rm = TRUE)
  errsd     = apply(prederrmat, 2L, sd, na.rm = TRUE)

  object = list('error.mean'    = errmean, 
                'error.median'  = errmedian, 
                'error.sd'      = errsd, 
                'predict.error.matrix' = prederrmat)
  class(object) = 'enpls.od'
  return(object)

}

#' core function for enpls.od
#'
#' select the best ncomp with cross-validation and
#' use it to fit the complete training set again, 
#' then predict on the test set. scale = TRUE
#'
#' @return the error vector between predicted y and real y
#' 
#' @keywords internal

enpls.od.core = function(plsdf.sample, plsdf.remain, maxcomp) {

  plsr.cvfit = plsr(y ~ ., data = plsdf.sample, 
                    ncomp  = maxcomp, 
                    scale  = TRUE, 
                    method = 'simpls', 
                    validation = 'CV', segments = 5L)

  # choose best component number using adjusted CV
  cv.bestcomp = which.min(RMSEP(plsr.cvfit)[['val']][2L, 1L, -1L])

  plsr.fit = plsr(y ~ ., data = plsdf.sample, 
                  ncomp  = cv.bestcomp, 
                  scale  = TRUE, 
                  method = 'simpls', 
                  validation = 'none')

  pred = predict(plsr.fit, ncomp = cv.bestcomp, 
                 newdata = plsdf.remain[, !(colnames(plsdf.remain) %in% c('y'))])[, 1L, 1L]

  error = plsdf.remain[, 'y'] - pred
  names(error) = NULL

  return(error)

}
