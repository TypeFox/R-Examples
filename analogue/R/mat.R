###########################################################################
##                                                                       ##
## mat - function to perform the modern analogue technique for           ##
##       environmental reconstruction                                    ##
##                                                                       ##
## Created       : 17-Apr-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.0-1                                                 ##
## Last modified : 17-Apr-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
##                                                                       ##
###########################################################################
## x = training data, y = env var of interest
mat <- function(x, ...) UseMethod("mat")

mat.default <- function(x, y,
                        method = c("euclidean", "SQeuclidean", "chord",
                        "SQchord", "bray", "chi.square", "SQchi.square",
                        "information", "chi.distance", "manhattan",
                        "kendall", "gower", "alt.gower", "mixed"),
                        kmax, ...) {
    ##dims <- dim(x) # the numbers of samples / species
    site.nams <- rownames(x) # store sample names for later
    .call <- match.call()
    ## need to reset due to method dispatch
    .call[[1]] <- as.name("mat")
    if(missing(method))
        method <- "euclidean"
    method <- match.arg(method)
    dis <- distance(x, method = method, ...) # calculate the distances
    ## new speed-ups might leave dimnames on dis
    x <- as.matrix(x) # convert to matrix for speed (?)
    nams <- dimnames(x)
    dimnames(x) <- NULL # clear the dimnames for speed (?)
    ## insure sample under test is not chosen as analogue for itself
    diag(dis) <- NA
    ## process the kmax
    if(missing(kmax)) {
        kmax <- nrow(dis) - 1
    }
    if(kmax > (nrow(dis) - 1)) {
        kmax <- nrow(dis) - 1
    }
    if(kmax < 1) {
        kmax <- 1
    }
    ## drop = FALSE in next calls as we now make sure sample cannot be
    ## chosen as analogue for itself
    Wmeans <- apply(dis, 2, cumWmean, y, drop = FALSE, kmax = kmax) # Estimated values
    means <- apply(dis, 2, cummean, y, drop = FALSE, kmax = kmax)
    minDC <- apply(dis, 2, minDij, drop = FALSE) # minimum Dij per sample
    Werror <- sweep(Wmeans, 2, y, "-") # residuals for Wmeans
    error <- sweep(means, 2, y, "-") # residuals for mean
    WRMSE <- sqrt(rowMeans(Werror^2))
    k.w <- which.min(WRMSE)
    RMSE <- sqrt(rowMeans(error^2))
    k <- which.min(RMSE)
    Wbias <- rowMeans(Werror)
    bias <- rowMeans(error)
    Wmax.bias <- apply(Werror, 1, maxBias, y) # maximum bias
    max.bias <- apply(error, 1, maxBias, y)
    r2.mean <- apply(means, 1, function(x, y) {cor(x, y)^2}, y) # r.squared
    r2.Wmean <- apply(Wmeans, 1, function(x, y) {cor(x, y)^2}, y)
    ## re-apply samples names and n. closest
    colnames(Wmeans) <- colnames(means) <- site.nams
    colnames(Werror) <- colnames(error) <- site.nams
    rownames(Wmeans) <- rownames(means) <-
        rownames(Werror) <- rownames(error) <- seq_len(kmax)
    dimnames(x) <- nams
    ## return results
    retval <- structure(list(standard = list(est = means, resid = error,
                             rmsep = RMSE, avg.bias = bias, max.bias = max.bias,
                             r.squared = r2.mean, k = k, auto = TRUE),
                             weighted = list(est = Wmeans, resid = Werror,
                             rmsep = WRMSE, avg.bias = Wbias, max.bias = Wmax.bias,
                             r.squared = r2.Wmean, k = k.w, auto = TRUE),
                             Dij = dis,
                             orig.x = x,
                             orig.y = y,
                             call = .call,
                             method = method),
                        class = "mat")
    attr(retval, "method") <- method
    retval
}

mat.formula <- function(formula, data, subset, na.action,
                        method = c("euclidean", "SQeuclidean", "chord",
                          "SQchord", "bray", "chi.square", "SQchi.square",
                          "information", "chi.distance", "manhattan",
                          "kendall", "gower", "alt.gower", "mixed"),
                        model = FALSE, ...) {
  if(missing(method))
    method <- "euclidean"
  ## the function call
  .call <- match.call()
  ## need to reset due to method dispatch
  .call[[1]] <- as.name("mat")
  ## keep only the arguments which should go into the model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  ## drop the intercept
  attr(attr(mf, "terms"), "intercept") <- 0
  ## 1) allow model.frame to update the terms object before saving it.
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  res <- mat.default(x, y, method = method, ...)
  res$na.action <- attr(mf, "na.action")
  res$call <- .call
  if(model) {
    res$terms <- mt
    res$model <- mf
  }
  return(res)
}

print.mat <- function(x, k = 10,
                      digits = min(3, getOption("digits") - 4),
                      ...)
  {
    ##if(is.null(k))
    ##  k <- k(x)
    cat("\n")
    writeLines(strwrap("Modern Analogue Technique", prefix = "\t"))
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    tbl <- cbind(x$standard$rmsep[1:k], x$standard$r.squared[1:k],
                 x$standard$avg.bias[1:k], x$standard$max.bias[1:k])
    tbl.w <- cbind(x$weighted$rmsep[1:k], x$weighted$r.squared[1:k],
                   x$weighted$avg.bias[1:k], x$weighted$max.bias[1:k])
    tbl <- as.matrix(format(tbl, digits = digits))
    tbl.w <- as.matrix(format(tbl.w, digits = digits))
    tbl <- cbind(as.integer(1:k), tbl)
    tbl.w <- cbind(as.integer(1:k), tbl.w)
    rownames(tbl) <- rownames(tbl.w) <- rep("", nrow(tbl))
    colnames(tbl) <- colnames(tbl.w) <- c("k",
                                          "RMSEP","R2","Avg Bias","Max Bias")
    cat("\nPercentiles of the dissimilarities for the training set:\n\n")
    print(quantile(x$Dij[lower.tri(x$Dij)],
                   probs = c(0.01, 0.02, 0.05, 0.1, 0.2)),
          digits = digits)
    cat("\nInferences based on the mean of k-closest analogues:\n\n")
    print(tbl, quote = FALSE, right = TRUE)
    cat("\nInferences based on the weighted mean of k-closest analogues:\n\n")
    print(tbl.w, quote = FALSE, right = TRUE)
    cat("\n")
    invisible(x)
  }
