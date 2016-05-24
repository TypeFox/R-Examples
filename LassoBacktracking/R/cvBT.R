#' Cross-validation for \code{LassoBT}
#' 
#' Perform k-fold cross-validation potentially multiple times on permuted version of the data.
#' 
#' @inheritParams LassoBT
#' @param nfolds number of folds. Default is 5.
#' @param nperms the number of permuted datasets to apply k-folds corss-validation to. Default is 1 so
#'  we carry out vanilla cross-validation.
#' @param mc.cores the number of cores to use. Only applicable when not in Windows as it uses
#'  the \pkg{parallel} package to parallelise the computations.
#' @param ... other arguments that can be passed to \code{LassoBT}. 
#' @return A list with components as below.
#' \describe{
#'   \item{\code{lambda}}{the sequence of \code{lambda} values used}
#'   \item{\code{cvm}}{a matrix of error estimates (with squared error loss). The rows correspond
#'    to different \code{lambda} values whilst the columns correspond to different iterations}
#'   \item{\code{BT_fit}}{a "\code{BT}" object from a fit to the full data.}
#'   \item{\code{cv_opt}}{a two component vector giving the cross-validation optimal \code{lambda} index
#'    and iteration}
#'   \item{\code{cv_opt_err}}{the minimal cross-validation error.}
#' }
#' @examples
#' x <- matrix(rnorm(100*250), 100, 250)
#' y <- x[, 1] + x[, 2] - x[, 1]*x[, 2] + x[, 3] + rnorm(100)
#' out <- cvLassoBT(x, y, iter_max=10, nperms=2)
#' @export
#' @importFrom parallel mclapply
cvLassoBT <- function(x, y, lambda=NULL, nlambda = 100L, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 0.0001),
                  nfolds=5L, nperms=1L, mc.cores=1L, ...) {
  nfolds <- as.integer(nfolds)
  nperms <- as.integer(nperms)
  if (nfolds < 2L) stop("nfolds must be at least 2")
  tot_reps <- nfolds * nperms
  
  if (!is.matrix(x)) stop("x should be a matrix with two or more columns")
  np <- dim(x)
  if (is.null(np) | (np[2] < 1L)) 
    stop("x should be a matrix at least one non-constant column")
  nobs <- n <- nrow(x)
  nvars <- ncol(x)
  # needed for lambda.min.ratio
  
  # fit on all data first
  out_full_fit <- LassoBT(x, y, lambda=lambda, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, ...)
  
  lambda <- out_full_fit$lambda
    
  foldid1 <- rep(0L:(nfolds-1L), length=n) # folds go 0, 1,..., nfolds-1
  foldid <- sapply(1L:nperms, function(i) foldid1[sample.int(n)])
  rm(foldid1)
  
  rep_fun <- function(i) {
    # apply on folds
    cur_fold_ind <- i %% nfolds
    cur_perm <- i %/% nfolds + 1L
    cur_fold <- foldid[, cur_perm] == cur_fold_ind 
    
    out <- LassoBT(x[!cur_fold, , drop=FALSE], y[!cur_fold], lambda=lambda, ...)
    out <- apply((y[cur_fold] - predict.BT(out, newx=x[cur_fold, , drop=FALSE], s=NULL, iter=NULL, type="response"))^2,
                 c(2, 3), sum)
    return(out) # should be a nlambda by max_iter matrix
  }
  
  out <- parallel::mclapply(0L:(tot_reps-1L), rep_fun, mc.cores=mc.cores)
  Max_iter <- max(sapply(out, ncol))
  out_av <- matrix(0, nrow=nlambda, ncol=Max_iter)
  for (j in 1L:tot_reps) {
    cur_max_iter <- ncol(out[[j]])
    out_av[, 1L:cur_max_iter] <- out_av[, 1L:cur_max_iter] + out[[j]]
    # pad out remaining with last col
    for (k in seq_len(Max_iter - cur_max_iter)) {
      out_av[, k + cur_max_iter] <- out_av[, k + cur_max_iter] + out[[j]][, cur_max_iter]
    }
  }
  out_av <- out_av / (n*nperms)
  opt_err <- min(out_av)
  ## out_av must be a matrix here
  opt <- as.integer(which(out_av == opt_err, arr.ind=TRUE)[1L, ]) # take the smallest number of iterations
  #opt[1L] <- lambda[opt[1L]]
  names(opt) <- c("lambda_index", "iteration")
  out <- list("lambda"=lambda,
              "cvm"=out_av,
              "BT_fit"=out_full_fit,
              "cv_opt"=opt,
              "cv_opt_err"=opt_err)
  return(out)
}
