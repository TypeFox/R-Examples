#' Robust Signature Selection for Survival Outcomes
#'
#' Find a robust signature, i.e. a set of features, using averaged and shrucken generalized linear models.
#' Subsamples are taken to fit models, via \eqn{\ell_1}-penalized Cox regression (lasso) or preconditioned lasso (prlasso) algorithm.
#'
#' @param surv [\code{Surv}]\cr
#'   Survival object, see \code{\link[survival]{Surv}}.
#' @param X [\code{data.frame}]\cr
#'   Data frame or matrix or matrix of input data (rows: examples, columns: features). Columns must have names assigned.
#' @param model [\code{character(1)}]\cr
#'   Model to use. One of \cr
#'   "rs.prlasso" (preconditioned lasso with robust selection), \cr
#'   "rs.lasso" (penalized Cox regression with robust selection), \cr
#'   "prlasso" (preconditioned lasso), or \cr
#'   "lasso" (penalized Cox regression)
#' @param n.rep [\code{integer}]\cr
#'   The number in replicates to be used for model aggregation. A large enough number is suggested.
#' @param plapply [\code{function}]\cr
#'   Function used for internal parallelization.
#'   Default is \code{\link[parallel]{mclapply}} for multi-core parallel execution. Change it to \code{\link[base]{lapply}} for single-core execution.
#' @param sd.filter [\code{list}]\cr
#'   Pre-filter features by their standard deviation, by one of the options specified:\cr
#' topk: no. of features to be selected with largest standard devations, or\cr
#' quant: the min percentile in standard deviations of features to be selected. 
#' @param verbose [\code{logical}]\cr
#'   Controls message output.
#' @return Object of class \dQuote{rsig}; a list consisting of
#'   \item{model}{model specified by the user}
#'   \item{sd.filter}{sd.filter object}
#'   \item{beta}{coefficient vector}
#'   \item{intercept}{intercept}
#' @seealso \code{\link{predict.rsig}}, \code{\link{rsig.eval}}, \code{\link{rsig.all}}
#' @examples
#' # An example adapted from glmnet package
#' 
#' set.seed(11011)
#' n = 300
#' p = 10
#' nz = 3
#' X = matrix(rnorm(n*p),n,p,dimnames=list(NULL,seq_len(p)))
#' beta = rnorm(nz)
#' f = X[,seq_len(nz)] %*% beta
#' h = exp(f) / 365.25
#' t = rexp(n,h)
#' tcens = rbinom(n=n,prob=.3,size=1) # censoring indicator
#' S = Surv(t, 1-tcens)
#' 
#' fit = rsig(S, X, "rs.prlasso", n.rep=2)
#' pred = predict(fit, X)
#' perf = rsig.eval(pred, S, X)
#' @export
rsig = function(surv, X, model, n.rep=10L, plapply = mclapply, sd.filter=NULL, verbose=TRUE) {
  #TBD: clinical = data.frame(), 
  ### Argument checks
  checkArg(surv, "Surv")
  checkArg(X, c("matrix", "data.frame"))
#  checkArg(clinical, c("matrix", "data.frame"))
#   use.clinical = any(dim(clinical))
  if (nrow(surv) != nrow(X))
    stopf("%s and %s must have the same number of observations", sQuote("surv"), sQuote("X"))
#   if (use.clinical && nrow(clinical) != nrow(X))
#     stopf("%s and %s must have the same number of observations", sQuote("clinical"), sQuote("X"))
  model = match.arg(model, c("rs.prlasso", "rs.lasso", "prlasso", "lasso")) #TBD: "sd", "rank", "clinical"
  checkArg(plapply, "function")

  if(!is.null(sd.filter)) {
    # Filter features by sd
    checkArg(sd.filter, "list")
    tmp = filterSD(X, sd.filter)
    X = X[, tmp, drop=FALSE]
    
    if(verbose) cat(sprintf('No. of features after sd.filter: %d\n', dim(X)[2]))
  }
  
  ### Conversions
  X = data.matrix(X)
#   if (model %in% c("rs.prlasso", "rs.lasso", "prlasso", "lasso", "sd", "rank")) {
#     if (use.clinical)
#       clinical = data.matrix(clinical)
#   }
#   if (use.clinical) {
# #     is.clinical = c(rep.int(TRUE, ncol(clinical)), rep.int(FALSE, ncol(X)))
#     X = cbind(clinical, X)
#   }

  # FIXME we might need to set mc.cores=1 on windows 

  ### Fit the model
  if(substr(model, 1L, 3L) == "rs.") {
    if(verbose)
      catf("Robust selection: %s (with %d inner replicates)\n", model, n.rep)
    
    inner.model = substr(model, 4L, nchar(model))
    fit = run_rs(surv, X, model = inner.model, n.rep=n.rep,  plapply = plapply)
  } else {
    if(verbose)
      catf("Running: %s (no robust selection)\n", model)
    
    learner = get(sprintf("fit_%s", model), mode = "function", envir = environment(rsig))
    fit = learner(surv, X)
  }
  
  # This should be a list, right? Yes. Members of the list "fit" is to be copied to another list "slots".
  slots = c(list(model=model, sd.filter=sd.filter), fit)

  class(slots) = "rsig"
  return(slots)
}


run_rs = function(surv, X, cv.folds=10L, n.rep, model, n.thres=20L, plapply)
{
  subsets = balanced.subset.sample(n.rep, surv[,2])
  learner = get(sprintf("fit_%s", model), mode = "function", envir = environment(rsig))
  
  fits <- plapply(subsets, function(s) {  
    train = s$train.idx
    return(learner(surv[train,], X[train,,drop=FALSE], cv.folds=cv.folds, sparse.output=TRUE))
  })
  
  ### Aggregation
  # SL: cbind2 doesn't work with more than 2 elements in the list
  fit.mean = rowMeans(do.call(cBind, extractSubList(fits, "beta")))
  names(fit.mean) = rownames(fits[[1L]]$beta)
  intercept.mean = mean(extractSubList(fits, "intercept"))
  
  ### Shrinking
  # ML: something like this? 
  # SL: Not exactly the same. My code try to get n.thres values with equal gaps between min and max vals in abs(fit.means)
  # SL: Let's stop changing this part. The structure of output here affects many other things...
#   ind = abs(fit.mean) > sqrt(.Machine$double.eps) & rank(abs(fit.mean)) >= (length(fit.mean) - n.thres)
#   thres.list = sort(abs(fit.mean[ind]), decreasing=TRUE)
  ind = abs(fit.mean) > sqrt(.Machine$double.eps)
  sorted = sort(abs(fit.mean[ind]), decreasing=TRUE)
  if(length(sorted) > n.thres) {
    thres.list = sorted[seq(1L, length(sorted), by=(length(sorted)-1)/n.thres)]
  } else {
    thres.list = sorted
  }
  thres.nnz = sapply( thres.list, function(t) { length(which(abs(fit.mean)>=t)) })
  
  idx = which(thres.nnz >= 5L) # min. # of features for superpc.cv is 5.
  thres.list = thres.list[idx]
  thres.nnz = thres.nnz[idx]
  names(thres.list) = thres.nnz
#   print(thres.list)
  
  # CV to find the optimal threshold
  list.ll = numeric(length(thres.list))
  for(i in 1:length(thres.list)) {
    thres = thres.list[i]
    feat.idx = which(abs(fit.mean) >= thres)
    beta = fit.mean[feat.idx]

    ll = numeric(n.rep)
    for(j in seq_len(n.rep)) {
      idx.tt = subsets[[j]]$test.idx
      
      surv.tt = surv[idx.tt]
      X.tt = X[idx.tt, feat.idx, drop=FALSE]
      
      pred = X.tt %*% beta + intercept.mean
      ll[j] = concordance.index(x=pred, surv.time = surv.tt[,1], surv.event = surv.tt[,2], method="noether", na.rm=TRUE)$c.index
    }
    if(length(which(is.na(ll)))>0) {
      cat('Warning: ', length(which(is.na(ll))), ' C values were NA in inner replicate ', i, '\n')
    }
    list.ll[i] = mean(ll, na.rm=TRUE)
  }
  
  best.thres = thres.list[which.max(list.ll)]
  best.feats = which(abs(fit.mean) >= best.thres)
  
  return(list(beta=fit.mean[best.feats], intercept=intercept.mean))
}

