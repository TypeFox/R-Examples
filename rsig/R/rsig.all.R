#' Robust Signature Selection for Survival Outcomes with Estimation of Selection Probabilities of Features
#' 
#' Fit a specified model using subsamples and evaluate its performance on out-of-subsample data.
#'
#' @param surv [\code{Surv}]\cr
#'   Survival object, see \code{\link[survival]{Surv}}.
#' @param X [\code{data.frame}]\cr
#'   Data frame or matrix or matrix of input data (rows: examples, columns: features).
#' @param model [\code{character(1)}]\cr
#'   Model to use. One of \cr
#'   "rs.prlasso" (preconditioned lasso with robust selection), \cr
#'   "rs.lasso" (penalized Cox regression with robust selection), \cr
#'   "prlasso" (preconditioned lasso), or \cr
#'   "lasso" (penalized Cox regression)
#' @param n.rep.out [\code{integer}]\cr
#'   The number of replicates to be used to estimate selection probability of features (outer subsampling)
#' @param n.rep.in [\code{integer}]\cr
#'   The number of replicates to be used for model aggregation (inner subsampling)
#' @param plapply [\code{function}]\cr
#'   Function used for internal parallelization.
#'   Default is \code{\link[parallel]{mclapply}} for multi-core parallel execution.
#' @param sd.filter [\code{list}]\cr
#'   Pre-filter features by their standard deviation, by one of the options specified:\cr
#' topk: no. of features to be selected with largest standard devations.\cr
#' quant: the min percentile in standard deviations of features to be selected. 
#' @return Object of class \dQuote{list}.
#'   \item{selection.frequency}{a named vector of selected features with their estimated selection frequencies amongst n.rep.out replicates.}
#'   \item{perf}{performance measured on out-of-sample data in n.rep.out replicates}
#' @seealso \code{\link{rsig}}
#' @export
rsig.all = function(surv, X, model, n.rep.out=10L, n.rep.in=10L, plapply = mclapply, sd.filter=NULL) {
  checkArg(surv, "Surv")
  checkArg(X, c("data.frame", "matrix"))
  checkArg(model, choices = c("rs.prlasso", "rs.lasso", "prlasso", "lasso"))
  n.rep.out = convertInteger(n.rep.out); checkArg(n.rep.out, "integer", len=1L, lower=1L, na.ok=FALSE)
  n.rep.in = convertInteger(n.rep.in); checkArg(n.rep.in, "integer", len=1L, lower=1L, na.ok=FALSE)
  checkArg(plapply, "function", formals = c("X", "FUN"))
#   checkArg(sd.filter, "logical", len=1L, na.ok=FALSE)
#   checkArg(sd.filter, "list")
  
  # such info messages should be printed with message()
  messagef("Robust selection with %s (with %d outer, %d inner replicates)\n", model, n.rep.out, n.rep.in)
  
  # Do not overwrite user options!
  # options("mc.cores"=ceiling(detectCores()*.65))
  
  subsets = balanced.subset.sample(n.rep.out, surv[,2L])
  
  list.f = list.perf = vector("list", n.rep.out)
  
  out <- plapply(subsets, function(s) {    
    train.idx = s$train.idx
    test.idx = s$test.idx
    
    cat('.', append=TRUE); flush.console()
    fit = rsig(surv[train.idx,,drop=FALSE], X[train.idx,,drop=FALSE], model=model, 
               n.rep=n.rep.in, plapply = lapply, sd.filter=sd.filter, verbose=FALSE)            # serial excution here
    
    pred = predict(fit, X[test.idx,,drop=FALSE])
    
    cat('o', append=TRUE); flush.console()
    
    perf = rsig.eval(pred, surv[test.idx,], X[test.idx,])
    return(list(f=names(fit$beta), perf=perf))
  })
  cat('\n')
#   print(out)
  
  list.f = extractSubList(out, "f")
  list.perf = extractSubList(out, "perf")
  list.perf = matrix(unlist(list.perf), nrow=ncol(list.perf), ncol=nrow(list.perf), byrow=TRUE, dimnames=list(NULL, rownames(list.perf)))
#   list.perf = t(sapply(out$perf, unlist))
#   all.f = unique(unlist(list.f))
#   freq.f = rep(0, length(all.f)); names(freq.f) = all.f
#   for(i in seq_along(list.f)) {
#     freq.f[ list.f[[i]] ] <- freq.f[ list.f[[i]] ] + 1    
#   }
  freq.f = table(unlist(list.f))
  freq.f = sort(freq.f, decreasing=TRUE)

  return(list(selection.freq = freq.f / length(list.f), perf=list.perf))
}
