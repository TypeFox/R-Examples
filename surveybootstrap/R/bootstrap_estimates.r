
#####################################################
##' bootstrap.estimates
##'
##' this function contains the core of the rescaled bootstrap
##' method for estimating uncertainty in our estimates
##' it should be designed so that it can be passed in to
##' estimation functions as an argument\cr
##' OR\cr
##' \cr
##' @section TODO:
##' \itemize{
##'   \item{estimator.fn/bootstrap.fn and summary.fn are treated differently
##'         (one expects characters, one expects an actual fn. fix!)}
##'   \item{write description block, including estimator.fn, bootstrap.fn,
##'         summary.fn, more?}
##'}
##'
##' @param survey.data the dataset to use
##' @param survey.design a formula describing the design of the survey
##'                      (see below - TODO)
##' @param estimator.fn name of a function which, given a dataset like
##'                     survey data and arguments in \code{...},
##'                     will produce an estimate of interest
##' @param bootstrap.fn name of the method to be used to take
##'                     bootstrap resamples; see below
##' @param num.reps the number of bootstrap replication samples to draw
##' @param weights weights to use in estimation (or NULL, if none)
##' @param summary.fn (optional) name of a function which, given the set of estimates
##'                   produced by estimator.fn, summarizes them. if not specified, all of
##'                   the estimates are returned in a list
##' @param parallel if TRUE, use the plyr library's .parallel argument to
##'                 produce bootstrap resamples and estimates in parallel
##' @param paropts if not NULL, additional arguments to pass along to the
##'                parallelization routine
##' @param verbose if TRUE, produce lots of feedback about what is going on
##' @param ... additional arguments which will be passed on to the estimator fn
##' @return if no summary.fn is specified, then return the list of estimates
##'         produced by estimator.fn; if summary.fn is specified, then return
##'         its output
##' @export
##' @examples
##' \donttest{
##' # code goes here
##' 
##' }
bootstrap.estimates <- function(survey.data,
                                survey.design,
                                bootstrap.fn,
                                estimator.fn,
                                num.reps,
                                weights=NULL,
                                ...,
                                summary.fn=NULL,
                                verbose=TRUE,
                                parallel=FALSE,
                                paropts=NULL)
{

  ## get the weights
  weights <- get.weights(survey.data, weights)

  ## build up a single call to obtain an actual bootstrap
  ## replicate; we'll call this once for each one...
  boot.call <- match.call()

  boot.arg.idx <- match(c("survey.data",
                          "survey.design",
                          "num.reps",
                          "parallel",
                          "paropts"),
                        names(boot.call),
                        0L)
  boot.call <- boot.call[c(1,boot.arg.idx)]
  
  boot.call[[1]] <- get.fn(bootstrap.fn)

  ## also build up a call to obtain an estimate from the data
  est.call <- match.call(expand.dots=TRUE)

  ## these are the args we *won't* use when we call the estimator
  ## (ie, we use them here or in the bootstrap fn instead)
  est.arg.idx <- match(c("survey.design",
                         "estimator.fn",
                         "num.reps",
                         ##"weights",
                         "summary.fn",
                         "bootstrap.fn",
                         "parallel",
                         "paropts",
                         ## don't use survey.data b/c we're going to pass
                         ## in a bootstrap resample each time instead...
                         "survey.data"),
                       names(est.call),
                       0L)
  est.call <- est.call[-est.arg.idx]
  est.call[[1]] <- get.fn(estimator.fn)

  ## get the bootstrap samples
  boot.idx <- eval(boot.call, parent.frame())

  ## produce our estimate for each one
  res <- plyr::llply(boot.idx,

               function(this.rep) {

                 ## use the resampled indices to construct
                 ## a full resampled dataset
                 tmpdat <- survey.data[this.rep$index,]
                 tmpweights <- weights[this.rep$index]

                 ## apply the weight.scale to the estimation weights
                 tmpweights <- tmpweights * this.rep$weight.scale

                 ## add the information about which rows in the individual dataset
                 ## these resamples come from as an attribute
                 attr(tmpdat, "resampled.rows.orig.idx") <- this.rep$index

                 est.call[["survey.data"]] <- tmpdat
                 est.call[["weights"]] <- tmpweights

                 ## call estimator.fn to produce an estimate from
                 ## the bootstrap-resampled dataset
                 this.est <- eval(est.call, parent.frame(2))

                 return(this.est)
               },
               .parallel=parallel,
               .paropts=paropts)

  ## if the user specified a summary function, use it
  if (! is.null(summary.fn)) {

    this.sf <- get.fn(summary.fn)
    res <- do.call(this.sf, list(res))

  }

  return(res)

}

