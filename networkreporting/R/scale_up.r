#####################################################
## scale_up.r
##
## the network scale-up method for estimating the
## size of hidden populations
##
## TODO -- handle NAs in code below (this is partially done)
## TODO -- make handling of column vars, etc, uniformly handled
##         when passed in as args. (ie, use get.var; now some fns
##         still don't use get.var)
## TODO -- perhaps an easy to use interface for nsum estimates
##         with bootstrap (and then use this in nsum.internal.validation)
## TODO -- when using defaults (for example, taking
##         popn size info from dataframe attributes,
##         should we print out a message to the screen?
##         or perhaps have a (default) verbose mode?)
## TODO -- think about code to get 45q15 from these data...

#####################################################
##' total.degree
##'
##' estimate the total degree of the population network
##' from sample degrees
##'
##' this computes the weighted sum of the respondents'
##' estimated degrees.\cr
##''
##' TODO -- for now, it doesn't worry about missing values
##' OR about differences between the frame and the universe
##'
##' @param survey.data the dataframe with survey results
##' @param d.hat.vals the name or index of the column that contains
##'                  each respondent's estimated degree
##' @param weights if not NULL, weights to use in computing the estimate. this
##'                should be the name of the column in the survey.data which has
##'                the variable with the appropriate weights. these weights
##'                should be construted so that, eg, the mean of the degrees is
##'                estimated as (1/n) * \\sum_i {w_i * d_i}
##' @param missing if "ignore", then proceed with the analysis without
##'                doing anything about missing values. if "complete.obs"
##'                then only use rows that have no missingness for the
##'                computations (listwise deletion). care
##'                must be taken in using this second option
##' @return the estimated total degree
##' @export
total.degree.estimator <- function(survey.data,
                                   d.hat.vals="d",
                                   weights=NULL,
                                   missing="ignore")
{

  ## get the weights;
  ## weights will default to 1 for everyone, unless the user specified
  ## a weights variable
  weights <- surveybootstrap:::get.weights(survey.data, weights)

  ## get the estimated degrees
  d.hat.vals <- surveybootstrap:::get.var(survey.data, d.hat.vals)

  if (missing == 'complete.obs') {
    touse <- which(! is.na(d.hat.vals))
    res <- sum(d.hat.vals[touse]*weights[touse])
  } else if (missing == 'ignore') {
    res <- sum(d.hat.vals*weights)
  } else {
    stop("error in specifying procedure for handling missing values in total.degree.estimator. invalid option.\n")
  }

  return(res)

}


#####################################################
##' nsum.estimator
##'
##' compute network scale-up (nsum) estimate of the
##' hidden population's size. if the degree ratio
##' and information transmission rate are both 1
##' (the defaults), this is the Killworth estimator.
##'
##' TODO -- cite Killworth estimator, our methods paper
##' TODO -- add refs to deg ratio and tx rate stuff...
##'
##' @param survey.data the dataframe with survey results
##' @param d.hat.vals the name or index of the column that contains
##'                  each respondent's estimated degree
##' @param y.vals the name or index of the column that contains
##'              the count of hidden popn members known
##' @param total.popn.size NULL, NA, or a size
##' @param weights if not NULL, weights to use in computing the estimate. this
##'                should be the name of the column in the survey.data which has
##'                the variable with the appropriate weights. these weights
##'                should be construted so that, eg, the mean of the degrees is
##'                estimated as (1/n) * \\sum_i {w_i * d_i}
##' @param deg.ratio the degree ratio, \\frac{\\bar{d_T}}{\\bar{d}}; defaults to 1
##' @param tx.rate the information transmission rate; defaults to 1
##' @param killworth.se if not NA, return the Killworth et al estimate of
##                 the standard error (not generally recommended)
##' @param missing if "ignore", then proceed with the analysis without
##'                doing anything about missing values. if "complete.obs"
##'                then only use rows that have no missingness for the
##'                computations (listwise deletion). care
##'                must be taken in using this second option
##' @param verbose if TRUE, print messages to the screen
##' @param ... extra parameters to pass on to the bootstrap fn, if applicable
##' @return the nsum estimate of the hidden population's size (as a prevalence or
##'         an absolute number, depending on total.popn.size)
##' @export
nsum.estimator <- function(survey.data,
                           d.hat.vals="d",
                           y.vals="y",
                           total.popn.size=NULL,
                           deg.ratio=1,
                           tx.rate=1,
                           weights=NULL,
                           killworth.se=FALSE,
                           missing="ignore",
                           verbose=FALSE,
                           ...)
{

  if (! missing %in% c("ignore", "complete.obs")) {
    stop("error in specifying procedure for handling missing values in nsum.estimator. invalid option.\n")
  }

  total.popn.size <- parse.total.popn.size(total.popn.size,
                                           survey.data,
                                           verbose=verbose)

  ## get the weights;
  ## weights will default to 1 for everyone, unless the user specified
  ## a weights variable
  weights <- surveybootstrap:::get.weights(survey.data, weights)

  raw.d.hat.vals <- surveybootstrap:::get.var(survey.data, d.hat.vals)

  raw.y.vals <- surveybootstrap:::get.var(survey.data, y.vals)

  #### compute the actual estimates
  y.vals <- raw.y.vals * weights
  d.hat.vals <- raw.d.hat.vals * weights

  ## figure out if we have to only use non-missing entries
  touse.idx <- 1:length(y.vals)
  if (missing == "complete.obs") {
    touse.idx <- which( (! is.na(y.vals) & ! is.na(d.hat.vals)))

    notused <- length(y.vals) - length(touse.idx)

    if (verbose & notused > 0) {
      cat(str_c("missing=='complete.obs', so dropping ",
                notused, " rows with missing data\n"))
    }

  }

  ## NB: for now, this will return NA if either the degrees or
  ##     the y_i's has any NAs
  res <- sum(y.vals[touse.idx])/sum(d.hat.vals[touse.idx])

  if (deg.ratio != 1) {
    res <- res * (1/deg.ratio)
  }

  if (tx.rate != 1) {
    res <- res * (1/tx.rate)
  }

  if (! is.na(total.popn.size)) {
    res <- res * total.popn.size
  }

  toret <- list(estimate=res,
                tot.connections=sum(y.vals[touse.idx]) *
                (1/deg.ratio) *
                (1/tx.rate),
                sum.d.hat=sum(d.hat.vals[touse.idx]))

  ## not recommended, but interesting in some cases:
  ## the killworth estimate for the standard error
  if (killworth.se) {

      ## NOTE: this is not really defined for the case of
      ## non-trivial degree ratio or transmission rate
      ## we'll use unadjusted proportion in all cases
      #p.hat <- toret$tot.connections / toret$sum.d.hat
      p.hat <- sum(y.vals[touse.idx])/sum(d.hat.vals[touse.idx])
      p.hat.raw <- sum(raw.y.vals[touse.idx])/sum(raw.d.hat.vals[touse.idx])

      ## see Killworth et al, 1998 (Evaluation Review)
      kse <- sqrt((p.hat * (1-p.hat))/toret$sum.d.hat)
      kse.wgtdenom <- sqrt((p.hat.raw * (1-p.hat.raw))/sum(raw.d.hat.vals[touse.idx]))

      if (! is.na(total.popn.size)) {
          kse <- kse * total.popn.size
          kse.wgtdenom <- kse.wgtdenom * total.popn.size
      }

      toret$killworth.se <- kse
      toret$killworth.se.wgtdenom <- kse.wgtdenom
  }

  return(toret)

}

#####################################################
##' nsum.internal.validation
##'
##' use a hold-one-out method to estimate the predictive
##' accuracy of the network scale-up estimator on the
##' known populations
##'
##' given a set of estimated degrees, responses to a group
##' of ARD questions, and the total size of the populations
##' that the ARD questions ask about, this function estimates
##' the accuracy of the network scale-up method by dropping
##' each known population in turn, using the non-dropped
##' populations to compute the degree and an estimate of the
##' size of the known population, and comparing the result
##' to the actual size of the known population
##'
##' TODO -- document bootstrap ci option better
##' TODO -- add example of usage to the comments...\cr
##' TODO -- make amenable to parallelization
##'
##' @param survey.data the dataframe with the survey results
##' @param known.popns if not NULL, a vector whose entries are the size of the known
##'                    populations, and whose names are the variable names in the dataset
##'                    corresponding to each one. if NULL, then assume that the survey.data
##'                    dataframe has an attribute called 'known.popns' containing this vector.
##' @param total.popn.size the size of the entire population. if NA,
##'                        this function works with proportions;
##'                        if NULL, it looks for the 'total.popn.size' attribute of the
##'                        dataset \code{survey.data};
##'                        if not NULL or NA, it
##'                        works with absolute numbers (ie, the proportions * total popn size)
##' @param degrees if not NULL, then the name or index of the column in the datset
##'                containing the degree estimates. if NULL, then use the known population
##'                method to estimate the degrees (see \code{\link{kp.degree.estimator}})
##' @param missing if "ignore", then proceed with the analysis without
##'                doing anything about missing values. if "complete.obs"
##'                then only use rows that have no missingness for the
##'                computations (listwise deletion). care
##'                must be taken in using this second option
##' @param kp.method if TRUE, then we're using known population method estimates of the
##'                  degrees. this means we have to recompute the degrees each time we
##'                  hold out a known subgroup. if the degrees come from another estimator,
##'                  like the summation method, then we don't need to do that since we
##'                  don't use the ARD questions in coming up with the degree estimate.
##' @param weights if not NULL, weights to use in computing the estimate. this
##'                should be the name of the column in the survey.data which has
##'                the variable with the appropriate weights. these weights
##'                should be construted so that, eg, the mean of the degrees is
##'                estimated as (1/n) * \\sum_i {w_i * d_i}
##' @param killworth.se if TRUE, return the Killworth et al estimate of the standard error
##' @param return.plot if TRUE, make and return a ggplot2 plot object
##' @param verbose if TRUE, report more detailed information about what's going on
##' @param bootstrap if TRUE, use \code{bootstrap.estimates} to take bootstrap resamples
##'                  in order to obtain intervals around each estimate. in this case,
##'                  you are expected to also pass in at least \code{bootstrap.fn},
##'                  \code{survey.design}, and \code{num.reps}
##' @param ... additional arguments, which are passed on to \code{bootstrap.estimates}
##'            if \code{bootstrap} is TRUE
##' @return a list with a dataset containing the subpopn-specific estimates, as well as
##'         several summaries of the accuracy of those estimates, including
##'         mae (mean absolute error), mse (mean squared error),
##'         rmse (root mean squared error), and are (average relative error)
##' @export
nsum.internal.validation <- function(survey.data,
                                     known.popns=NULL,
                                     total.popn.size=NULL,
                                     degrees=NULL,
                                     missing="ignore",
                                     kp.method=FALSE,
                                     weights=NULL,
                                     killworth.se=FALSE,
                                     return.plot=FALSE,
                                     verbose=FALSE,
                                     bootstrap=FALSE,
                                     ...)
{

  if (! missing %in% c("ignore", "complete.obs")) {
    stop("error in specifying procedure for handling missing values in nsum.internal.validation. invalid option.\n")
  }

  if (is.null(known.popns)) {
    known.popns <- attr(survey.data, "known.popns")
  }

  if (bootstrap) {

    boot.call <- match.call()

    ## for use below, in building up the call to the estimator fn
    boot.expected.args <- c("bootstrap.fn",
                            "survey.design",
                            "num.reps")
    boot.other.args <- c("parallel", "paropts")

    ## be sure that the expected args have been passed in
    chk.idx <- match(boot.expected.args, names(boot.call))
    if(any(is.na(chk.idx))) {
      stop("bootstrap was specified, but you are missing one of the required arguments")
    }
  }

  total.popn.size <- parse.total.popn.size(total.popn.size,
                                           survey.data,
                                           verbose=verbose)

  ## go through each known population...
  res.all <- plyr::llply(names(known.popns),

               function(this.kp) {

                 surveybootstrap:::vcat(verbose, "staring known popn: ", this.kp)

                 known.size <- known.popns[this.kp]

                 ## if we're using the known population method, then
                 ## we need to recompute the degrees without using the
                 ## holdout variable
                 if( kp.method ) {

                   kp.minus <- known.popns[-match(this.kp,names(known.popns))]
                   deg.minus <- kp.degree.estimator(survey.data=survey.data,
                                                    known.popns=kp.minus,
                                                    total.popn.size=total.popn.size,
                                                    missing=missing,
                                                    verbose=verbose)

                   thisdat <- survey.data
                   thisdat$deg.minus <- deg.minus

                 } else {

                   ##TODO -- handle NAs better in the near future...
                   thisdat <- survey.data
                   thisdat$deg.minus <- thisdat[,degrees]

                 }

                 degsum <- total.degree.estimator(thisdat,
                                                 d.hat.vals="deg.minus",
                                                 missing=missing)

                 ## note that this attribute could be NULL and we're
                 ## still OK here...
                 tps <- attr(survey.data, 'total.popn.size')
                 attr(thisdat, 'total.popn.size') <- tps

                 ## stick the y values and weights into the actual
                 ## dataframe; this is useful because if we take
                 ## bootstrap samples below, we only have to worry about
                 ## the dataset and not other, parallel, vectors
                 thisdat$y.val.minus <- surveybootstrap:::get.var(survey.data, this.kp)
                 thisdat$weights.minus <- surveybootstrap:::get.weights(survey.data, weights)

                 ## build up a call to nsum.estimator
                 ## (we do this rather than just calling directly
                 ##  because it's useful to have the args from this
                 ##  call when we do the bootstrap, below)
                 est.call <- call("nsum.estimator",
                                  survey.data=thisdat,
                                  d.hat.vals="deg.minus",
                                  y.vals="y.val.minus",
                                  total.popn.size=total.popn.size,
                                  weights="weights.minus",
                                  missing=missing,
                                  verbose=verbose,
                                  killworth.se=killworth.se)

                 nsum.holdout.res <- eval(est.call)
                 nsum.holdout.est <- nsum.holdout.res$estimate
                 nsum.holdout.sum.d.hat <- nsum.holdout.res$sum.d.hat

                 ##nsum.holdout.kse <- NULL
                 nsum.holdout.kse <- NA
                 nsum.holdout.kse.wgtdenom <- NA
                 if (killworth.se) {
                     nsum.holdout.kse <- nsum.holdout.res$killworth.se
                     nsum.holdout.kse.wgtdenom <- nsum.holdout.res$killworth.se.wgtdenom
                 }

                 boot.res <- NULL

                 if (bootstrap) {

                   ## if we're using the bootstrap, build up
                   ## a call to bootstrap.estimates
                   ## (note that at the top of this fn, we check
                   ##  to be sure that the minimal expected args
                   ##  for bootstrap.estimates are present if
                   ##  bootstrap==TRUE)
                   boot.expected.args <- c("bootstrap.fn",
                                           "survey.design",
                                           "num.reps")
                   boot.other.args <- c("parallel", "paropts")
                   boot.arg.idx <- match(c(boot.expected.args,
                                           boot.other.args),
                                         names(boot.call),
                                         0L)

                   boot.args <- as.list(boot.call)[boot.arg.idx]
                   boot.args[["estimator.fn"]] <- "nsum.estimator"
                   boot.args <- c(boot.args, as.list(est.call[-1]))
                   boot.call <- as.call(c(as.name("bootstrap.estimates"),
                                          boot.args))
                   boot.res <- eval(boot.call)

                   boot.res.ests <- plyr::laply(boot.res, function(x) { x$estimate })
                   boot.res.d.hat.sum <- plyr::laply(boot.res,
                                               function(x) { x$sum.d.hat })
                 } else {
                   boot.res.ests <- NULL
                   boot.res.d.hat.sum <- NULL
                 }

                 return(list(data=data.frame(name=this.kp,
                                             nsum.holdout.est=nsum.holdout.est,
                                             known.size=known.size,
                                             d.hat.sum=degsum,
                                             killworth.se=nsum.holdout.kse,
                                             killworth.se.wgtdenom=nsum.holdout.kse.wgtdenom),
                             boot=list(name=this.kp,
                                       values=boot.res.ests,
                                       d.hat.sum=boot.res.d.hat.sum)))

               })

  res <- plyr::ldply(res.all,
               function(x) { x$data })

  res.boot <- NULL
  res.boot.summ <- NULL
  res.boot.d <- NULL

  if (bootstrap) {
    res.boot <- plyr::llply(res.all,
                      function(x) {
                        x$boot$values
                      })
    res.boot.d <- plyr::llply(res.all,
                        function(x) {
                          x$boot$d.hat.sum
                        })
    names(res.boot) <- plyr::llply(res.all, function(x) { x$boot$name })
    names(res.boot.d) <- plyr::llply(res.all, function(x) { x$boot$name })

  }

  ## now compute summaries of the prediction errors
  errs <- estimate.error(estimate=res$nsum.holdout.est,
                         truth=res$known.size)
  res <- cbind(res,errs)

  mse <- mean(res$sqerr)
  rmse <- sqrt(mse)
  are <- mean(abs(res$relerr))
  mae <- mean(res$abserr)

  if(return.plot) {

    ## to placate R CMD CHECK
    name <- NULL
    known.size <- NULL
    nsum.holdout.est <-  NULL

    iv.plot <- ggplot(res) +
               geom_text(aes(x=known.size, y=nsum.holdout.est, label=name),
                         ##color=alpha("red", 0.3), size=3) +
                         ##color="black",
                         size=3) +
               coord_equal(ratio=1) +
               xlim(with(res,
                         range(known.size,nsum.holdout.est,na.rm=TRUE))) +
               ylim(with(res,
                         range(known.size,nsum.holdout.est,na.rm=TRUE))) +
               geom_abline(intercept=0, slope=1) +
               ggtitle(paste("Hold-out estimate versus known popn size")) +
               xlab("known population size") +
               ylab("hold-out NSUM estimate\nof popn size")
  } else {
    iv.plot <- NULL
  }

  return(list(results=res,
              mse=mse,
              rmse=rmse,
              are=are,
              mae=mae,
              plot=iv.plot,
              boot=res.boot,
              boot.d.hat=res.boot.d))

}

#####################################################
##' plot_meanties_truth
##'
##' plot the relationship between the mean number of
##' ties in the survey dataset and the true popn sizes
##'
##' TODO - more in-depth description of this function
##'
##' @param survey.data the dataframe with the survey results
##' @param known.popns if not NULL, a vector whose entries are the
##' size of the known
##' populations, and whose names are the variable names in the dataset
##' corresponding to each one. if NULL, then assume that the survey.data
##' dataframe has an attribute called 'known.popns' containing this vector.
##' @param weights if not NULL, weights to use in computing the estimate. this
##'                should be the name of the column in the survey.data which
##'               has the variable with the appropriate weights. these weights
##'                should be construted so that, eg, the mean of the degrees
##'                is estimated as (1/n) * \\sum_i {w_i * d_i}
##' @return a ggplot2 object with the relationship plot
##' @export
plot_meanties_truth <- function(survey.data, weights=NULL, known.popns=NULL)
{

  if (is.null(known.popns)) {
    known.popns <- attr(survey.data, "known.popns")
  }

  ## weights will default to 1 for everyone, unless the user specified
  ## a weights variable
  weights <- surveybootstrap:::get.weights(survey.data, weights)

  ard.q <- subset(survey.data, select=names(known.popns))

  ard.means <- colMeans(ard.q*weights, na.rm=TRUE)

  res <- cbind(ard.means,
               known.popns[match(names(known.popns),names(ard.means))])
  colnames(res)[2] <- "truth"
  res <- data.frame(res)
  res$name <- rownames(res)

  ## to placate R CMD CHECK
  truth <- NULL
  name <- NULL

  resplot <- ggplot(res) +
             geom_text(aes(x=truth, y=ard.means,label=name), size=4) +
             xlab("known population size") +
             ylab("average response to ARD question")

  return(list(plot=resplot, data=res))
}
