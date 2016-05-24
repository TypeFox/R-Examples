
#####################################################
##' srs.bootstrap.sample
##'
##' given a survey dataset and a description of the survey
##' design (ie, which combination of vars determines primary sampling
##' units, and which combination of vars determines strata), take
##' a bunch of bootstrap samples under a simple random sampling
##' (with repetition) scheme
##'
##' @param survey.data the dataset to use
##' @param num.reps the number of bootstrap replication samples to draw
##' @param parallel if TRUE, use parallelization (via \code{plyr})
##' @param paropts an optional list of arguments passed on to \code{plyr} to control
##'        details of parallelization
##' @param ... ignored, but useful because it allows params like
##\code{survey.design},
##' which are used in other bootstrap designs, to be passed in without error
##' @return a list with \code{num.reps} entries. each entry is a dataset which has
##' at least the variables \code{index} (the row index of the original dataset that
##' was resampled) and \code{weight.scale} (the factor by which to multiply the
##' sampling weights in the original dataset).
##'
##' @export
srs.bootstrap.sample <- function(survey.data,
                                 num.reps=1,
                                 parallel=FALSE,
                                 paropts=NULL,
                                 ...)
{

  survey.data$.internal_id <- 1:nrow(survey.data)

  res <- plyr::llply(1:num.reps,
               function(rep.idx) {

                 these.samples <- sample(1:nrow(survey.data),
                                         nrow(survey.data),
                                         replace=TRUE)

                 this.rep <- data.frame(index=these.samples,
                                        weight.scale=1)

                 return(this.rep)
               },
               .parallel=parallel,
               .paropts=paropts)

  return(res)
}
