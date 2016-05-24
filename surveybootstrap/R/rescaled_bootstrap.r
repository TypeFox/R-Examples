
#####################################################
##' rescaled.bootstrap.sample
##'
##' C++ version: given a survey dataset and a description of the survey
##' design (ie, which combination of vars determines primary sampling
##' units, and which combination of vars determines strata), take
##' a bunch of bootstrap samples for the rescaled bootstrap estimator
##' (see, eg, Rust and Rao 1996).
##'
##' Note that we assume that the formula uniquely specifies PSUs.
##' This will always be true if the PSUs were selected without replacement.
##' If they were selected with replacement, then it will be necessary
##' to make each realization of a given PSU in the sample a unique id.
##' Bottom line: the code below assumes that all observations within
##' each PSU (as identified by the design formula) are from the same draw
##' of the PSU.
##'
##' The rescaled bootstrap technique works by adjusting the
##' estimation weights based on the number of times each
##' row is included in the resamples. If a row is never selected,
##' it is still included in the returned results, but its weight
##' will be set to 0. It is therefore important to use estimators
##' that make use of the estimation weights on the resampled
##' datasets.
##'
##' We always take m_i = n_i - 1, according to the advice presented
##' in Rao and Wu (1988) and Rust and Rao (1996).
##'
##' @param survey.data the dataset to use
##' @param survey.design a formula describing the design of the survey (see below - TODO)
##' @param num.reps the number of bootstrap replication samples to draw
##' @param parallel if TRUE, use parallelization (via \code{plyr})
##' @param paropts an optional list of arguments passed on to \code{plyr} to control
##'        details of parallelization
##' @return a list with \code{num.reps} entries. each entry is a dataset which
##' has at least the variables \code{index} (the row index of the original
##' dataset that was resampled) and \code{weight.scale}
##' (the factor by which to multiply the sampling weights
##' in the original dataset).
##' @details \code{survey.design} is a formula of the form\cr
##'    weight ~ psu_vars + strata(strata_vars),
##' where weight is the variable with the survey weights and psu
##' is the variable denoting the primary sampling unit
##' @export
rescaled.bootstrap.sample <- function(survey.data,
                                      survey.design,
                                      parallel=FALSE,
                                      paropts=NULL,
                                      num.reps=1)
{

  survey.data$.internal_id <- 1:nrow(survey.data)

  design <- parse_design(survey.design)

  ## drop the "~" at the start of the formula
  psu.vars <- design$psu.formula[c(-1)][[1]]

  ## in the special case where there are no PSU vars, treat each row as
  ## its own PSU
  if (length(psu.vars)==1 & psu.vars=="1") {
    psu.vars <- as.name(".internal_id")
  }

  ## create a single variable with an id number for each PSU
  ## (we need this to use the C++ code, below)
  survey.data$.cluster_id <- group_indices_(survey.data, .dots=all.vars(psu.vars))

  ## if no strata are specified, enclose the entire survey all in
  ## one stratum
  if (is.null(design$strata.formula)) {
    strata <- list(survey.data)
  } else {
    strata <- plyr::dlply(survey.data, design$strata.formula, identity)
  }

  ## get num.reps bootstrap resamples within each stratum,
  ## according to the rescaled bootstrap scheme
  ## (see, eg, Rust and Rao 1996)

  ## this llply call returns a list, with one entry for each stratum
  ## each stratum's entry contains a list with the bootstrap resamples
  ## (see the note for the inner llply call below)
  bs <- plyr::llply(strata,
              function(stratum.data) {

                ## (this part is written in c++)
                res <- resample_stratum(stratum.data$.cluster_id,
                                        num.reps)

                colnames(res) <- paste0("rep.", 1:ncol(res))
                res <- cbind("index"=stratum.data$.internal_id,
                             res)
                
                return(res)
              })

  ## bs: list, one entry for each stratum
  ## each stratum's entry is a matrix.
  ## first column of the matrix is called 'index',
  ## which is the row number for the observation in the
  ## original dataset; there is one remaining column for each
  ## bootstrap resample. the entries of each column are the factors
  ## by which the original weights should be scaled

  bs.all <- do.call("rbind", bs)

  res <- plyr::alply(bs.all[,-1],
               2,
               function(this_col) {
                   return(data.frame(index=bs.all[,1],
                                     weight.scale=this_col))
               })

  return(res)

}

#####################################################
##' rescaled.bootstrap.sample.pureR
##'
##' (this is the pure R version; it has been supplanted by
##'  \code{rescaled.bootstrap.sample}, which is partially written in C++)
##' 
##' given a survey dataset and a description of the survey
##' design (ie, which combination of vars determines primary sampling
##' units, and which combination of vars determines strata), take
##' a bunch of bootstrap samples for the rescaled bootstrap estimator
##' (see, eg, Rust and Rao 1996).
##'
##' Note that we assume that the formula uniquely specifies PSUs.
##' This will always be true if the PSUs were selected without replacement.
##' If they were selected with replacement, then it will be necessary
##' to make each realization of a given PSU in the sample a unique id.
##' Bottom line: the code below assumes that all observations within
##' each PSU (as identified by the design formula) are from the same draw
##' of the PSU.
##'
##' The rescaled bootstrap technique works by adjusting the
##' estimation weights based on the number of times each
##' row is included in the resamples. If a row is never selected,
##' it is still included in the returned results, but its weight
##' will be set to 0. It is therefore important to use estimators
##' that make use of the estimation weights on the resampled
##' datasets.
##'
##' We always take m_i = n_i - 1, according to the advice presented
##' in Rao and Wu (1988) and Rust and Rao (1996).
##'
##' @param survey.data the dataset to use
##' @param survey.design a formula describing the design of the survey (see below - TODO)
##' @param num.reps the number of bootstrap replication samples to draw
##' @param parallel if TRUE, use parallelization (via \code{plyr})
##' @param paropts an optional list of arguments passed on to \code{plyr} to control
##'        details of parallelization
##' @return a list with \code{num.reps} entries. each entry is a dataset which
##' has at least the variables \code{index} (the row index of the original
##' dataset that was resampled) and \code{weight.scale}
##' (the factor by which to multiply the sampling weights
##' in the original dataset).
##' @details \code{survey.design} is a formula of the form\cr
##'    weight ~ psu_vars + strata(strata_vars),
##' where weight is the variable with the survey weights and psu
##' is the variable denoting the primary sampling unit
##' @export
rescaled.bootstrap.sample.pureR <- function(survey.data,
                                      survey.design,
                                      parallel=FALSE,
                                      paropts=NULL,
                                      num.reps=1)
{

  survey.data$.internal_id <- 1:nrow(survey.data)

  design <- parse_design(survey.design)

  ## drop the "~" at the start of the formula
  psu.vars <- design$psu.formula[c(-1)][[1]]

  ## in the special case where there are no PSU vars, treat each row as
  ## its own PSU
  if (length(psu.vars)==1 & psu.vars=="1") {
    psu.vars <- as.name(".internal_id")
  }

  ## if no strata are specified, enclose the entire survey all in
  ## one stratum
  if (is.null(design$strata.formula)) {
    strata <- list(survey.data)
  } else {
    strata <- plyr::dlply(survey.data, design$strata.formula, identity)
  }

  ## get num.reps bootstrap resamples within each stratum,
  ## according to the rescaled bootstrap scheme
  ## (see, eg, Rust and Rao 1996)

  ## this llply call returns a list, with one entry for each stratum
  ## each stratum's entry contains a list with the bootstrap resamples
  ## (see the note for the inner llply call below)
  bs <- plyr::llply(strata,
              function(stratum.data) {

                ## figure out how many PSUs we have in our sample
                psu.count <- count(stratum.data,
                                   psu.vars)

                n.h <- nrow(psu.count)

                ## take m_h = n_h - 1, which is what the literature
                ## most commonly recommends
                m.h <- n.h - 1

                ## this llply call returns a list, with one entry for each bootstrap rep.
                ## each list entry has a data frame with the same number of rows as
                ## stratum.data, 
                ## and with colums for 
                ## the survey design variables, the .internal_id,
                ## r.hi, and weight.scale
                resamples <- plyr::llply(1:num.reps,
                                   ## for each bootstrap rep, this function returns
                                   function(rep) {

                                     ## sample m.h PSUs with replacement
                                     these.psu.samples <- sample(1:nrow(psu.count),
                                                                 m.h,
                                                                 replace=TRUE)

                                     ## r.hi is the number of times PSU i in stratum
                                     ## h was chosen in our resample
                                     r.hi <- count(data.frame(psu.row=these.psu.samples))

                                     r.hi$weight.scale <- r.hi$freq * (n.h / m.h)

                                     psu.count$freq <- NULL

                                     psu.count$r.hi <- 0
                                     psu.count$weight.scale <- 0

                                     psu.count[ r.hi$psu.row, "r.hi" ] <- r.hi$freq

                                     ## this is the factor by which we
                                     ## need to multiply
                                     ## the sampling weights
                                     ## (again, see, eg, Rust and Rao 1996, pg 292)
                                     psu.count[ r.hi$psu.row,
                                                "weight.scale" ] <- r.hi$weight.scale

                                     this.resample <- merge(stratum.data[,c(all.vars(psu.vars),
                                                                            ".internal_id")],
                                                            psu.count,
                                                            by=all.vars(psu.vars))
                                     this.resample$.internal_id.1 <- NULL

                                     return(this.resample)
                                   },
                                   .parallel=parallel,
                                   .paropts=paropts)

                return(resamples)
              })

  # now reassemble each stratum...
  res <- plyr::llply(1:num.reps,
               function(rep.idx) {
                 this.rep <- plyr::ldply(bs,
                                   function(this.stratum.samples) {
                                     return(this.stratum.samples[[rep.idx]])
                                   })

                 this.rep <- plyr::rename(this.rep,
                                    c(".internal_id"="index"))

                 return(this.rep)
               })

  return(res)
}

