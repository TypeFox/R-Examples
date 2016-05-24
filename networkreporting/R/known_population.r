#####################################################
## known_population.r
##
## the known population method for estimating respondents'
## degrees
##

#####################################################
##' Average personal network size estimates using known population method
##'
##' If given \code{attribute.names}, then this function produces
##' estimated average network sizes given by the groups that are defined by 
##' all combinations of the attributes; otherwise, it estimates the
##' average personal network size for the entire frame population.
##' 
##' @section Technical note:
##' The estimated average degree is 
##' \eqn{(\sum y_{F_\alpha, A} / N_A) \times
##'      N_F / N_{F_\alpha}}
##' here, we estimate \eqn{N_F / N_{F_\alpha}} by dividing the
##' total of all respondents' weights by the sum of
##' the weights for respondents in each cell \eqn{\alpha}.
##'
##' @section TODO:
##' \itemize{
##' \item{ handle case where attribute.names is NULL (should compute overall average)}
##' \item{ handle missing values }
##' \item{ integrate the individual-level estimator above, kp.degree.estimator}
##' \item{ finish documentation for NSE version }
##' \item{ make unit tests }
##' \item{ think about how to elegantly add options for dbar_(P,Q) vs dbar_(Q,P)}
##' }
##'
##' @param resp.data the dataframe that has the survey responses
##' @param known.populations  the names of the columns in \code{resp.data}
##'          that have respondents' reports about connections to known populations
##' @param attribute.names the names of the columns in \code{resp.data} that
##'                        determine the subgroups for which average degree is estimated
##' @param weights weights to use in computing the estimate
##' @param total.kp.size the size of the probe alters; i.e., the sum of the known population
##'                      sizes. if NULL, then this is set to 1
##' @param alter.popn.size the size of the population of alters; this is most
##'        often the frame population, which is the default if nothing else is
##'        specified; the size of the frame population is taken to be the sum
##'        of the weights over all of resp.data
##' @return the estimated average degree for respondents in each
##'         of the categories given by \code{attribute.names}
##'
##' @rdname kp.estimator
##' @export
kp.estimator_ <- function(resp.data, 
                          known.populations,
                          attribute.names,
                          weights,
                          total.kp.size=NULL,
                          alter.popn.size=NULL) {

  wdat <- select_(resp.data, .dots=weights)
  kpdat <- select_(resp.data, .dots=known.populations)
  adat <- select_(resp.data, .dots=attribute.names)

  alter.popn.size <- ifelse(is.null(alter.popn.size) ||
                            is.null(lazy_eval(alter.popn.size)),
                            sum(wdat[,1]),
                            lazy_eval(alter.popn.size))

  total.kp.size <- ifelse(is.null(total.kp.size),
                          1,
                          lazy_eval(total.kp.size))

  # get individual-level sums for known population connections
  # (NB: might want to eventually allow for this to be adjusted by
  #      some kind of response / reporting model?)
  kptot <- data_frame(kptot=rowSums(kpdat))

  df <- bind_cols(kptot, wdat, adat)

  atnames <- paste(colnames(adat))

  # now aggregate by attributes
  agg <- report.aggregator_(resp.data=df,
                            attribute.names=atnames,
                            qoi='kptot',
                            weights=weights,
                            qoi.name='y.kp')

  tograb <- lapply(c(colnames(adat),
                     'sum.y.kp', 'wgt.total.y.kp', 'num.obs.y.kp'),
                   as.symbol)

  ## to placate R CMD CHECK
  sum.y.kp <- NULL
  sum.y.kp.over.kptot <- NULL
  wgt.total.y.kp <- NULL

  res <- select_(agg, .dots=tograb) %>%
         dplyr::mutate(sum.y.kp.over.kptot = sum.y.kp / total.kp.size,
                ## here, we estimate N_F / N_{F_\alpha} by dividing the
                ## total of all respondents' weights by the sum of
                ## the weights for respondents in each cell \alpha
                ## (see technical note in documentation)
                dbar.Fcell.F = sum.y.kp.over.kptot * 
                               (alter.popn.size / wgt.total.y.kp))

  return(res)

}

#####################################################
##' @rdname kp.estimator
kp.estimator <- function(resp.data, 
                         known.populations,
                         attribute.names,
                         weights,
                         total.kp.size=1,
                         alter.popn.size=NULL) {

    kp.estimator_(resp.data,
                  known.populations=lazy(known.populations),
                  attribute.names=lazy(attribute.names, resp.data),
                  weights=lazy(weights, resp.data),
                  total.kp.size=lazy(total.kp.size),
                  alter.popn.size=lazy(alter.popn.size))

}



#####################################################
##' Individual personal network size estimates using the known population method
##'
##' In most situations, the known population method will be
##' used to estimate the average personal network size;
##' this can be done with \code{\link{kp.estimator_}}. If, instead, you wish
##' to estimate the personal network size of each individual
##' respondent, then you can use this function.
##'
##' Note that this is not making inference about any larger population;
##' it estimates a property of each individual respondent. So the
##' sampling weights are not used here.
##'
##' @section TODO:
##' \itemize{
##' \item{ handle missing values! }
##' \item{ make unit tests }
##' }
##'
##' @param resp.data the respondent (survey) data
##' @param known.populations the names of the known populations
##' @param total.kp.size the sum of the sizes of all of the known populations
##' @param alter.popn.size the size of the population respondents
##'        are reporting about connections to; typically this will
##'        be the frame population, so \code{alter.popn.size} should
##'        be the size of the frame population, N.F
##' @return a data frame with an estimate of each individual respondent's personal
##'         network size
##' @rdname kp.individual.estimator
##' @export
kp.individual.estimator <- function(resp.data, 
                                    known.populations,
                                    total.kp.size=1,
                                    alter.popn.size) {

    resp.data$.rowid <- 1:nrow(resp.data)
    resp.data$.noweight <- 1

    res <- kp.estimator_(resp.data,
                         known.populations=lazy(known.populations),
                         attribute.names=".rowid",
                         weights=".noweight",
                         total.kp.size=lazy(total.kp.size),
                         alter.popn.size=lazy(alter.popn.size))

}

#####################################################
##' @rdname kp.individual.estimator
##' @export 
kp.individual.estimator_ <- function(resp.data, 
                                     known.populations,
                                     total.kp.size=1,
                                     alter.popn.size) {

    resp.data$.rowid <- 1:nrow(resp.data)
    resp.data$.noweight <- 1

    res <- kp.estimator_(resp.data,
                         known.populations=known.populations,
                         attribute.names=".rowid",
                         weights=".noweight",
                         total.kp.size=total.kp.size,
                         alter.popn.size=alter.popn.size)

}


#####################################################
##' kp.degree.estimator (DEPRECATED)
##'
##' see \code{\link{kp.individual.estimator}} instead.
##'
##' compute an estimate of the respondents' degrees using
##' the known population method\cr
##'
##' note that this function does not take survey weights, since
##' these estimates are not for total degree, but just for the
##' individual degree of each respondent
##'
##' @param survey.data the dataframe with the survey results
##' @param known.popns if not NULL, a vector whose entries are the size of the known
##'                    populations, and whose names are the variable names in the dataset
##'                    corresponding to each one. if NULL, then assume that the survey.data
##'                    dataframe has an attribute called 'known.popns' containing this vector.
##' @param total.popn.size the size of the entire population. if NULL,
##'                        this function returns proportions; if not NULL, it
##'                        returns absolute numbers (ie, the proportions * total popn size)
##' @param missing if "ignore", then proceed with the analysis without
##'                doing anything about missing values. if "complete.obs"
##'                then, for each row, use only the known populations
##'                that have no missingness for the
##'                computations. care must be taken in using this second option
##' @param verbose if TRUE, print messages to the screen
##' @return a vector with an estimate of the degree for each row
##'         in survey.data. if missing=="ignore", then the degree for rows that have
##'         missingness in the 'how many X' questions will be set
##'         to NA
##' @export
kp.degree.estimator <- function(survey.data,
                                known.popns=NULL,
                                total.popn.size=NULL,
                                missing="ignore",
                                verbose=FALSE)
{

  ## TODO - PLAN: this should be a special wrapper around the
  ## kp.estimator_ function which produces estimates for the
  ## individual-level degree
  warning("kp.degree.estimator will be deprecated soon!")

  #### TODO
  #### this is returning a matrix now (see, eg, the rwanda mortality analysis)
  #### it should return a vector

  if (! missing %in% c("ignore", "complete.obs")) {
    stop("error in specifying procedure for handling missing values in kp.degree.estimator. invalid option.\n")
  }


  if (is.null(known.popns)) {
    known.popns <- attr(survey.data, "known.popns")
  }

  total.popn.size <- parse.total.popn.size(total.popn.size,
                                           survey.data,
                                           verbose=verbose)

  if (missing == "complete.obs") {

    ## use the modified estimator: for each row, use the known
    ## population estimator, taking populations whose responses we
    ## don't have (ie are missing) out of the numerator and
    ## denominator

    kp.dat <- survey.data[, names(known.popns)]

    ## mask for missing values: this matrix has the same shape
    ## as kp.dat, but its entries are 1 if the entry in kp.dat is
    ## observed and 0 if missing
    miss.mask <- data.matrix(as.data.frame(llply(kp.dat,
                                                 function(x) {
                                                   as.numeric(!is.na(x))
                                                 })))

    ## use the missing mask to get the sum of the N_k's for each
    ## indiviudal respondent
    ind.overall.known <- miss.mask %*% known.popns

    ind.tot.known <- (rowSums(kp.dat, na.rm=TRUE))

    res <- ind.tot.known/ind.overall.known

  } else {

    tot.known <- (rowSums(subset(survey.data,
                                 select=names(known.popns))))

    overall.known <- sum(known.popns)

    res <- tot.known/overall.known

  }

  if (! is.na(total.popn.size)) {
    res <- res * total.popn.size
  }

  return(res)

}



