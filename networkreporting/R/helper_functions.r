#######################################################################
## helper_functions.R
##
## various helper functions for use in the networkreporting
## package
##

##########################################################################
##' turn a dataframe into a known population vector
##'
##' \code{df.to.kpvec} takes a dataframe which has a column with
##' known population names, and a column with known population
##' totals, and turns it into a known population vector. if the
##' names of the survey variables corresponding to each known population
##' are available, they can be passed in as well
##'
##' @param kp.data the known population dataset
##' @param kp.var the column of \code{kp.data} that has known population names;
##'               either a column name, a column index, or a vector of values
##' @param kp.value the column of \code{kp.data} that has known population sizes;
##'               either a column name, a column index, or a vector of value
##' @return a vector whose entries have the known population values and whose
##' names have the corresponding \code{kp.var} value
##' @export
##' @seealso \link{add.kp}
##' @examples \dontrun{
##'   ## see example in add.kp
##' }
##'
df.to.kpvec <- function(kp.data,
                        kp.var,
                        kp.value) {

  vals <- surveybootstrap:::get.var(kp.data, kp.value)
  var <- surveybootstrap:::get.var(kp.data, kp.var)

  kp.vec <- vals
  names(kp.vec) <- var

  return(kp.vec)

}

##########################################################################
##' attach known populations to a dataframe
##'
##' @description
##' take a known population vector (see \code{\link{df.to.kpvec}}) and
##' associate it with a survey dataframe. this makes it more convenient
##' to use some of the \code{networksampling} package's function
##'
##' @details
##' The \code{total.popn.size} parameter is interpreted as follows:
##' \itemize{
##' \item NA if total.popn.size is NA then work with proportions
##' \item NULL if total.popn.size is NULL (nothing passed in), then
##'          assume that there's a total.popn.size attribute
##'          associated with the dataset we're using
##' \item numerical value if an actual total.popn.size was passed in,
##'        use that value
##' }
##' 
##'
##' @param survey.data the survey dataframe
##' @param kp.vec the known population vector
##' @return the survey dataframe with the known population vector
##' attached as an attribute
##' @param total.pop.size (optional) the total population size to use (see below)
##' @export
##' @seealso \link{df.to.kpvec}
##' @examples \dontrun{
##'
##'   # if kp.dat is a dataframe with columns 'kp' with known popn names
##'   # and 'total.size' with the total size,
##'   # and my.survey is the dataframe with survey responses
##'
##'   my.kp.vec <- df.to.kpvec(kp.data, kp.var='kp', kp.value='total.size')
##'   my.survey <- add.kp(my.survey, my.kp.vec)
##'
##'   # now we can call estimator functions like
##'   # kp.degree.estimator without having to specify known
##'   # populations each time
##' }
add.kp <- function(survey.data, kp.vec, total.pop.size=NULL) {
  attr(survey.data, "known.popns") <- kp.vec

  if (! is.null(total.pop.size)) {
      attr(survey.data, "total.popn.size") <- total.pop.size
  }
    
  return(survey.data)
}


##########################################################################
##' topcode a vector of numerical values
##'
##' this function topcodes one vector; it's used by the \code{topcode}
##' function to topcode a set of columns in a data frame
##'
##' @param x the vector of values to topcode
##' @param max the maximum value; all values > max are recoded to max
##' @param to.na a vector of values to recode to NA (this happens before topcoding)
##' @param ignore a vector of values to leave unchanged
##' @return the topcoded vector
##' @export
##' @examples \dontrun{
##'    ## TODO write example
##' }
topcode.var <- function(x, max, to.na=NULL, ignore=NA) {

  if (! is.numeric(x)) {
    stop("You can only topcode a numeric vector.")
  }

  if (! is.null(to.na)) {
    x[x %in% to.na] <- NA
  }

  x[(x > max) & (! x %in% ignore)] <- max

  return(x)
}

##########################################################################
##' topcode a group of variables
##'
##' this function uses \code{topcode.var} to topcode a set of variables.
##' it's useful for topcoding a whole set of aggregated relational data
##' ("how many X are you connected to?") questions in the same way.
##'
##' @param survey.data  the dataset with the survey responses
##' @param vars a vector with the names or indices of the columns in the
##'             dataframe that are to be topcoded
##' @param max the maximum value; all values > max are recoded to max
##' @param ignore a vector of values to leave unchanged
##' @param to.na a vector of values to recode to NA (this happens before topcoding)
##' @return the topcoded vector
##' @export
##' @examples \dontrun{
##'    data(hh.survey) # example data included with the package
##'    example.survey <- topcode.data(example.survey,
##'                                   vars=known.popn.vars,
##'                                   max=30)
##' }
topcode.data <- function(survey.data, vars, max, to.na=NULL, ignore=NA) {

  ## TODO -- eventually check that vars are found in the columns of survey.data

  survey.data[,vars] <- plyr::colwise(topcode.var)(survey.data[,vars,drop=FALSE],
                                             max=max,
                                             to.na=to.na,
                                             ignore=ignore)

  return(survey.data)

}


##########################################################################
##' handle the total.popn.size argument in a uniform way across
##' several functions
##'
##' handle the total.popn.size argument in a uniform way across
##' several functions, including
##' \code{\link{kp.degree.estimator}},
##' \code{\link{nsum.internal.validation}}, and
##' \code{\link{nsum.estimator}}.
##'
##' The result depends upon the value that was passed in:
##' \itemize{
##' \item NA if total.popn.size is NA then work with proportions
##' \item NULL if total.popn.size is NULL (nothing passed in), then
##'          assume that there's a total.popn.size attribute
##'          associated with the dataset we're using
##' \item numerical value if an actual total.popn.size was passed in,
##        use that value
##' }
##'
##' @param total.popn.size value to parse
##' @param survey.data the dataframe we're analyzing, which may or may not
##'                    have an attribute called 'total.popn.size'
##' @param verbose if TRUE, print messages to the screen
##' @return the parsed total population size
##' @keywords internal
parse.total.popn.size <- function(total.popn.size, survey.data, verbose=FALSE) {

  ## this is a little complicated.
  ## (see also kp.degree.estimator)
  if ((! is.null(total.popn.size)) && is.na(total.popn.size)) {

    surveybootstrap:::vcat(verbose, "working in proportions\n")

  } else if (is.null(total.popn.size)) {

    total.popn.size <- attr(survey.data, "total.popn.size")

    if(! is.numeric(total.popn.size)) {
      stop("error - no suitable attribute 'total.popn.size' for dataframe.\n")
    } else {
      surveybootstrap:::vcat(verbose, "using dataframe's attribute for total population size.\n")
    }

  } else {
    surveybootstrap:::vcat(verbose, "working in absolute numbers\n")
  }

  return(total.popn.size)
}

##########################################################################
##' given an estimated subpopn size or prevalence and the correct value,
##' produce some measurements of how close the esimate is
##'
##' @param estimate the estimate
##' @param truth the correct answer
##' @return a vector whose entries have various summaries of fit
##' @export
##' @examples \dontrun{
##'    ## TODO add example
##' }
estimate.error <- function(estimate, truth) {

  err <- estimate - truth
  abserr <- abs(err)
  sqerr <- err^2
  relerr <- abserr/truth

  return(cbind(err=err,abserr=abserr,sqerr=sqerr,relerr=relerr))

}
