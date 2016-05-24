# --------------------------------------------
# Generic function for Expectations
# --------------------------------------------

#' Get expectations of a randomization list
#' 
#' Generates a matrix of the expectations of the included patients in the 
#' clinical trial.
#' 
#' @param randSeq object of the class randSeq.
#' @param issue object of the class issue (optional).
#' @param endp object of the class endpoint (optional).
#'
#' @details
#' It is assumed that the expectations of the included patients in a clinical trial
#' can be influenced in three different ways:
#' \itemize{
#' \item The strength of selection bias and the guessing strategy
#' of the investigator (see \code{\link{selBias}}).
#' \item The strength of a linear time trend, which is described by an object
#'  of the class \code{\link{chronBias}}.
#' \item The expectations of the investigated treatement groups can be different 
#' (see e.g. \code{\link{normEndp}}).
#' }
#'
#' @examples 
#' myPar <- bsdPar(10, 2)
#' M <- genSeq(myPar, 2)
#' cs <- selBias("CS", 2, "sim")
#' endp <- normEndp(mu = c(2, 2), sigma = c(1, 1))
#' getExpectation(M, cs, endp)
#' 
#' @name getExpectation
NULL


#' @rdname getExpectation
#' 
#' @export
setGeneric("getExpectation", function(randSeq, issue, endp) standardGeneric("getExpectation"))


