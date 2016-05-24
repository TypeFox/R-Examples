################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


#' Formula Response for Recurrent Event Data
#' 
#' \code{Survr} is an S3 class that represents
#' formula response for recurrent event data
#' modeled by methods based on counts and rate function.
#' The last letter 'r' in 'Survr' represents 'rate'.
#'
#' This is a similar function to \code{Survr} in package
#' \pkg{survrec} but with a better embedded checking procedure for
#' recurrent event data modeled by methods
#' based on counts and rate function.
#' The checking rules include that
#' \itemize{
#'     \item Identification of each subject cannot be missing.
#'     \item Event indicator must be coded as 0 (censored) or 1 (event).
#'     \item Event time and censoring time cannot be missing.
#'     \item Each subject must have one and only one censoring time.
#'     \item Event time cannot not be later than censoring time.
#' }
#'  
#' @param ID Identificator of each subject. 
#' @param time Time of reccurence event or censoring.
#' @param event The status indicator, 0 = censored, 1 = event. 
#' @aliases Survr
#' @seealso \code{\link{rateReg}} for model fitting.
#' @export
Survr <- function (ID, time, event) {
    inpDat <- data.frame(ID, time, event)
    dat <- check_Survr(inpDat)
    outDat <- with(dat, as.matrix(cbind(ID, time, event)))
    attr(outDat, "ID") <- attr(dat, "ID")
    oldClass(outDat) <- "Survr"
    invisible(outDat)
}


#' An S4 Class to Represent a Fitted Model
#' 
#' \code{rateReg-class} is an S4 class that represents a fitted model. 
#' \code{\link{rateReg}} produces objects of this class.  
#' See ``Slots'' for details.
#' 
#' @slot call Function call.
#' @slot formula Formula.
#' @slot nObs A positive integer
#' @slot knots A numeric vector.
#' @slot boundaryKnots A numeric vector.
#' @slot degree A nonnegative integer.
#' @slot df A list of nonnegative numeric vectors.
#' @slot estimates A list.
#' @slot control A list.
#' @slot start A list.
#' @slot na.action A length-one character vector.
#' @slot xlevels A list.
#' @slot contrasts A list.
#' @slot convergCode A nonnegative integer.
#' @slot logL A numeric value.
#' @slot fisher A numeric matrix.
#' @aliases rateReg-class
#' @seealso \code{\link{rateReg}} for details of slots.
#' @export
setClass(Class = "rateReg", 
         slots = c(call = "call", 
                   formula = "formula",
                   nObs = "integer",
                   knots = "numeric",
                   boundaryKnots = "numeric",
                   degree = "integer",
                   df = "list",
                   estimates = "list",
                   control = "list",
                   start = "list",
                   na.action = "character",
                   xlevels = "list",
                   contrasts = "list",
                   convergCode = "integer",
                   logL = "numeric",
                   fisher = "matrix"),
         validity = function (object) {
             ## check on nObs
             if (object@nObs <= 0) {
                 return("Number of observations must be a positive integer.")
             }
             ## check on knots
             if (length(object@knots) > 0) { # if there exists internal knots
                 if (min(object@knots) < min(object@boundaryKnots) ||
                     max(object@knots) > max(object@boundaryKnots)) {
                     return(paste("Internal knots must all lie in the", 
                                  "coverage of boundary knots."))
                 }
             }
             ## check on degree
             if (object@degree < 0) {
                 return("Degree of spline bases must be a nonnegative integer.")
             }
             ## check on df
             dfVec <- do.call("c", object@df)
             dfValid <- is.integer(dfVec) && all(dfVec >= 0)
             if (! dfValid) {
                 return("Degree of freedom must be nonnegative integers.")
             }
             ## else return
             TRUE
         })


#' An S4 Class to Represent Summary of a Fitted Model
#' 
#' \code{summaryRateReg-class} is an S4 class with selective slots 
#' of \code{rateReg-class} object.  See ``Slots'' for details.  
#' \code{\link{summary,rateReg-method}} produces objects of this class. 
#'  
#' @slot call Function call.
#' @slot knots A numeric vector.
#' @slot boundaryKnots A numeric vector.
#' @slot covarCoef A numeric matrix.
#' @slot frailtyPar A numeric matrix.
#' @slot degree A nonnegative integer.
#' @slot baseRateCoef A numeric matrix.
#' @slot logL A numeric value.
#' @aliases summaryRateReg-class
#' @seealso \code{\link{summary,rateReg-method}} for details of slots.
#' @export
setClass(Class = "summaryRateReg", 
         slots = c(call = "call", 
                   knots = "numeric",
                   boundaryKnots = "numeric",
                   covarCoef = "matrix",
                   frailtyPar = "matrix",
                   degree = "integer",
                   baseRateCoef = "matrix",
                   logL = "numeric"),
         validity = function (object) {
             ## check on knots
             if (length(object@knots) > 0) { # if there exists internal knots
                 if (min(object@knots) < min(object@boundaryKnots) ||
                     max(object@knots) > max(object@boundaryKnots)) {
                     return(paste("Internal knots must all lie in the", 
                                  "coverage of boundary knots."))
                 }
             }
             ## check on degree
             if (object@degree < 0) {
                 return("Degree of spline bases must be a nonnegative integer.")
             }
             ## else return
             TRUE
         })


#' An S4 Class to Represent Sample MCF
#' 
#' An S4 class that represents sample mean cumulative function (MCF) from data.
#' \code{\link{mcf}} produces objects of this class.  
#'
#' @slot call Function call
#' @slot formula Formula.
#' @slot na.action A length-one character vector.
#' @slot level A numeric value.
#' @slot MCF A data frame.
#' @slot multiGroup A logical value.
#' @aliases sampleMcf-class
#' @seealso \code{\link{mcf,formula-method}} for details of slots.
#' @importFrom methods setClass
#' @export
setClass(Class = "sampleMcf", 
         slots = c(call = "call",
                   formula = "formula",
                   na.action = "character",
                   level = "numeric",
                   MCF = "data.frame", 
                   multiGroup = "logical"))


#' An S4 Class to Respresent Estimated MCF from a Fitted Model
#' 
#' An S4 class that represents estimated mean cumulative function (MCF) 
#' from Models.
#' \code{\link{mcf}} produces objects of this class.  
#' 
#' @slot call Function call.
#' @slot formula Formula.
#' @slot knots A numeric vector.
#' @slot degree A nonnegative integer.
#' @slot boundaryKnots A numeric vector.
#' @slot newdata A numeric matrix.
#' @slot MCF A data frame.
#' @slot level A numeric value between 0 and 1.
#' @slot na.action A length-one character vector.
#' @slot control A list.
#' @slot multiGroup A logical value. 
#' @aliases rateRegMcf-class
#' @seealso \code{\link{mcf,rateReg-method}} for details of slots.
#' @export
setClass(Class = "rateRegMcf", 
         slots = c(call = "call",
                   formula = "formula",
                   knots = "numeric",
                   degree = "integer",
                   boundaryKnots = "numeric",
                   newdata = "matrix",
                   MCF = "data.frame",
                   level = "numeric", 
                   na.action = "character",
                   control = "list", 
                   multiGroup = "logical"),
         validity = function (object) {
             ## check on knots
             if (length(object@knots) > 0) { # if there exists internal knots
                 if (min(object@knots) < min(object@boundaryKnots) ||
                     max(object@knots) > max(object@boundaryKnots)) {
                     return(paste("Internal knots must all lie in the", 
                                  "coverage of boundary knots."))
                 }
             }
             ## check on degree
             if (object@degree < 0) {
                 return("Degree of spline bases must be a nonnegative integer.")
             }
             ## check on level
             if (object@level <= 0 || object@level >= 1) {
                 return("Confidence level mush be between 0 and 1.")
             }
             ## else return
             TRUE
         })

