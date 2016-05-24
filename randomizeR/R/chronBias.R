#' @include doublyT.R
#' @include getStat.R
#' @include testDec.R
#' @include getExpectation.R
#' @include endpoint.R
NULL

###############################################
# --------------------------------------------#
# Class chronBias                             #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the chronBias class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateChronBias <- function(object) {
  errors <- character()
  lengthType <- length(object@type)
  if (lengthType != 1) {
    msg <- paste("type is length ", lengthType, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  type <- object@type[1]
  if (!(type %in% c("linT", "stepT", "logT"))) {
    msg <- paste("(First) Argument of type is ", type, ". Should be in linT, 
                 logT, or stepT.", sep = "")
    errors <- c(errors, msg)
  }
  
  lengthMethod <- length(object@method)
  if (lengthMethod != 1) {
    msg <- paste("method is length ", lengthMethod, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  method <- object@method[1]
  if (!(method %in% c("exact", "sim"))) {
    msg <- paste("(First) Argument of method is ", method, ". Should be exact or sim."
                 , sep = "")
    errors <- c(errors, msg)
  }
  
  lengthAlpha<- length(object@alpha)
  if (lengthAlpha != 1) {
    msg <- paste("alpha is length ", lengthAlpha, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  alpha <- object@alpha[1]
  if (!(alpha >= 0 && alpha <= 1 )) {
    msg <- paste("(First) Argument of alpha is ", alpha, ". Should be in [0,1]."
                 , sep = "")
    errors <- c(errors, msg)
  }
   
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for chronBias
# --------------------------------------------

# Randomization paramters generic
setClass("chronBias",
         slots = c("type" = "character", "theta" = "numeric",
                   method = "character", alpha = "numeric"),
         validity = validateChronBias)

# --------------------------------------------
# Constructor function for chronBias
# --------------------------------------------

#' Representing chronological bias
#' 
#' Represents the issue of chronological bias in a clinical trial.
#' 
#' @inheritParams overview
#' @param method character string, should be one of \code{"sim"} or \code{"exact"}, see Description.
#' @param type character string, should be one of "\code{linT}", "\code{logT}", or "\code{stepT}", 
#' see Details.
#' @param alpha  significance level
#'
#'
#' @details
#' Chronological bias can be an issue in the design of a clinical trial. The 
#' \code{chronBias} function is a constructor function
#' for an S4 object of the class \code{chronBias} representing the issue of
#' chronological bias, s.a. time trends, in a clinical trial. It supports two possible modes,
#' \code{method="sim"} and \code{method="exact"}, and three different types of trend. 
#'
#' If \code{method="sim"}, the object represents the simulated type-I-error rate given 
#'  the level \code{alpha}, the selection effect \code{eta} and the biasing 
#'  strategy \code{type}. When calling \code{assess} for a \code{selBias} object 
#'  with \code{method="sim"}, one test decision is computed for each sequence of
#' \code{randSeq}. The type-I-error rate (power) is the proportion of falsely
#' (correctly) rejected null hypotheses.
#' 
#' If \code{method="exact"}, the object represents the exact type-I-error proabability 
#'  given the level \code{alpha}, the selection effect \code{eta} and the 
#'  biasing strategy \code{type}. When calling \code{assess} for a \code{selBias} 
#'  object with \code{method="exact"}, the exact \emph{p}-value of each 
#'  randomization sequence is computed. So far, this is only supported for
#'  normal endpoints. Then the type-I-error probability is
#'  the sum of the corresponding quantiles of the doubly noncentral t-distribution.
#'
#' \subsection{Types of chronological bias}{
#' \describe{
#' 	 \item{\code{type = "linT"}}{
#'    Represents linear time trend. Linear time trend means that the expected response
#' 	  of the patients increases evenly by \code{theta} with 
#'    every patient included in the study, until reaching \code{N theta} after \code{N} patients.
#'    Linear time trend may occur as a result of gradually relaxing in- or exlusion criteria 
#'    throughout the trial.
#'    It can be presented by the formula: \deqn{f(i) = i  \theta}{f(i) = i \theta}
#'   }
#' 	 \item{\code{type = "logT"}}{
#'    Represents logistic time trend. Logistic time trend means that the expected response
#' 	  of the patients increases logistically in the patient index by \code{theta} with 
#'    every patient included in the study, until reaching \code{log(N) theta} after \code{N} patients.
#'    Logistic time trend may occur as a result of a learning curve, i.e. in a surgical trial.
#'    It can be presented by the formula:
#'   \deqn{\log(i) \theta}{f(i) = log(i/N) \theta}
#'   }
#' 	 \item{\code{type = "stepT"}}{
#'    Represents step trend. Step trend means that the expected response of the patients increases
#'    by \code{theta} after a given point (\code{"saltus"}) in the allocation process.
#'    Step trend may occur if a new device is used after the point \eqn{c} = \code{"saltus"}, or if 
#'    the medical personal changes after after this point.
#'    Step time trend can be presented by the formula:
#'   \deqn{f(i) = 1_{c \leq i \leq N} \theta}{f(i) = 1_{c \le i\le N} \theta }
#'   }
#' }
#' }
#'
#' @return
#' \code{S4} object of class \code{chronBias}, a formal representation of the
#' issue of chronological bias in a clinical trial.
#'
#' @references
#' G. K. Rosenkranz (2011) The impact of randomization on the analysis of
#' clinical trials. \emph{Statistics in Medicine}, \strong{30}, 3475-87. 
#' 
#' M. Tamm and R.-D. Hilgers (2014) Chronological bias in randomized clinical 
#' trials under different types of unobserved time trends. 
#' \emph{Methods of Information in Medicine}, \strong{53}, 501-10. 
#' 
#' @family issues
#' 
#'
#' @export
chronBias <- function(type, theta, method, saltus, alpha = 0.05) {
  if(missing(saltus) && type %in% c("linT", "logT")) {
    new("chronBias", type = type, theta = theta, method = method,  alpha = alpha)
  } else {
    new("chronBiasStepT", type = type, theta = theta, method = method, saltus = saltus, alpha = alpha)
  }  
}


# --------------------------------------------
# Methods for chronBias
# --------------------------------------------

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "chronBias", endp = "missing"),
          function(randSeq, issue, endp) stop("Need an object of endpoint class."))

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "chronBias", endp = "normEndp"),
          function(randSeq, issue, endp) {
            validObject(randSeq); validObject(issue); validObject(endp)
            if (issue@method == "sim") {
              D <- data.frame(testDec = testDec(randSeq, issue, endp))
              colnames(D) <- paste("testDec(", issue@type, ")", sep = "")
              D
            } else {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("P(rej)(", issue@type, ")", sep = "")
              D
            }
          }
)


#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "chronBias", endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2, randSeq@K == length(endp@mu))
            validObject(randSeq); validObject(issue); validObject(endp)
            n <- N(randSeq)
            # linear time trend
            if (issue@type == "linT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                1:n * issue@theta
              }))
            }
			# logarithmic time trend			
			else if (issue@type == "logT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                log((1:n)/n) * issue@theta
              }))  
            } 
			# step time trend
			else if (issue@type == "stepT") {
              stopifnot(randSeq@N > issue@saltus)
              issue <- t(apply(randSeq@M, 1, function(x) {
                c(rep(0, issue@saltus - 1), rep(issue@theta, n - issue@saltus + 1))
              }))  
            }

            issue[randSeq@M == 0] <- issue[randSeq@M == 0] + endp@mu[1]
            issue[randSeq@M == 1] <- issue[randSeq@M == 1] + endp@mu[2]
            issue
          }
)


#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "chronBias", endp = "missing"),
          function(randSeq, issue) {
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(issue)
            n <- N(randSeq)
            # linear time trend
            if (issue@type == "linT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                1:n * issue@theta
              }))
            # logarithmic time trend			
            } else if (issue@type == "logT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                log((1:n)/n) * issue@theta
              })) 
            # step time trend
            } else if (issue@type == "stepT") {
              stopifnot(randSeq@N > issue@saltus, is(issue,"chronBiasStepT"))
              issue <- t(apply(randSeq@M, 1, function(x) {
                c(rep(0, issue@saltus - 1), rep(issue@theta, n - issue@saltus + 1))
              })) 
            }
            issue
          }
)

