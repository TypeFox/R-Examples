#' @include getStat.R
NULL

###############################################
# --------------------------------------------#
# Class power                                 #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the power class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validatePower <- function(object) {
  errors <- character()
  
  lengthMethod <- length(object@method)
  if (lengthMethod != 1) {
    msg <- paste("Method is length ", lengthMethod, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  method <- object@method[1]
  if (!(method %in% c("exact", "sim"))) {
    msg <- paste("(First) Argument of method is ", method, ". Should be \"exact\" or \"sim\"."
                 , sep = "")
    errors <- c(errors, msg)
  }
  
  lengthAlpha<- length(object@alpha)
  if (lengthAlpha != 1) {
    msg <- paste("Alpha is length ", lengthAlpha, ". Should be 1.", sep = "")
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
# Class definition for power
# --------------------------------------------

# The power class
setClass("power",
         slots = c("d" = "numeric", method = "character", alpha = "numeric"),
         validity = validatePower)


# --------------------------------------------
# Constructor function for power
# --------------------------------------------

#' Representing the power 
#' 
#' @description 
#' Represents the expected power of the individual randomization sequences.
#'
#' @family issues
#' 
#' @inheritParams overview 
#' @param method character string, should be one of \code{"sim"} or \code{"exact"}, see Description.
#' @param alpha significance level.
#'
#' @details
#' The attained power of an individual randomization sequence can be an issue
#' in the design of a clinical trial. The power of a randomization sequence is
#' is computed dependent on the effect size \code{d} and the difference in 
#' group sizes in the end if. 
#' 
#' If \code{method="sim"}, the object represents the simulated power of an
#' individual randomization sequence. When calling \code{assess} for a 
#' \code{power} object with \code{method="sim"}, one test decision is computed
#' for each randomization sequence of \code{randSeq}. The power is the 
#' proportion of falsely (correctly) rejected null hypotheses.
#' 
#' If \code{method="exact"}, the object represents the exact power of an 
#' individual randomization sequence. When calling \code{assess} for a 
#' \code{power} object with \code{method="exact"}, the exact \emph{p}-value
#' of each randomization sequence is computed. So far, this is only supported
#' for normal endpoints. Then the power is the sum of the corresponding 
#' quantiles of the noncentral t-distribution.
#' 
#' @return
#' \code{S4} object of class \code{power}, a formal representation of the
#' issue of power in a clinical trial.
#'
#' @export
setPower <- function(d, method, alpha = 0.05) {
  new("power", d = d, method = method, alpha = alpha)
}

# --------------------------------------------
# Generic function for power
# --------------------------------------------

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "power", 
                                      endp = "normEndp"), 
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == length(endp@mu))
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(normEndp); validObject(issue)
            expectation <- matrix(numeric(0), ncol = ncol(randSeq@M), 
                                  nrow = nrow(randSeq@M))
            for(i in 0:1) {
              expectation[randSeq@M == i] <- endp@mu[i+1]
            }  
            
            expectation[randSeq@M == 1] <- expectation[randSeq@M == 1] + issue@d

            expectation
          }
)

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "power", endp = "normEndp"),
          function(randSeq, issue, endp) {
            validObject(randSeq); validObject(issue); validObject(endp)
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("power", "(", issue@method, ")", sep = "")
              D
          }
)



