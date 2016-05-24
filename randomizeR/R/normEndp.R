#' @include getExpectation.R
NULL

###############################################
# --------------------------------------------#
# Class normEndp                              #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the endpoint class
#
# @inheritParams overview
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateNormEndp <- function(object) {
  errors <- character() 
  if (!(all(object@sigma > 0))) {
    msg <- ("The standard deviation must be positive.")
    errors <- c(errors, msg)
  }
  
  if (!(length(object@sigma) == length(object@mu))) {
    msg <- ("sigma and mu should have same length.")
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for normEndp
# --------------------------------------------

# Representation of the normal endpoints
setClass("normEndp", 
         slots = c(mu = "numeric", sigma = "numeric"),
         validity = validateNormEndp)



# --------------------------------------------
# Constructor function for normEndp
# --------------------------------------------

#' Representation of normally distributed endpoints
#' 
#' Represents normally distributed endpoints in clinical trials.
#'
#' @inheritParams overview
#'
#' @details
#' The \code{normEnd} function is a constructor function
#' for an S4 object of the class \code{normEnd} representing 
#' a normally distributed endpoint in a clinical trial.
#' In conjunction with the assess function, normal endpoints
#' admit the calculation of the exact type-I-error probability and power.
#'
#' @family endopoint types
#'
#' @seealso Compute exact or simulated type-I-error: \code{\link{assess}}.
#' 
#' @export
normEndp <- function(mu, sigma) {
  new("normEndp", mu = mu, sigma = sigma)
}


# --------------------------------------------
# Generic function for normEndp
# --------------------------------------------

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "missing", 
                                      endp = "normEndp"), 
          function(randSeq, endp) {
            stopifnot(randSeq@K == length(endp@mu))
            validObject(randSeq); validObject(normEndp)
            expectation <- matrix(numeric(0), ncol = ncol(randSeq@M), 
                      nrow = nrow(randSeq@M))
            for(i in 0:(randSeq@K-1)) {
              expectation[randSeq@M == i] <- endp@mu[i+1]
            }  
            expectation
          }
)

