#' @include getExpectation.R
#' @include getStat.R
#' @include issue.R
#' @include randSeq.R
#' @include util.R
#' @include endpoint.R
NULL

###############################################
# --------------------------------------------#
# Class combinedBias                          #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Class definition for combinedBias
# --------------------------------------------

# The combinedBias class
setClass("combinedBias", slots = c(eta = "numeric", typeSB = "character",
                   theta = "numeric", typeCB = "character",
                   method = "character", alpha = "numeric"))
setClass("combinedBiasStepTrend", slots = c(saltus = "numeric"),
         contains = "combinedBias")


#' Combined additive bias criterion
#'
#' This class combines a \code{selBias} object and a \code{chronBias} object
#' to a new object. In the analysis within the new object
#' the two types of bias are treated as additive.
#' effect.
#' 
#' @param selBias object of class \code{selBias}
#' @param chronBias object of class \code{chronBias}
#'
#' @family issues
#'
#' @examples
#' chronBias <- chronBias(type="linT", theta=1, method="sim")
#' selBias <- selBias(type="CS", eta=1, method="sim")
#' combineBias(selBias, chronBias)
#'
#' @export
combineBias <- function(selBias, chronBias) {
  stopifnot(is(selBias, "selBias"), is(chronBias, "chronBias"))
  if (selBias@alpha != chronBias@alpha) {
    warning("Parameter alpha taken from object selBias.")
  }
  if (selBias@method != chronBias@method) {
    warning("Parameter method taken from object selBias.")
  }
  if (!.hasSlot(chronBias, "saltus")) {
    new("combinedBias", eta = selBias@eta, typeSB = selBias@type, theta = chronBias@theta,
        method = selBias@method, alpha = selBias@alpha, typeCB = chronBias@type)
  } else {
    new("combinedBiasStepTrend", eta = selBias@eta, typeSB = selBias@type,
        theta = chronBias@theta,  method = selBias@method, 
        alpha = selBias@alpha, saltus = chronBias@saltus, typeCB = chronBias@type)
  }  
}


# --------------------------------------------
# Class union of bias and combinedBias
# --------------------------------------------

setClassUnion("bias", c("combinedBias", "combinedBiasStepTrend"))
setClassUnion("issue", c("combinedBias", "combinedBiasStepTrend")) 


# --------------------------------------------
# Methods for combinedBias
# --------------------------------------------

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "combinedBias",
                                      endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2, randSeq@K == length(endp@mu))
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            expectationSB <- getExpectation(randSeq, selBias)
            expectationCB <- getExpectation(randSeq, chronBias, endp)
            expectationSB + expectationCB
})

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "combinedBiasStepTrend",
                                      endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2, randSeq@K == length(endp@mu))
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@saltus,
                                   issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            expectationSB <- getExpectation(randSeq, selBias)
            expectationCB <- getExpectation(randSeq, chronBias, endp)
            expectationSB + expectationCB
})

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "combinedBias",
                               endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(validObject(randSeq), validObject(endp))
            if (issue@method == "sim") {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("testDec", " ", issue@method, "(combined)", sep = "")
              D
            } else {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("rejection prob.", " ", issue@method, "(combined)",
                                   sep = "")
              D
            }
          }
)



# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "combinedBiasStepTrend",
                               endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(validObject(randSeq), validObject(endp))
            if (issue@method == "sim") {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("testDec", " ", issue@method, "(combined)", sep = "")
              D
            } else {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("rejection prob.", " ", issue@method, "(combined)", sep = "")
              D
            }
          }
)
