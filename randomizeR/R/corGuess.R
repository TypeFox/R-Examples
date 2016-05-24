#' @include getStat.R
NULL

###############################################
# --------------------------------------------#
# Class corGuess                              #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the selBias class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validatecG <- function(object) {
  errors <- character()
  lengthType <- length(object@type)
  if (lengthType != 1) {
    msg <- paste("Type is length ", lengthType, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  type <- object@type[1]
  if (!(type %in% c("CS", "DS"))) {
    msg <- paste("(First) Argument of type is ", type, ". Should be CS or DS."
                 , sep = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for corGuess
# --------------------------------------------

# The corGuess class
setClass("corGuess",
         slots = c("type" = "character"),
         validity = validatecG)


# --------------------------------------------
# Constructor function for corGuess
# --------------------------------------------

#' Representing the expected number of correct guesses
#' 
#' @description 
#' Represents the expected number of correct guesses of randomization sequences.
#'
#' @family issues
#' 
#' @param type character string, should be one of \code{"CS"} or \code{"DS"},
#' see Details.
#'
#' @details
#' Selection bias can be an issue in the design of a clinical trial. The 
#' expected number of correct guesses is one measure for selection bias.
#' The \code{corGuess} function is a constructor function
#' for an S4 object of the class \code{corGuess} representing the issue of
#' correct guesses in a clinical trial. The parameter \code{type} takes the 
#' following values:
#' \describe{
#'   \item{\code{"CS"}}{refers to "convergence strategy", i.e. the investigator
#'   predicts the treatment which has hitherto occured less often.}
#'   \item{\code{"DS"}}{refers to "divergence strategy", i.e. the investigator
#'   predicts the treatment which has hitherto occured more often.}
#' }
#' 
#' @return
#' \code{S4} object of class \code{corGuess}, a formal representation of the
#' issue of correct guesses in a clinical trial.
#'
#' @references
#' D. Blackwell and J.L. Hodges Jr. (1957) Design for the control of 
#' selection bias. \emph{Annals of Mathematical Statistics}, \strong{25}, 449-60. 
#'
#' @export
corGuess <- function(type) new("corGuess", type = type)


# --------------------------------------------
# Methods for corGuess
# --------------------------------------------

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "corGuess", endp = "missing"),
          function(randSeq, issue) {
            validObject(randSeq); validObject(issue)
            G <- getCorGuesses(randSeq, issue)
            rL <- getRandList(randSeq)
            res <- (apply(G == rL, 1, sum) + apply(G == "nG", 1, sum)/2)/randSeq@N
            D <- data.frame(res)
            colnames(D) <- paste("propCG(", issue@type, ")", sep = "")
            D
          }
)

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "corGuess", endp = "endpoint"),
          function(randSeq, issue, endp) {
            validObject(randSeq); validObject(issue)
            G <- getCorGuesses(randSeq, issue)
            rL <- getRandList(randSeq)
            res <- (apply(G == rL, 1, sum) + apply(G == "nG", 1, sum)/2)/randSeq@N
            D <- data.frame(res)
            colnames(D) <- paste("propCG(", issue@type, ")", sep = "")
            D
          }
)


# --------------------------------------------
# Helper functions for corGuess
# --------------------------------------------

#' Matrix of the guesses of the investigator
#'
#' Calculates the guesses of the investigator of a randomization list following
#' the specified guessing strategy.
#'
#' @param randSeq object of the class randSeq.
#' @param guessing object of the class corGuess.
#'
#' @examples 
#' myPar <- bsdPar(10, 2)
#' M <- genSeq(myPar, 2)
#' type <- corGuess("CS")
#' getCorGuesses(M, type)
#' 
#' @return
#' Matrix of the guesses of the investigator following the specified guessing
#'  strategy. No guess is abbreviated with \code{"nG"}.
#'
#' @export
getCorGuesses <- function(randSeq, guessing) {
    stopifnot(is(guessing, "corGuess"), randSeq@K == 2,
            is(randSeq, "randSeq"), all(randSeq@groups != "nG"),
             validObject(randSeq), validObject(guessing))
    # Guessing matrix
    G <- matrix(numeric(0), ncol = dim(randSeq@M)[2], nrow = dim(randSeq@M)[1])
    # convergence strategey
    if (guessing@type == "CS") {
      t(apply(randSeq@M, 1, function(x) {
        sImb <- sign(cumsum(2*x - 1))
        sImb[sImb == 0] <- "nG"
        sImb[sImb == 1] <- randSeq@groups[1]
        sImb[sImb == -1] <- randSeq@groups[2]
        sImb <- sImb[-length(sImb)]
        sImb <- c("nG", sImb)
        return(sImb)
      }))
    } else if (guessing@type == "DS") {
      t(apply(randSeq@M, 1, function(x) {
        sImb <- sign(cumsum(2*x - 1))
        sImb[sImb == 0] <- "nG"
        sImb[sImb == 1] <- randSeq@groups[2]
        sImb[sImb == -1] <- randSeq@groups[1]
        sImb <- sImb[-length(sImb)]
        sImb <- c("nG", sImb)
        return(sImb)
      }))          
    }    
}

