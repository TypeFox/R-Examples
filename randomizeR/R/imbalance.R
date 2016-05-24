#' @include randSeq.R
#' @include endpoint.R
#' @include getStat.R
#' @include normEndp.R
NULL

###############################################
# --------------------------------------------#
# Class imbal                                 #
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
validateImbalance <- function(object) {
  errors <- character()
  lengthType <- length(object@type)
  if (lengthType != 1) {
    msg <- paste("type is length ", lengthType, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  type <- object@type[1]
  if (!(type %in% c("imb", "absImb", "loss", "maxImb"))) {
    msg <- paste("(First) Argument of type is ", type, ". Should be imb, absImb, or loss."
                 , sep = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for imbal
# --------------------------------------------

# Imbalance of allocation assignments in a clinical trial
#
# The \code{imbal} class provides a formal representation
# of the imbalance in a clinical trial. The balancing
# properties of a randomization procedure is one criterion
# for its assessment.
# 
# @description \code{imbal} is an S4 class inheriting 
# from \code{issue}. 
#
# @references
# A.C. Atkinson (2014) Selecting a biased coin design. \emph{Statistical Science}, \strong{29}, Vol. 1, 144-163.
#
setClass("imbal",
         slots = c("type" = "character"),
         validity = validateImbalance)


# --------------------------------------------
# Constructor function for imbal
# --------------------------------------------

#' Representing the allocation imbalance
#' 
#' @description
#' Represents the imbalance of the treatment assignments 
#' of patients in a clinical trial.
#' 
#' @family issues
#' 
#' @param type character string, should be one of \code{"imb"}, \code{"absImb"}, 
#' \code{"loss"}, or \code{"maxImb"}, see Details.
#'
#' @details 
#' Balance of the treatment assignment of patients can
#' be an issue in the design of a clinical trial. The \code{imbal} function is
#' a constructor function for an S4 object of class \code{imbal} representing
#' the issue of imbalance of a clinical trial. The parameter \code{type} can 
#' take the following values: 
#' The \code{type}
#' \describe{
#'   \item{\code{"imb"}}{the final imbalance, i.e. difference in group sizes at the end of a trial}
#'   \item{\code{"absImb"}}{the absolute value of the final imbalance}
#'	 \item{\code{"loss"}}{the loss in power estimation, i.e. \code{imb^2/N}}
#'   \item{\code{"maxImb"}}{the maximal attained imbalance during the trial}
#' }
#'
#' @return
#' \code{S4} object of class \code{imbal}, a formal represenation of the issue of
#' imbalance in a clinical trial.
#'
#' @references
#' A.C. Atkinson (2014) Selecting a biased coin design. \emph{Statistical Science}, \strong{29}, Vol. 1, 144-163. 
#' 
#' @export
imbal <- function(type) new("imbal", type = type)


# --------------------------------------------
# Methods for imbal
# --------------------------------------------

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "imbal", endp = "missing"),
          function(randSeq, issue) {
            validObject(randSeq); validObject(issue)
            D <- data.frame(imbal = imbSeq(randSeq, issue))
            colnames(D) <- issue@type
            D
          }
)

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "imbal", endp = "endpoint"),
          function(randSeq, issue, endp) {
            validObject(randSeq); validObject(issue)
            D <- data.frame(imbal = imbSeq(randSeq, issue))
            colnames(D) <- issue@type
            D
          }
)

# --------------------------------------------
# Helper functions for getStat and imbal
# --------------------------------------------

# Workhouse function for getStat for imbal
#
# Calculates the imbalance of a randomization list.
#
# @param randSeq object of the class randSeq.
# @param imbal object of the class imbal.
#
# @examples 
# myPar <- bsdPar(10, 2)
# M <- genSeq(myPar, 2)
# imb <- imbal("imb")
# imbSeq(M, imb)
#
# @name imbSeq
# 
# @return
# vector of the achieved imbalances of the sequences
imbSeq <- function(randSeq, imbal) {
  stopifnot(is(imbal, "imbal"), randSeq@K == 2,
            is(randSeq, "randSeq"))
  # Guessing matrix
  G <- matrix(numeric(0), ncol = dim(randSeq@M)[2], nrow = dim(randSeq@M)[1])
  # final imbal of the sequeces
  if (imbal@type == "imb") {
    apply(randSeq@M, 1, function(x) {
      imb <- sum(2*x - 1)
      return(imb)
    })
  } else if (imbal@type == "absImb") {
    apply(randSeq@M, 1, function(x) {
      imb <- abs(sum(2*x - 1))
      return(imb)
    })     
  } else if (imbal@type == "loss") {
    apply(randSeq@M, 1, function(x) {
      imb <- (sum(2*x - 1))^2/ randSeq@N
      return(imb)
    })      
  } else if (imbal@type == "maxImb") {
    apply(randSeq@M, 1, function(x) {
      imb <- max(abs(cumsum(2*x - 1)))
      return(imb)
    })      
  }       
}

