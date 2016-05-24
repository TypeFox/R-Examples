#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class crSeq                                 #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for crSeq
# --------------------------------------------

# Representation of sequences for Complete Randomization
# 
# @description This set of classes provides functionality of storing 
# randomization sequences of Complete Randomization along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial
# @slot M matrix containing randomisation sequences of length \code{N} in its
# rows. 
setClass("crSeq", contains = "randSeq")

# --------------------------------------------
# Class definition for rCrSeq
# --------------------------------------------
# Representation of sequences for Complete Randomization
# 
# @description This set of classes provides functionality of storing random 
# randomization sequences of Complete Randomization along with the parameters 
# representing the design.
setClass("rCrSeq", contains = c("rRandSeq","crSeq"))

# --------------------------------------------
# Methods for crSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "crSeq"),
          function(obj) {
            validObject(obj)
            # every sequence equiprobable, generate same probability
            r <- nrow(obj@M)
            N <- N(obj)
            if(obj@K == 2) rep(2^{-N}, r)
            else "Only supported for K=2."
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "crSeq"),
          function(obj) {
            "CR"
          }
)

