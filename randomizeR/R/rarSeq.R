#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class rarSeq                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for rarSeq
# --------------------------------------------

# Representation of sequences for Random Allocation Rule
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Random Allocation Rule along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
setClass("rarSeq", contains = "randSeq")


# --------------------------------------------
# Class definition for rRarSeq
# --------------------------------------------

# Representation of sequences for Random Allocation Rule
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of Random Allocation Rule along with the parameters 
# representing the design.
setClass("rRarSeq", contains = c("rRandSeq", "rarSeq"))


# --------------------------------------------
# Methods for rarSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "rarSeq"),
          function(obj) {
            if(obj@K == 2) {
             # every sequence equally probable, generate same probability
             r <- nrow(obj@M)
             N <- N(obj)
             rep(1/choose(N, N/2), r)
           }
           else "Only supported for K=2."
          } 
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "rarSeq"),
          function(obj) {
            "RAR"
          }
)
