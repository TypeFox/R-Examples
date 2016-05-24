#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class bsdSeq                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for bsdSeq
# --------------------------------------------

# Representation of sequences for the Big Stick Design
# 
# @description This set of classes provides functionality of storing randomization
# sequences of the Big Stick Design along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot mti The maximum tolerated imbalance during the trial.
# @slot M matrix containing randomisation sequences of length \code{N} in its
# rows.
setClass("bsdSeq", slots = c(mti = "numeric"), contains = "randSeq")

# --------------------------------------------
# Class definition for rBsdSeq
# --------------------------------------------

# Representation of random sequences from the Big Stick Design
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of the Big Stick Design along with the parameters 
# representing the design.
setClass("rBsdSeq", contains = c("rRandSeq", "bsdSeq"))

# --------------------------------------------
# Methods for bsdSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "bsdSeq"),
          function(obj) {
            if(obj@K == 2) {
              # compute number of non-deterministic allocations nd
              # cut out last element, due to no deterministic allocation possible
              # after last allocation
              # add one for firs element which is always random
              apply(obj@M, 1, function(x, MTI) {
                nd <- sum(abs(cumsum(2*x - 1))[-length(x)] < mti(obj)) + 1
                0.5^nd
              }, MTI = mti(obj))
            }
            else "Only supported for K=2."
          }  
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "bsdSeq"),
          function(obj) {
            paste("BSD(", obj@mti, ")", sep = "")
          }
)

