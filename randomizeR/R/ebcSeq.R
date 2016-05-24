#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class ebcSeq                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for ebcSeq
# --------------------------------------------

# Representation of sequences for Efron's Biased Coin Design
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Efron's Biased Coin Design along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot p success probability of the biased coin.
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
setClass("ebcSeq", slots=c(p = "numeric"), contains = "randSeq")


# --------------------------------------------
# Class definition for rEbcSeq
# --------------------------------------------

# Representation of sequences for Efron's Biased Coin Design
# 
# @description This set of classes provides functionality of storing random 
# randomization sequences of Efron's Biased Coin Design along with the parameters 
# representing the design.
setClass("rEbcSeq", contains = c("rRandSeq", "ebcSeq"))

# --------------------------------------------
# Methods for ebcSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "ebcSeq"),
          function(obj) {
            if(obj@K == 2) {
            apply(obj@M, 1, function(x, p){
              rw <- abs(c(0, cumsum(2*x - 1)))
              # hitting zero imbalance
              origin <- sum(rw[-length(rw)] == 0) 
              # reducing the imbalance
              together <- sum(rw[-length(rw)] > rw[-1]) 
              # increasing the imbalance
              apart <- sum((rw[-length(rw)] < rw[-1])*(rw[-length(rw)] > 0))
              0.5^origin * p^together * (1 - p)^apart
            }, p = coin(obj))
          } 
          else "Only supported for K=2."
          }  
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "ebcSeq"),
          function(obj) {
            paste("EBC(", round(obj@p, digits = 2), ")", sep = "")
          }
)

